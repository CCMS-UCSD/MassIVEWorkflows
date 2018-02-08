#!/usr/bin/python

import sys
import getopt
import os
import json
import xmltodict
import subprocess
import requests
import ming_proteosafe_library
import ming_fileio_library
import time
import uuid
import shutil
from collections import defaultdict

SERVER_BASE_URL = "proteomics2.ucsd.edu"
PATH_TO_DATASETS = "/data/massive"

def get_mzTab_parameter(mztab_parameter, scratch_dir, path_to_front_end_tasks):
    print mztab_parameter
    target_tabname = os.path.basename(mztab_parameter["tab_name"])
    target_task_id = mztab_parameter["task"]

    target_task_object = ming_proteosafe_library.get_task_information(SERVER_BASE_URL, target_task_id)

    target_task_user = target_task_object["user"]
    target_task_status = target_task_object["status"]
    target_task_workflow = target_task_object["workflow"]

    print(target_task_workflow)

    target_task_mztab_folder = "mzTab"

    target_spectrum_mztab_file = ""

    # if target_task_workflow == "MASSIVE" or target_task_workflow == "MASSIVE-COMPLETE":
    #     #target_task_mztab_folder = "result"
    #     dataset_information = ming_proteosafe_library.get_dataset_information(target_task_id)
    #     print(json.dumps(dataset_information, indent = 4))
    #     target_spectrum_mztab_file = os.path.join(PATH_TO_DATASETS, dataset_information["dataset_id"], mztab_parameter["tab_name"])
    # elif target_task_workflow.find("ADD-MASSIVE-REANALYSIS") != -1:
    #     target_task_mztab_folder = "result"
    #     target_spectrum_mztab_file = os.path.join(PATH_TO_DATASETS, dataset_information["dataset_id"], mztab_parameter["tab_name"])
    # elif target_task_workflow.find("CONVERT-TSV") != -1:
    #     target_task_mztab_folder = "result"
    #     target_spectrum_mztab_file = os.path.join(path_to_front_end_tasks, target_task_user, target_task_id, target_task_mztab_folder, target_tabname)
    # else:
    #     target_spectrum_mztab_file = os.path.join(path_to_front_end_tasks, target_task_user, target_task_id, target_task_mztab_folder, target_tabname)

    #Using File Descriptor
    file_descriptor_path = mztab_parameter["tab_file_descriptor"]
    if file_descriptor_path[0] == "f":
        target_spectrum_mztab_file = os.path.join("/data/ccms-data/uploads", file_descriptor_path[2:])
    if file_descriptor_path[0] == "d":
        target_spectrum_mztab_file = os.path.join("/data/ccms-data/uploads", file_descriptor_path[2:])
    if file_descriptor_path[0] == "u":
        target_spectrum_mztab_file = os.path.join("/data/ccms-data/tasks", file_descriptor_path[2:])



    #copy mztab to location
    print(target_spectrum_mztab_file)
    mztab_working_copy = os.path.join(scratch_dir, str(uuid.uuid4()) + "-" + target_tabname)
    shutil.copy(target_spectrum_mztab_file, mztab_working_copy)

    #lets only keep the PSH and PSM lines
    scratch_file = mztab_working_copy + "-onlypsm"
    filtered_file = open(scratch_file, "w")

    for line in open(target_spectrum_mztab_file):
        prefix = line[:3]
        if prefix == "PSM" or prefix == "PSH":
            filtered_file.write(line)

    filtered_file.close()

    #open as a tsv file
    row_count, table_data = ming_fileio_library.parse_table_with_headers(scratch_file)

    all_peptides = table_data["sequence"]

    return all_peptides


def parse_task_list(all_tabs, scratch_dir, task_id, user, path_to_front_end_tasks):
    all_peptides = []
    for tab_parameter in all_tabs:
        tab_peptides = get_mzTab_parameter(tab_parameter, scratch_dir, path_to_front_end_tasks)

        all_peptides += tab_peptides

    all_peptides = list(set(all_peptides))

    return all_peptides


def usage():
    print "<input parameters> <output peptide list>"


def main():
    usage()

    input_parameters = sys.argv[1]
    output_peptide_list = sys.argv[2]
    path_to_front_end_tasks = sys.argv[3]
    path_to_scratch_space = "./scratch"

    ming_fileio_library.make_sure_path_exists(path_to_scratch_space)

    #Parsing the parameters
    parsed_parameters = ming_proteosafe_library.parse_xml_file(open(input_parameters))

    #Execute Whatever you need here
    left_parameters = parsed_parameters["left_objects"][0]
    right_parameters = parsed_parameters["right_objects"][0]
    task_id = parsed_parameters["task"][0]
    user = parsed_parameters["user"][0]

    left_peptides = parse_task_list(json.loads(left_parameters), path_to_scratch_space, task_id, user, path_to_front_end_tasks)
    right_peptides = parse_task_list(json.loads(right_parameters), path_to_scratch_space, task_id, user, path_to_front_end_tasks)

    left_peptides_map = defaultdict(lambda: 0)
    for peptide in left_peptides:
        left_peptides_map[peptide] += 1

    right_peptides_map = defaultdict(lambda: 0)
    for peptide in right_peptides:
        right_peptides_map[peptide] += 1

    #Dummy output
    output_file_map = {}
    output_file_map["left_peptides"] = left_peptides_map
    output_file_map["right_peptides"] = right_peptides_map
    open(output_peptide_list, "w").write(json.dumps(output_file_map))



if __name__ == "__main__":
    main()
