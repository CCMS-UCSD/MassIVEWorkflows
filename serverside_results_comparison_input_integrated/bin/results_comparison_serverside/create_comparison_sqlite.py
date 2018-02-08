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
from collections import defaultdict

SERVER_BASE_URL = "proteomics2.ucsd.edu"
SERVER_URL = "http://" + SERVER_BASE_URL + "/ProteoSAFe"

workflow_server_url_map = {}
workflow_server_url_map["MASSIVE-COMPLETE"] = "massive.ucsd.edu"
workflow_server_url_map["MASSIVE"] = "massive.ucsd.edu"
workflow_server_url_map["ADD-MASSIVE-REANALYSIS"] = "massive.ucsd.edu"
workflow_server_url_map["CONVERT-TSV"] = "massive.ucsd.edu"
workflow_server_url_map["MSGFDB"] = "proteomics2.ucsd.edu"
workflow_server_url_map["MSGF_PLUS"] = "proteomics2.ucsd.edu"
workflow_server_url_map["MODA"] = "proteomics2.ucsd.edu"
workflow_server_url_map["MSPLIT_NEW"] = "proteomics2.ucsd.edu"
workflow_server_url_map["MSPLITNEW_MSGFDB_MODA_2PASS"] = "proteomics2.ucsd.edu"
workflow_server_url_map["MSPLITNEW_MSGFPLUS_MODA_2PASS"] = "proteomics2.ucsd.edu"

def usage():
    print "<input parameters> <output log file> <output folder> <path to front end tasks>"

#Lets also keep hitting this until we get a 200
def get_mztab_results_data(task_id, view_name, tab_file_descriptor, workflow):
    url = SERVER_URL + "/result_json.jsp"

    if workflow in workflow_server_url_map:
        url = "http://" + workflow_server_url_map[workflow] + "/ProteoSAFe" + "/result_json.jsp"

    payload = {'task': task_id, 'view': view_name, 'file' : tab_file_descriptor}

    print url
    print payload
    r = requests.get(url, params=payload)
    status_code = r.status_code

    while(status_code != 200):
        print("retry", payload)
        r = requests.get(url, params=payload)
        status_code = r.status_code

    return json.loads(r.text)

def wait_for_sqlite_file(path, retry_count, sleep_time):
    if os.path.isfile(path) == False:
        if retry_count == 0:
            return -1
        time.sleep(sleep_time)
        return wait_for_sqlite_file(path, retry_count - 1, sleep_time)
    else:
        return 0


#This takes in a tab file
#Copies the remote sqlite file locally
#Keeps only necessary columns
#Adds extra columns
#Creates a dump
#Dumps data into merged sqlite file
def process_tab_fileinto_unified_sqlite(view_name, task_parameters, compare_partition, merged_db_filename, execute_file, remove_columns_sqlite_filename, path_to_sqlite3, path_to_front_end_tasks, path_to_scratch_space):
    target_task_id = task_parameters["task"]
    target_tabname = os.path.basename(task_parameters["tab_name"])
    full_target_tabname = task_parameters["tab_name"]
    tabname_file_descriptor = task_parameters["tab_file_descriptor"]

    target_task_object = ming_proteosafe_library.get_task_information(SERVER_BASE_URL, target_task_id)
    target_task_user = target_task_object["user"]
    target_task_status = target_task_object["status"]
    target_task_workflow = target_task_object["workflow"]


    if target_task_status != "DONE":
        print("Task: %s is not DONE" % (target_task_id) )
        exit(1)

    #target_spectrum_sqlite_file = os.path.join(path_to_front_end_tasks, target_task_user, target_task_id, "sqlite", view_name + "-main_" + target_tabname.split(".")[0] + ".db")
    #local_spectrum_sqlite_file = os.path.join(path_to_scratch_space, target_task_id + "-" +  view_name + "-main_" + target_tabname.split(".")[0] + ".db")

    target_spectrum_sqlite_file = ""
    local_spectrum_sqlite_file = ""

    if os.path.isfile(target_spectrum_sqlite_file) == False:
        result_json = get_mztab_results_data(target_task_id, view_name, tabname_file_descriptor, target_task_workflow)
        print(result_json)
        target_spectrum_sqlite_file = os.path.join(path_to_front_end_tasks, target_task_user, target_task_id, "sqlite", result_json["blockData"]["file"])
        local_spectrum_sqlite_file = os.path.join(path_to_scratch_space, target_task_id + "-" + result_json["blockData"]["file"])

        time.sleep(30)

        if wait_for_sqlite_file(target_spectrum_sqlite_file, 10, 30) == -1:
            print(target_spectrum_sqlite_file)
            print "SQLITE FILE MISSING SPECTRUM"
            exit(1)

    execute_file.write("cp " + target_spectrum_sqlite_file + " " + local_spectrum_sqlite_file + "\n")

    if remove_columns_sqlite_filename != None:
        execute_file.write("cat " + remove_columns_sqlite_filename + " | " + path_to_sqlite3 + " " + local_spectrum_sqlite_file + "\n")

    execute_file.write("echo \"alter table Result add column ComparePartition VARCHAR(16);\" | " + path_to_sqlite3 + " "  + local_spectrum_sqlite_file + "\n");
    #Writing number in partition
    execute_file.write("echo \"update Result SET ComparePartition=" + "\\\"" + compare_partition + "\\\"" + ";\" | " + path_to_sqlite3 + " "  + local_spectrum_sqlite_file + "\n");

    #Adding Task Identifier
    execute_file.write("echo \"alter table Result add column taskid TEXT;\" | " + path_to_sqlite3 + " "  + local_spectrum_sqlite_file + "\n");
    #Writing Task
    execute_file.write("echo \"update Result SET taskid=\\\"" + target_task_id + "\\\";\" | " + path_to_sqlite3 + " "  + local_spectrum_sqlite_file + "\n");

    #Adding Tab Identifier
    execute_file.write("echo \"alter table Result add column tabid TEXT;\" | " + path_to_sqlite3 + " "  + local_spectrum_sqlite_file + "\n");
    #Writing Tab
    execute_file.write("echo \"update Result SET tabid=\\\"" + target_tabname + "\\\";\" | " + path_to_sqlite3 + " "  + local_spectrum_sqlite_file + "\n");

    #Creating the db dump
    local_dump_file = local_spectrum_sqlite_file + ".dump";
    execute_file.write("" + path_to_sqlite3 + " " + local_spectrum_sqlite_file + " .dump > " + local_dump_file + "\n");

    #Merging it
    execute_file.write("cat " + local_dump_file + " | " + path_to_sqlite3 + " " + merged_db_filename + "\n");

def add_peptide_mapping_table(peptide_mapping, sqlite_file, scratch_folder, execute_file, path_to_sqlite3):
    dump_file_location = os.path.join(scratch_folder, "adding_peptide_mappings.dump")

    dump_file = open(dump_file_location, "w")

    #Writing dump file
    create_lines = []
    create_lines.append("PRAGMA foreign_keys=OFF;")
    create_lines.append("BEGIN TRANSACTION;")
    create_lines.append("CREATE TABLE PeptideMapping")
    create_lines.append("('sequence' TEXT,'proteins' TEXT);")
    dump_file.write("\n".join(create_lines))
    dump_file.write("\n")

    for peptide in peptide_mapping:
        peptide_mapping_object = peptide_mapping[peptide]
        proteins_string = ";".join(peptide_mapping_object["proteins"])
        insert_statement = "INSERT INTO \"PeptideMapping\" VALUES('%s', '%s');" % (peptide, proteins_string)
        dump_file.write(insert_statement + "\n")

    dump_file.write("COMMIT;")
    dump_file.close()

    execute_file.write("cat " + dump_file_location + " | " + path_to_sqlite3 + " " + sqlite_file + "\n");

def create_comparison_sqlite_spectrum_level(left_comparison, right_comparison, task_id, user, path_to_scratch_space, output_folder, path_to_front_end_tasks, path_to_sqlite3, peptide_remapping_map, ignore_fileprefix_extension):
    left_compare_obj = json.loads(left_comparison)
    right_compare_obj = json.loads(right_comparison)

    #Ultimate output files
    spectrum_level_merged_db_filename = os.path.join(output_folder, "spectrum_level.db")

    execute_filename = os.path.join(path_to_scratch_space, "execute.sh")
    execute_file = open(execute_filename, "w")

    #SQLITE Code to remove columns
    remove_columns_sqlite_filename = os.path.join(path_to_scratch_space, "remove_columns.sqlite")
    remove_columns_sqlite_file = open(remove_columns_sqlite_filename, "w")
    remove_columns_sqlite_file.write("""BEGIN TRANSACTION;
            CREATE TEMPORARY TABLE t1_backup(sequence, sequence_li, accession, modifications, nativeID, [#SpecFile], modified_sequence, modified_sequence_li, nativeID_index, nativeID_scan, charge, baseFilename);
            INSERT INTO t1_backup SELECT sequence, sequence_li, accession, modifications, nativeID, [#SpecFile], modified_sequence, modified_sequence_li, nativeID_index, nativeID_scan, charge, baseFilename FROM Result;
            DROP TABLE Result;
            CREATE TABLE Result(sequence, sequence_li, accession, modifications, nativeID, [#SpecFile], modified_sequence, modified_sequence_li, nativeID_index, nativeID_scan, charge, baseFilename);
            INSERT INTO Result SELECT sequence, sequence_li, accession, modifications, nativeID, [#SpecFile], modified_sequence, modified_sequence_li, nativeID_index, nativeID_scan, charge, baseFilename FROM t1_backup;
            DROP TABLE t1_backup;
            COMMIT;""")
    remove_columns_sqlite_file.close()

    left_tasks = json.loads(left_comparison)
    right_tasks = json.loads(right_comparison)

    for task in left_tasks:
        process_tab_fileinto_unified_sqlite("group_by_spectrum", task, "A", spectrum_level_merged_db_filename, execute_file, remove_columns_sqlite_filename, path_to_sqlite3, path_to_front_end_tasks, path_to_scratch_space)

    for task in right_tasks:
        process_tab_fileinto_unified_sqlite("group_by_spectrum", task, "B", spectrum_level_merged_db_filename, execute_file, remove_columns_sqlite_filename, path_to_sqlite3, path_to_front_end_tasks, path_to_scratch_space)

    #Reducing redundancy of the PSMs by grouping by file, scan, and unmodified sequence
    reduce_site_localization_redundancy_remapping_queries = "echo \"CREATE TABLE Result_Filtered AS SELECT * FROM Result  GROUP BY  sequence_li, nativeID, #SpecFile, ComparePartition;\" | " + path_to_sqlite3 + " " + spectrum_level_merged_db_filename + "\n"
    reduce_site_localization_redundancy_remapping_queries += "echo \"DROP TABLE Result;\" | " + path_to_sqlite3 + " " + spectrum_level_merged_db_filename + "\n"
    reduce_site_localization_redundancy_remapping_queries += "echo \"ALTER TABLE Result_Filtered RENAME TO Result;\" | " + path_to_sqlite3 + " " + spectrum_level_merged_db_filename + "\n"

    execute_file.write(reduce_site_localization_redundancy_remapping_queries)


    #Creating a table in the database that includes all the mappings from peptide to protein
    add_peptide_mapping_table(peptide_remapping_map, spectrum_level_merged_db_filename, path_to_scratch_space, execute_file, path_to_sqlite3)

    precompute_remapping_queries = "echo \"ALTER TABLE Result ADD remapped_accession TEXT;\" | " + path_to_sqlite3 + " " + spectrum_level_merged_db_filename + "\n"
    precompute_remapping_queries += "echo \"CREATE INDEX peptide_mapping_sequence_index ON PeptideMapping (sequence);\" | " + path_to_sqlite3 + " " + spectrum_level_merged_db_filename  + "\n"
    precompute_remapping_queries += "echo \"CREATE INDEX result_sequence_index ON Result (sequence);\" | " + path_to_sqlite3 + " " + spectrum_level_merged_db_filename  + "\n"
    precompute_remapping_queries += "echo \"UPDATE Result SET remapped_accession=(SELECT PeptideMapping.proteins FROM PeptideMapping WHERE PeptideMapping.sequence=Result.sequence);\" | " + path_to_sqlite3 + " " + spectrum_level_merged_db_filename  + "\n"





    precompute_mismatch_query = ""
    precompute_match_query = ""
    precompute_first_unique_query = ""
    precompute_second_unique_query = ""
    if ignore_fileprefix_extension:
        precompute_mismatch_query = "echo \"CREATE TABLE mismatched_rows AS SELECT * FROM Result as first JOIN Result as second ON (first.nativeID = second.nativeID) AND (first.baseFilename = second.baseFilename) AND (first.ComparePartition <> second.ComparePartition) AND (first.sequence_li <> second.sequence_li) WHERE first.ComparePartition='A';\" | " + path_to_sqlite3 + " " + spectrum_level_merged_db_filename;
        precompute_match_query = "echo \"CREATE TABLE matched_rows AS SELECT * FROM Result as first JOIN Result as second ON (first.nativeID = second.nativeID) AND (first.baseFilename = second.baseFilename) AND (first.ComparePartition <> second.ComparePartition) AND (first.sequence_li = second.sequence_li) WHERE first.ComparePartition='A';\" | " + path_to_sqlite3 + " " + spectrum_level_merged_db_filename;
        precompute_first_unique_query = "echo \"CREATE TABLE unique_first AS SELECT * FROM Result as first LEFT JOIN Result as second ON (first.nativeID = second.nativeID) AND (first.baseFilename = second.baseFilename) AND (first.ComparePartition <> second.ComparePartition) WHERE first.ComparePartition='A' AND second.[#SpecFile] is NULL;\" | " + path_to_sqlite3 + " " + spectrum_level_merged_db_filename;
        precompute_second_unique_query = "echo \"CREATE TABLE unique_second AS SELECT * FROM Result as first LEFT JOIN Result as second ON (first.nativeID = second.nativeID) AND (first.baseFilename = second.baseFilename) AND (first.ComparePartition <> second.ComparePartition) WHERE first.ComparePartition='B' AND second.[#SpecFile] is NULL;\" | " + path_to_sqlite3 + " " + spectrum_level_merged_db_filename;
    else:
        precompute_mismatch_query = "echo \"CREATE TABLE mismatched_rows AS SELECT * FROM Result as first JOIN Result as second ON (first.nativeID = second.nativeID) AND (first.[#SpecFile] = second.[#SpecFile]) AND (first.ComparePartition <> second.ComparePartition) AND (first.sequence_li <> second.sequence_li) WHERE first.ComparePartition='A';\" | " + path_to_sqlite3 + " " + spectrum_level_merged_db_filename;
        precompute_match_query = "echo \"CREATE TABLE matched_rows AS SELECT * FROM Result as first JOIN Result as second ON (first.nativeID = second.nativeID) AND (first.[#SpecFile] = second.[#SpecFile]) AND (first.ComparePartition <> second.ComparePartition) AND (first.sequence_li = second.sequence_li) WHERE first.ComparePartition='A';\" | " + path_to_sqlite3 + " " + spectrum_level_merged_db_filename;
        precompute_first_unique_query = "echo \"CREATE TABLE unique_first AS SELECT * FROM Result as first LEFT JOIN Result as second ON (first.nativeID = second.nativeID) AND (first.[#SpecFile] = second.[#SpecFile]) AND (first.ComparePartition <> second.ComparePartition) WHERE first.ComparePartition='A' AND second.[#SpecFile] is NULL;\" | " + path_to_sqlite3 + " " + spectrum_level_merged_db_filename;
        precompute_second_unique_query = "echo \"CREATE TABLE unique_second AS SELECT * FROM Result as first LEFT JOIN Result as second ON (first.nativeID = second.nativeID) AND (first.[#SpecFile] = second.[#SpecFile]) AND (first.ComparePartition <> second.ComparePartition) WHERE first.ComparePartition='B' AND second.[#SpecFile] is NULL;\" | " + path_to_sqlite3 + " " + spectrum_level_merged_db_filename;


    precompute_matching_peptide_sequence_in_mismatched_rows = "echo \"ALTER TABLE mismatched_rows ADD matching_sequence INT;\" | " + path_to_sqlite3 + " " + spectrum_level_merged_db_filename + "\n"
    precompute_matching_peptide_sequence_in_mismatched_rows += "echo \"UPDATE mismatched_rows SET matching_sequence=CASE WHEN sequence_li=\`sequence_li:1\` THEN 1 ELSE 0 END;\" | " + path_to_sqlite3 + " " + spectrum_level_merged_db_filename + "\n"


    execute_file.write(precompute_remapping_queries)
    execute_file.write(precompute_mismatch_query + "\n")
    execute_file.write(precompute_match_query + "\n")
    execute_file.write(precompute_first_unique_query + "\n")
    execute_file.write(precompute_second_unique_query + "\n")
    execute_file.write(precompute_matching_peptide_sequence_in_mismatched_rows)

    execute_file.close()

    #Running the execute scripts
    subprocess.call("sh " + execute_filename, shell=True)

def create_comparison_sqlite_peptide_level(left_comparison, right_comparison, task_id, user, path_to_scratch_space, output_folder, path_to_front_end_tasks, path_to_sqlite3, peptide_remapping_map):
    left_compare_obj = json.loads(left_comparison)
    right_compare_obj = json.loads(right_comparison)

    #Ultimate output files
    peptide_level_merged_db_filename = os.path.join(output_folder, "peptide_level.db")

    execute_filename = os.path.join(path_to_scratch_space, "execute.sh")
    execute_file = open(execute_filename, "w")

    #SQLITE Code to remove columns
    remove_columns_sqlite_filename = os.path.join(path_to_scratch_space, "remove_columns.sqlite")
    remove_columns_sqlite_file = open(remove_columns_sqlite_filename, "w")
    remove_columns_sqlite_file.write("""BEGIN TRANSACTION;
            CREATE TEMPORARY TABLE t1_backup(sequence, sequence_li, accession, modifications, modified_sequence, modified_sequence_li, charge);
            INSERT INTO t1_backup SELECT sequence, sequence_li, accession, modifications, modified_sequence, modified_sequence_li, charge FROM Result;
            DROP TABLE Result;
            CREATE TABLE Result(sequence, sequence_li, accession, modifications, modified_sequence, modified_sequence_li, charge);
            INSERT INTO Result SELECT sequence, sequence_li, accession, modifications, modified_sequence, modified_sequence_li, charge FROM t1_backup;
            DROP TABLE t1_backup;
            COMMIT;""")
    remove_columns_sqlite_file.close()

    left_tasks = json.loads(left_comparison)
    right_tasks = json.loads(right_comparison)

    for task in left_tasks:
        process_tab_fileinto_unified_sqlite("group_by_peptide_derived", task, "A", peptide_level_merged_db_filename, execute_file, remove_columns_sqlite_filename, path_to_sqlite3, path_to_front_end_tasks, path_to_scratch_space)


    for task in right_tasks:
        process_tab_fileinto_unified_sqlite("group_by_peptide_derived", task, "B", peptide_level_merged_db_filename, execute_file, remove_columns_sqlite_filename, path_to_sqlite3, path_to_front_end_tasks, path_to_scratch_space)


    #Creating a table in the database that includes all the mappings from peptide to protein
    add_peptide_mapping_table(peptide_remapping_map, peptide_level_merged_db_filename, path_to_scratch_space, execute_file, path_to_sqlite3)

    precompute_remapping_queries = "echo \"ALTER TABLE Result ADD remapped_accession TEXT;\" | " + path_to_sqlite3 + " " + peptide_level_merged_db_filename + "\n"
    precompute_remapping_queries += "echo \"CREATE INDEX peptide_mapping_sequence_index ON PeptideMapping (sequence);\" | " + path_to_sqlite3 + " " + peptide_level_merged_db_filename  + "\n"
    precompute_remapping_queries += "echo \"CREATE INDEX result_sequence_index ON Result (sequence);\" | " + path_to_sqlite3 + " " + peptide_level_merged_db_filename  + "\n"
    precompute_remapping_queries += "echo \"UPDATE Result SET remapped_accession=(SELECT PeptideMapping.proteins FROM PeptideMapping WHERE PeptideMapping.sequence=Result.sequence);\" | " + path_to_sqlite3 + " " + peptide_level_merged_db_filename  + "\n"


    #Making things unique, and then moving tables back.
    get_unique_peptides = "echo \"CREATE TABLE Result_Unique AS SELECT * FROM Result GROUP BY modified_sequence,ComparePartition;\" | " + path_to_sqlite3 + " " + peptide_level_merged_db_filename  + "\n"
    get_unique_peptides += "echo \"DROP TABLE Result;\" | " + path_to_sqlite3 + " " + peptide_level_merged_db_filename  + "\n"
    get_unique_peptides += "echo \"ALTER TABLE Result_Unique RENAME TO Result;\" | " + path_to_sqlite3 + " " + peptide_level_merged_db_filename  + "\n"


    precompute_match_query = "echo \"CREATE TABLE matched_rows AS SELECT * FROM Result as first JOIN Result as second ON (first.ComparePartition <> second.ComparePartition) AND (first.modified_sequence_li = second.modified_sequence_li) WHERE first.ComparePartition='A';\" | " + path_to_sqlite3 + " " + peptide_level_merged_db_filename;
    precompute_first_unique_query = "echo \"CREATE TABLE unique_first AS SELECT * FROM Result as first LEFT JOIN Result as second ON (first.modified_sequence = second.modified_sequence) AND (first.ComparePartition <> second.ComparePartition) WHERE first.ComparePartition='A' AND second.modified_sequence is NULL;\" | " + path_to_sqlite3 + " " + peptide_level_merged_db_filename;
    precompute_second_unique_query = "echo \"CREATE TABLE unique_second AS SELECT * FROM Result as first LEFT JOIN Result as second ON (first.modified_sequence = second.modified_sequence) AND (first.ComparePartition <> second.ComparePartition) WHERE first.ComparePartition='B' AND second.modified_sequence is NULL;\" | " + path_to_sqlite3 + " " + peptide_level_merged_db_filename;

    execute_file.write(precompute_remapping_queries)
    execute_file.write(get_unique_peptides + "\n")
    execute_file.write(precompute_match_query + "\n")
    execute_file.write(precompute_first_unique_query + "\n")
    execute_file.write(precompute_second_unique_query + "\n")

    remove_columns_sqlite_file.close()

    execute_file.close()

    #Running the execute scripts
    subprocess.call("sh " + execute_filename, shell=True)

    return


def create_comparison_sqlite_protein_level(left_comparison, right_comparison, task_id, user, path_to_scratch_space, output_folder, path_to_front_end_tasks, path_to_sqlite3):
    left_compare_obj = json.loads(left_comparison)
    right_compare_obj = json.loads(right_comparison)

    #Ultimate output files
    protein_level_merged_db_filename = os.path.join(output_folder, "protein_level.db")


    execute_filename = os.path.join(path_to_scratch_space, "execute.sh")
    execute_file = open(execute_filename, "w")

    left_tasks = json.loads(left_comparison)
    right_tasks = json.loads(right_comparison)

    for task in left_tasks:
        process_tab_fileinto_unified_sqlite("group_by_protein", task, "A", protein_level_merged_db_filename, execute_file, None, path_to_sqlite3, path_to_front_end_tasks, path_to_scratch_space)


    for task in right_tasks:
        process_tab_fileinto_unified_sqlite("group_by_protein", task, "B", protein_level_merged_db_filename, execute_file, None, path_to_sqlite3, path_to_front_end_tasks, path_to_scratch_space)

    precompute_match_query = "echo \"CREATE TABLE matched_rows AS SELECT * FROM Result as first JOIN Result as second ON (first.ComparePartition <> second.ComparePartition) AND (first.accession = second.accession) WHERE first.ComparePartition='A';\" | " + path_to_sqlite3 + " " + protein_level_merged_db_filename;
    precompute_first_unique_query = "echo \"CREATE TABLE unique_first AS SELECT * FROM Result as first LEFT JOIN Result as second ON (first.accession = second.accession) AND (first.ComparePartition <> second.ComparePartition) WHERE first.ComparePartition='A' AND second.accession is NULL;\" | " + path_to_sqlite3 + " " + protein_level_merged_db_filename;
    precompute_second_unique_query = "echo \"CREATE TABLE unique_second AS SELECT * FROM Result as first LEFT JOIN Result as second ON (first.accession = second.accession) AND (first.ComparePartition <> second.ComparePartition) WHERE first.ComparePartition='B' AND second.accession is NULL;\" | " + path_to_sqlite3 + " " + protein_level_merged_db_filename;

    execute_file.write(precompute_match_query + "\n")
    execute_file.write(precompute_first_unique_query + "\n")
    execute_file.write(precompute_second_unique_query + "\n")


    execute_file.close()

    #Running the execute scripts
    subprocess.call("sh " + execute_filename, shell=True)

    return


def create_comparison_sqlite_remapped_protein_level(input_peptide_mapping, path_to_protein_remapped_scratch_space, output_result_folder, path_to_sqlite3):
    protein_level_merged_db_filename = os.path.join(output_result_folder, "protein_remapped_level.db")
    execute_filename = os.path.join(path_to_protein_remapped_scratch_space, "execute.sh")
    execute_file = open(execute_filename, "w")

    #Write some stuff
    A_proteins = defaultdict(int)
    B_proteins = defaultdict(int)

    for peptide in input_peptide_mapping:
        peptide_mapping_object = input_peptide_mapping[peptide]
        if peptide_mapping_object["A_peptide"] != 0:
            for protein in peptide_mapping_object["proteins"]:
                A_proteins[protein] += peptide_mapping_object["A_peptide"]

        if peptide_mapping_object["B_peptide"] != 0:
            for protein in peptide_mapping_object["proteins"]:
                B_proteins[protein] += peptide_mapping_object["B_peptide"]

    #Adding tables to database
    dump_filename = os.path.join(path_to_protein_remapped_scratch_space, "remapped_protein_dump.dump")
    dump_file = open(dump_filename, "w")

    create_lines = []
    create_lines.append("PRAGMA foreign_keys=OFF;")
    create_lines.append("BEGIN TRANSACTION;")
    create_lines.append("CREATE TABLE Result")
    create_lines.append("('accession' TEXT,'Hits' INTEGER,ComparePartition VARCHAR(16));")
    dump_file.write("\n".join(create_lines))
    dump_file.write("\n")

    for protein in A_proteins:
        insert_statement = "INSERT INTO \"Result\" VALUES('%s', %d, '%s');" % (protein, A_proteins[protein], "A")
        dump_file.write(insert_statement + "\n")

    for protein in B_proteins:
        insert_statement = "INSERT INTO \"Result\" VALUES('%s', %d, '%s');" % (protein, B_proteins[protein], "B")
        dump_file.write(insert_statement + "\n")

    dump_file.write("COMMIT;")

    dump_file.close()

    execute_file.write("cat " + dump_filename + " | " + path_to_sqlite3 + " " + protein_level_merged_db_filename + "\n");

    precompute_match_query = "echo \"CREATE TABLE matched_rows AS SELECT * FROM Result as first JOIN Result as second ON (first.ComparePartition <> second.ComparePartition) AND (first.accession = second.accession) WHERE first.ComparePartition='A';\" | " + path_to_sqlite3 + " " + protein_level_merged_db_filename;
    precompute_first_unique_query = "echo \"CREATE TABLE unique_first AS SELECT * FROM Result as first LEFT JOIN Result as second ON (first.accession = second.accession) AND (first.ComparePartition <> second.ComparePartition) WHERE first.ComparePartition='A' AND second.accession is NULL;\" | " + path_to_sqlite3 + " " + protein_level_merged_db_filename;
    precompute_second_unique_query = "echo \"CREATE TABLE unique_second AS SELECT * FROM Result as first LEFT JOIN Result as second ON (first.accession = second.accession) AND (first.ComparePartition <> second.ComparePartition) WHERE first.ComparePartition='B' AND second.accession is NULL;\" | " + path_to_sqlite3 + " " + protein_level_merged_db_filename;

    execute_file.write(precompute_match_query + "\n")
    execute_file.write(precompute_first_unique_query + "\n")
    execute_file.write(precompute_second_unique_query + "\n")

    execute_file.close()

    #Running the execute scripts
    subprocess.call("sh " + execute_filename, shell=True)

def main():
    usage()

    input_parameters = sys.argv[1]
    output_log_filename = sys.argv[2]
    output_result_folder = sys.argv[3]
    path_to_front_end_tasks = sys.argv[4]
    path_to_scratch_space = "./scratch"
    path_to_peptide_scratch_space = "./peptide_scratch"
    path_to_protein_scratch_space = "./protein_scratch"
    path_to_protein_remapped_scratch_space = "./protein_remapped_scratch"
    path_to_sqlite3 = sys.argv[5]
    input_peptide_mapping = sys.argv[6]

    ming_fileio_library.make_sure_path_exists(path_to_protein_remapped_scratch_space)
    ming_fileio_library.make_sure_path_exists(path_to_scratch_space)
    ming_fileio_library.make_sure_path_exists(path_to_peptide_scratch_space)
    ming_fileio_library.make_sure_path_exists(path_to_protein_scratch_space)

    #Parsing the parameters
    parsed_parameters = ming_proteosafe_library.parse_xml_file(open(input_parameters))

    #Parsing the peptide mappings
    peptide_mappings = json.loads(open(input_peptide_mapping).read())

    #Execute Whatever you need here
    left_parameters = parsed_parameters["left_objects"][0]
    right_parameters = parsed_parameters["right_objects"][0]
    task_id = parsed_parameters["task"][0]
    user = parsed_parameters["user"][0]

    ignore_fileprefix_extension = False
    try:
        if parsed_parameters["ignore_path_prefix_extension"][0] == "Yes":
            ignore_fileprefix_extension = True
    except:
        ignore_fileprefix_extension = False


    create_comparison_sqlite_spectrum_level(left_parameters, right_parameters, task_id, user, path_to_scratch_space, output_result_folder, path_to_front_end_tasks, path_to_sqlite3, peptide_mappings, ignore_fileprefix_extension)
    create_comparison_sqlite_peptide_level(left_parameters, right_parameters, task_id, user, path_to_peptide_scratch_space, output_result_folder, path_to_front_end_tasks, path_to_sqlite3, peptide_mappings)
    create_comparison_sqlite_protein_level(left_parameters, right_parameters, task_id, user, path_to_protein_scratch_space, output_result_folder, path_to_front_end_tasks, path_to_sqlite3)
    create_comparison_sqlite_remapped_protein_level(peptide_mappings, path_to_protein_remapped_scratch_space, output_result_folder, path_to_sqlite3)

    #Output Log File
    output_log_file = open(output_log_filename, "w")
    spectrum_level_refer_url = "/ProteoSAFe/result_comparison/results_comparison_server_view.jsp?task=" + task_id + "&compare_type=spectrum_mismatch"
    output_log_file.write('<div style="text-align: center">')
    output_log_file.write('<a href="' + spectrum_level_refer_url + '">View Comparison</a>')
    output_log_file.write('</div>')
    output_log_file.write('<meta http-equiv="refresh" content="0; url=' + spectrum_level_refer_url + '" />');
    output_log_file.close()


if __name__ == "__main__":
    main()
