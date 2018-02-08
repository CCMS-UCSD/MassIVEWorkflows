#!/usr/bin/python


import sys
import getopt
import os
import json
import ming_fileio_library

def usage():
    print "<input folder> <output merged file> <format: json/tsv>"


def main():
    input_intermediate_folder = sys.argv[1]
    output_filename = sys.argv[2]
    output_format = sys.argv[3]

    merged_data = {}

    #Creating a command line for each partition
    all_intermediate_files = ming_fileio_library.list_files_in_dir(input_intermediate_folder)
    for parallel_output_filename in all_intermediate_files:
        #parallel_output_filename = os.path.join(scratch_folder, str(i) + ".json")
        json_parition = json.loads(open(parallel_output_filename, "r").read())

        for key in json_parition:
            merged_data[key] = json_parition[key]

    if output_format == "json":
        open(output_filename, "w").write(json.dumps(merged_data))
    if output_format == "tsv":
        ming_fileio_library.write_dictionary_table_data(merged_data, output_filename)

    exit(0)


if __name__ == "__main__":
    main()
