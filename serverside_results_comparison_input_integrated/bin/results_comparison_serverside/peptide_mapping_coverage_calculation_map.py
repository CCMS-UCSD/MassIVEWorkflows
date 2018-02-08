#!/usr/bin/python


import sys
import getopt
import os
import json
import ming_protein_library
import ming_fileio_library
import cPickle as pickle

def usage():
    print "<input fasta> <input peptides> <input parameters file> <output results folder>"

def main():
    print sys.argv

    input_fasta_filename = sys.argv[1]
    input_peptide_list_filename = sys.argv[2]
    input_parameters_filename = sys.argv[3]
    output_results_folder = sys.argv[4]

    partition_total = 1
    partition_of_node = 0

    params_map = json.loads(open(input_parameters_filename).read())
    partition_total = params_map["total_paritions"]
    partition_of_node = params_map["node_partition"]

    output_coverage_filename = os.path.join(output_results_folder, str(partition_of_node) + ".json")

    proteome = ming_protein_library.parse_fasta_proteome_file(input_fasta_filename)

    peptide_data = json.loads(open(input_peptide_list_filename).read())

    left_peptides = peptide_data["left_peptides"]
    right_peptides = peptide_data["right_peptides"]
    peptide_list = left_peptides.keys() + right_peptides.keys()
    peptide_list = list(set(peptide_list))
    peptide_list = peptide_list[partition_of_node::partition_total]

    peptide_to_protein_map = {}

    peptide_mapping = proteome.get_peptides_mapped_to_proteins_efficient(peptide_list)

    for peptide in peptide_list:
        proteins_list =  peptide_mapping[peptide]
        peptide_information_map = {}
        peptide_information_map["proteins"] = proteins_list
        if peptide in left_peptides:
            peptide_information_map["A_peptide"] = left_peptides[peptide]
        else:
            peptide_information_map["A_peptide"] = 0

        if peptide in right_peptides:
            peptide_information_map["B_peptide"] = right_peptides[peptide]
        else:
            peptide_information_map["B_peptide"] = 0

        peptide_to_protein_map[peptide] = peptide_information_map

    open(output_coverage_filename, "w").write(json.dumps(peptide_to_protein_map))

    exit(0)


if __name__ == "__main__":
    main()
