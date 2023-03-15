import subprocess
from Query import Query
import os
import sys
import json
import glob

queries_working_set = {}


# all queries is a dict with structure:
#   {
#       databaseA:
#           {
#               query1:[QueryObj, QueryObj, QueryObj],
#               query2:[QueryObj]
#           }
#       databaseB:
#           {
#               query1:[QueryObj, QueryObj, QueryObj],
#               query2:[QueryObj]
#           }
#   }   etc.


def find_filterfasta_output_files(dirname):
    fasta_dict = {}
    for subdir, dirs, files in os.walk(dirname):
        for one_file in files:
            filename = os.path.join(subdir, one_file)
            if filename.find(".fasta") != -1:
                # fasta_dict = {
                #    "upper": "",
                #    "lower": "",
                #    "wildtype": ""
                # }
                # save the first part of the filename as allele
                rfile = filename.split("/")[1]
                allele = rfile.split(".")[0]
                store_and_construct_sequences(filename, allele)
                ##with open(filename, "r") as fastafile:
                #   for line in fastafile:
                #       if line[0] == ">":
                #           position = line.split("_")[-1]      # upper, lower, or wildtype
                #       else:
                #           fasta_dict[position.strip()] = line.strip()
                #   queries_working_set[allele].upper_sequence = fasta_dict["upper"]
                #   queries_working_set[allele].lower_sequence = fasta_dict["lower"]
                #   queries_working_set[allele].wildtype_sequence = fasta_dict["wildtype"]
                #   queries_working_set[allele].__build_insertion_sequence__()

    return


def get_data_from_json(filename):
    '''
    Parse a JSON file and return a JSON object (key:value format)
    '''
    json_file = open(filename)
    data_from_json = json.load(json_file)
    json_file.close()

    return data_from_json


def create_query_struct(filename):
    data_from_json = get_data_from_json(filename)
    for genome in data_from_json:
        for allele in data_from_json[genome]:
            for i in range(0, len(data_from_json[genome][allele])):
                if data_from_json[genome][allele][i]["best_query"]:
                    queries_working_set[allele] = Query(None, None, None)
                    queries_working_set[allele].__Query_from_JSON__(data_from_json[genome][allele][i])


def store_and_construct_sequences(filename, allele):
    fasta_dict = {
        "upper": "",
        "lower": "",
        "wildtype": ""
    }
    with open(filename, "r") as fastafile:
        for line in fastafile:
            if line[0] == ">":
                position = line.split("_")[-1]  # upper, lower, or wildtype
            else:
                fasta_dict[position.strip()] = line.strip()
        queries_working_set[allele].upper_sequence = fasta_dict["upper"]
        queries_working_set[allele].lower_sequence = fasta_dict["lower"]
        queries_working_set[allele].wildtype_sequence = fasta_dict["wildtype"]
        queries_working_set[allele].__build_insertion_sequence__()

    # TODO: figire out how to make this happen


def write_to_primer3_input_files(sequence, query, strand, somestring):
    '''
    Call f(x) to write Primer3 input to a file with relevant data and
    parameters.
    '''
    path = "data/primer3/input/" + somestring + "/" + somestring
    if strand == "1":
        # 3dsgg as right primer file
        support.write_primer3_input(path + "_plus_B_left_primers.txt", query, sequence, False, "left", "b")
        # gfp3utr as left primer file
        support.write_primer3_input(path + "_plus_A_right_primers.txt", query, sequence, False, "right", "a")
    elif strand == "-1":
        # 3dsgg as left primer file
        support.write_primer3_input(path + "_minus_B_right_primers.txt", query, sequence, False, "right", "b")
        # gfp3utr as right primer file
        support.write_primer3_input(path + "_minus_A_left_primers.txt", query, sequence, False, "left", "a")


def main():
    create_query_struct("AllBlastData04.json")
    find_filterfasta_output_files("filterfasta_files")
    queries_working_set["R78G02"].__print_query__()


if __name__ == "__main__":
    main()