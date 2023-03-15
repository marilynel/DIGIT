import subprocess
from Query import Query
import os
import sys
import json
import glob

queries_working_set = {}

def find_filterfasta_output_files(dirname):
    fasta_dict = {}
    for subdir, dirs, files in os.walk(dirname):
        for one_file in files:
            filename = os.path.join(subdir, one_file)
            if filename.find(".fasta") != -1:
                rfile = filename.split("/")[1]
                allele = rfile.split(".")[0]
                store_and_construct_sequences(filename, allele)
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
                position = line.split("_")[-1]      # upper, lower, or wildtype
            else:
                fasta_dict[position.strip()] = line.strip()
        queries_working_set[allele].upper_sequence = fasta_dict["upper"]
        queries_working_set[allele].lower_sequence = fasta_dict["lower"]
        queries_working_set[allele].wildtype_sequence = fasta_dict["wildtype"]
        queries_working_set[allele].__build_insertion_sequence__()




def main():
    create_query_struct("AllBlastData05.json")
    find_filterfasta_output_files("filterfasta_files")
    with open("primer3_input.txt", "w") as p3file:
        for allele in queries_working_set:
            p3file.write(queries_working_set[allele].__create_primer3_input__())



    # TODO:
    # pass this file to primer3 and direct results somehwere



if __name__ == "__main__":
    main()