import subprocess
from Query import Query
from Sequences import Sequences
from Primer3Object import Primer3Object
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
                position = line.split("_")[-1]  # upper, lower, or wildtype
            else:
                fasta_dict[position.strip()] = line.strip()
        # TODO: may change the below to be handled in Query method
        queries_working_set[allele].upper_sequence = fasta_dict["upper"]
        queries_working_set[allele].lower_sequence = fasta_dict["lower"]
        queries_working_set[allele].wildtype_sequence = fasta_dict["wildtype"]
        queries_working_set[allele].__build_insertion_sequence__()


def queries_to_json(filename):
    '''
    Converts dictionary to JSON object and writes it to a file.
    '''
    json_object = {}
    for q in queries_working_set:
        if queries_working_set[q].query not in json_object:
            json_object[queries_working_set[q].query] = []
        json_object[queries_working_set[q].query] = dict(queries_working_set[q])
    new_json_object = json.dumps(json_object, indent=4)

    with open(filename, "w") as outfile:
        outfile.write(new_json_object)

    return


def run_primer3():
    subprocess.run(
        [
            "SGE_Batch",
            "-c",
            "primer3_core Primer3Input08.txt > Primer3Output.txt",
            "-q",
            "bpp",
            "-P",
            "8",
            "-r",
            "sge.primer3"
        ]
    )


def main():
    create_query_struct("AllBlastData08.json")
    find_filterfasta_output_files("filterfasta_files")
    with open("Primer3Input.txt", "w") as p3file:
        for allele in queries_working_set:
            leftP3 = Primer3Object("generic", "left", queries_working_set[allele])
            # p3file.write(queries_working_set[allele].__create_primer3_input__())
            rightP3 = Primer3Object("generic", "right", queries_working_set[allele])
            inputStr = leftP3.input_str + rightP3.input_str
            p3file.write(inputStr)

    queries_to_json("WorkingQueries.json")

    if sys.argv[1] == "true":
        run_primer3()


if __name__ == "__main__":
    main()