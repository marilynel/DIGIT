import subprocess
from Query import Query
import os
import sys
import json
import glob

queries_working_set = {}


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
    for allele in data_from_json:
        if allele not in queries_working_set:
            queries_working_set[allele] = Query(None, None, None)
            queries_working_set[allele].__Query_from_JSON__(data_from_json[allele])


def read_primer3_output(filename):
    with open(filename, "r") as p3file:
        output_dict = {}  # create a dictionary for each output item
        for line in p3file:
            if line[0] != "=":  # go through the file until reaching a break mark ("=")
                lineitems = line.split("=")  # key:value pairs for OUTPUTLABEL:data
                output_dict[lineitems[0].strip()] = lineitems[1].strip()
            elif line[0] == "=":  # end of a section,time to process data
                allele = output_dict["SEQUENCE_ID"][:-1]

                # print(line)
                # print(allele)
                # print(output_dict.keys())

                queries_working_set[allele].__update_Query_primer3__(output_dict)
                output_dict.clear()
                # print(output_dict)
            else:
                print("ERROR")


# TODO: is this generalizable? similar function appears in another part of program
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


def main():
    # get all the primer data
    create_query_struct("WorkingQueriesWithPrimers03.json")
    # read_primer3_output("Primer3Output08.txt")

    # make an input filefor primer3
    with open("VerifyPrimer3Input.txt", "w") as p3file:
        for allele in queries_working_set:
            p3file.write(queries_working_set[allele].__create_primer3_input__())

    queries_to_json("")


if __name__ == "__main__":
    main()