import subprocess
from Query import Query
import os
import sys
import json
import glob

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
    '''
    This function goes to the specified directory, locates relevant (.fasta) 
    files
    '''
    for subdir, dirs, files in os.walk(dirname):
        for one_file in files:
            filename = os.path.join(subdir, one_file)
            if filename.find(".fasta") != -1:
                pass
                # process fasta files

    return


def get_data_from_json(filename):
    '''
    Parse a JSON file and return a JSON object (key:value format)
    '''
    json_file = open(filename)
    data_from_json = json.load(json_file)
    json_file.close()
    return data_from_json


def main():
    all_queries = get_data_from_json("AllBlastData03.json")
    find_filterfasta_output_files("filterfasta_files")


if __name__ == "__main__":
    main()
