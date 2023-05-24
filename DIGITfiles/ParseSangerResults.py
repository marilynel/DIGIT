import sys
import os

from Query import Query
from Utils import *

sangerQueriesAll = {
    "A188v1": {},
    "B73v5": {},
    "W22v2": {}
}


####################################################################


# sanger_blast_output/good_matches/
def lookup_in_file(db, all_query_objs, flanking_file_name, list_name):
    '''
    Lookup query data in file corresponding to best genome for that query
    '''
    # print(flanking_file_name)
    with open(
            "order_" + flanking_file_name + "/blastResults/good_matches/" + db + "_vs_sanger_seq_" + flanking_file_name +
            "_good_matches.txt", "r+") as parsed_file:
        first_line = parsed_file.readline()
        for line in parsed_file:
            line_items = line.split(',')
            line_items[-1].strip()
            # print(list_name)
            if line_items[0] in list_name:
                query_data = line_items

                col_headings = ("query acc.ver,subject acc.ver,percent identity"
                                ",alignment length,mismatches,gap opens,q. start,q. end,s. sta"
                                "rt,s. end,evalue,bit score,database,number of hits,percent di"
                                "fference to next value,strand,q. start status,bit score statu"
                                "s\n")
                col_headings_list = col_headings.split(',')

                query_obj = {
                    col_headings_list[i]: query_data[i] for i in range(len(col_headings_list))}
                # print(query_obj)
                all_query_objs.append(query_obj)
    return


def get_best_genome(best_genome_dict, flanking_file_name):
    '''
    Find the best fit genome for a query based on Step 2 results, and direct to
    look up relevant data on the query
    '''
    list_a, list_b, list_w = [], [], []

    for query in best_genome_dict:
        # B73 is preferred if B, BW, AB, or ABW
        if "B73v5" in best_genome_dict[query]["best_genomes"]:
            list_b.append(query)
        # W22 is preferred if W, AW
        elif "W22v2" in best_genome_dict[query]["best_genomes"]:
            list_w.append(query)
            # A188 only for A
        elif "A188v1" in best_genome_dict[query]["best_genomes"]:
            list_a.append(query)

    all_query_objs = []
    lookup_in_file("A188v1", all_query_objs, flanking_file_name, list_a)
    lookup_in_file("B73v5", all_query_objs, flanking_file_name, list_b)
    lookup_in_file("W22v2", all_query_objs, flanking_file_name, list_w)
    return all_query_objs


def main():
    flankseq = "SangerSeq_Order24639"

    # 1. parse blast output, save as queries in struct
    findBlastOutputFiles(f"../DIGITfiles/BlastOutputSangerSeqs/{flankseq}", sangerQueriesAll)

    # 2. get best queries to compare

    print(sangerQueriesAll)


if __name__ == '__main__':
    main()