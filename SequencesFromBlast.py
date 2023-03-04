import subprocess
from Query import Query
import os
import sys
import json

all_queries = {
    "A188v1": {},
    "B73v5": {},
    "W22v2": {}
}


def find_blast_output_files(dirname):
    '''
    This function goes to the specified directory, locates relevant (.tab)
    files, and sends the file names to function process_blast_output() for
    parsing.
    '''
    for subdir, dirs, files in os.walk(dirname):
        for one_file in files:
            filename = os.path.join(subdir, one_file)
            if filename.find(".tab") != -1:
                process_blast_output(filename)  # , all_queries)

    # return all_queries
    return


def process_blast_output(filename):  # , all_queries):
    '''
    Parse blast output files line by line, gathering data that will be used to
    create Query objects, which will be added to the all_queries dictionary.
    Functions get_num_hits() and get_genome() are called to help parse hash (#)
    lines.
    '''
    genome = ""
    num_hits = -1
    with open(filename, "r+") as blastfile:
        for line in blastfile:
            if line[0] == '#':
                # hash lines may hold data that will be necessary for Query
                # object initialization. genome and num_hits will remain
                # unchanged until after all of the hits for that query are read
                genome = get_genome(line, genome)
                num_hits = get_num_hits(line, num_hits)
            else:
                new_query = Query(line, genome, num_hits)
                new_query.__set_values__()
                if new_query.query not in all_queries[genome]:
                    # Add that query to all_queries[genome] as a key, value
                    # being an empty list
                    all_queries[genome][new_query.query] = []
                # Place Query object in all_queries object
                all_queries[genome][new_query.query].append(new_query)
    return


def get_genome(line, genome):
    '''
    Parse a line for string indicating the genome name. If line does not have
    that string, return genome as original value.
    '''
    blast_data = line.split(' ')
    if blast_data[1] == "Database:":
        return blast_data[2].split('/')[-1].strip()
    return genome


def get_num_hits(line, num_hits):
    '''
    Parse a line for string indicating the number of hits for that query in a
    genome. If line does not have that string, return num_hits as original
    value.
    '''
    blast_data = line.split(' ')
    if blast_data[2] == "hits":
        return int(blast_data[1])
    return num_hits


def best_queries_per_genome():  # all_queries):
    '''
    Create dictionary identifying the best BLAST hit for each query, within
    each genome.
    '''
    best_query_per_genome = {
        "query_ids": [],
        "A188v1": {},
        "B73v5": {},
        "W22v2": {}
    }

    for genome in all_queries:
        for allele in all_queries[genome]:
            if allele not in best_query_per_genome["query_ids"]:
                best_query_per_genome["query_ids"].append(allele)
            if allele not in best_query_per_genome[genome]:
                # best_query_per_genome[genome][allele] = get_best_query(genome, allele)#all_queries[genome][allele])
                get_best_query(genome, allele)
    # write_to_best_queries_file(best_query_per_genome)

    return best_query_per_genome


def allele_in_genome_write(genome_dict, allele):
    if allele in genome_dict:
        return genome_dict[allele].__make_best_genome_string__()
    return "none,none,none,"


def write_to_best_queries_file(best_query_per_genome):
    with open("best_queries_by_genome.csv", "w+") as newfile:
        newfile.write(
            "query,A188_bit_score,A188_num_hits,A188_qstart_status,B73_bit_sc" +
            "ore,B73_num_hits,B73_qstart_status,W22_bit_score,W22_num_hits,W2" +
            "2_qstart_status,best_genomes\n"
        )
        for allele in best_query_per_genome["query_ids"]:
            newfile.write(
                allele + "," +
                allele_in_genome_write(best_query_per_genome["A188v1"], allele) +
                allele_in_genome_write(best_query_per_genome["B73v5"], allele) +
                allele_in_genome_write(best_query_per_genome["W22v2"], allele) +
                "\n"
            )

            # TODO:
            # compare accross genomes
            # does this happen elsewhere???? --> yes
            # call allele in gehome comp before hand??
            # change these damn function names!
    return


def get_best_query(genome, query):  # list_of_hits):
    '''
    Sort hits to identify best hit for a query. Find the percent difference
    between the best and second best hits.
    Returns:
    best_query              Query object
    '''
    # best_query = list_of_hits[0]
    best_query = all_queries[genome][query][0]
    second_best_query = None
    # for i in range(1, len(list_of_hits)):
    for i in range(1, len(all_queries[genome][query])):
        # if list_of_hits[i].bit_score > best_query.bit_score:
        if all_queries[genome][query][i].bit_score > best_query.bit_score:
            second_best_query = best_query
            # best_query = list_of_hits[i]
            best_query = all_queries[genome][query][i]
        # elif list_of_hits[i].bit_score < best_query.bit_score:
        elif all_queries[genome][query][i].bit_score < best_query.bit_score:
            # if not second_best_query or list_of_hits[i].bit_score > second_best_query.bit_score:
            if not second_best_query or all_queries[genome][query][i].bit_score > second_best_query.bit_score:
                # second_best_query = list_of_hits[i]
                second_best_query = all_queries[genome][query][i]
    # TODO: come back and make this work
    # if second_best_query:
    #    best_query.diff = abs(
    #        ((best_query.bit_score - second_best_query.bit_score) / (best_query.bit_score + second_best_query.bit_score)) / 2) * 100

    # best_query.best_for_genome = True
    for i in range(0, len(all_queries[genome][query])):
        if all_queries[genome][query][i] == best_query:
            all_queries[genome][query][i].best_for_genome = True

    return  # best_query


def allele_in_genome_comp(genome_dict, allele):
    if allele in genome_dict:
        return genome_dict[allele].bit_score
    return -1


def best_queries(best_query_per_genome):
    '''
    working_query_set = []
    for query in best_query_per_genome["query_ids"]:
        best_bit_scores = [
            allele_in_genome_comp(best_query_per_genome["A188v1"],query),
            allele_in_genome_comp(best_query_per_genome["B73v5"],query),
            allele_in_genome_comp(best_query_per_genome["W22v2"],query),
        ]

        if best_query_per_genome["B73v5"][query].bit_score == max(best_bit_scores):
            working_query_set.append(best_query_per_genome["B73v5"][query])
        elif best_query_per_genome["W22v2"][query].bit_score == max(best_bit_scores):
            working_query_set.append(best_query_per_genome["W22v2"][query])
        elif best_query_per_genome["A188v1"][query].bit_score == max(best_bit_scores):
            working_query_set.append(best_query_per_genome["A188v1"][query])
    return working_query_set
    '''
    pass


def queries_to_json(filename):  # data, filename):
    '''
    Converts dictionary to JSON object and writes it to a file.
    '''
    json_object = {
        "A188v1": {},
        "B73v5": {},
        "W22v2": {}
    }

    for database in all_queries:
        for query in all_queries[database]:
            if query not in json_object[database]:
                json_object[database][query] = []
            for hit in all_queries[database][query]:
                json_object[database][query] = dict(hit)

    new_json_object = json.dumps(json_object, indent=4)
    with open(filename, "w") as outfile:
        outfile.write(new_json_object)

    return


if __name__ == '__main__':
    # all_queries =
    find_blast_output_files(sys.argv[1])
    # queries_to_json("AllBlastDataBefore.json")
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

    best_query_per_genome = best_queries_per_genome()  # all_queries)
    # print("after best_query_per_genome")
    # print("before working query set")
    # working_query_set = best_queries(best_query_per_genome)
    # print("after working query set")
    # print(working_query_set)
    # queries_to_json(working_query_set, "outfile.json")
    # for query in working_query_set:
    #    print(dict(query))
    queries_to_json("AllBlastData.json")

    ''' PERL STUFF
    for q in working_query_set:
        my_cmd = "perl practice.pl " + q.query
        pipe = os.popen(my_cmd)
        for line in pipe:
            print(line.rstrip())
        print("ok")
    '''

    # filterfasta Usage:
    # filterfasta.pl --match '$chromosome' --start '$start' --end '$end' '$db'/Zm*
