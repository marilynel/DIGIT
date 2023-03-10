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
                process_blast_output(filename)

    return


def process_blast_output(filename):
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


def set_best_for_genome():
    '''
    Create dictionary identifying the best BLAST hit for each query, within
    each genome.
    '''
    list_of_queries = []

    for genome in all_queries:
        for allele in all_queries[genome]:
            if allele not in list_of_queries:
                list_of_queries.append(allele)
            # if allele not in best_query_per_genome[genome]:
            # if all_queries[genome][allele]
            get_best_query(genome, allele)
    write_to_best_queries_file(list_of_queries)

    return list_of_queries


def allele_in_genome_write(genome, allele, best_bit_score):
    '''
    Find the specific hit that was specified as the best for that query ID and
    that genome (Query.best_for_genome == True) and return a string reporting
    data for the best_queries_by_genome.csv outfile.
    '''
    if allele in all_queries[genome]:
        for i in range(0, len(all_queries[genome][allele])):
            if all_queries[genome][allele][i].best_for_genome == True:
                best_bit_score[genome] = all_queries[genome][allele][i].bit_score
                return all_queries[genome][allele][i].__make_best_genome_string__()
    return "none,none,none,"


def best_genomes_write(best_bit_score):
    '''
    Find the best of the best hits, return strinigied list for writing to file.
    '''
    genome_list = []
    best_score = max(best_bit_score.values())
    for genome in best_bit_score:
        if best_bit_score[genome] == best_score:
            genome_list.append(genome)
        else:
            best_bit_score[genome] = 0
    return str(genome_list)


def working_query_selection(genome, allele):
    '''
    set ID for best_query --> indicates that this will belong to working set
    '''
    for i in range(0, len(all_queries[genome][allele])):
        if all_queries[genome][allele][i].best_for_genome == True:
            all_queries[genome][allele][i].best_query = True
            return


def pick_genome(allele, best_bit_score):
    if best_bit_score["B73v5"] != 0:
        working_query_selection("B73v5", allele)
        return
    elif best_bit_score["W22v2"] != 0:
        working_query_selection("W22v2", allele)
        return
    elif best_bit_score["A188v1"] != 0:
        working_query_selection("A188v1", allele)
    else:
        print("you really messed something up")
        exit()


def write_to_best_queries_file(list_of_queries):
    with open("best_queries_by_genome.csv", "w+") as newfile:
        newfile.write(
            "query,A188_bit_score,A188_num_hits,A188_qstart_status,B73_bit_sc" +
            "ore,B73_num_hits,B73_qstart_status,W22_bit_score,W22_num_hits,W2" +
            "2_qstart_status,best_genomes\n"
        )
        for allele in list_of_queries:
            best_bit_score = {}
            newfile.write(
                allele + "," +
                allele_in_genome_write("A188v1", allele, best_bit_score) +
                allele_in_genome_write("B73v5", allele, best_bit_score) +
                allele_in_genome_write("W22v2", allele, best_bit_score) +
                best_genomes_write(best_bit_score) +
                "\n"
            )
            pick_genome(allele, best_bit_score)
    return


def get_best_query(genome, query):
    '''
    Sort hits to identify best hit for a query. Find the percent difference
    between the best and second best hits.
    Returns:
    best_query              Query object
    '''
    best_query = all_queries[genome][query][0]
    second_best_query = None
    for i in range(1, len(all_queries[genome][query])):
        if all_queries[genome][query][i].bit_score > best_query.bit_score:
            second_best_query = best_query
            best_query = all_queries[genome][query][i]
        elif all_queries[genome][query][i].bit_score < best_query.bit_score:
            if not second_best_query or all_queries[genome][query][i].bit_score > second_best_query.bit_score:
                second_best_query = all_queries[genome][query][i]
    # TODO: come back and make this work
    # if second_best_query:
    #    best_query.diff = abs(
    #        ((best_query.bit_score - second_best_query.bit_score) / (best_query.bit_score + second_best_query.bit_score)) / 2) * 100

    for i in range(0, len(all_queries[genome][query])):
        if all_queries[genome][query][i] == best_query:
            all_queries[genome][query][i].best_for_genome = True

    return


def allele_in_genome_comp(genome, allele):
    if allele in all_queries[genome]:
        for i in range(0, len(all_queries[genome][allele])):
            if all_queries[genome][allele][i].best_for_genome == True:
                return all_queries[genome][allele][i].bit_score
    return -1


def set_best_query(list_of_queries):
    for query in list_of_queries:
        best_bit_scores = [
            allele_in_genome_comp("A188v1", query),
            allele_in_genome_comp("B73v5", query),
            allele_in_genome_comp("W22v2", query),
        ]


def queries_to_json(filename):
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
    find_blast_output_files(sys.argv[1])
    list_of_queries = set_best_for_genome()
    set_best_query(list_of_queries)
    queries_to_json("AllBlastData.json")

    # PERL STUFF

    # TODO: this doesn't work yet, I'm still figuring out why
    
    for gen in all_queries:
        for q in all_queries[gen]:
            for i in range(0, len(all_queries[gen][q])):
                if all_queries[gen][q][i].best_query == True:
                    # $1 = query
                    # $2 = chromosome
                    # $3 = wildtype start
                    # $4 = wildtype end
                    # $5 = upper start
                    # $6 = upper end
                    # $7 = lower start
                    # $8 = lower end
                    # $9 = genome
                    subprocess.run(
                        [
                            "SGE_Batch",
                            "-c",
                            f"./filterfasta.sh {all_queries[gen][q][i].query} {all_queries[gen][q][i].chromosome} {all_queries[gen][q][i].wildtype_coordinates[0]} {all_queries[gen][q][i].wildtype_coordinates[1]} {all_queries[gen][q][i].upper_coordinates[0]} {all_queries[gen][q][i].upper_coordinates[1]} {all_queries[gen][q][i].lower_coordinates[0]} {all_queries[gen][q][i].lower_coordinates[1]} {gen[:-2]}",
                            # "filterfasta.sh",
                            # all_queries[gen][q][i].query,
                            # all_queries[gen][q][i].chromosome,
                            # str(all_queries[gen][q][i].wildtype_coordinates[1]),
                            # str(all_queries[gen][q][i].upper_coordinates[0]),
                            # str(all_queries[gen][q][i].upper_coordinates[1]),
                            # str(all_queries[gen][q][i].lower_coordinates[0]),
                            # str(all_queries[gen][q][i].lower_coordinates[1]),
                            # gen[:-2],
                            "-q",
                            "bpp",
                            "-P",
                            "8",
                            "-r",
                            f"sge.{all_queries[gen][q][i].query}"
                        ]
                    )
                    exit()

    # for genome in all_queries:
    #    for query_list in all_queries[genome]:
    #        for i in range(0,len(all_queries[genome][query_list])):
    #            if all_queries[genome][query_list][i].best_query == True:
    #                try:
    #                    subprocess.run(["practice.sh", all_queries[genome][query_list][i].query])
    #                except:
    #                    print("nope")
    #                    exit(-1)

    # filterfasta Usage:
    # filterfasta.pl --match '$chromosome' --start '$start' --end '$end' '$db'/Zm*
'''
Input:
go through queries in all_queries[genome] (for each genome)
    check if Query.best_query == True
        call shell one or three times? which would be better?
        filterfasta.pl --match '$chromosome' --start '$start' --end '$end' '$db'/Zm*


Upper
chromosome          Query.chromosome            string
start               Query.upper_coordinates[0]  int? may need to convert to string
end                 Query.upper_coordinates[1]  int? may need to convert to string
database            Query.genome

Lower
start               Query.lower_coordinates[0]
end                 Query.lower_coordinates[1]

WT
start               Query.wildtype_coordinates[0]
end                 Query.wildtype_coordinates[1]

can I go straight to SGE_Batch from the subprocess?

could I do one fasta file? just to make life a little easier
    -> that might be weird with multiple sge jobs
    -> maybe instead a new fasta file for each query?

    wt_name=${query}_in_${db}_wildtype_${strand}_${chromosome}_${wt_start}_${wt_end}
    echo ">$query,$chromosome,$strand,$db,wildtype,s_start placeholder" >> data/filterfasta/filterfasta_output/$temp/$query/$wt_name.fasta
    SGE_Batch -c 'part1/filterfasta.pl --match '$chromosome' --start '$wt_start' --end '$wt_end' resources/part0/reference_genomes/'$db'/Zm* >> data/filterfasta/filterfasta_output/'$temp'/'$query'/'$wt_name'.fasta' -q bpp -P 8 -r sge.${temp}_wild_${query}_${now}



'''