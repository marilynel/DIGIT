from query import Query
import os
import sys
# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

def read_blast_output(dirname):
    all_queries = {
        "A188v1": {},
        "B73v5": {},
        "W22v2": {}
    }

    for subdir, dirs, files in os.walk(dirname):
        for one_file in files:
            filename = os.path.join(subdir, one_file)
            if filename.find(".tab") != -1:
                process_blast_output_line(filename, all_queries)

    return all_queries


'''def process_data_line(line, db_id, num_hits, all_queries):
    new_query = Query(line, db_id, num_hits)
    print(new_query.print_query())
    if new_query.query not in all_queries[db_id]:
        all_queries[db_id][new_query.query] = []
    all_queries[db_id][new_query.query].append(new_query)'''


'''#def process_hash_line(line, db_id, num_hits):
def process_hash_line(line):
    blast_data = line.split(' ')
    if blast_data[1] == "Database:":
        db_id = blast_data[2].split('/')
        db_id = db_id[-1].strip()
    if blast_data[2] == "hits":
        num_hits = int(blast_data[1])'''


def get_db_id(line, db_id):
    blast_data = line.split(' ')
    if blast_data[1] == "Database:":
        return blast_data[2].split('/')[-1].strip()
    return db_id


def get_num_hits(line, num_hits):
    blast_data = line.split(' ')
    if blast_data[2] == "hits":
        return int(blast_data[1])
    return num_hits


def process_blast_output_line(filename, all_queries):
    db_id = ""
    num_hits = -1
    with open(filename, "r+") as blastfile:
        for line in blastfile:
            if line[0] == '#':
                db_id = get_db_id(line, db_id)
                num_hits = get_num_hits(line, num_hits)
            else:
                new_query = Query(line, db_id, num_hits)
                if new_query.query not in all_queries[db_id]:
                    all_queries[db_id][new_query.query] = []
                all_queries[db_id][new_query.query].append(new_query)




def best_queries_per_genome(all_queries):
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
                best_query_per_genome[genome][allele] = get_best_query(all_queries[genome][allele])
    write_to_best_queries_file(best_query_per_genome)
    return best_query_per_genome


def allele_in_genome_write(genome_dict,allele):
    if allele in genome_dict:
        return genome_dict[allele].make_best_genome_string()
    return "none,none,none,"


def write_to_best_queries_file(best_query_per_genome):
    with open("best_queries_by_genome.csv", "w+") as newfile:
        newfile.write("query,A188_bit_score,A188_num_hits,A188_qstart_status,B73_bit_score,B73_num_hits,B73_qstart_stat"
            "us,W22_bit_score,W22_num_hits,W22_qstart_status,best_genomes\n")
        for allele in best_query_per_genome["query_ids"]:
            newfile.write(allele + "," + allele_in_genome_write(best_query_per_genome["A188v1"],allele) + allele_in_genome_write(best_query_per_genome["B73v5"],allele) + allele_in_genome_write(best_query_per_genome["W22v2"],allele) + "\n")
    return





def get_best_query(queries_by_genome):
    best_query = queries_by_genome[0]
    second_best_query = None
    for i in range(1, len(queries_by_genome)):
        if queries_by_genome[i].bit_score > best_query.bit_score:
            second_best_query = best_query
            best_query = queries_by_genome[i]
        elif queries_by_genome[i].bit_score < best_query.bit_score:
            if not second_best_query or queries_by_genome[i].bit_score > second_best_query.bit_score:
                second_best_query = queries_by_genome[i]
    if second_best_query:
        best_query.diff = abs(
            ((best_query.bit_score - second_best_query.bit_score) / (best_query.bit_score + second_best_query.bit_score)) / 2) * 100

    return best_query

def allele_in_genome_comp(genome_dict, allele):
    if allele in genome_dict:
        return genome_dict[allele].bit_score
    return -1



def best_queries(best_query_per_genome):
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



        # Press the green button in the gutter to run the script.
if __name__ == '__main__':
    all_queries = read_blast_output(sys.argv[1])
    best_query_per_genome = best_queries_per_genome(all_queries)
    working_query_set = best_queries(best_query_per_genome)

    # TODO:
    # send working query set through filter fasta
    # or....create my own filterfasta????


# See PyCharm help at https://www.jetbrains.com/help/pycharm/
