from DIGITfiles.Query import Query
import json

queriesWorkingSet = {}


def getDataFromJSON(filename):
    '''
    Parse a JSON file and return a JSON object (key:value format)
    '''
    # TODO: can i pulls this out and put it elsewhere? I keep reusing it
    jsonFile = open(filename)
    dataFromJSON = json.load(jsonFile)
    jsonFile.close()

    return dataFromJSON


def createQueryStruct(filename):
    dataFromJSON = getDataFromJSON(filename)
    for allele in dataFromJSON:
        if allele not in queriesWorkingSet:
            queriesWorkingSet[allele] = Query(None, None, None)
            queriesWorkingSet[allele].__QueryFromJSON__(dataFromJSON[allele])


def readPrimer3Output(filename):
    with open(filename, "r") as p3file:
        outputDict = {}  # create a dictionary for each output item
        for line in p3file:
            if line[0] != "=":  # go through the file until reaching a break mark ("=")
                lineItems = line.split("=")  # key:value pairs for OUTPUTLABEL:data
                outputDict[lineItems[0].strip()] = lineItems[1].strip()
            elif line[0] == "=":  # end of a section,time to process data
                allele = outputDict["SEQUENCE_ID"][:-1]

                queriesWorkingSet[allele].__updateQueryPrimer3__(outputDict)
                outputDict.clear()
                # print(output_dict)
            else:
                print("ERROR")


# TODO: is this generalizable? similar function appears in another part of program
def queriesToJSON(filename):
    '''
    Converts dictionary to JSON object and writes it to a file.
    '''
    jsonObject = {}
    for q in queriesWorkingSet:
        if queriesWorkingSet[q].query not in jsonObject:
            jsonObject[queriesWorkingSet[q].query] = []
        jsonObject[queriesWorkingSet[q].query] = dict(queriesWorkingSet[q])
    newJSONobject = json.dumps(jsonObject, indent=4)
    with open(filename, "w") as outfile:
        outfile.write(newJSONobject)

    return


def main():
    # get all the primer data
    createQueryStruct("WorkingQueriesWithPrimers03.json")
    # read_primer3_output("Primer3Output08.txt")

    # make an input filefor primer3
    with open("VerifyPrimer3Input.txt", "w") as p3file:
        for allele in queriesWorkingSet:
            p3file.write(queriesWorkingSet[allele].__createPrimer3Input__())

    queriesToJSON("")


if __name__ == "__main__":
    main()



    import sys
import os
parent_dir_name = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(parent_dir_name + "/part1")
import support
#from part1.support import Primer, read_csv, get_datetime

def wt_verification_length(flankseq):
    filepath = os.listdir("data/primer3/verification_output/")
    for afile in filepath:
        if afile.find(flankseq) != -1:
            wt_dict = {}
            wt_len, query = "", ""
            new_filepath = os.path.join("data/primer3/verification_output", afile)
            with open(new_filepath, "r+") as wt_file:
                for line in wt_file:
                    lineparts = line.split("=")
                    if lineparts[0] == "SEQUENCE_ID":
                        query = lineparts[1].strip()
                    if lineparts[0] == "PRIMER_PAIR_0_PRODUCT_SIZE":
                        wt_len = lineparts[1].strip()
                    if lineparts[0] == "":
                        wt_dict[query] = wt_len
    return wt_dict

def set_wt_var(all_primers, wt_product_lens):
    # for each primer object in a bucket
    for primer in all_primers:
        # is the value primer.query a key in wt_product_lens?
        if primer.query in wt_product_lens:
            # if so, set productWT_verify_size to the value at that key
            primer.productWT_verify_size = wt_product_lens[primer.query]


def read_primer_files(list_of_files):
    bucketlist = [[],[],[],[],[]]
    for i in range(0,len(list_of_files)):
        with open(list_of_files[i], "r+") as infile:
            # reaad csv and save parts as a new Primer object
            infile.readline() # header

            for line in infile:
                #new_primer = Primer()
                lineparts = line.split(',')
                #new_primer.primer_name = lineparts[0]
                new_primer = support.Primer(lineparts[1])
                new_primer.primer_name = lineparts[0]
                new_primer.sequence = lineparts[2]
                new_primer.product_size = lineparts[3]
                new_primer.tm = lineparts[5]
                new_primer.primer_pair_penalty = lineparts[6]
                new_primer.primer_left_penalty = lineparts[7]
                new_primer.primer_right_penalty = lineparts[8].strip()

                bucketlist[i].append(new_primer)
    return bucketlist

# new_primer_names is a list of strings
def append_to_completed_list(completed_primers):
    with open("resources/part3/completed_primers.txt", "a+") as primer_file:
        for query in completed_primers:
            primer_file.write(query + "\n")


def get_unver_parsed_primer3_files(path_to_unver_output):
    list_of_files = []

    for subdir, dirs, files in os.walk(sys.argv[1]):
        for one_file in files:
            filename = os.path.join(subdir, one_file)
            list_of_files.append(filename)

    return list_of_files


def create_primer_dataset_driver(path_to_parsed_unver_output):
    flankseq = path_to_parsed_unver_output.split('/')[3]

    list_unverified_files = get_unver_parsed_primer3_files(path_to_parsed_unver_output)

    # parse verification primer3 results
    wt_product_lens = wt_verification_length(flankseq)
    bucketlist = read_primer_files(list_unverified_files)

    # set wt length in primer obj
    # run set_wt_var for each bucket
    for bucket in bucketlist:
        set_wt_var(bucket, wt_product_lens)

    # list of lists; each inner list is a line of data from csv
    # read_csv returns a list of lines (strings) in the csv
    lookup_data = support.read_csv("data/lookup_data/lookup_" + flankseq + ".csv", True)
    dict_lookup_data = {}
    for line in lookup_data:
        if line[0] == "\n":
            continue
        else:
            dict_lookup_data[line[0]] = {
                "database":line[12],
                "bit_score":float(line[11]),
                "num_hits":int(line[13])
                }

    completed_primers = []
    for bucket in bucketlist:
        for primer in bucket:
            if primer.query in dict_lookup_data:

                primer.database = dict_lookup_data[primer.query]["database"]
                primer.bit_score = dict_lookup_data[primer.query]["bit_score"]
                primer.num_hits = dict_lookup_data[primer.query]["num_hits"]

    time = support.get_datetime()
    for i in range(0, len(bucketlist)):
        with open("data/primer3/parsed_output_verified/complete_" + flankseq + "_set" + str(i+1) + "_primers_" + time + ".csv", "w+") as outfile:
            outfile.write("primer_name,query,database,sequence,product_size,productWT_verify_size,tm,primer_pair_penalty,primer_left_penalty,primer_right_penalty,bit_score,num_hits\n")
            for primer in bucketlist[i]:
                outfile.write(primer.csv_format_full_data_verified())
                completed_primers.append(primer.primer_name)
    append_to_completed_list(completed_primers)


def main():
    create_primer_dataset_driver(sys.argv[1])
    # sys.argv[1] is the path to the unverified csv data from the first round
    # with primer3 --> data/primer3/parsed_output_unverified/<flankseq>/

    #list_of_files = []

    #for subdir, dirs, files in os.walk(sys.argv[1]):
        #for one_file in files:
            #filename = os.path.join(subdir, one_file)
            #list_of_files.append(filename)
    '''
    flankseq = sys.argv[1].split('/')[3]

    list_unverified_files = get_unver_parsed_primer3_files(sys.argv[1])

    # parse verification primer3 results
    wt_product_lens = wt_verification_length(flankseq)
    bucketlist = read_primer_files(list_unverified_files)

    # set wt length in primer obj
    # run set_wt_var for each bucket
    for bucket in bucketlist:
        set_wt_var(bucket, wt_product_lens)

    # list of lists; each inner list is a line of data from csv
    # read_csv returns a list of lines (strings) in the csv
    lookup_data = support.read_csv("data/lookup_data/lookup_FlankingSequences_SingleGenomicSequences.csv", True)
    dict_lookup_data = {}
    for line in lookup_data:
        dict_lookup_data[line[0]] = {
            "database":line[12],
            "bit_score":float(line[11]),
            "num_hits":int(line[13])
            }

    completed_primers = []
    for bucket in bucketlist:
        for primer in bucket:
            if primer.query in dict_lookup_data:

                primer.database = dict_lookup_data[primer.query]["database"]
                primer.bit_score = dict_lookup_data[primer.query]["bit_score"]
                primer.num_hits = dict_lookup_data[primer.query]["num_hits"]

    time = support.get_datetime()
    for i in range(0, len(bucketlist)):
        with open("data/primer3/complete_FlankingSequences_SingleGenomicSequences_set" + str(i+1) + "_primers_" + time + ".csv", "w+") as outfile:
            outfile.write("primer_name,query,database,sequence,product_size,productWT_verify_size,tm,primer_pair_penalty,primer_left_penalty,primer_right_penalty,bit_score,num_hits\n")
            for primer in bucketlist[i]:
                outfile.write(primer.csv_format_full_data_verified())    
                completed_primers.append(primer.primer_name)
    append_to_completed_list(completed_primers)
    '''
if __name__ == '__main__':
    main()