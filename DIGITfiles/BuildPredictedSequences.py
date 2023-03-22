'''
BuildPredictedSequences.py
Marilyn Leary 2023

WARNING: make sure all of your Part 1 jobs are complete before running this program!
Use "qstat" on the command line to see the status of your jobs.

This program has three main functions:
-   creating the predicted DsGFP insertion sequences for each allele
-   creates the Primer3 input file using those insertions
-   calls Primer3 on the SGE queue system

Part 1 will have had to be completed in order for this program to work. The resultant
fasta files will be in directory DIGITfiles/Fastafiles/. If they cannot be found there,
check the corresponding sge folders for error messages.

Output:
-   WorkingQueries_<flankseq>.json
-   Primer3Input_<flankseq>.txt
-   Primer3Output_<flankseq>.txt
'''


import subprocess
from Query import Query
from Sequences import Sequences
from Primer3Object import Primer3Object
import os
import sys
import json
import glob

queriesWorkingSet = {}


def findFilterfastaOutputFiles(dirname):
    # TODO: I think the error is occuring here
    #fasta_dict = {}
    for subdir, dirs, files in os.walk(dirname):
        for one_file in files:
            filename = os.path.join(subdir, one_file)
            if filename.find(".fasta") != -1:
                rfile = filename.split("/")[1]
                allele = rfile.split(".")[0]
                # a little worried this won't work, but there may be more fasta files than are actually needed
                if allele in queriesWorkingSet:
                    storeAndConstructSequences(filename, allele)
    return


def getDataFromJSON(filename):
    '''
    Parse a JSON file and return a JSON object (key:value format)
    '''
    jsonFile = open(filename)
    dataFromJSON = json.load(jsonFile)
    jsonFile.close()

    return dataFromJSON


def createQueryStruct(filename):
    '''
    Initialize Query objects from JSON file and add to workingQueriesSet.
    '''
    dataFromJSON = getDataFromJSON(filename)
    for genome in dataFromJSON:
        for allele in dataFromJSON[genome]:
            for i in range(0, len(dataFromJSON[genome][allele])):
                if dataFromJSON[genome][allele][i]["bestHitForAllele"]:
                    queriesWorkingSet[allele] = Query(None, None, None)
                    queriesWorkingSet[allele].__QueryFromJSON__(dataFromJSON[genome][allele][i])


def storeAndConstructSequences(filename, allele):
    # TODO: issues here; see above note
    fastaDict = {
        "upper": "",
        "lower": "",
        "wildtype": ""
    }
    print(fastaDict)
    with open(filename, "r") as fastafile:
        for line in fastafile:
            if line[0] == ">":
                position = line.split("_")[-1]  # upper, lower, or wildtype
            else:
                fastaDict[position.strip()] = line.strip()
        # TODO: may change the below to be handled in Query method\
        # TODO: above is done, just needs to be tested
        #print(allele)
        queriesWorkingSet[allele].__assignSeqsToQuery__(fastaDict["upper"], fastaDict["lower"], fastaDict["wildtype"])
        #queriesWorkingSet[allele].upperSequence = fastaDict["upper"]
        #queriesWorkingSet[allele].lowerSequence = fastaDict["lower"]
        #queriesWorkingSet[allele].wildtypeSequence = fastaDict["wildtype"]
        #queriesWorkingSet[allele].wildtypeSequence = fastaDict["wildtype"]
        #queriesWorkingSet[allele].__buildInsertionSequence__()


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


def runPrimer3(inputFile):
    subprocess.run(
        [
            "SGE_Batch",
            "-c",
            "primer3_core " + inputFile + " > DIGITfiles/Primer3Output.txt",
            "-q",
            "bpp",
            "-P",
            "8",
            "-r",
            "sge.primer3"
        ]
    )

def createPrimer3Input(flankseq):
    with open("DIGITfiles/Primer3Input" + flankseq + ".txt", "w") as p3file:
        for allele in queriesWorkingSet:
            leftP3 = Primer3Object("generic", "left", queriesWorkingSet[allele])
            rightP3 = Primer3Object("generic", "right", queriesWorkingSet[allele])
            inputStr = leftP3.input_str + rightP3.input_str
            p3file.write(inputStr)


def main():
    createQueryStruct("DIGITfiles/" + sys.argv[1])
    findFilterfastaOutputFiles("DIGITfiles/filterfasta_files")
    createPrimer3Input(sys.argv[1][12:])
    queriesToJSON("DIGITfiles/WorkingQueries.json")
    runPrimer3("DIGITfiles/Primer3Input" + sys.argv[1][12:] + ".txt")


if __name__ == "__main__":
    main()
