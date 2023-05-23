import subprocess
from Query import Query
from FilterFasta import filterFasta
from Sequences import Sequences
from Primer3Object import Primer3Object
from PrimerIDs import PrimerIDs
from Utils import *
import os
import sys
import json
import glob
import time

# allQueries is a dict with structure:
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


allQueries = {
    "A188v1": {},
    "B73v5": {},
    "W22v2": {}
}

queriesWorkingSet = {}


def findBlastOutputFiles(dirname):
    '''
    This function locates the Blast output file for the user selected set of flanking sequences.
    '''
    for subdir, dirs, files in os.walk(dirname):
        for oneFile in files:
            filename = os.path.join(subdir, oneFile)
            if filename.find(".tab") != -1:
                processBlastOutput(filename)

    return


def processBlastOutput(filename):
    '''
    Parse blast output files line by line, gathering data that will be used to create Query
    objects, which will be added
    to the allQueries dataset. Functions get_num_hits() and get_genome() are called to help parse
    hash (#) lines.
    '''
    genome = ""
    num_hits = -1
    with open(filename, "r+") as blastfile:
        for line in blastfile:
            if line[0] == '#':
                genome = getGenome(line, genome)
                num_hits = getNumHits(line, num_hits)
            else:
                newQuery = Query(line, genome, num_hits)
                newQuery.__setValues__()
                if newQuery.query not in allQueries[genome]:
                    allQueries[genome][newQuery.query] = []
                allQueries[genome][newQuery.query].append(newQuery)
    return


def getGenome(line, genome):
    '''
    Parse a line for string indicating the genome name. If line does not have that string,
    return genome as original
    value.
    '''
    blastData = line.split(' ')
    if blastData[1] == "Database:":
        return blastData[2].split('/')[-1].strip()
    return genome


def getNumHits(line, numHits):
    '''
    Parse a line for string indicating the number of hits for that query in a genome. If line
    does not have that string,
    return num_hits as original value.
    '''
    blastData = line.split(' ')
    if blastData[2] == "hits":
        return int(blastData[1])
    return numHits


def setBestForGenome():
    '''
    Create dictionary identifying the best BLAST hit for each query, within each genome.
    '''
    listQueries = []

    for genome in allQueries:
        for allele in allQueries[genome]:
            if allele not in listQueries:
                listQueries.append(allele)
            getBestQuery(genome, allele)

    return listQueries


def makeStrBestInGenome(genome, allele, best_bit_score):
    '''
    Find the specific hit that was specified as the best for that query ID and that genome (
    Query.best_for_genome ==
    True) and return a string reporting data for the best_queries_by_genome.csv outfile.
    '''
    if allele in allQueries[genome]:
        for i in range(0, len(allQueries[genome][allele])):
            if allQueries[genome][allele][i].bestAlleleForGenome == True:
                best_bit_score[genome] = allQueries[genome][allele][i].bitScore
                return allQueries[genome][allele][i].__makeBestGenomeString__()

    return "none,none,none,"


def bestGenomesWrite(bestBitScore):
    '''
    Find the best of the best hits, return strinigied list for writing to file.
    '''
    genomeList = []
    bestScore = max(bestBitScore.values())
    for genome in bestBitScore:
        if bestBitScore[genome] == bestScore:
            genomeList.append(genome)
        else:
            bestBitScore[genome] = 0
    return str(genomeList)


def workingQuerySelection(genome, allele):
    '''
    set ID for best_query --> indicates that this will belong to working set
    '''
    for i in range(0, len(allQueries[genome][allele])):
        if allQueries[genome][allele][i].bestAlleleForGenome == True:
            allQueries[genome][allele][i].bestHitForAllele = True
            return


def pickGenome(allele, bestBitScore):
    '''
    TODO: why did I write this?
    '''
    if bestBitScore["B73v5"] != 0:
        workingQuerySelection("B73v5", allele)
        return
    elif bestBitScore["W22v2"] != 0:
        workingQuerySelection("W22v2", allele)
        return
    elif bestBitScore["A188v1"] != 0:
        workingQuerySelection("A188v1", allele)
    else:
        print("you really messed something up")


def writeToBestQueriesFile(listQueries, flankseq):
    '''
    Write to file the best bit scores per genome per query, for manual review if needed.
    '''
    with open("DIGIToutput/" + flankseq + "/DataStats/BestQueriesByGenome_" + flankseq + ".csv",
              "w+") as newfile:
        newfile.write(
            f"query,A188_bit_score,A188_num_hits,A188_qstart_status,B73_bit_score,B73_num_hits,"
            f"B73_qstart_status,W22_" +
            f"bit_score,W22_num_hits,W22_qstart_status,best_genomes\n"
        )
        for allele in listQueries:
            bestBitScores = {
                "A188v1": 0,
                "B73v5": 0,
                "W22v2": 0
            }
            newfile.write(
                allele + "," +
                makeStrBestInGenome("A188v1", allele, bestBitScores) +
                makeStrBestInGenome("B73v5", allele, bestBitScores) +
                makeStrBestInGenome("W22v2", allele, bestBitScores) +
                bestGenomesWrite(bestBitScores) +
                "\n"
            )
            pickGenome(allele, bestBitScores)
    return


def getBestQuery(genome, query):
    '''
    Sort hits to identify best hit for a query. Find the percent difference between the best and
    second best hits.
    '''
    bestQuery = allQueries[genome][query][0]  # Query
    secondBestQuery = None  # Query
    for i in range(1, len(allQueries[genome][
                              query])):  # iterate through Query objects in allQueries[genome][
        # allele]
        ### NEW: changed > to >= make sure it works before deleting this comment!!!! 5/22 ###
        ### NOTE: this does change the output from before. yikes. still valid but the results are
        # inconsistant...going to put the = on the "less than" comparison
        ### note again: that fixed it.
        if allQueries[genome][query][
            i].bitScore > bestQuery.bitScore:  # if the new bit score is better than bestQuery's
            # bit score
            secondBestQuery = bestQuery  # redefine secondBestQuery as Query in bestQuery
            bestQuery = allQueries[genome][query][i]  # redefine bestQuery as the current Query

        elif allQueries[genome][query][
            i].bitScore <= bestQuery.bitScore:  # if the new bit score is with than bestQuery's
            # bit score
            ###  NEW: changed > to >= make sure it works before deleting this comment!!!! 5/22 ###
            if not secondBestQuery or allQueries[genome][query][
                i].bitScore >= secondBestQuery.bitScore:  # if either (there is no
                # secondBestQuery) OR (new query's bit score is better than the secondBestWuery
                # bit score)
                secondBestQuery = allQueries[genome][query][
                    i]  # redefine second best query as the new query
    # TODO: come back and make this work
    perDiff = -1
    if secondBestQuery:  # if there is a second-best query
        perDiff = abs(
            (
                    ((bestQuery.bitScore - secondBestQuery.bitScore) / (
                                bestQuery.bitScore + secondBestQuery.bitScore)) / 2
            ) * 100
        )

    for i in range(0, len(allQueries[genome][query])):
        if allQueries[genome][query][i] == bestQuery:
            allQueries[genome][query][i].bestAlleleForGenome = True
            allQueries[genome][query][i].percentDiff = perDiff

    return


def allQueriesToJSON(filename):
    '''
    Converts dictionary to JSON object and writes it to a file.
    '''
    jsonObject = {
        "A188v1": {},
        "B73v5": {},
        "W22v2": {}
    }

    for database in allQueries:
        for query in allQueries[database]:
            if query not in jsonObject[database]:
                jsonObject[database][query] = []
            for hit in allQueries[database][query]:
                jsonObject[database][query].append(dict(hit))
    newJSONbject = json.dumps(jsonObject, indent=4)

    with open(filename, "w") as outfile:
        outfile.write(newJSONbject)

    return


def buildQueriesWorkingSet():
    '''
    Iterates through entire dataset, and creates a new smaller working dataset with just the
    "best" hit for each query.
    '''
    for gen in allQueries:
        for q in allQueries[gen]:
            for i in range(0, len(allQueries[gen][q])):
                if allQueries[gen][q][i].bestHitForAllele:
                    queriesWorkingSet[allQueries[gen][q][i].query] = allQueries[gen][q][i]


# TODO: is there really no way to make this prettier??
def buildCoordinateSetsForFilterFasta():
    '''
    This function creates a dataset from queriesWorkingSet that can be passed to FilterFasta.py
    in the following format:

        coordinates = {
            alleleName0 : [
                chromosomeID0, [startingPosition0, endingPosition0]
            ],
            alleleName1 : [
                chromosomeID1, [startingPosition1, endingPosition1]
            ],
            etc.
        }

    '''
    coorA, coorB, coorW = {}, {}, {}
    for q in queriesWorkingSet:
        if queriesWorkingSet[q].genome.strip() == "A188v1":
            coorA[queriesWorkingSet[q].query + "_wt"] = [
                queriesWorkingSet[q].chromosome, queriesWorkingSet[q].wildtypeCoordinates
            ]
            coorA[queriesWorkingSet[q].query + "_up"] = [
                queriesWorkingSet[q].chromosome, queriesWorkingSet[q].upperCoordinates
            ]
            coorA[queriesWorkingSet[q].query + "_lo"] = [
                queriesWorkingSet[q].chromosome, queriesWorkingSet[q].lowerCoordinates
            ]
        elif queriesWorkingSet[q].genome.strip() == "B73v5":
            coorB[queriesWorkingSet[q].query + "_wt"] = [
                queriesWorkingSet[q].chromosome, queriesWorkingSet[q].wildtypeCoordinates
            ]
            coorB[queriesWorkingSet[q].query + "_up"] = [
                queriesWorkingSet[q].chromosome, queriesWorkingSet[q].upperCoordinates
            ]
            coorB[queriesWorkingSet[q].query + "_lo"] = [
                queriesWorkingSet[q].chromosome, queriesWorkingSet[q].lowerCoordinates
            ]
        elif queriesWorkingSet[q].genome.strip() == "W22v2":
            coorW[queriesWorkingSet[q].query + "_wt"] = [
                queriesWorkingSet[q].chromosome, queriesWorkingSet[q].wildtypeCoordinates
            ]
            coorW[queriesWorkingSet[q].query + "_up"] = [
                queriesWorkingSet[q].chromosome, queriesWorkingSet[q].upperCoordinates
            ]
            coorW[queriesWorkingSet[q].query + "_lo"] = [
                queriesWorkingSet[q].chromosome, queriesWorkingSet[q].lowerCoordinates
            ]

    return coorA, coorB, coorW


def runFilterfasta():
    coorA, coorB, coorW = buildCoordinateSetsForFilterFasta()

    seqDataA = filterFasta("DIGITfiles/Genomes/A188v1/Zm-A188-REFERENCE-KSU-1.0.fa", coorA)
    seqDataB = filterFasta("DIGITfiles/Genomes/B73v5/Zm-B73-REFERENCE-NAM-5.0.fa", coorB)
    seqDataW = filterFasta("DIGITfiles/Genomes/W22v2/Zm-W22-REFERENCE-NRGENE-2.0.fa", coorW)

    parseFilterfastaData(seqDataA)
    parseFilterfastaData(seqDataB)
    parseFilterfastaData(seqDataW)


def parseFilterfastaData(seqData):
    for allele in seqData:
        upper, lower, wildtype = "", "", ""
        if allele.find("wt") != -1:
            queriesWorkingSet[allele[:-3]].wildtypeSequence = seqData[allele]
        if allele.find("up") != -1:
            queriesWorkingSet[allele[:-3]].upperSequence = seqData[allele]
        if allele.find("lo") != -1:
            queriesWorkingSet[allele[:-3]].lowerSequence = seqData[allele]
    for allele in queriesWorkingSet:
        if queriesWorkingSet[allele].upperSequence and queriesWorkingSet[allele].lowerSequence:
            queriesWorkingSet[allele].__buildInsertionSequence__()
        else:
            queriesWorkingSet[allele].insertionSequence = "__failed__"


# TODO: move to utils?
def makeDirectories(flankseq):
    necessaryDirs = [
        f"DIGIToutput/{flankseq}",
        f"DIGIToutput/{flankseq}/DataSets",
        f"DIGIToutput/{flankseq}/DataStats",
        f"DIGIToutput/{flankseq}/Primer3Files",
        f"DIGIToutput/{flankseq}/Primer3Files/Input",
        f"DIGIToutput/{flankseq}/Primer3Files/Output",
        f"DIGIToutput/{flankseq}/WTVerification",
        f"DIGIToutput/{flankseq}/WTVerification/Input",
        f"DIGIToutput/{flankseq}/WTVerification/Output"
    ]
    for dir in necessaryDirs:
        if not os.path.exists(dir):
            os.makedirs(dir)


def main():
    flankseq = sys.argv[1]
    makeDirectories(flankseq)

    findBlastOutputFiles(f"DIGITfiles/BlastOutput/{flankseq}")
    listQueries = setBestForGenome()
    writeToBestQueriesFile(listQueries, flankseq)

    allQueriesToJSON(f"DIGIToutput/" + flankseq + "/DataSets/AllBlastData_" + flankseq + ".json")
    buildQueriesWorkingSet()

    runFilterfasta()

    setPrimerIDsObj = PrimerIDs()
    setPrimerIDsObj.__setPrimerIDs__(queriesWorkingSet)

    inputFile = writeP3InputFile(
        "DIGIToutput/" + flankseq + "/Primer3Files/Input/Primer3Input_" + flankseq + ".txt",
        queriesWorkingSet,
        flankseq,
        "generic"
    )

    runPrimer3(
        inputFile,
        flankseq,
        "DIGIToutput/" + flankseq + "/Primer3Files/Output/Primer3Output_" + flankseq + ".txt"
    )

    queriesToJSON(f"DIGIToutput/{flankseq}/DataSets/WorkingSet_{flankseq}.json", queriesWorkingSet)


if __name__ == '__main__':
    main()
