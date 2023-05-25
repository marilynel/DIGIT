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


def main():
    flankseq = sys.argv[1]

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

    makeDirectories(necessaryDirs)

    findBlastOutputFiles(f"DIGITfiles/BlastOutput/{flankseq}", allQueries)

    listQueries = setBestForGenome(allQueries)

    writeToBestQueriesFile(
        listQueries,
        flankseq,
        allQueries,
        "DIGIToutput/" + flankseq + "/DataStats/BestQueriesByGenome_" + flankseq + ".csv"
    )

    allQueriesToJSON(f"DIGIToutput/" + flankseq + "/DataSets/AllBlastData_" + flankseq + ".json",
                     allQueries)

    buildWorkingSetFromAllQueries(allQueries, queriesWorkingSet)

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
