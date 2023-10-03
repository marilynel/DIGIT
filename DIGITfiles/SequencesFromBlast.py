import sys

from QueriesWorkingSet import QueriesWorkingSet
from Query import Query
from FilterFasta import filterFasta
from Sequences import Sequences
from Primer3Object import Primer3Object
from PrimerIDs import PrimerIDs
from Utils import *

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


def main():
    flankseq = sys.argv[1]
    b73only = int(sys.argv[2])
    outputDir = "DIGIToutput/FlankingSequences"

    queriesWorkingSet = QueriesWorkingSet()

    newFlankseq = flankseq
    genomeSpec = ""
    if b73only == 2:
        genomeSpec = "B73"
        newFlankseq += "B73Only"

    necessaryDirs = [
        f"{outputDir}",
        f"{outputDir}/{newFlankseq}",
        f"{outputDir}/{newFlankseq}/QueryData",
        f"{outputDir}/{newFlankseq}/QueryData/JSONfiles",
        f"{outputDir}/{newFlankseq}/QueryData/CSVfiles",
        f"{outputDir}/{newFlankseq}/QueryData/GFFfiles",
        f"{outputDir}/{newFlankseq}/QueryData/Primer3",
        f"{outputDir}/{newFlankseq}/QueryData/Primer3/FindingPrimers",
        f"{outputDir}/{newFlankseq}/QueryData/Primer3/VerifyingPrimers"
    ]

    makeDirectories(necessaryDirs)

    findBlastOutputFiles(f"{outputDir}/{flankseq}/BlastOutput/FlankingSequenceBlast", allQueries,
                         genomeSpec)

    listQueries = setBestForGenome(allQueries)

    writeToBestQueriesFile(listQueries, flankseq, allQueries,
                           f"{outputDir}/{newFlankseq}/QueryData/CSVfiles/BestQueriesByGenome_"
                           f"{newFlankseq}.csv")

    allQueriesToJSON(
        f"{outputDir}/{newFlankseq}/QueryData/JSONfiles/AllBlastData_{newFlankseq}.json",
        allQueries)

    queriesWorkingSet = buildWorkingSetFromAllQueries(allQueries)

    runFilterfasta(queriesWorkingSet)

    queriesWorkingSet.__setPrimerIDs__()

    inputFile = f"{outputDir}/{newFlankseq}/QueryData/Primer3/FindingPrimers/Primer3Input_" \
                f"{newFlankseq}.txt"

    queriesWorkingSet.__writeP3InputFile__(inputFile, "generic")

    runPrimer3(inputFile, newFlankseq,
               f"{outputDir}/{newFlankseq}/QueryData/Primer3/FindingPrimers/Primer3Output_{
    newFlankseq}.txt")

    queriesWorkingSet.__printToJson__(f"{outputDir}/{newFlankseq}/QueryData", newFlankseq)


if __name__ == '__main__':
    main()
