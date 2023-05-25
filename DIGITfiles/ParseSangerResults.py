import sys
import os

from Query import Query
from Utils import *

sangerQueriesAll = {
    "A188v1": {},
    "B73v5": {},
    "W22v2": {}
}

bestSangerQueries = {}
queriesWorkingSet = {}


def main():
    flankseq = "SangerSeq_Order24639"

    necessaryDirs = [
        f"DIGIToutput/{flankseq}",
        f"DIGIToutput/{flankseq}/DataSets",
        f"DIGIToutput/{flankseq}/DataStats",
    ]

    makeDirectories(necessaryDirs)

    # 1. parse blast output, save as queries in struct
    findBlastOutputFiles(f"DIGITfiles/BlastOutputSangerSeqs/{flankseq}", sangerQueriesAll)

    # 2. get best queries to compare
    sangerList = setBestForGenome(sangerQueriesAll)

    writeToBestQueriesFile(
        sangerList,
        flankseq,
        sangerQueriesAll,
        "DIGIToutput/" + flankseq + "/DataStats/BestQueriesByGenome_" + flankseq + ".csv"
    )

    # 3. Save as json for reference
    allQueriesToJSON(f"DIGIToutput/" + flankseq + "/DataSets/AllBlastData_" + flankseq + ".json",
                     sangerQueriesAll)

    # 4. Pull from Working Sets of Queries
    workingSetFilename = "DIGIToutput/FlankingSequences_SingleGenomicSequences/DataSets" \
                         "/WorkingSet_FlankingSequences_SingleGenomicSequences.json"
    createQueryStruct(workingSetFilename, queriesWorkingSet)

    # 5. you forgot to make a snger working set!
    buildWorkingSetFromAllQueries(sangerQueriesAll, bestSangerQueries)

    # 6. compare sangers and originals, mote differences

    for query in bestSangerQueries:
        if query in queriesWorkingSet:
            print(
                f"{bestSangerQueries[query].query}\t{bestSangerQueries[query].genome}\t{queriesWorkingSet[query].genome}")
        else:
            print(f"{bestSangerQueries[query].query} not in queries working set")


'''
TODO: 
some of these queries are  in a different set of alleles...how to account for that?
none yet, but eventually there will be an imperfect match allele, need to account for that
'''

if __name__ == '__main__':
    main()