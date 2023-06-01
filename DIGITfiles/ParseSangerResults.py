import sys
import os
from collections import ChainMap

from Query import Query
from Utils import *

sangerQueriesAll = {
    "A188v1": {},
    "B73v5": {},
    "W22v2": {}
}

bestSangerQueries = {}
queriesWorkingSet = {}


def sameQuery(queryS, queryOG):
    if queryS.query != queryOG.query:
        print("how did this get called?")
        return None

    if (queryS.genome != queryOG.genome) or (queryS.chromosome != queryOG.chromosome) or (
            queryS.sStart != queryOG.sStart):
        print(f"Sanger:\t\t{queryS.genome}\t{queryS.chromosome}\t{queryS.sStart}")
        print(f"Original:\t{queryOG.genome}\t{queryOG.chromosome}\t{queryOG.sStart}")

        return False
    else:
        return True


def buildFullWorkingSet():
    listOfSets = []
    for subdir, dirs, files in os.walk("DIGIToutput"):
        if subdir.find("DataSets") != -1:
            # print(f"subdir: {subdir}")
            # print(f"dir: {dirs}")
            for file in files:
                if file.find("Working") != -1:
                    # print(f"file: {file}")
                    newSet = {}
                    createQueryStruct(subdir + "/" + file, newSet)
                    listOfSets.append(newSet)
    return dict(ChainMap(*listOfSets))


def main():
    flankseq = sys.argv[1]

    necessaryDirs = [
        f"DIGIToutput/{flankseq}",
        f"DIGIToutput/{flankseq}/DataSets",
        f"DIGIToutput/{flankseq}/DataStats",
        f"DIGIToutput/{flankseq}/Comparisons"
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
    # workingSetFilename = "DIGIToutput/FlankingSequences_SingleGenomicSequences/DataSets
    # /WorkingSet_FlankingSequences_SingleGenomicSequences.json"
    # createQueryStruct(workingSetFilename, queriesWorkingSet)
    queriesWorkingSet = buildFullWorkingSet()

    # 5. you forgot to make a snger working set!
    buildWorkingSetFromAllQueries(sangerQueriesAll, bestSangerQueries)
    queriesToJSON("DIGIToutput/" + flankseq + "/DataSets/WorkingSet_" + flankseq + ".json",
                  bestSangerQueries)

    # 6. compare sangers and originals, mote differences
    # query,sanger_order,sanger_database,sanger_chr,sanger_strand,sanger_s_start,sanger_s_end,sanger_align,sanger_q_start,sanger_q_end,sanger_bit_score,sanger_num_hits,sanger_per_identity,sanger_mismatches,sanger_evalue,sanger_per_diff,sanger_q_start_status,sanger_bit_score_status,sanger_gap,original_database,original_chr,original_strand,original_s_start,original_s_end,original_align,original_q_start,original_q_end,original_bit_score,original_num_hits,original_per_identity,original_mismatches,original_evalue,original_per_diff,original_q_start_status,original_gapalt_database,alt_chr,alt_strand,alt_s_start,alt_s_end,alt_align,alt_q_start,alt_q_end,alt_bit_score,alt_num_hits,alt_per_identity,alt_mismatches,alt_evalue,alt_per_diff,alt_q_start_status,alt_gap

    data = []
    badData = []

    for query in bestSangerQueries:
        if query in queriesWorkingSet:
            if sameQuery(bestSangerQueries[query], queriesWorkingSet[query]):
                data.append(query)
            else:
                badData.append(query)
        else:
            print(f"{bestSangerQueries[query].query} not in queries working set")

    print(f"Good Queries: {data}")
    with open("DIGIToutput/" + flankseq + "/Comparisons/Matching_" + flankseq + ".csv",
              "w+") as matchingFile:
        matchingFile.write(f"Query,RefGenome,Chromosome,SeqStart,SeqEndSanger,SeqEndOriginal\n")
        for query in data:
            matchingFile.write(
                f"{bestSangerQueries[query].query},{bestSangerQueries[query].genome}," +
                f"{bestSangerQueries[query].chromosome},{bestSangerQueries[query].sStart}," +
                f"{bestSangerQueries[query].sEnd},{queriesWorkingSet[query].sEnd}\n"
            )
    print(f"Bad Queries: {badData}")

    # Get data for Bad Queries


'''
TODO: 
some of these queries are  in a different set of alleles...how to account for that?
none yet, but eventually there will be an imperfect match allele, need to account for that
'''

if __name__ == '__main__':
    main()