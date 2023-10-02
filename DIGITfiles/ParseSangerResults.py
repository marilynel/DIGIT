
'''
Sanger Sequence file name structure:

28247_A01_R91F11_DsGFP3_UTR_A01_008.seq

28247       -> order number (should be the same for all .seq files in a folder)
A01         -> position in well plate
R91F11      -> allele/query/R number
DsGFP3_UTR  -> matching primer
'''

import sys
import os
from collections import ChainMap

from Utils import *

sangerQueriesAll = {
    "A188v1": {},
    "B73v5": {},
    "W22v2": {}
}




def sameQuery(queryS, queryOG):
    '''
    Given a Query object from the sanger set and a query object from the original set, check if the metadata of the
    Query is the same for genome, database, or sequence start position. Print if anythin is mismatched and return False.
    Return True if the Query objects match.
    '''
    # if queryS.query != queryOG.query:
    #    print("how did this get called?")
    #    return None

    if (queryS.genome != queryOG.genome) or (queryS.chromosome != queryOG.chromosome) or \
            (queryS.sStart != queryOG.sStart):
        # print(f"{queryS.query}")
        # print(f"Sanger:\t\t{queryS.genome}\t{queryS.chromosome}\t{queryS.sStart}")
        # print(f"Original:\t{queryOG.genome}\t{queryOG.chromosome}\t{queryOG.sStart}")
        return False
    return True


# TODO: file str
def buildFullWorkingSet():
    '''
    Collect all query data from existing working sets and save to a dictionary. Return dictionary.
    '''
    listOfSets = []
    for subdir, dirs, files in os.walk("DIGIToutput/SangerSequences"):
        if subdir.find("DataSets") != -1:
            for file in files:
                if file.find("Working") != -1 and file.find("Sanger") == -1 and file.find(
                        "json") != -1:
                    newSet = {}
                    createQueryStruct(subdir + "/" + file, newSet)
                    listOfSets.append(newSet)
    return dict(ChainMap(*listOfSets))


def sortSangerQueries(originalData, bestSangerQueries):
    data = []
    badData = []

    for query in bestSangerQueries:

        queryOG = query.split("_")[0]

        if queryOG in originalData:
            if sameQuery(bestSangerQueries[query], originalData[queryOG]):
                # print("add to data")
                data.append(query)
            else:
                # print("add to badData")
                badData.append(query)
        else:
            i = 1
            # print(f"{bestSangerQueries[query].query} not in queries working set")
    return data, badData


def writeSangerComparisonFile(filename, dataList, originalData, bestSangerQueries):
    with open(filename, "w+") as f:
        f.write(
            f"Query,RefGenomeSanger,ChromosomeSanger,RefGenomeOriginal,ChromosomeOriginal,SeqStartSanger,SeqEndSanger" +
            f",SeqStartOriginal,SeqEndOriginal,EValueSanger,EValueOriginal,BitScoreSanger,BitScoreOriginal,NumHitsSan" +
            f"ger,NumHitsOriginal,StrandSanger,StrandOriginal,QStartStatusSanger,QStartStatusOriginal\n")
        for query in dataList:
            queryOG = query.split("_")[0]
            f.write(
                f"{bestSangerQueries[query].query},{bestSangerQueries[query].genome}," +
                f"{bestSangerQueries[query].chromosome},{originalData[queryOG].genome}," +
                f"{originalData[queryOG].chromosome},{bestSangerQueries[query].sStart},{bestSangerQueries[query].sEnd}," +
                f"{originalData[queryOG].sStart},{originalData[queryOG].sEnd},{bestSangerQueries[query].eValue}," +
                f"{originalData[queryOG].eValue},{bestSangerQueries[query].bitScore},{originalData[queryOG].bitScore}," +
                f"{bestSangerQueries[query].numHits},{originalData[queryOG].numHits},{bestSangerQueries[query].strand}," +
                f"{originalData[queryOG].strand},{bestSangerQueries[query].qStartStatus}," +
                f"{originalData[queryOG].qStartStatus}\n"
            )


def writeToFile(dataList, flankseq, genomeDetail, seqType, originalData, bestSangers):
    if len(dataList) > 0:
        writeSangerComparisonFile(
            f"DIGIToutput/SangerSequences/{flankseq}/CSVfiles/{genomeDetail}_{flankseq}_{seqType}.csv",
            dataList,
            originalData,
            bestSangers
        )


def main():
    flankseq = sys.argv[1]

    necessaryDirs = [
        f"DIGIToutput/SangerSequences",
        f"DIGIToutput/SangerSequences/{flankseq}",
        f"DIGIToutput/SangerSequences/{flankseq}/JSONfiles",
        f"DIGIToutput/SangerSequences/{flankseq}/CSVfiles"
    ]

    makeDirectories(necessaryDirs)

    # 1. parse blast output, save as queries in struct
    findBlastOutputFiles(f"DIGIToutput/SangerSequences/{flankseq}/BlastOutput", sangerQueriesAll,
                         flankseq)

    matchingSeqs, nonMatchingSeqs = {}, {}

    for genome in sangerQueriesAll:
        if genome not in matchingSeqs:
            matchingSeqs[genome] = {}
        if genome not in nonMatchingSeqs:
            nonMatchingSeqs[genome] = {}
        for sangerSeq in sangerQueriesAll[genome]:
            if sangerSeq.find("_ImperfectMatchDsGfp") == -1:
                matchingSeqs[genome][sangerSeq] = sangerQueriesAll[genome][sangerSeq]
            else:
                nonMatchingSeqs[genome][sangerSeq] = sangerQueriesAll[genome][sangerSeq]

    # 2. get best queries to compare
    sangerListMatching = setBestForGenome(matchingSeqs)
    sangerListNonMatching = setBestForGenome(nonMatchingSeqs)

    writeToBestQueriesFile(
        sangerListMatching,
        flankseq,
        sangerQueriesAll,
        f"DIGIToutput/SangerSequences/{flankseq}/CSVfiles/BestQueriesByGenome_{flankseq}.csv"
    )

    writeToBestQueriesFile(
        sangerListNonMatching,
        flankseq,
        sangerQueriesAll,
        f"DIGIToutput/SangerSequences/{flankseq}/CSVfiles/BestQueriesByGenome_{flankseq}_ImperfectSeq.csv"
    )

    # 3. Save as json for reference
    allQueriesToJSON(
        f"DIGIToutput/SangerSequences/{flankseq}/JSONfiles/AllBlastData_{flankseq}.json",
        sangerQueriesAll)

    # 4. Pull from Working Sets of Queries
    originalData = buildFullWorkingSet()

    bestSangerQueriesMatching = {}
    bestSangerQueriesNonMatching = {}
    # 5. you forgot to make a snger working set!
    buildWorkingSetFromAllQueries(matchingSeqs, bestSangerQueriesMatching)
    buildWorkingSetFromAllQueries(nonMatchingSeqs, bestSangerQueriesNonMatching)

    # queriesToJSON(filepath, flankseq, queriesWorkingSet)
    queriesToJSON(
        f"DIGIToutput/SangerSequences/{flankseq}",
        flankseq,
        bestSangerQueriesMatching
    )
    queriesToJSON(
        f"DIGIToutput/SangerSequences/{flankseq}",
        flankseq,
        bestSangerQueriesNonMatching
    )

    # 6. compare sangers and originals, mote differences
    # query,sanger_order,sanger_database,sanger_chr,sanger_strand,sanger_s_start,sanger_s_end,sanger_align,sanger_q_start,sanger_q_end,sanger_bit_score,sanger_num_hits,sanger_per_identity,sanger_mismatches,sanger_evalue,sanger_per_diff,sanger_q_start_status,sanger_bit_score_status,sanger_gap,original_database,original_chr,original_strand,original_s_start,original_s_end,original_align,original_q_start,original_q_end,original_bit_score,original_num_hits,original_per_identity,original_mismatches,original_evalue,original_per_diff,original_q_start_status,original_gapalt_database,alt_chr,alt_strand,alt_s_start,alt_s_end,alt_align,alt_q_start,alt_q_end,alt_bit_score,alt_num_hits,alt_per_identity,alt_mismatches,alt_evalue,alt_per_diff,alt_q_start_status,alt_gap

    dataMatching, badDataMatching = sortSangerQueries(originalData, bestSangerQueriesMatching)
    dataNonMatching, badDataNonMatching = sortSangerQueries(originalData,
                                                            bestSangerQueriesNonMatching)

    # data includes sanger seq IDs that were found in original set
    # badData includes sanger seq IDs that were not found in any original set
    # NOTE: this does not mean it isn't there. I need to check for partial ID matches asthere may be a suffix on the ID
    # name, or the original might be in the two seqs set
    # I think this can be done in sameQuery()

    writeToFile(dataMatching, flankseq, "SameGenome", "MatchingSeq", originalData,
                bestSangerQueriesMatching)
    writeToFile(badDataMatching, flankseq, "DiffGenome", "MatchingSeq", originalData,
                bestSangerQueriesMatching)
    writeToFile(dataNonMatching, flankseq, "SameGenome", "NonMatchingSeq", originalData,
                bestSangerQueriesNonMatching)
    writeToFile(badDataNonMatching, flankseq, "DiffGenome", "NonMatchingSeq", originalData,
                bestSangerQueriesNonMatching)


'''
TODO: 
some of these queries are  in a different set of alleles...how to account for that?
none yet, but eventually there will be an imperfect match allele, need to account for that
'''

if __name__ == '__main__':
    main()