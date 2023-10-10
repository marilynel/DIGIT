'''
Sanger Sequence file name structure:

28247_A01_R91F11_DsGFP3_UTR_A01_008.seq

28247       -> order number (should be the same for all .seq files in a folder)
A01         -> position in well plate
R91F11      -> allele/query/R number
DsGFP3_UTR  -> matching primer
'''

import sys

from Utils import *
from SangerSeq import *


sangerQueriesAll = {
    "A188v1": {},
    "B73v5": {},
    "W22v2": {}
}


def main():
    ordernum = sys.argv[1]
    necessaryDirs = [
        f"DIGIToutput/SangerSequences",
        f"DIGIToutput/SangerSequences/{ordernum}",
        f"DIGIToutput/SangerSequences/{ordernum}/JSONfiles",
        f"DIGIToutput/SangerSequences/{ordernum}/GFFfiles",
        f"DIGIToutput/SangerSequences/{ordernum}/CSVfiles"
    ]

    makeDirectories(necessaryDirs)
    outSangerDir = f"DIGIToutput/SangerSequences/{ordernum}"

    findBlastOutputFiles(f"{outSangerDir}/BlastOutput", sangerQueriesAll, ordernum)

    allQueriesToJSON(f"{outSangerDir}/JSONfiles/AllBlastData_{ordernum}.json", sangerQueriesAll)

    ### NEW

    for genome in sangerQueriesAll:
        for allele in sangerQueriesAll[genome]:
            for i in range(len(sangerQueriesAll[genome][allele])):
                if sangerQueriesAll[genome][allele][i].numHits == 0:
                    print(f"{allele} does not appear in {genome}")



    ### END NEW


    matchingSeqs, nonMatchingSeqs = splitSangerQueriesAllIntoMatchingAndNonMatchingSets(sangerQueriesAll)

    # Get best queries to compare
    sangerListMatching = setBestForGenome(matchingSeqs)
    sangerListNonMatching = setBestForGenome(nonMatchingSeqs)

    if len(sangerListMatching) > 0:
        writeToBestQueriesFile(sangerListMatching, ordernum, sangerQueriesAll,
                               f"{outSangerDir}/CSVfiles/BestQueriesByGenome_{ordernum}.csv")

    if len(sangerListNonMatching) > 0:
        writeToBestQueriesFile(sangerListNonMatching, ordernum, sangerQueriesAll,
                               f"{outSangerDir}/CSVfiles/BestQueriesByGenome_{ordernum}_ImperfectSeq.csv")

    # Make working sets of sanger queries
    workingSetSangerQueriesMatching = buildWorkingSetFromAllQueries(matchingSeqs)
    workingSetSangerQueriesNonMatching = buildWorkingSetFromAllQueries(nonMatchingSeqs)

    # TODO: how to set this so it only happens if there is stuff in the object? don't want to make a file with no data
    if workingSetSangerQueriesMatching.workingSet:
        workingSetSangerQueriesMatching.__printToJson__(f"{outSangerDir}", ordernum)
    if workingSetSangerQueriesNonMatching.workingSet:
        workingSetSangerQueriesNonMatching.__printToJson__(f"{outSangerDir}", f"{ordernum}_NonMatching")

    # Pull from Working Sets of Queries
    originalData = buildFullWorkingSet(False)
    # Sort into 6 categories:
    #   dataMatching                    sequence and genome match
    #   badDataMatching                 sequence matches, genome does not
    #   queryDataInTwoSeqsMatching      sequence matches, but query is part of the "two seqs" set of flankseqs
    #   dataNonMatching                 sequences doesn't match, but genome matches
    #   badDataNonMatching              neither sequence nor genome matche
    #   queryDataInTwoSeqsNonMatching   sequence doesn't match, query is in "two seqs" subset
    dataMatching, badDataMatching, queryDataInTwoSeqsMatching = \
        sortSangerQueries(originalData, workingSetSangerQueriesMatching)
    dataNonMatching, badDataNonMatching, queryDataInTwoSeqsNonMatching = \
        sortSangerQueries(originalData, workingSetSangerQueriesNonMatching)

    writeToFile(dataMatching, ordernum, "SameGenome", "MatchingSeq", originalData, workingSetSangerQueriesMatching)
    writeToFile(badDataMatching, ordernum, "DiffGenome", "MatchingSeq", originalData, workingSetSangerQueriesMatching)
    writeToFile(dataNonMatching, ordernum, "SameGenome", "NonMatchingSeq", originalData,
                workingSetSangerQueriesNonMatching)
    writeToFile(badDataNonMatching, ordernum, "DiffGenome", "NonMatchingSeq", originalData,
                workingSetSangerQueriesNonMatching)

    if len(queryDataInTwoSeqsMatching) > 0:
        writeSangerComparisonFileTwoSeqs(f"{outSangerDir}/CSVfiles/TwoSeqsSet_{ordernum}_MatchingSeq.csv",
                                         queryDataInTwoSeqsMatching, originalData, workingSetSangerQueriesMatching)

    if len(queryDataInTwoSeqsNonMatching) > 0 :
        writeSangerComparisonFileTwoSeqs(f"{outSangerDir}/CSVfiles/TwoSeqsSet_{ordernum}_NonMatching.csv",
                                         queryDataInTwoSeqsNonMatching, originalData,
                                         workingSetSangerQueriesNonMatching)


if __name__ == '__main__':
    main()