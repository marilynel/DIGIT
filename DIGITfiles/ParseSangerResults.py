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
from BlastOutput import *


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

    # Find output of blasting sanger sequences and store in allSangers;
    # # TODO: note if any have 0 hits in genome
    allSangers = BlastOutput(ordernum, "SangerSequences")
    allSangers.__parseBlastOutputFiles__("")
    allSangers.__findBestForGenome__()
    allSangers.__pickGenome__()
    allSangers.__allBlastOutputDataToJson__()


    ### read in sangerjson
    sangerSeqs = SangerSeqSet(ordernum)
    sangerSeqs.__readFromJson__()

    sangerSeqs.__addQueriesToSangerObjs__(allSangers)
    sangerSeqs.__setBestQueryToSanger__()

    #noteNoHitSangers(allSangers)


### HERE
    # Split allSangers into two dicts with the same internal structure as allSangers
    #matchingSeqs, nonMatchingSeqs = splitallSangersIntoMatchingAndNonMatchingSets(allSangers)

    # Compare and find the "best" hits for each sanger sequence query; save results to files
    #sangerListMatching = setBestForGenome(matchingSeqs)
    #sangerListNonMatching = setBestForGenome(nonMatchingSeqs)

    #if len(sangerListMatching) > 0:
    #    writeToBestQueriesFile(sangerListMatching, ordernum, allSangers,
    #                           f"{outSangerDir}/CSVfiles/BestQueriesByGenome_{ordernum}.csv")

    #if len(sangerListNonMatching) > 0:
    #    writeToBestQueriesFile(sangerListNonMatching, ordernum, allSangers,
    #                           f"{outSangerDir}/CSVfiles/BestQueriesByGenome_{ordernum}_ImperfectSeq.csv")

    # Make working sets of sanger queries and save to files
    #workingSetSangerQueriesMatching = buildWorkingSetFromAllQueries(matchingSeqs)
    #workingSetSangerQueriesNonMatching = buildWorkingSetFromAllQueries(nonMatchingSeqs)


    # Pull from Working Sets of Queries
    originalData = buildFullWorkingSet(False)
    sangerSeqs.__setOriginalQuery__(originalData)


    print("Query\tOrder\tDsGFPmismatch\tNumBlastHits\tQuery\tGenome\tNumHits")
    for s in sangerSeqs.sangerWorkingSet:
        #"{self.queryName}\t{self.orderNum}\t{self.dsgfpMismatch}\t{len(self.possibleQueryMatches)}"
        print(sangerSeqs.sangerWorkingSet[s].__dataLine__())




    #TODO: write this stuff to file!!!!

    exit()
    # Sort into 6 categories:
    #   dataMatching                    sequence and genome match
    #   badDataMatching                 sequence matches, genome does not
    #   queryDataInTwoSeqsMatching      sequence matches, but query is part of the "two seqs" set of flankseqs
    #   dataNonMatching                 sequences doesn't match, but genome matches
    #   badDataNonMatching              neither sequence nor genome matche
    #   queryDataInTwoSeqsNonMatching   sequence doesn't match, query is in "two seqs" subset

    # sortSangerQueries() returns 3 lists of Query objects
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