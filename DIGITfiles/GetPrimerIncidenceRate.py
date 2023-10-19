'''
GetPrimerIncidenceRate.py creates a csv file out of number of times each primer occurs in each of the genomes, A188,
B73, and W22.

Written by: Marilyn Leary 2023

    Input:
    flankseq (sys.argv[1])  String representing the Flanking Sequence Set of interest. Should match with a subfolder of
                            DIGIToutput/FlankingSequences

    Output:
    Three Blast output files with path and naming convention:
        DIGIToutput/FlankingSequences/<flankseq>/PrimerData/<flankseq>_PrimerIncidenceRate.csv
'''

import sys

from Utils import *
from Query import *
from Primer import *


def main():
    allPrimers = {}
    flankseq = sys.argv[1]
    genomeSpec = "PrimerBlast.tab"
    dirname = f"DIGIToutput/FlankingSequences/{flankseq}/BlastOutput/PrimerBlast"

    findBlastOutputFiles(dirname, allPrimers, genomeSpec)
    # Create dict of Primer objects; key=query (R num), val=list of Primer objects
    primerDict = makePrimerDict(allPrimers)

    workingQueriesJsonFile = f"DIGIToutput/FlankingSequences/{flankseq}/QueryData/JSONfiles/WorkingSet_{flankseq}.json"
    queriesWorkingSet = QueriesWorkingSet()
    queriesWorkingSet.__createQueryStructFromJson__(workingQueriesJsonFile)

    # Check if primers are in the right spot
    setPrimerFoundNearFlankSeq(primerDict, queriesWorkingSet)
    writePrimerIncidenceDataToFile(primerDict, flankseq)


if __name__ == '__main__':
    main()