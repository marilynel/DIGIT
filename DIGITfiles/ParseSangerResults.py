'''
Parses results of Blasting Sanger sequences against the selected maize genomes. Identifies the "best" hit for each
Sanger query (according to the same prioritization scheme as the flanking sequences, see QueriesWorkingSet.py) and
compares that result to the original flanking sequence data, if available.

Written by: Marilyn Leary 2023
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
        f"DIGIToutput/SangerSequences/{ordernum}/QueryData",
        f"DIGIToutput/SangerSequences/{ordernum}/QueryData/JSONfiles",
        f"DIGIToutput/SangerSequences/{ordernum}/QueryData/GFFfiles",
        f"DIGIToutput/SangerSequences/{ordernum}/QueryData/CSVfiles"
    ]

    makeDirectories(necessaryDirs)
    originalData = QueriesWorkingSet()
    originalData.__buildCompleteWorkingSet__("")

    queryList = originalData.workingSet.keys()

    allSangers = BlastOutput(ordernum, "SangerSequences")

    # TODO: enable setting flag to true for b73only
    allSangers.__fillBlastOutputObject__(False, False)

    sangerSeqs = SangerSeqSet(ordernum)
    sangerSeqs.__readFromJson__()

    sangerSeqs.__addQueriesToSangerObjs__(allSangers)
    sangerSeqs.__setBestQueryToSanger__()
    sangerSeqs.__setOriginalQuery__(originalData)

    sangerSeqs.__printDataFile__()
    sangerSeqs.__printToGff__(f"DIGIToutput/SangerSequences/{ordernum}/QueryData/GFFfiles/SangerSet_{ordernum}.gff")


if __name__ == '__main__':
    main()