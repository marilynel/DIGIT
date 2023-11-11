'''
BlastSangerOutput.py:
    1.  Creates a fasta file out of the .seq files returned from sanger sequencing.
    2.  Blasts the fasta file against the three reference genomes, A188, B73, and W22.
    3.  Builds preliminary Sanger Working SetD

Written by: Marilyn Leary 2023

Input:
    orderDir (sys.argv[1])  Path to directory containing sanger sequencing .seq files. Should be a subdirectory within
                            PutSangerOutputFilesHere/.

Output:
    Three Blast output files with path and naming convention:
        DIGIToutput/SangerSequences/<orderNumber>/BlastOutput/<referenceGenome_vs_orderNumber>.tab
'''

import sys

from Utils import *
from SangerSeq import *


def main():
    orderDir = sys.argv[1]
    order = orderDir.split("/")[-1]

    sangerSeqs = SangerSeqSet(order)
    sangerSeqs.__appendSangerSeqObjsToWorkingSet__(orderDir)

    originalData = QueriesWorkingSet()
    originalData.__buildCompleteWorkingSet__("")
    sangerSeqs.__setOriginalQuery__(originalData)

    fastaFile = sangerSeqs.__prepFasta__()

    newDirs = [f"DIGIToutput/SangerSequences",
               f"DIGIToutput/SangerSequences/{order}",
               f"DIGIToutput/SangerSequences/{order}/BlastOutput",
               f"DIGIToutput/SangerSequences/{order}/QueryData",
               f"DIGIToutput/SangerSequences/{order}/QueryData/JSONfiles"]

    makeDirectories(newDirs)
    sangerSeqs.__printToJson__()

    callBlastScript(fastaFile, f"DIGIToutput/SangerSequences/{order}/BlastOutput", order)



if __name__ == '__main__':
    main()