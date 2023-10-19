'''
BlastSangerOutput.py creates a fasta file out of the .seq files returned from sanger sequencing and blasts those
sequences against the three reference genomes, A188, B73, and W22.

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
from SangerSeqWorkingSet import *


def main():
    orderDir = sys.argv[1]

    order = orderDir.split("/")[-1]

    sangerSeqSet = SangerSeqWorkingSet(order)
    sangerSeqSet.__appendSangerSeqObjsToWorkingSet__(orderDir)
    fastaFile = sangerSeqSet.__prepFasta__()

    newDirs = [f"DIGIToutput/SangerSequences",
               f"DIGIToutput/SangerSequences/{order}",
               f"DIGIToutput/SangerSequences/{order}/BlastOutput",
               f"DIGIToutput/SangerSequences/{order}/JSONfiles"]
    makeDirectories(newDirs)

    callBlastScript(fastaFile, f"DIGIToutput/SangerSequences/{order}/BlastOutput", order)
    sangerSeqSet.__printToJson__()



if __name__ == '__main__':
    main()