'''
SequencesFromBlast.py:
    1.  Parses flanking sequence blast output files (.tab format only).
    2.  Aggregates data for all hits (occurances of the blasted sequence anywhere in the
    reference genome).
    3.  Identifies the "best" hit for each query according to the following criteria:
        -   First, the hit with the best bit score over all is chosen
        -   If bit scores are equivalent across genomes, the best score is chosen based on
        reference genome where B73 is
            preferred, then W22, then A188
        -   If there are multiple hits within the same preferred genome with the same bit score,
        one hit is chosen at
            random to be the best hit
    4.  Produces working set of queries (QueriesWorkingSet object) based on #3, where each allele
    only appears one time.
    5.  Composes and submits Primer3 input files.
'''

import sys

from BlastOutput import *
from QueriesWorkingSet import QueriesWorkingSet
from Query import Query
from FilterFasta import filterFasta
from Sequences import Sequences
from Primer3Object import Primer3Object
from PrimerIDs import PrimerIDs
from Utils import *


def main():
    flankseq = sys.argv[1]
    b73only = int(sys.argv[2])
    outputDir = "DIGIToutput/FlankingSequences"

    queriesWorkingSet = QueriesWorkingSet()

    newFlankseq = flankseq
    genomeSpec = ""
    x = False

    if b73only == 2:
        genomeSpec = "B73"
        newFlankseq += "B73Only"
        x = True

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

    # Read Blast output data into allQueries BlastOutput Object
    allQueries = BlastOutput(newFlankseq, "FlankingSequences")
    allQueries.__fillBlastOutputObject__(x, False)

    # Identify "best" hits for each allele and use to create a QueriesWorkingSet Object,
    # queriesWorkingSet
    queriesWorkingSet.__buildWorkingSetFromBlastOutputObj__(allQueries)

    queriesWorkingSet.__getOriginalFlankingSequences__(flankseq)

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
