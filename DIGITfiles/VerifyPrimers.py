'''
VerifyPrimers.py:
    1.  Parses Primer3 verification output and stores relevant data in the QueriesWorkingSet for that set.
    2.  Creates a fasta file of the potential primers and gives user the option to blast primers against genomes.
    3.  Sorts primers based on Primer3 verification and produces two files, one with successful primers and one with
        failed primers.
'''

import json
import os
import sys

from QueriesWorkingSet import QueriesWorkingSet
from Query import Query
from Utils import *


def main():
    flankseq = sys.argv[1]
    outputDir = f"DIGIToutput/FlankingSequences/{flankseq}"
    makeDirectories([f"{outputDir}/PrimerData", f"{outputDir}/BlastOutput/PrimerBlast"])

    workingQueriesJsonFile = f"{outputDir}/QueryData/JSONfiles/WorkingSet_{flankseq}.json"

    # Get Query Working Set from JSON file
    queriesWorkingSet = QueriesWorkingSet()
    queriesWorkingSet.__createQueryStructFromJson__(workingQueriesJsonFile)

    # Read P3 verification output & update primers
    readPrimer3Output(f"{outputDir}/QueryData/Primer3/VerifyingPrimers/Primer3VerificationOutput_{flankseq}.txt",
                      "check_primers", queriesWorkingSet)

    primerFastaFile = f"{outputDir}/PrimerData/{flankseq}_PrimerFastaFile.fasta"
    queriesWorkingSet.__makeFastaFromPrimers__(primerFastaFile)

    print("Would you like to Blast the primers against all three genomes in order to find their incidence rates and location? y/n")
    ans = input()

    if ans == "y":
        callPrimerBlastScript(primerFastaFile, f"{outputDir}/BlastOutput/PrimerBlast")

    queriesWithBadPrimers = queriesWorkingSet.__createBadQueryStruct__()

    queriesWorkingSet.__printToJson__(f"{outputDir}/QueryData", flankseq)

    queriesWorkingSet.__makePrimerDataFile__(f"{outputDir}/PrimerData/PrimerResults_{flankseq}.csv")

    if queriesWithBadPrimers.workingSet:
        queriesWithBadPrimers.__makePrimerDataFile__(f"{outputDir}/PrimerData/FailedPrimers_{flankseq}.csv")


if __name__ == "__main__":
    main()
