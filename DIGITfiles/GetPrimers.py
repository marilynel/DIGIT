'''
GetPrimers.py
    1. parses the initial Primer3 output file for primer information
    2. saves primer data to the WorkingSet JSON document
    3. composes a Primer3 input file to be used for verification
    4. calls primer3_core with new Primer3 input file

Written by: Marilyn Leary 2023

Input:
flankseq            Name of working folder, should be the same as the filename given to the original fasta file of
                    flanking sequences.

Output:
WorkingSet_flankseq .json               Updated JSON document with parsed Primer3 data.
Primer3VerificationInput_flankseq.txt   Input file to pass to primer3_core; contains wildtype sequence information as
                                        well as the parsed Primer3 data.
Primer3VerificationOutput_flankseq.txt  Output file produced by Primer3.

'''

#import json
import sys
#import subprocess

from QueriesWorkingSet import QueriesWorkingSet
#from Query import Query
#from Sequences import Sequences
#from Primer3Object import Primer3Object
from Utils import *


def main():
    flankseq = sys.argv[1]
    outputDir = f"DIGIToutput/FlankingSequences/{flankseq}"
    queriesWorkingSet = QueriesWorkingSet()

    workingQueriesJsonFile = f"{outputDir}/QueryData/JSONfiles/WorkingSet_{flankseq}.json"

    queriesWorkingSet.__createQueryStructFromJson__(workingQueriesJsonFile)

    p3OutFile = f"{outputDir}/QueryData/Primer3/FindingPrimers/Primer3Output_{flankseq}.txt"

    readPrimer3Output(p3OutFile, "generic", queriesWorkingSet)

    queriesWorkingSet.__printToJson__(f"{outputDir}/QueryData", flankseq)

    verfInput = f"{outputDir}/QueryData/Primer3/VerifyingPrimers/Primer3VerificationInput_{flankseq}.txt"
    queriesWorkingSet.__writeP3InputFile__(verfInput, "check_primers")

    runPrimer3(verfInput, flankseq,
               f"{outputDir}/QueryData/Primer3/VerifyingPrimers/Primer3VerificationOutput_{flankseq}.txt")


if __name__ == "__main__":
    main()