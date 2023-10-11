'''
GetPrimers.py
    1. parses the initial Primer3 output file for primer information
    2. saves primer data to the WorkingSet JSON document
    3. composes a Primer3 input file to be used for verification
    4. calls primer3_core with new Primer3 input file for verification


Written by: Marilyn Leary 2023

    Input:
    flankseq (sys.argv[1])  String representing the Flanking Sequence Set of interest. Should match with a subfolder of
                            DIGIToutput/FlankingSequences

    Output:
    Working set files in CSV, JSON, and GFF formats:
        DIGIToutput/FlankingSequences/<flankseq>/QueryData/JSONfiles/WorkingSet_<flankseq>.json
        DIGIToutput/FlankingSequences/<flankseq>/QueryData/CSVfiles/WorkingSet_<flankseq>.csv
        DIGIToutput/FlankingSequences/<flankseq>/QueryData/GFFfiles/WorkingSet_<flankseq>.gff

    Primer3 Verification input and output:
        DIGIToutput/FlankingSequences/<flankseq>/QueryData/Primer3/VerifyingPrimers/Primer3VerificationInput_flankseq.txt
        DIGIToutput/FlankingSequences/<flankseq>/QueryData/Primer3/VerifyingPrimers/Primer3VerificationOutput_flankseq.txt

'''

import sys

from QueriesWorkingSet import QueriesWorkingSet
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