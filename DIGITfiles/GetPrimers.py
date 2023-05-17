from Query import Query
from Sequences import Sequences
from Primer3Object import Primer3Object
from Utils import *
import json
import sys
import subprocess

queriesWorkingSet = {}


# def runPrimer3(inputFile, flankseq):
#    '''
#    Submit wildtype verification (as inputFile) to primer3_core and redirect output to a similarly stored file.
#    '''
#    with open(
#        "DIGIToutput/" + flankseq + "/WTVerification/Output/Primer3VerificationOutput_" + flankseq + ".txt", "w") as outfile:
#        subprocess.run(
#            [
#                "primer3_core",
#                inputFile
##            ],
#            stdout=outfile
#        )


def main():
    flankseq = sys.argv[1]

    workingQueriesFile = "DIGIToutput/" + flankseq + "/DataSets/WorkingSet_" + flankseq + ".json"
    primer3OutputFile = "DIGIToutput/" + flankseq + "/Primer3Files/Output/Primer3Output_" + flankseq + ".txt"

    createQueryStruct(workingQueriesFile, queriesWorkingSet)
    readPrimer3Output(primer3OutputFile, "generic", queriesWorkingSet)
    # originally created new file, now just rewriting old file
    # queriesToJSON("DIGIToutput/" + flankseq + "/DataSets/WorkingSet_" + flankseq +
    # "_withPrimers.json", queriesWorkingSet)
    queriesToJSON(workingQueriesFile, queriesWorkingSet)

    verfInput = writeP3InputFile(
        "DIGIToutput/" + flankseq + "/WTVerification/Input/Primer3VerificationInput_" + flankseq + ".txt",
        queriesWorkingSet,
        flankseq,
        "check_primers"
    )

    runPrimer3(
        verfInput,
        flankseq,
        "DIGIToutput/" + flankseq + "/WTVerification/Output/Primer3VerificationOutput_" + flankseq + ".txt"
    )


if __name__ == "__main__":
    main()