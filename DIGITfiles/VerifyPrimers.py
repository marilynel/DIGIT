from Query import Query
from Utils import *
import json
import os
import sys

queriesWorkingSet = {}
'''
def composeOutput(filename):
    primerList = []
    for query in queriesWorkingSet:
        leftPrimerMatch, rightPrimerMatch = "", ""
        if queriesWorkingSet[query].sideMatchGFP3UTR == "right" and queriesWorkingSet[query].sideMatch3DsgG == "left":
            leftPrimerMatch = "GFP3UTR"
            rightPrimerMatch = "3DsgG"
        elif queriesWorkingSet[query].sideMatchGFP3UTR == "left" and queriesWorkingSet[query].sideMatch3DsgG == "right":
            leftPrimerMatch = "3DsgG"
            rightPrimerMatch = "GFP3UTR"
        else:
            leftPrimerMatch = "FAIL"
            rightPrimerMatch = "FAIL"
        primerList.append(
            f"{queriesWorkingSet[query].primerNameLeft},{queriesWorkingSet[query].query}," + 
            f"{queriesWorkingSet[query].genome},{queriesWorkingSet[query].primerSequenceLeft}," + 
            f"{leftPrimerMatch}," +
            f"{queriesWorkingSet[query].primerLeftProductSize},{queriesWorkingSet[query].primerPairProductSize}," +
            f"{queriesWorkingSet[query].tmLeft},{queriesWorkingSet[query].primerPenaltyLeft}," +
            f"{queriesWorkingSet[query].primerPairPenalty},{queriesWorkingSet[query].bitScore}," +
            f"{queriesWorkingSet[query].numHits},{queriesWorkingSet[query].qStartStatus}\n"
        )
        primerList.append(
            f"{queriesWorkingSet[query].primerNameRight},{queriesWorkingSet[query].query}," + 
            f"{queriesWorkingSet[query].genome},{queriesWorkingSet[query].primerSequenceRight}," + 
            f"{rightPrimerMatch}," +
            f"{queriesWorkingSet[query].primerRightProductSize},{queriesWorkingSet[query].primerPairProductSize}," +
            f"{queriesWorkingSet[query].tmRight},{queriesWorkingSet[query].primerPenaltyRight}," +
            f"{queriesWorkingSet[query].primerPairPenalty},{queriesWorkingSet[query].bitScore}," +
            f"{queriesWorkingSet[query].numHits},{queriesWorkingSet[query].qStartStatus}\n"
        )

    primerList.sort()
    with open(filename, "w") as primerFile:
        primerFile.write(
            "PrimerID,Allele,ReferenceGenome,PrimerSequence,MatchingDsGFPPrimer,ExpectedProductSizeWithDsPrimer,Expec" +
            "tedWTProductSize,TM,PrimerPenaltyWithDSPrimer,PrimerPenaltyWT,BlastBitScore,TotalBlastHitsForRefGenome,Q" +
            "ueryStartStatus\n"
        )
        for item in primerList:
            primerFile.write(item)

'''


def createBadPrimerQueryStruct():
    queriesWithBadPrimers = {}
    # with open("testBadPrimerOutput.csv", "w+") as outfile:
    #    outfile.write("Query,RightPrimer,LeftPrimer\n")
    for q in queriesWorkingSet:
        if queriesWorkingSet[q].primerSequenceRight == "FAIL" or queriesWorkingSet[
            q].primerSequenceLeft == "FAIL":
            queriesWithBadPrimers[q] = queriesWorkingSet[q]
            # outfile.write(f"{q},{queriesWorkingSet[q].primerSequenceRight},{queriesWorkingSet[
            # q].primerSequenceLeft}\n")
    return queriesWithBadPrimers


def main():
    # get all the primer data
    flankseq = sys.argv[1]

    # Convert JSON file to query struct, function in Utils
    createQueryStruct("DIGIToutput/" + flankseq + "/DataSets/WorkingSet_" + flankseq + ".json",
                      queriesWorkingSet)

    # read P3 verification output & update primers
    readPrimer3Output(
        "DIGIToutput/" + flankseq + "/WTVerification/Output/Primer3VerificationOutput_" + flankseq + ".txt",
        "check_primers",
        queriesWorkingSet
    )

    ### START NEW ###

    '''
    TODO:
    -   make a fasta of the primers to blast

    '''
    primerDataPath = "DIGIToutput/" + flankseq + "/PrimerData"
    if not os.path.exists(primerDataPath):
        os.makedirs(primerDataPath)

    with open(primerDataPath + "/testPrimerFastaFile.fasta", "w+") as fastafile:
        for q in queriesWorkingSet:
            if queriesWorkingSet[q].primerSequenceLeft != "FAIL":
                fastafile.write(
                    f">{queriesWorkingSet[q].primerNameLeft}_ogRefGen_{queriesWorkingSet[q].genome}\n")
                fastafile.write(f"{queriesWorkingSet[q].primerSequenceLeft}\n")
            if queriesWorkingSet[q].primerSequenceLeft != "FAIL":
                fastafile.write(
                    f">{queriesWorkingSet[q].primerNameLeft}_ogRefGen_{queriesWorkingSet[q].genome}\n")
                fastafile.write(f"{queriesWorkingSet[q].primerSequenceLeft}\n")

    callBlastScript(
        primerDataPath + "/testPrimerFastaFile.fasta",
        primerDataPath,
        "testPrimerFastaFile"
    )

    ### END NEW ###

    queriesWithBadPrimers = createBadPrimerQueryStruct()

    queriesToJSON("DIGIToutput/" + flankseq + "/DataSets/WorkingSet_" + flankseq + ".json",
                  queriesWorkingSet)

    composePrimerListOutput(
        "DIGIToutput/" + flankseq + "/DataSets/PrimerResults_" + flankseq + ".csv",
        queriesWorkingSet)
    composePrimerListOutput(
        "DIGIToutput/" + flankseq + "/DataSets/FailedPrimers_" + flankseq + ".csv",
        queriesWithBadPrimers)


if __name__ == "__main__":
    main()


