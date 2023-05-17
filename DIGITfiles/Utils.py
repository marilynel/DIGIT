import json
import subprocess
import time
import sys
import os

# sys.path.append(os.path.dirname(os.path.realpath(__file__)))

from Query import Query
from Primer3Object import Primer3Object


def callBlastScript(pathToFastaInput, pathToDesitinationDir, fastaContents):
    now = time.time()
    subprocess.run(
        [
            "sh",
            f"./DIGITfiles/RunBlastInitial.sh",
            pathToFastaInput,
            pathToDesitinationDir,
            fastaContents
        ]
    )
    print(
        f"SGE may take a while to run. To check the status of your jobs, enter 'qstat' in the command line. Runni" +
        f"ng time may vary wildly.\n"
    )


def readPrimer3Output(filename, task, queriesWorkingSet):
    with open(filename, "r") as p3file:
        # outputDict contains key:value pairs for each line of a single P3 response. Remember there may be  multiple
        # responses for each P3 output file.
        outputDict = {}
        for line in p3file:
            # A single = indicates the end of a P3 response. Next line will be a new response.
            if line[0] != "=":
                lineItems = line.split("=")  # key:value pairs for OUTPUTLABEL:data
                outputDict[lineItems[0].strip()] = lineItems[1].strip()
            elif line[0] == "=":
                # process data
                allele = outputDict["SEQUENCE_ID"].split("_")[1]

                # Initialize a P3 object using the output values
                newP3Obj = Primer3Object(task)
                newP3Obj.__initValsFromP3Output__(outputDict)

                # Pass P3 object to Query method to update Query object
                queriesWorkingSet[allele].__updateQueryWithP3Output__(newP3Obj)

                outputDict.clear()
            else:
                print("ERROR")


def createQueryStruct(filename, queriesWorkingSet):
    dataFromJSON = getDataFromJSON(filename)
    for allele in dataFromJSON:
        if allele not in queriesWorkingSet:
            queriesWorkingSet[allele] = Query(None, None, None)
            queriesWorkingSet[allele].__QueryFromJSON__(dataFromJSON[allele])


def getDataFromJSON(filename):
    '''
    Parse a JSON file and return a JSON object (key:value format)
    '''
    jsonFile = open(filename)
    dataFromJSON = json.load(jsonFile)
    jsonFile.close()

    return dataFromJSON


def writeP3InputFile(inputFilepath, queriesWorkingSet, flankseq, task):
    with open(inputFilepath, "w") as p3file:
        for allele in queriesWorkingSet:
            p3obj = Primer3Object(task)
            p3obj.__initValsFromQueryObj__(queriesWorkingSet[allele])

            # TODO: refactor this so I'm handling it inside p3obj
            if task == "generic":
                p3obj.__buildInputStrGeneric__()
            else:
                p3obj.__buildInputStrValidate__()
            p3file.write(p3obj.inputStr)

    if task == "generic":
        return "DIGIToutput/" + flankseq + "/Primer3Files/Input/Primer3Input_" + flankseq + ".txt"
    else:
        return "DIGIToutput/" + flankseq + "/WTVerification/Input/Primer3VerificationInput_" + flankseq + ".txt"


def queriesToJSON(filename, queriesWorkingSet):
    '''
    Converts dictionary to JSON object and writes it to a file.
    '''
    jsonObject = {}
    for q in queriesWorkingSet:
        if queriesWorkingSet[q].query not in jsonObject:
            jsonObject[queriesWorkingSet[q].query] = []
        jsonObject[queriesWorkingSet[q].query] = dict(queriesWorkingSet[q])
    newJSONobject = json.dumps(jsonObject, indent=4)

    with open(filename, "w") as outfile:
        outfile.write(newJSONobject)

    return


def checkIfPrimerHasID():
    # mr017b_R153F05
    # mr018a_R153F05
    listPreviousPrimers = []
    with open("DIGITfiles/compListPreviousPrimers", "r") as prevFile:
        for line in prevFile:
            listPreviousPrimers.append(line.strip())
    return listPreviousPrimers


def runPrimer3(inputFile, flankseq, outputFilePath):
    with open(outputFilePath, "w") as outfile:
        subprocess.run(
            [
                "primer3_core",
                inputFile
            ],
            stdout=outfile, text=True
        )


def composePrimerListOutput(filename, queriesWorkingSet):
    primerList = []
    for query in queriesWorkingSet:
        leftPrimerMatch, rightPrimerMatch = "", ""
        if queriesWorkingSet[query].sideMatchGFP3UTR == "right" and queriesWorkingSet[
            query].sideMatch3DsgG == "left":
            leftPrimerMatch = "GFP3UTR"
            rightPrimerMatch = "3DsgG"
        elif queriesWorkingSet[query].sideMatchGFP3UTR == "left" and queriesWorkingSet[
            query].sideMatch3DsgG == "right":
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

