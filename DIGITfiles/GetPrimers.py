from Query import Query
from Sequences import Sequences
from Primer3Object import Primer3Object
import json
import sys
import subprocess

queriesWorkingSet = {}


# TODO: add parsing to primer3object?

def getDataFromJSON(filename):
    print("running getDataFromJSON")
    '''
    Parse a JSON file and return a JSON object (key:value format)
    '''
    jsonFile = open(filename)
    dataFromJSON = json.load(jsonFile)
    jsonFile.close()

    return dataFromJSON


def createQueryStruct(filename):
    print("running createqueryStruct")
    dataFromJSON = getDataFromJSON(filename)
    for allele in dataFromJSON:
        if allele not in queriesWorkingSet:
            queriesWorkingSet[allele] = Query(None, None, None)
            queriesWorkingSet[allele].__QueryFromJSON__(dataFromJSON[allele])


def readPrimer3Output(filename):
    print("running readPrimer3Outpu")
    with open(filename, "r") as p3file:
        outputDict = {}  # create a dictionary for each output item
        for line in p3file:
            if line[0] != "=":  # go through the file until reaching a break mark ("=")
                # print("ok")
                lineItems = line.split("=")  # key:value pairs for OUTPUTLABEL:data
                outputDict[lineItems[0].strip()] = lineItems[1].strip()
            elif line[0] == "=":  # end of a section,time to process data
                # print("--")
                allele = outputDict["SEQUENCE_ID"][:-1]
                # print(outputDict)
                queriesWorkingSet[allele].__updateQueryPrimer3__(outputDict)
                outputDict.clear()
            else:
                print("ERROR")


# TODO: is this generalizable? similar function appears in another part of program
def queriesToJSON(filename):
    '''
    Converts dictionary to JSON object and writes it to a file.
    '''
    print("running queriesToJSON")
    jsonObject = {}

    for q in queriesWorkingSet:
        if queriesWorkingSet[q].query not in jsonObject:
            jsonObject[queriesWorkingSet[q].query] = []
        jsonObject[queriesWorkingSet[q].query] = dict(queriesWorkingSet[q])
    newJSONobject = json.dumps(jsonObject, indent=4)
    with open(filename, "w") as outfile:
        outfile.write(newJSONobject)

    return


def createPrimer3VerificationInput(flankseq):
    print("running createPrimer3VerificationINput")
    with open(
            "DIGIToutput/" + flankseq + "/WTVerification/Primer3VerificationInput_" + flankseq + ".txt",
            "w") as p3file:
        for allele in queriesWorkingSet:
            leftP3 = Primer3Object("check_primers", "left", queriesWorkingSet[allele])
            rightP3 = Primer3Object("check_primers", "right", queriesWorkingSet[allele])
            newInputStr = leftP3.inputStr + rightP3.inputStr
            p3file.write(newInputStr)
    return "DIGIToutput/" + flankseq + "/WTVerification/Primer3VerificationInput_" + flankseq + ".txt"


def runPrimer3(inputFile, flankseq):
    print("running runPrimer3")
    # now = int(time.time())
    with open(
            "DIGIToutput/" + flankseq + "/WTVerification/Primer3VerificationOutput_" + flankseq + ".txt",
            "w") as outfile:
        subprocess.run(
            [
                "primer3_core",
                inputFile
            ],
            stdout=outfile  # ,text=True
        )


def main():
    print("running main")
    flankseq = sys.argv[1]

    workingQueriesFile = "DIGIToutput/" + flankseq + "/DataSets/WorkingSet_" + flankseq + ".json"
    primer3OutputFile = "DIGIToutput/" + flankseq + "/Primer3Files/Output/Primer3Output_" + flankseq + ".txt"

    createQueryStruct(workingQueriesFile)
    readPrimer3Output(primer3OutputFile)
    queriesToJSON(
        "DIGIToutput/" + flankseq + "/DataSets/WorkingSet_" + flankseq + "_withPrimers.json")

    # make verification input???
    verfInput = createPrimer3VerificationInput(flankseq)
    print(verfInput)
    runPrimer3(verfInput, flankseq)


if __name__ == "__main__":
    main()