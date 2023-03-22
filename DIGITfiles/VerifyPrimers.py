from DIGITfiles.Query import Query
import json

queriesWorkingSet = {}


def getDataFromJSON(filename):
    '''
    Parse a JSON file and return a JSON object (key:value format)
    '''
    # TODO: can i pulls this out and put it elsewhere? I keep reusing it
    jsonFile = open(filename)
    dataFromJSON = json.load(jsonFile)
    jsonFile.close()

    return dataFromJSON


def createQueryStruct(filename):
    dataFromJSON = getDataFromJSON(filename)
    for allele in dataFromJSON:
        if allele not in queriesWorkingSet:
            queriesWorkingSet[allele] = Query(None, None, None)
            queriesWorkingSet[allele].__QueryFromJSON__(dataFromJSON[allele])


def readPrimer3Output(filename):
    with open(filename, "r") as p3file:
        outputDict = {}  # create a dictionary for each output item
        for line in p3file:
            if line[0] != "=":  # go through the file until reaching a break mark ("=")
                lineItems = line.split("=")  # key:value pairs for OUTPUTLABEL:data
                outputDict[lineItems[0].strip()] = lineItems[1].strip()
            elif line[0] == "=":  # end of a section,time to process data
                allele = outputDict["SEQUENCE_ID"][:-1]

                queriesWorkingSet[allele].__updateQueryPrimer3__(outputDict)
                outputDict.clear()
                # print(output_dict)
            else:
                print("ERROR")


# TODO: is this generalizable? similar function appears in another part of program
def queriesToJSON(filename):
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


def main():
    # get all the primer data
    createQueryStruct("WorkingQueriesWithPrimers03.json")
    # read_primer3_output("Primer3Output08.txt")

    # make an input filefor primer3
    with open("VerifyPrimer3Input.txt", "w") as p3file:
        for allele in queriesWorkingSet:
            p3file.write(queriesWorkingSet[allele].__createPrimer3Input__())

    queriesToJSON("")


if __name__ == "__main__":
    main()