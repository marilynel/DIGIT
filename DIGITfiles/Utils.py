import json

from Query import Query
from Primer3Object import Primer3Object


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
        return "DIGIToutput/" + flankseq + "/WTVerification/Input/Primer3VerificationInput_" + \
               flankseq + ".txt"


def workingQueriesToJSON(filename, queriesWorkingSet):
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