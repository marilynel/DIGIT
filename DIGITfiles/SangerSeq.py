import re
from difflib import SequenceMatcher

import json
# import os

# import subprocess

from Query import *


class SangerSeq:
    def __init__(self, filename):
        self.queryName = ""
        self.metadata = ""
        if filename:
            self.queryName = re.findall(r"R\d{1,4}[A-Z]*\d{1,4}", filename)[0]
            self.metadata = filename.split(".")[0]
        self.orderNum = ""
        self.primerMatch = None
        self.rawSequence = ""
        self.sequence = ""
        self.cleanSequence = ""
        self.dsgfpMismatch = False
        self.sangerQueryObj = None
        self.originalQueryObj = None
        self.dsgfpEnd = "TTTTACCGACCGTTCCCGACCGTTTTCATCCCTA"
        self.possibleQueryMatches = []
        self.foundInA188 = False
        self.foundInB73 = False
        self.foundInW22 = False

    def __dataLine__(self):
        if self.originalQueryObj:
            return f"{self.queryName}\t{self.sangerQueryObj.genome}\t{self.originalQueryObj.genome}\t{self.__compareSangerAndOriginal__()}"
        else:
            return f"{self.queryName}\t{self.sangerQueryObj.genome}\tnone\t{self.__compareSangerAndOriginal__()}"

    def __compareSangerAndOriginal__(self):
        # return (f"{type(self.sangerQueryObj)}\t{type(self.originalQueryObj)}")
        if not self.originalQueryObj:
            return False
        if (self.sangerQueryObj.genome != self.originalQueryObj.genome) or (
                self.sangerQueryObj.chromosome != self.originalQueryObj.chromosome) or (
                self.sangerQueryObj.sStart != self.originalQueryObj.sStart):
            return False
        return True

    def __findBestHit__(self):
        for i in range(len(self.possibleQueryMatches)):
            if self.possibleQueryMatches[i].bestHitForAllele:
                self.sangerQueryObj = self.possibleQueryMatches[i]
                break

    def __readSeqFile__(self, filename, order):
        # 28329_C02_R179H03_5-Ds-3_C02_006.seq
        self.orderNum = order
        # self.query =  re.findall (r"R\d{1,4}[A-Z]*\d{1,4}", filename)
        if filename.find("DsGFP_3UTR") != -1:
            self.primerMatch = "DsGFP_3UTR"
        if filename.find("5-Ds-3") != -1:
            self.primerMatch = "5-Ds-3"
        with open(filename, "r") as seqFile:
            firstLine = seqFile.readline()
            if firstLine[0] == ">":
                self.metadata = firstLine.strip()
            else:
                self.rawSequence += firstLine.strip()
            for line in seqFile:
                self.rawSequence += line.strip()

        self.__removeDsgfp__()
        self.__cleanConsecutiveNs__()

    def __printToFasta__(self):
        return f">{self.queryName}\n{self.cleanSequence}\n"

    def __removeDsgfp__(self):
        dsIdx = 0
        # fastaFormatStr += sequence[:endIdx] + "\n"
        if self.rawSequence.find(self.dsgfpEnd) != -1:
            dsIdx = self.rawSequence.find(self.dsgfpEnd)
            #       raw seq     index
            # return sequence, sequence.find(dsgfpEnd), False
        else:
            self.dsgfpMismatch = True
            try:
                seqMatchObj = SequenceMatcher(lambda x: x == "ACGT", self.rawSequence,
                                              self.dsgfpEnd)
                matchObj = seqMatchObj.find_longest_match(0, len(self.rawSequence), 0,
                                                          len(self.dsgfpEnd))
                dsIdx = matchObj.a
            except:
                dsIdx = None
        self.sequence = self.rawSequence[dsIdx + len(self.dsgfpEnd):-22]

    def __cleanConsecutiveNs__(self):
        previousBasePairs = self.sequence[0]
        for i in range(0, len(self.sequence)):
            if previousBasePairs == "N" and self.sequence[i] == "N":
                self.cleanSequence = self.sequence[:i - 1]
                return
            else:
                previousBasePairs = self.sequence[i]
        self.cleanSequence = self.sequence

    def __addToPossibleQueries__(self, queryObj):
        if queryObj.query == self.querName:
            self.possibleQueryMatches.append(queryObj)

    def __SangerSeqFromJSON__(self, jsonObject):
        self.queryName = jsonObject["queryName"]
        self.metadata = jsonObject["metadata"]
        self.orderNum = jsonObject["orderNum"]
        self.primerMatch = jsonObject["primerMatch"]
        self.rawSequence = jsonObject["rawSequence"]
        self.sequence = jsonObject["sequence"]
        self.cleanSequence = jsonObject["cleanSequence"]
        self.dsgfpMismatch = jsonObject["dsgfpMismatch"]
        self.sangerQueryObj = jsonObject["sangerQueryObj"]
        self.originalQueryObj = jsonObject["originalQueryObj"]
        self.dsgfpEnd = jsonObject["dsgfpEnd"]
        self.possibleQueryMatches = jsonObject["possibleQueryMatches"]

    def __toJSON__(self):
        json0bj = json.dumps(self, default=lambda o: o.__dict__, sort_keys=True, indent=4)
        return json0bj

    def __iter__(self):
        iters = dict((x, y) for x, y in SangerSeq.__dict__.items() if x[:2] != '__')
        iters.update(self.__dict__)
        for x, y in iters.items():
            yield x, y


########################################################################################################################


class SangerSeqSet:

    def __init__(self, order):
        self.order = order
        self.sangerWorkingSet = {}

    def __setBestQueryToSanger__(self):
        for s in self.sangerWorkingSet:
            self.sangerWorkingSet[s].__findBestHit__()

    def __appendSangerSeqObjsToWorkingSet__(self, orderDir):
        order = orderDir.split("/")[-1][-5:]
        # print(order)
        sangerSeqs = []
        for subdir, dirs, files in os.walk(orderDir):
            for oneFile in files:
                filename = os.path.join(subdir, oneFile)
                if filename.find(".seq") != -1:
                    newSangerObj = SangerSeq(filename)
                    newSangerObj.__readSeqFile__(filename, order)
                    # TODO: conditional: what if there are two of these queireis?
                    self.sangerWorkingSet[newSangerObj.queryName] = newSangerObj

    def __prepFasta__(self):
        filename = f"PutSangerOutputFilesHere/{self.order}/{self.order}.fasta"
        with open(filename, "w") as outfile:
            for q in self.sangerWorkingSet:
                outfile.write(self.sangerWorkingSet[q].__printToFasta__())
        return filename

    def __printToJson__(self):
        jsonfile = f"DIGIToutput/SangerSequences/{self.order}/JSONfiles/SangerSet_{self.order}.json"
        jsonObject = {}
        for s in self.sangerWorkingSet:
            if self.sangerWorkingSet[s].queryName not in jsonObject:
                jsonObject[self.sangerWorkingSet[s].queryName] = []
            jsonObject[self.sangerWorkingSet[s].queryName] = dict(self.sangerWorkingSet[s])
        newJSONobject = json.dumps(jsonObject, indent=4)

        with open(jsonfile, "w") as outfile:
            outfile.write(newJSONobject)

    def __readFromJson__(self):
        sangerJsonFile = f"DIGIToutput/SangerSequences/{self.order}/JSONfiles/SangerSet_{self.order}.json"
        jsonFile = open(sangerJsonFile)
        dataFromJSON = json.load(jsonFile)
        jsonFile.close()

        # Read old sangerseq from json
        for sangerSeq in dataFromJSON:
            newSangerSeq = SangerSeq(None)
            newSangerSeq.__SangerSeqFromJSON__(dataFromJSON[sangerSeq])
            if newSangerSeq.queryName not in self.sangerWorkingSet:
                self.sangerWorkingSet[newSangerSeq.queryName] = None
            self.sangerWorkingSet[newSangerSeq.queryName] = newSangerSeq
        ## end todo

    def __addQueriesToSangerObjs__(self, blastOutputObj):
        # add the Query objects taht represent blast output hits to the correct item in sangerObjDict

        # for genome in sangerQueriesAll.hits:
        for genome in blastOutputObj.hits:
            #    for query in sangerQueriesAll.hits[genome]:
            for allele in blastOutputObj.hits[genome]:
                #        if query in sangerObjDict:
                if allele in self.sangerWorkingSet:
                    #            sangerObjDict[query].possibleQueryMatches += sangerQueriesAll.hits[genome][query]
                    self.sangerWorkingSet[allele].possibleQueryMatches += \
                    blastOutputObj.hits[genome][allele]
        # writeToBestQueriesFile(listQ, ordernum, sangerQueriesAll,
        # f"DIGIToutput/SangerSequences/{ordernum}/CSVfiles/BestSangersByGenome_{ordernum}.csv")

    def __setOriginalQuery__(self, originalData):
        # queryData, badQueryData, queryDataInTwoSeqs = [], [], []

        for query in self.sangerWorkingSet:
            # queryOG = query.split("_")[0]
            if query in originalData.workingSet:
                self.sangerWorkingSet[query].originalQueryObj = originalData.workingSet[query]

#            # If the query does not appear in originalData, check if it appears as part of the twoSeq datasets
#           else:
#              queryA = queryOG + "a"
#             queryB = queryOG + "b"
#
#           if (queryA in originalData.workingSet) or (queryB in originalData.workingSet):
#              queryDataInTwoSeqs.append(query)
#
#        # If it still doesn't appear, its not there
#       else:
#          print(f"{query} does not exist in original dataset")
#
# return queryData, badQueryData, queryDataInTwoSeqs
