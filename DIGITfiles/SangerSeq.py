import re
from difflib import SequenceMatcher
import json
import os

from Query import *
from QueriesWorkingSet import *
from GFF import *


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

        self.dsgfpMatches = True
        self.dsgfpPresent = True
        self.blastable = True
        self.foundOriginalQuery = True
        self.bestSangerMatchesBestOriginal = True
        self.sangerSequenceFoundInGenome = True

        self.sangerQueryObj = None
        self.originalQueryObj = None
        self.dsgfpEnd = "TTTTACCGACCGTTCCCGACCGTTTTCATCCCTA"
        self.possibleQueryMatches = []

    def __setOriginalQueryObj__(self, ogDict):
        # {'R46C04a': <Query.Query object at 0x7f28cce9a438>, 'R46C04b': <Query.Query object at 0x7f28cce9a470>}
        # No potential matches at all
        if len(ogDict) == 0:
            # print(f"no matches found for {self.queryName}")
            emptyQ = Query("", "", 0)
            emptyQ.query = self.queryName
            emptyQ.numHits = 0
            self.originalQueryObj = emptyQ
            self.foundOriginalQuery = False

        # Potential match in single sequence set
        if len(ogDict) == 1:
            # print(f"{self.queryName} is in single seq set")
            self.originalQueryObj = ogDict[self.queryName]

        if len(ogDict) > 1:
            x = True
            for query in ogDict:
                # print(query)
                if self.cleanSequence[:20] == ogDict[query].flankingSequence[:20]:
                    x = False
                    # print(f"{query} is in two seq set")
                    self.queryName = ogDict[query].query
                    self.originalQueryObj = ogDict[query]
            if x:
                # print(f"{query} is not identifiable in two seq set")
                emptyQ = Query("", "", 0)
                emptyQ.query = query
                self.originalQueryObj = emptyQ
                self.foundOriginalQuery = False

    def __dataLine__(self):
        # For printing sanger data to CSV
        return f"{self.queryName},{self.dsgfpMatches},{self.dsgfpPresent},{self.blastable}," + \
               f"{self.foundOriginalQuery},{self.bestSangerMatchesBestOriginal},{self.sangerSequenceFoundInGenome},{self.sangerQueryObj.genome}," + \
               f"{self.sangerQueryObj.chromosome},{self.originalQueryObj.genome}," + \
               f"{self.originalQueryObj.chromosome},{self.sangerQueryObj.sStart},{self.sangerQueryObj.sEnd}," + \
               f"{self.originalQueryObj.sStart},{self.originalQueryObj.sEnd},{self.sangerQueryObj.eValue}," + \
               f"{self.originalQueryObj.eValue},{self.sangerQueryObj.bitScore},{self.originalQueryObj.bitScore}," + \
               f"{self.sangerQueryObj.qStartStatus},{self.originalQueryObj.qStartStatus}\n"
        # return ""

    def __findBestHit__(self):
        if len(self.possibleQueryMatches) == 0:
            self.sangerQueryObj = Query("", "", 0)
            self.bestSangerMatchesBestOriginal = False
            self.sangerSequenceFoundInGenome = False
        for i in range(len(self.possibleQueryMatches)):
            if self.possibleQueryMatches[i].bestHitForAllele:
                self.sangerQueryObj = self.possibleQueryMatches[i]
                self.__validate__()
                break

    def __validate__(self):
        if self.sangerQueryObj.genome != self.originalQueryObj.genome:
            self.bestSangerMatchesBestOriginal = False
        if self.sangerQueryObj.chromosome != self.originalQueryObj.chromosome:
            self.bestSangerMatchesBestOriginal = False
        if self.sangerQueryObj.sStart != self.originalQueryObj.sStart:
            self.bestSangerMatchesBestOriginal = False

    def __readSeqFile__(self, filename, order):
        # 28329_C02_R179H03_5-Ds-3_C02_006.seq
        self.orderNum = order
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
        if len(self.sequence) > 0:
            self.__cleanConsecutiveNs__()

        if len(self.cleanSequence) < 22 or not self.dsgfpPresent:
            self.blastable = False

    def __printToFasta__(self):
        if self.blastable:
            return f">{self.queryName}\n{self.cleanSequence}\n"
        pass

    def __removeDsgfp__(self):
        dsIdx = 0
        if self.rawSequence.find(self.dsgfpEnd) != -1:
            dsIdx = self.rawSequence.find(self.dsgfpEnd)
            self.dsgfpMatches = True

        else:
            self.dsgfpMatches = False
            seqMatchObj = SequenceMatcher(lambda x: x == "ACGT", self.rawSequence, self.dsgfpEnd)
            matchObj = seqMatchObj.find_longest_match(0, len(self.rawSequence), 0,
                                                      len(self.dsgfpEnd))
            if matchObj.size >= len(self.dsgfpEnd) / 2:
                dsIdx = matchObj.a
            else:
                dsIdx = -1
                self.sequence = ""
                self.cleanSequence = ""
                self.dsgfpPresent = False
                return
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
        if queryObj.query == self.queryName:
            self.possibleQueryMatches.append(queryObj)

    def __SangerSeqFromJSON__(self, jsonObject):

        self.queryName = jsonObject["queryName"]
        self.metadata = jsonObject["metadata"]
        self.orderNum = jsonObject["orderNum"]
        self.primerMatch = jsonObject["primerMatch"]
        self.rawSequence = jsonObject["rawSequence"]
        self.sequence = jsonObject["sequence"]
        self.cleanSequence = jsonObject["cleanSequence"]
        self.dsgfpMatches = jsonObject["dsgfpMatches"]
        self.dsgfpPresent = jsonObject["dsgfpPresent"]
        self.blastable = jsonObject["blastable"]
        self.foundOriginalQuery = jsonObject["foundOriginalQuery"]
        self.bestSangerMatchesBestOriginal = jsonObject["bestSangerMatchesBestOriginal"]
        self.sangerSequenceFoundInGenome = jsonObject["sangerSequenceFoundInGenome"]
        # self.sangerQueryObj = jsonObject["sangerQueryObj"]
        # self.originalQueryObj = jsonObject["originalQueryObj"]
        self.dsgfpEnd = jsonObject["dsgfpEnd"]
        self.possibleQueryMatches = jsonObject["possibleQueryMatches"]

        if jsonObject["originalQueryObj"]:
            self.originalQueryObj = Query("", "", 0)
            self.originalQueryObj.__QueryFromJSON__(jsonObject["originalQueryObj"])
        else:
            self.originalQueryObj = None
        if jsonObject["sangerQueryObj"]:
            self.sangerQueryObj = Query("", "", 0)
            self.sangerQueryObj.__QueryFromJSON__(jsonObject["sangerQueryObj"])
        else:
            self.sangerQueryObj = None

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

    def __printDataFile__(self):
        filename = f"DIGIToutput/SangerSequences/{self.order}/QueryData/CSVfiles/SangerVsOriginal_{self.order}.csv"
        with open(filename, "w") as f:
            f.write(
                "Query,DsGFPMatch,DsGFPPresent,AbleToBlastSequence,OriginalQueryFound,BestSangerMatchesBestOrigin" +
                f"al,SangerSequenceFoundInAnyGenome,SangerGenome,SangerChromosome,OriginalGenome,OriginalChromosome,SangerStart,SangerEnd,Origin" +
                f"alStart,OriginalEnd,SangerEValue,OriginalEValue,SangerBitScore,OriginalBitScore,SangerQStartSta" +
                f"tus,OriginalQStartStatus\n")
            for s in self.sangerWorkingSet:
                if s:
                    f.write(self.sangerWorkingSet[s].__dataLine__())
                else:
                    f.write(f"{s}_is_missing\n")

    def __setBestQueryToSanger__(self):
        for s in self.sangerWorkingSet:
            self.sangerWorkingSet[s].__findBestHit__()

    def __appendSangerSeqObjsToWorkingSet__(self, orderDir):
        order = orderDir.split("/")[-1][-5:]
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
                if self.sangerWorkingSet[q].__printToFasta__():
                    outfile.write(self.sangerWorkingSet[q].__printToFasta__())
        return filename

    def __printToJson__(self):
        jsonfile = f"DIGIToutput/SangerSequences/{self.order}/QueryData/JSONfiles/SangerSet_{self.order}.json"
        jsonObject = {}
        for s in self.sangerWorkingSet:
            if self.sangerWorkingSet[s].queryName not in jsonObject:
                jsonObject[self.sangerWorkingSet[s].queryName] = []
            jsonObject[self.sangerWorkingSet[s].queryName] = dict(self.sangerWorkingSet[s])
        # newJSONobject = json.dumps(jsonObject, indent=4)
        newJSONobject = json.dumps(jsonObject, indent=4, default=lambda o: o.__dict__)

        with open(jsonfile, "w") as outfile:
            outfile.write(newJSONobject)

    def __printToGff__(self, filename):
        # Save working set in GFF file format. Used internally in class.
        qwsFormat = QueriesWorkingSet()
        for s in self.sangerWorkingSet:
            if self.sangerWorkingSet[s].queryName not in qwsFormat.workingSet:
                qwsFormat.workingSet[s] = None
            qwsFormat.workingSet[s] = self.sangerWorkingSet[s].sangerQueryObj

        newGffFile = GffFile()
        newGffFile.__formatFromQueriesWorkingSet__(qwsFormat.workingSet, "sanger_sequence")
        newGffFile.__writeToGffFile__(filename)

    def __readFromJson__(self):
        sangerJsonFile = f"DIGIToutput/SangerSequences/{self.order}/QueryData/JSONfiles/SangerSet_{self.order}.json"
        jsonFile = open(sangerJsonFile)
        dataFromJSON = json.load(jsonFile)
        jsonFile.close()

        for sangerSeq in dataFromJSON:
            newSangerSeq = SangerSeq(None)
            newSangerSeq.__SangerSeqFromJSON__(dataFromJSON[sangerSeq])
            if newSangerSeq.queryName not in self.sangerWorkingSet:
                self.sangerWorkingSet[newSangerSeq.queryName] = None
            self.sangerWorkingSet[newSangerSeq.queryName] = newSangerSeq

    def __addQueriesToSangerObjs__(self, blastOutputObj):
        for genome in blastOutputObj.hits:
            for allele in blastOutputObj.hits[genome]:
                if allele in self.sangerWorkingSet:
                    self.sangerWorkingSet[allele].possibleQueryMatches += \
                    blastOutputObj.hits[genome][allele]

    def __setOriginalQuery__(self, originalData):
        for query in self.sangerWorkingSet:
            querya = query + "a"
            queryb = query + "b"
            ogDict = {}
            if query in originalData.workingSet:
                ogDict[query] = originalData.workingSet[query]
            if querya in originalData.workingSet:
                ogDict[querya] = originalData.workingSet[querya]
            if queryb in originalData.workingSet:
                ogDict[queryb] = originalData.workingSet[queryb]
            # print(len(ogDict))
            self.sangerWorkingSet[query].__setOriginalQueryObj__(ogDict)



