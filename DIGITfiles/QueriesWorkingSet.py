import json
import logging

from GFF import *
from Query import *
from PrimerIDs import *
from Primer3Object import *


class QueriesWorkingSet:
    def __init__(self):
        self.workingSet = {}

    def __addToWorkingSet__(self, queryObj):
        # for adding a single Query object to workingSet as value (key is string Query.query)
        if queryObj.query not in self.workingSet:
            self.workingSet[queryObj.query] = None
            self.workingSet[queryObj.query] = queryObj
        else:
            logging.info(f"\t{queryObj.query} already exists in working set")

    def __getQuery__(self, query):
        if query in self.workingSet:
            return self.workingSet[query]
        return None

    def __createQueryStructFromJson__(self, filename):
        # Build workingSet from JSON file
        jsonFile = open(filename)
        dataFromJSON = json.load(jsonFile)
        jsonFile.close()
        for allele in dataFromJSON:
            if allele not in self.workingSet:
                self.workingSet[allele] = Query(None, None, None)
                self.workingSet[allele].__QueryFromJSON__(dataFromJSON[allele])
        self.__findOrderedPrimers__()

    def __updateQueryWithPrimer3Data__(self, outputDict, task):
        allele = outputDict["SEQUENCE_ID"].split("_")[-1]
        newP3Obj = Primer3Object(task)
        newP3Obj.__initValsFromP3Output__(outputDict)
        self.workingSet[allele].__updateQueryWithP3Output__(newP3Obj)

    def __lookupQuery__(self, query):
        # Specifically access Query object for Lookup task
        if query in self.workingSet:
            self.workingSet[query].__lookupPrint__()
        else:
            print(f"{query} is not in working set")

    def __printToJson__(self, filepath, flankseq):
        # Write workingSet to a JSON file
        filename = f"{filepath}/JSONfiles/WorkingSet_{flankseq}.json"
        jsonObject = {}
        for q in self.workingSet:
            if self.workingSet[q].query not in jsonObject:
                jsonObject[self.workingSet[q].query] = []
            jsonObject[self.workingSet[q].query] = dict(self.workingSet[q])
        newJSONobject = json.dumps(jsonObject, indent=4)
        with open(filename, "w") as outfile:
            outfile.write(newJSONobject)
        self.__printToCsv__(f"{filepath}/CSVfiles/WorkingSet_{flankseq}.csv")
        self.__printToGff__(f"{filepath}/GFFfiles/WorkingSet_{flankseq}.gff")

    def __printToCsv__(self, filename):
        csvName = filename.split(".")[0] + ".csv"
        with open(csvName, "w") as csvFile:
            csvFile.write(
                f"Query,Genome,BYearFamily,Chromosome,PercentIdentity,AlignmentLength,Mismatches,"
                f"GapOpens,q" +
                f"Start,qEnd,sStart,sEnd,eValue,BitScore,NumHitsForGenome,"
                f"PercentDifferenceToNextHit,Strand" +
                f",qStartStatus,PrimerNameLeft,PrimerSequenceLeft,PrimerLeftProductSize,"
                f"PrimerPenaltyLeft,P" +
                f"rimerLeftOrdered,PrimerLeftSangered,PrimerNameRight,PrimerSequenceRight,"
                f"PrimerRightProduc" +
                f"tSize,PrimerPenaltyRight,PrimerRightOrdered,PrimerRightSangered,"
                f"PrimerPairProductSize,Pri" +
                f"merPairPenalty\n")

            for q in self.workingSet:
                csvFile.write(self.workingSet[q].__workingSetCsvLine__())

    def __printToGff__(self, filename):
        gffName = filename.split(".")[0] + ".gff"
        newGffFile = GffFile("")
        newGffFile.__formatFromQueriesWorkingSet__(self.workingSet)
        newGffFile.__writeToGffFile__(gffName)

    def __getCoordinatesForFilterFasta__(self):
        coorA, coorB, coorW = {}, {}, {}
        for q in self.workingSet:
            if self.workingSet[q].genome.strip() == "A188v1":
                coorA[self.workingSet[q].query + "_wt"] = [self.workingSet[q].chromosome,
                                                           self.workingSet[q].wildtypeCoordinates]
                coorA[self.workingSet[q].query + "_up"] = [self.workingSet[q].chromosome,
                                                           self.workingSet[q].upperCoordinates]
                coorA[self.workingSet[q].query + "_lo"] = [self.workingSet[q].chromosome,
                                                           self.workingSet[q].lowerCoordinates]
            elif self.workingSet[q].genome.strip() == "B73v5":
                coorB[self.workingSet[q].query + "_wt"] = [self.workingSet[q].chromosome,
                                                           self.workingSet[q].wildtypeCoordinates]
                coorB[self.workingSet[q].query + "_up"] = [self.workingSet[q].chromosome,
                                                           self.workingSet[q].upperCoordinates]
                coorB[self.workingSet[q].query + "_lo"] = [self.workingSet[q].chromosome,
                                                           self.workingSet[q].lowerCoordinates]
            elif self.workingSet[q].genome.strip() == "W22v2":
                coorW[self.workingSet[q].query + "_wt"] = [self.workingSet[q].chromosome,
                                                           self.workingSet[q].wildtypeCoordinates]
                coorW[self.workingSet[q].query + "_up"] = [self.workingSet[q].chromosome,
                                                           self.workingSet[q].upperCoordinates]
                coorW[self.workingSet[q].query + "_lo"] = [self.workingSet[q].chromosome,
                                                           self.workingSet[q].lowerCoordinates]
        return coorA, coorB, coorW

    def __assignFilterFastaData__(self, seqData):
        for allele in seqData:
            upper, lower, wildtype = "", "", ""
            if allele.find("wt") != -1:
                self.workingSet[allele[:-3]].wildtypeSequence = seqData[allele]
            if allele.find("up") != -1:
                self.workingSet[allele[:-3]].upperSequence = seqData[allele]
            if allele.find("lo") != -1:
                self.workingSet[allele[:-3]].lowerSequence = seqData[allele]
        for allele in self.workingSet:
            if self.workingSet[allele].upperSequence and self.workingSet[allele].lowerSequence:
                self.workingSet[allele].__buildInsertionSequence__()
            else:
                self.workingSet[allele].insertionSequence = "__failed__"

    def __setPrimerIDs__(self):
        primerObject = PrimerIDs()
        for query in self.workingSet:
            leftNameExists, rightNameExists = False, False
            for primerName in primerObject.listOfPrimers:
                pts = primerName.split("_")
                # mr6091a, R98E09b
                noNum = pts[0][-1] + "_" + pts[1]
                # is left primer name in there? If so, rename left primer
                if noNum == self.workingSet[query].primerNameLeft:
                    self.workingSet[query].primerNameLeft = primerName
                    leftNameExists = True
                # is right primer in there? if so, rename right primer
                if noNum == self.workingSet[query].primerNameRight:
                    self.workingSet[query].primerNameRight = primerName
                    rightNameExists = True

            if not leftNameExists:
                self.workingSet[query].primerNameLeft = ("mr" + str(primerObject.numPrimers) +
                                                         self.workingSet[query].primerNameLeft)
                primerObject.listOfPrimers.append(self.workingSet[query].primerNameLeft)
                primerObject.numPrimers += 1

            if not rightNameExists:
                self.workingSet[query].primerNameRight = ("mr" + str(primerObject.numPrimers) +
                                                          self.workingSet[query].primerNameRight)
                primerObject.listOfPrimers.append(self.workingSet[query].primerNameRight)
                primerObject.numPrimers += 1
            primerObject.__rewriteRecordFile__()

    def __writeP3InputFile__(self, inputFilepath, task):
        '''
        Initializes a Primer3Object for each query in queriesWorkingSet, and creates the input
        file that will be used with
        Primer3. Return value is the path to the Primer3 input file (str).
        Appears in:
            GetPrimers.py
        '''
        with open(inputFilepath, "w") as p3file:
            for allele in self.workingSet:
                p3obj = Primer3Object(task)
                p3obj.__initValsFromQueryObj__(self.workingSet[allele])

                if task == "generic":
                    p3obj.__buildInputStrGeneric__()
                else:
                    p3obj.__buildInputStrValidate__()
                p3file.write(p3obj.inputStr)

    def __makeFastaFromPrimers__(self, filename):
        with open(filename, "w+") as fafile:
            for q in self.workingSet:
                if self.workingSet[q].primerSequenceLeft != "FAIL":
                    fafile.write(
                        f">{self.workingSet[q].primerNameLeft}_ogRefGen_"
                        f"{self.workingSet[q].genome}\n")
                    fafile.write(f"{self.workingSet[q].primerSequenceLeft}\n")
                if self.workingSet[q].primerSequenceRight != "FAIL":
                    fafile.write(
                        f">{self.workingSet[q].primerNameRight}_ogRefGen_"
                        f"{self.workingSet[q].genome}\n")
                    fafile.write(f"{self.workingSet[q].primerSequenceRight}\n")

    def __createBadQueryStruct__(self):
        queriesWithBadPrimers = QueriesWorkingSet()
        for q in self.workingSet:
            if self.workingSet[q].primerSequenceRight == "FAIL" or self.workingSet[
                q].primerSequenceLeft == "FAIL":
                queriesWithBadPrimers.__addToWorkingSet__(self.workingSet[q])
        return queriesWithBadPrimers

    def __makePrimerDataFile__(self, filename):
        listOfPrimerStrings = []
        for q in self.workingSet:
            listOfPrimerStrings.append(self.workingSet[q].__leftPrimerDataLine__())
            listOfPrimerStrings.append(self.workingSet[q].__rightPrimerDataLine__())
            listOfPrimerStrings.sort()

        with open(filename, "w") as primerFile:
            primerFile.write(
                "PrimerID,Allele,ReferenceGenome,PrimerSequence,MatchingDsGFPPrimer,"
                "ExpectedProductSizeWithD" +
                "sPrimer,ExpectedWTProductSize,TM,PrimerPenaltyWithDSPrimer,PrimerPenaltyWT,"
                "BlastBitScore,To" +
                "talBlastHitsForRefGenome,QueryStartStatus\n")
            for q in listOfPrimerStrings:
                primerFile.write(q)

    def __findOrderedPrimers__(self):
        orderedPrimersDict = {}
        with open("PutOrderedPrimersHere/OrderedPrimers", "r") as opfile:
            for line in opfile:
                primerName, primerSeq = line.split(",")
                orderedPrimersDict[primerName.strip()] = primerSeq.strip()

        for q in self.workingSet:
            if self.workingSet[q].primerNameLeft in orderedPrimersDict:
                self.workingSet[q].primerLeftOrdered = True
                if self.workingSet[q].primerSequenceLeft != orderedPrimersDict[
                    self.workingSet[q].primerNameLeft]:
                    print(
                        f"Error: Predicted left primer and ordered left primer do not match for "
                        f"{q}")

            if self.workingSet[q].primerNameRight in orderedPrimersDict:
                self.workingSet[q].primerRightOrdered = True
                if self.workingSet[q].primerSequenceRight != orderedPrimersDict[
                    self.workingSet[q].primerNameRight]:
                    print(
                        f"Error: Predicted right primer and ordered right primer do not match for {q}")


