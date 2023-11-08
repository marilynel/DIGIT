'''
QueriesWorkingSet is a class that stores the "best" data from a BlastOutput object. Queries are
stored as Query objects
in a dictionary struct. Each query should only appear once within one QueriesWorkingSet object.

Written by: Marilyn Leary
'''

import json
import logging
import os

from GFF import *
from Query import *
from PrimerIDs import *
from Primer3Object import *


class QueriesWorkingSet:
    def __init__(self):
        self.workingSet = {}

    def __buildCompleteWorkingSet__(self, b73spec):
        # Fills QueriesWorkingSet object with all previously found query data. Used to compare
        # sanger sequencing results
        # with original generated data. Appears in BlastSangerOutput and ParseSangerResults.
        for subdir, dirs, files in os.walk("DIGIToutput/FlankingSequences"):
            if subdir.find("JSON") != -1:
                for file in files:
                    if file.find("Working") != -1 and file.find("json") != -1 and file.find(
                            b73spec) != -1:
                        self.__createQueryStructFromJson__(subdir + "/" + file)

    def __buildWorkingSetFromBlastOutputObj__(self, blastOutputObj):
        # Does what the function name says. Used in SequencesFromBlast.
        print(f"calling __buildWorkingSetFromBlastOutputObj__")
        for gen in blastOutputObj.hits:
            for q in blastOutputObj.hits[gen]:
                for i in range(0, len(blastOutputObj.hits[gen][q])):
                    if blastOutputObj.hits[gen][q][i].bestHitForAllele:
                        self.__addToWorkingSet__(blastOutputObj.hits[gen][q][i])

    def __addToWorkingSet__(self, queryObj):
        # For adding a single Query object to workingSet as value (key is string Query.query).
        # Used internally in class.
        if queryObj.query not in self.workingSet:
            self.workingSet[queryObj.query] = None
            self.workingSet[queryObj.query] = queryObj
        else:
            logging.info(f"\t{queryObj.query} already exists in working set")

    def __getOriginalFlankingSequences__(self, flankseq):
        # Add flanking sequences to Query objects in working set from original fasta file that
        # was blasted. Used in
        # SequencesFromBlast.
        flankseqDict = {}
        with open(f"PutFlankingSequenceFilesHere/{flankseq}.fasta", "r") as fafile:
            newKey = ""
            for line in fafile:
                if line[0] == ">":
                    newKey = line[1:].strip()
                    flankseqDict[newKey] = ""
                else:
                    flankseqDict[newKey] = line.strip()
        for query in self.workingSet:
            if query in flankseqDict:
                self.workingSet[query].flankingSequence = flankseqDict[query]

    def __getQuery__(self, query):
        # Makes QueriesWorkingSet iterable.
        if query in self.workingSet:
            return self.workingSet[query]
        return None

    def __createQueryStructFromJson__(self, filename):
        # Build workingSet from JSON file. Used in GetPrimerIncidenceRate, GetPrimers,
        # VerifyPrimers, and internally in
        # this class.
        jsonFile = open(filename)
        dataFromJSON = json.load(jsonFile)
        jsonFile.close()
        for allele in dataFromJSON:
            if allele not in self.workingSet:
                self.workingSet[allele] = Query("", "", 0)
                self.workingSet[allele].__QueryFromJSON__(dataFromJSON[allele])

    def __updateQueryWithPrimer3Data__(self, outputDict, task):
        # Update Query objects in struct with parsed data from Primer3 initial output. Used in
        # Utils.
        allele = outputDict["SEQUENCE_ID"].split("_")[-1]
        newP3Obj = Primer3Object(task)
        newP3Obj.__initValsFromP3Output__(outputDict)
        self.workingSet[allele].__updateQueryWithP3Output__(newP3Obj)

    def __lookupQuery__(self, query):
        # Specifically access Query object for Lookup task
        # TODO: not currently in use
        if query in self.workingSet:
            self.workingSet[query].__lookupPrint__()
        else:
            print(f"{query} is not in working set")

    def __printToJson__(self, filepath, flankseq):
        # Write workingSet to a JSON file. Used in GetPrimers, SequencesFromBlast,
        # and VerifyPrimers.
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
        # Save working set in CSV file. Used internally in class.
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
        # Save working set in GFF file format. Used internally in class.
        gffName = filename.split(".")[0] + ".gff"
        newGffFile = GffFile()
        newGffFile.__formatFromQueriesWorkingSet__(self.workingSet, "nucleotide_match")
        newGffFile.__writeToGffFile__(gffName)

    def __getCoordinatesForFilterFasta__(self):
        # Finds and sorts location data for whole working set to send to Filterfasta. Used in Utils.
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
        # Assign sequence data returned from FilterFasta to the appropriate Query object. Used in
        # Utils.
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
        for p in primerObject.listOfPrimers:  # Ex: p = "mr5161a_R155D07b"
            pts = p.split("_")  # pts = ["mr5161a", "R155D07b"]
            qName = pts[1]  # qName = "R155D07b"
            num = pts[0][2:-1]  # num = "5161"
            side = pts[0][-1]  # side = "a"
            preNum = side + "_" + qName  # preNum = "a_R155D07b"
            if qName in self.workingSet:
                # Check if primer p matches primerNameLeft for Query q
                if preNum == self.workingSet[qName].primerNameLeft:
                    # If so, primerNameLeft has an assigned number already. Update workingSet
                    self.workingSet[qName].primerNameLeft = p
                # If its not on the left, check if primer p matches primerNameRight for Query q
                elif preNum == self.workingSet[qName].primerNameRight:
                    # If so, primerNameRIght has an assigned number already. Update workingSet
                    self.workingSet[qName].primerNameRight = p
                else:
                    # This condition should never be met
                    print(f"ERROR with {p}")

        # check working set for queries that do not have assoc numbers for primers already
        for q in self.workingSet:
            if self.workingSet[q].numHits == 0:
                continue
            if self.workingSet[q].primerNameLeft[0:2] != "mr":
                self.workingSet[q].primerNameLeft = f"mr{primerObject.numPrimers}" + \
                                                    self.workingSet[q].primerNameLeft
                primerObject.__append__(self.workingSet[q].primerNameLeft)
            if self.workingSet[q].primerNameRight[0:2] != "mr":
                self.workingSet[q].primerNameRight = f"mr{primerObject.numPrimers}" + \
                                                     self.workingSet[q].primerNameRight
                primerObject.__append__(self.workingSet[q].primerNameRight)
        primerObject.__rewriteRecordFile__()
        self.__havePrimersBeenOrderedOrTestedYet__()

    def __writeP3InputFile__(self, inputFilepath, task):
        # Initializes a Primer3Object for each query in queriesWorkingSet, and creates the input
        # file that will be used
        # with Primer3. Return value is the path to the Primer3 input file (str).   
        with open(inputFilepath, "w") as p3file:
            for allele in self.workingSet:
                if self.workingSet[allele].numHits != 0:
                    p3obj = Primer3Object(task)
                    p3obj.__initValsFromQueryObj__(self.workingSet[allele])

                    if task == "generic":
                        p3obj.__buildInputStrGeneric__()
                    else:
                        p3obj.__buildInputStrValidate__()
                    p3file.write(p3obj.inputStr)

    def __makeFastaFromPrimers__(self, filename):
        # Makes a fasta file out of the primers in a Query object.
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
        # Create a new QueriesWorkingSet object to hold failed primers to be saved in a different
        # location.
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
                f"PrimerID,Allele,ReferenceGenome,PrimerSequence,MatchingDsGFPPrimer,"
                f"ExpectedProductSize" +
                f"WithDsPrimer,ExpectedWTProductSize,TM,PrimerPenaltyWithDSPrimer,"
                f"PrimerPenaltyWT,BlastB" +
                f"itScore,TotalBlastHitsForRefGenome,QueryStartStatus,Ordered,Sangered\n")
            for q in listOfPrimerStrings:
                primerFile.write(q)

    def __havePrimersBeenOrderedOrTestedYet__(self):
        # Check if primers have been ordered or sangered already
        orderedPrimers = {}
        with open("PutOrderedPrimersHere/OrderedPrimers", "r") as opfile:
            for line in opfile:
                primerName, primerSequence = line.split(",")
                orderedPrimers[primerName.strip()] = primerSequence.strip()

        sangeredPrimers = []
        with open("PutOrderedPrimersHere/SucessfulPrimersFromSanger", "r") as spfile:
            for line in spfile:
                sangeredPrimers.append(line.strip())

        for q in self.workingSet:
            if self.workingSet[q].primerNameLeft in orderedPrimers:
                self.workingSet[q].primerLeftOrdered = True
                if self.workingSet[q].primerSequenceLeft != orderedPrimers[
                    self.workingSet[q].primerNameLeft]:
                    print(f"ERROR: {self.workingSet[q].primerNameLeft} predicted primer " +
                          f"{self.workingSet[q].primerSequenceLeft} does not match ordered primer "
                          f"" +
                          f"{orderedPrimers[self.workingSet[q].primerNameLeft]}")

            if self.workingSet[q].primerNameLeft in sangeredPrimers:
                self.workingSet[q].primerLeftSangered = True

            if self.workingSet[q].primerNameRight in orderedPrimers:
                self.workingSet[q].primerRightOrdered = True

                if self.workingSet[q].primerSequenceRight != orderedPrimers[
                    self.workingSet[q].primerNameRight]:
                    print(f"ERROR: {self.workingSet[q].primerNameRight} predicted primer " +
                          f"{self.workingSet[q].primerSequenceRight} does not match ordered "
                          f"primer " +
                          f"{orderedPrimers[self.workingSet[q].primerNameRight]}")

            if self.workingSet[q].primerNameRight in sangeredPrimers:
                self.workingSet[q].primerRightSangered = True