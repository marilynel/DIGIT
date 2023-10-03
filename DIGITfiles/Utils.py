# from collections import ChainMap
from difflib import SequenceMatcher

import json
import os
import re
import subprocess
# import time

from FilterFasta import *
from GFF import *
# from Primer3Object import Primer3Object
from Primer import Primer
from Query import Query
from QueriesWorkingSet import QueriesWorkingSet


### TODO: redo this what the hell
def thisIsDumb(query, originalData):
    '''
    Appears in:
        Utils.py
    '''
    if query in originalData.workingSet:
        return {
            "genome": originalData.workingSet[query].genome,
            "chromosome": originalData.workingSet[query].chromosome,
            "sStart": originalData.workingSet[query].sStart,
            "sEnd": originalData.workingSet[query].sEnd,
            "eValue": originalData.workingSet[query].eValue,
            "bitScore": originalData.workingSet[query].bitScore,
            "numHits": originalData.workingSet[query].numHits,
            "strand": originalData.workingSet[query].strand,
            "qStartStatus": originalData.workingSet[query].qStartStatus
        }
    return {"genome": None, "chromosome": None, "sStart": None, "sEnd": None, "eValue": None,
            "bitScore": None, "numHits":
                None, "strand": None, "qStartStatus": None}


# TODO: not yet in use
def addToTestedPrimersList(testedPrimers):
    with open("PutOrderedPrimersHere/SucessfulPrimersFromSanger", "r") as testedFile:
        for primer in testedPrimers:
            testedFile.write(primer)


def allQueriesToJSON(filename, allQueries):
    '''
    Converts dictionary to JSON object and writes it to a file.
    Appears in:
        ParseSangerResults.py
        SequencesFromBlast.py
    '''
    jsonObject = {"A188v1": {}, "B73v5": {}, "W22v2": {}}

    for database in allQueries:
        for query in allQueries[database]:
            if query not in jsonObject[database]:
                jsonObject[database][query] = []
            for hit in allQueries[database][query]:
                jsonObject[database][query].append(dict(hit))
    newJSONbject = json.dumps(jsonObject, indent=4)

    with open(filename, "w") as outfile:
        outfile.write(newJSONbject)
    return


def bestGenomesWrite(bestBitScore):
    '''
    Find the best of the best hits, return stringified list for writing to file.
    Appears in:
        Utils.py
    '''
    genomeList = []
    bestScore = max(bestBitScore.values())
    for genome in bestBitScore:
        if bestBitScore[genome] == bestScore:
            genomeList.append(genome)
        else:
            bestBitScore[genome] = 0
    return str(genomeList)


# TODO: delete? no longer in use
def buildCoordinateSetsForFilterFasta(queriesWorkingSet):
    pass


#    coorA, coorB, coorW = {}, {}, {}
#   for q in queriesWorkingSet:
#      if queriesWorkingSet[q].genome.strip() == "A188v1":
#         coorA[queriesWorkingSet[q].query + "_wt"] = [queriesWorkingSet[q].chromosome,
#                                                     queriesWorkingSet[q].wildtypeCoordinates]
#       coorA[queriesWorkingSet[q].query + "_up"] = [queriesWorkingSet[q].chromosome,
#                                                   queriesWorkingSet[q].upperCoordinates]
#     coorA[queriesWorkingSet[q].query + "_lo"] = [queriesWorkingSet[q].chromosome,
#                                                 queriesWorkingSet[q].lowerCoordinates]
#        elif queriesWorkingSet[q].genome.strip() == "B73v5":
#           coorB[queriesWorkingSet[q].query + "_wt"] = [queriesWorkingSet[q].chromosome,
#                                                       queriesWorkingSet[q].wildtypeCoordinates]
#         coorB[queriesWorkingSet[q].query + "_up"] = [queriesWorkingSet[q].chromosome,
#                                                     queriesWorkingSet[q].upperCoordinates]
#       coorB[queriesWorkingSet[q].query + "_lo"] = [queriesWorkingSet[q].chromosome,
#                                                   queriesWorkingSet[q].lowerCoordinates]
# elif queriesWorkingSet[q].genome.strip() == "W22v2":
#    coorW[queriesWorkingSet[q].query + "_wt"] = [queriesWorkingSet[q].chromosome,
#                                                queriesWorkingSet[q].wildtypeCoordinates]
#  coorW[queriesWorkingSet[q].query + "_up"] = [queriesWorkingSet[q].chromosome,
#                                              queriesWorkingSet[q].upperCoordinates]
# coorW[queriesWorkingSet[q].query + "_lo"] = [queriesWorkingSet[q].chromosome,
#                                            queriesWorkingSet[q].lowerCoordinates]
# return coorA, coorB, coorW


def buildFullWorkingSet(b73only):
    '''
    Collect all query data from existing working set JSON files and save to a dictionary. Return dictionary. Argument
    is a boolean that specifies which datasets to use.
    Appears in:
        Lookup.py
        ParseSangerResults.py
    '''
    completeWorkingSet = QueriesWorkingSet()
    for subdir, dirs, files in os.walk("DIGIToutput/FlankingSequences"):
        if subdir.find("JSON") != -1:
            for file in files:
                newSet = QueriesWorkingSet()
                if b73only:
                    if file.find("Working") != -1 and file.find("json") != -1 and file.find(
                            "B73") != -1:
                        newSet.__createQueryStructFromJson__(subdir + "/" + file)
                else:
                    if file.find("Working") != -1 and file.find("json") != -1 and file.find(
                            "B73") == -1:
                        newSet.__createQueryStructFromJson__(subdir + "/" + file)
                for query in newSet.workingSet:
                    completeWorkingSet.__addToWorkingSet__(newSet.workingSet[query])
    return completeWorkingSet


def buildWorkingSetFromAllQueries(allQueries):
    '''
    Iterates through entire dataset, and creates a new smaller working dataset with just the "best" hit for each query.
    Appears in:
        ParseSangerResults.py
        SequencesFromBlast.py
    '''
    queriesWorkingSet = QueriesWorkingSet()
    for gen in allQueries:
        for q in allQueries[gen]:
            for i in range(0, len(allQueries[gen][q])):
                if allQueries[gen][q][i].bestHitForAllele:
                    queriesWorkingSet.__addToWorkingSet__(allQueries[gen][q][i])
    return queriesWorkingSet


def callPrimerBlastScript(pathToFastaInput, pathToDesitinationDir):
    '''
    Calls script to blast primers against the A188, B73, and W22 genomes. Specialized script for short sequences.
    Appears in:
        VerifyPrimers.py
    '''
    subprocess.run(
        ["sh", f"./DIGITfiles/RunBlastPrimers.sh", pathToFastaInput, pathToDesitinationDir])

    print(
        f"SGE may take a while to run. To check the status of your jobs, enter 'qstat' in the command line. Running" +
        f" time may vary wildly.\n")


def callBlastScript(pathToFastaInput, pathToDesitinationDir, fastaContents):
    '''
    Calls script to blast the initial flanking sequences against the A188, B73, and W22 genomes.
    Appears in:
        DIGIT.py
        BlastSangerOutput.py
    '''
    subprocess.run(
        ["sh", f"./DIGITfiles/RunBlastInitial.sh", pathToFastaInput, pathToDesitinationDir,
         fastaContents])

    print(
        f"SGE may take a while to run. To check the status of your jobs, enter 'qstat' in the command line. Running" +
        f" time may vary wildly.\n")


def checkConsecutiveNs(sequence):
    '''
    Look for 2 N's in a row in a sequence. N indicates an unknown base pair. Returns an index to the consectutive N's if
    found or nonsense value (-2) if not.
    Appears in:
        Utils.py
    '''
    previousBasePairs = sequence[0]
    for i in range(0, len(sequence)):
        if previousBasePairs == "N" and sequence[i] == "N":
            return i - 1
        else:
            previousBasePairs = sequence[i]
    return -2


def checkCoordinates(primer, query):
    '''
    Appears in:
        Utils.py
    '''
    if matchChrAndGenome(primer, query):
        if query.strand == 1:
            if (primer.sStart >= query.wildtypeCoordinates[0]) and (
                    primer.sEnd <= query.wildtypeCoordinates[1]):
                return True
        elif query.strand == -1:
            if (primer.sStart <= query.wildtypeCoordinates[0]) and (
                    primer.sEnd >= query.wildtypeCoordinates[1]):
                return True
    return False


# TODO: remove if it runs ok 9/29
def checkIfPrimerHasID():
    '''
    Create and return a list of primers that have already been developed.
    TODO: do I need this in here????

    '''
    pass


#    listPreviousPrimers = []
#    with open("DIGITfiles/compListPreviousPrimers", "r") as prevFile:
##        for line in prevFile:
#           listPreviousPrimers.append(line.strip())
#   return listPreviousPrimers


def findBlastOutputFiles(dirname, allQueries, genomeSpec):
    '''
    This function locates the Blast output file for the user selected set of flanking sequences.
    Appears in:
        GetPrimerIncidenceRate.py
        ParseSangerResults.py
        SequencesFromBlast.py
    '''
    for subdir, dirs, files in os.walk(dirname):
        for oneFile in files:
            filepath = os.path.join(subdir, oneFile)
            if filepath.find(".tab") != -1 and filepath.find(genomeSpec) != -1:
                filename = filepath.split("/")[-1]
                genome = filename.split("_")[0]
                allQueries[genome] = processBlastOutput(filepath, genome)
    return


# TODO: what is this for and do I need it?
def findNumHitsPerPrimer(allPrimers):
    '''
    Appears in:
        none
    '''
    pass


#    numHitsDir = {}
#   for genome in allPrimers:
#      if genome not in numHitsDir:
#         numHitsDir[genome] = {}
#    for primer in allPrimers[genome]:
#       for hit in allPrimers[genome][primer]:
#          parts = hit.query.split("_")
#         primerName = parts[0] + "_" + parts[1]
#        if primerName not in numHitsDir[genome]:
#           numHitsDir[genome][primerName] = 0
#      numHitsDir[genome][primerName] = hit.numHits
# return numHitsDir


def getBestQuery(genome, query, allQueries):
    '''
    Sort hits to identify best hit for a query. Find the percent difference between the best and second best hits.
    Appears in:
        Utils.py
    '''
    bestQuery = allQueries[genome][query][0]
    secondBestQuery = None
    for i in range(1, len(allQueries[genome][query])):
        if allQueries[genome][query][i].bitScore > bestQuery.bitScore:
            secondBestQuery = bestQuery
            bestQuery = allQueries[genome][query][i]

        elif allQueries[genome][query][i].bitScore <= bestQuery.bitScore:
            if not secondBestQuery or allQueries[genome][query][
                i].bitScore >= secondBestQuery.bitScore:
                secondBestQuery = allQueries[genome][query][i]
    perDiff = -1
    if secondBestQuery:
        perDiff = abs((((
                                bestQuery.bitScore - secondBestQuery.bitScore) / (
                                    bestQuery.bitScore + secondBestQuery.bitScore)) / 2) * 100)

    for i in range(0, len(allQueries[genome][query])):
        if allQueries[genome][query][i] == bestQuery:
            allQueries[genome][query][i].bestAlleleForGenome = True
            allQueries[genome][query][i].percentDiff = perDiff

    return


# TODO: remove? not in use any more
def getDataFromJSON(filename):
    '''
    Parse a JSON file and return a JSON object (key:value format)
    Appears in:
        Utils.py
    '''
    pass


#    jsonFile = open(filename)
#   dataFromJSON = json.load(jsonFile)
#  jsonFile.close()

# return dataFromJSON


def getNumHits(line, numHits):
    '''
    Parse a line for string indicating the number of hits for that query in a genome. If line does not have that string,
    return num_hits as original value.
    Appears in:
        Utils.py
    '''
    blastData = line.split(' ')
    if blastData[2] == "hits":
        return int(blastData[1])
    return numHits


# TODO: not yet used?
def getOrderedAndTestedPrimers():
    orderedPrimers, testedPrimers = [], []
    with open("PutOrderedPrimersHere/OrderedPrimers", "r") as orderFile:
        for line in orderFile:
            primerName = line.split(",")[0]
            orderedPrimers.append(primerName)
    with open("PutOrderedPrimersHere/SucessfulPrimersFromSanger", "r") as testedFile:
        for line in testedFile:
            primerName = line.split(",")[0]
            testedPrimers.append(primerName)
    return orderedPrimers, testedPrimers


def getSequence(filename, dsgfpEnd):
    '''
    Read a .seq file and parse the data from it. Returns genomic sequence, the index where the last base pair of
    dsgfpEnd is found, and a boolean value where False indicates an exact match (not notExact) and True indicates a non-
    exact match (notExact).
    Appears in:
        Utils.py
    '''
    sequence = ""
    with open(filename, "r+") as sangerFile:
        for line in sangerFile:
            sequence += line.strip()

    if sequence.find(dsgfpEnd) != -1:
        return sequence, sequence.find(dsgfpEnd), False
    else:
        seqMatchObj = SequenceMatcher(lambda x: x == "ACGT", sequence, dsgfpEnd)
        matchObj = seqMatchObj.find_longest_match(0, len(sequence), 0, len(dsgfpEnd))
        return sequence, matchObj.a, True


def getSequenceID(metadata, notExact):
    '''
    Parse sequence ID (elsewhere known as query or allele) from the file name. Append string indicating its an imperfect
    match if relevant. Returns sequence ID.
    Appears in:
        Utils.py
    '''
    seqID = re.findall(r"R\d{1,4}[A-Z]*\d{1,4}", metadata)
    if notExact:
        return seqID[0] + "_ImperfectMatchDsGfp"
    return seqID[0]


def lookupFamilyAllele():
    '''
    Appears in:
        Utils.py
    '''
    familyDict = {}
    with open("DIGITfiles/YearB2022FlankingSequenceFamilyData", "r") as familyfile:
        familyfile.readline()
        for line in familyfile:
            allele, family = line.split("\t")
            if allele not in familyDict:
                familyDict[allele] = []
            familyDict[allele].append(int(family.strip()))
    return familyDict


def makeDirectories(necessaryDirs):
    '''
    Appears in:
        DIGIT.py
        BlastSangerOutput.py
        SequencesFromBlast.py
        Utils.py
        VerifyPrimers.py
    necessaryDirs is a list of strings
    '''
    for dir in necessaryDirs:
        if not os.path.exists(dir):
            os.makedirs(dir)


def makeFastaString(folder, metadata, sequence, setstring, notExact):
    '''
    Aggregates .seq data into a fasta-formatted string to be written to file. Returns that string.
    Appears in:
        Utils.py
    '''
    seqID = getSequenceID(metadata, notExact)
    fastaFormatStr = f">{seqID} {metadata}\n"
    endIdx = checkConsecutiveNs(sequence)
    if endIdx == -2:
        # There are no consecutive ends -> given sequence can be used in full as passed to this function
        fastaFormatStr += sequence + "\n"
    elif endIdx == -1 or endIdx == 0:
        fastaFormatStr += "ERROR\n"
    else:
        # sequence is further truncated before occurrance of consecutive Ns
        fastaFormatStr += sequence[:endIdx] + "\n"
    return fastaFormatStr


def makePrimerDict(allPrimers):
    '''
    Appears in:
        GetPrimerIncidenceRate.py
    '''
    primerDict = {"A188v1": {}, "B73v5": {}, "W22v2": {}}
    for genome in allPrimers:
        for primer in allPrimers[genome]:
            for i in range(0, len(allPrimers[genome][primer])):
                newPrimer = Primer(allPrimers[genome][primer][i])
                if newPrimer.query not in primerDict[genome]:
                    primerDict[genome][newPrimer.query] = []
                primerDict[genome][newPrimer.query].append(newPrimer)
    return primerDict


# TODO: delete this
def makePrimerList(qWorkSet):
    pass
    # primerList = []
    # for q in qWorkSet:
    #   leftPrimerMatch, rightPrimerMatch = "", ""

    #  if qWorkSet[q].sideMatchGFP3UTR == "right" and qWorkSet[q].sideMatch3DsgG == "left":
    #     leftPrimerMatch, rightPrimerMatch = "GFP3UTR", "3DsgG"

    # elif qWorkSet[q].sideMatchGFP3UTR == "left" and qWorkSet[q].sideMatch3DsgG == "right":
    #   leftPrimerMatch, rightPrimerMatch = "3DsgG", "GFP3UTR"


#        else:
#           leftPrimerMatch, rightPrimerMatch = "FAIL", "FAIL"

#      primerList.append(f"{qWorkSet[q].primerNameLeft},{qWorkSet[q].query},{qWorkSet[q].genome}," +
#                       f"{qWorkSet[q].primerSequenceLeft},{leftPrimerMatch},{qWorkSet[q].primerLeftProductSize}," +
#                      f"{qWorkSet[q].primerPairProductSize},{qWorkSet[q].tmLeft},{qWorkSet[q].primerPenaltyLeft}," +
#                     f"{qWorkSet[q].primerPairPenalty},{qWorkSet[q].bitScore},{qWorkSet[q].numHits}," +
#                    f"{qWorkSet[q].qStartStatus}\n")

# primerList.append(f"{qWorkSet[q].primerNameRight},{qWorkSet[q].query},{qWorkSet[q].genome}," +
#                  f"{qWorkSet[q].primerSequenceRight},{rightPrimerMatch}," +
#                 f"{qWorkSet[q].primerRightProductSize},{qWorkSet[q].primerPairProductSize}," +
#                f"{qWorkSet[q].tmRight},{qWorkSet[q].primerPenaltyRight},{qWorkSet[q].primerPairPenalty}," +
#               f"{qWorkSet[q].bitScore},{qWorkSet[q].numHits},{qWorkSet[q].qStartStatus}\n")

# primerList.sort()
# return primerList


def makeStrBestInGenome(genome, allele, bestBitScore, allQueries):
    '''
    Find the specific hit that was specified as the best for that query ID and that genome (Query.best_for_genome ==
    True) and return a string reporting data for the best_queries_by_genome.csv outfile.
    Appears in:
        Utils.py
    '''
    if allele in allQueries[genome]:
        for i in range(0, len(allQueries[genome][allele])):
            if allQueries[genome][allele][i].bestAlleleForGenome == True:
                bestBitScore[genome] = allQueries[genome][allele][i].bitScore
                return allQueries[genome][allele][i].__makeBestGenomeString__()
    return "none,none,none,"


def matchChrAndGenome(primer, query):
    '''
    Appears in:
        Utils.py
    '''
    if primer.chromosome != query.chromosome:
        return False
    if primer.genome != query.genome:
        return False
    return True


def pickGenome(allele, bestBitScore, allQueries):
    '''
    Selects best genome.
    Appears in:
        Utils.py
    '''
    if bestBitScore["B73v5"] != 0:
        workingQuerySelection("B73v5", allele, allQueries)
    elif bestBitScore["W22v2"] != 0:
        workingQuerySelection("W22v2", allele, allQueries)
    elif bestBitScore["A188v1"] != 0:
        workingQuerySelection("A188v1", allele, allQueries)
    else:
        print("you really messed something up")
    return


def prepFasta(orderDir):
    '''
    Create a fasta-formatted file from the .seq files contained in the given directory (orderDir). Returns path to
    resultant fasta file.
    Appears in:
        BlastSangerOutput.py
    '''
    # 34 BP tail at the end of the DsGFP insertion sequence
    dsgfpEnd = "TTTTACCGACCGTTCCCGACCGTTTTCATCCCTA"
    order = orderDir.split("/")[-1]
    processedSeqFilesContents = ""
    for subdir, dirs, files in os.walk(orderDir):
        for oneFile in files:
            filename = os.path.join(subdir, oneFile)
            if filename.find(".seq") != -1:
                sequence, idx, notExact = getSequence(filename, dsgfpEnd)
                # Note: sequence[idx+len(dsgfpEnd):-22] is the sequence from after the DsGFP insertion to 22 base pairs
                # prior to the end of the sanger sequence (22 bp tail cut off due to inaccuracy of sequencing).
                processedSeqFilesContents += makeFastaString(
                    orderDir, filename[:-4], sequence[idx + len(dsgfpEnd):-22], order, notExact)

    filename = f"PutSangerOutputFilesHere/{order}/{order}.fasta"
    with open(filename, "w") as outfile:
        outfile.write(processedSeqFilesContents)
    return filename


def processBlastOutput(filename, genome):
    '''
    Parse blast output files line by line, gathering data that will be used to create Query objects, which will be added
    to the allQueries dataset. Functions get_num_hits() and get_genome() are called to help parse hash (#) lines.
    Appears in:
        Utils.py
    '''
    bYearFamilies = lookupFamilyAllele()
    genomeDB = {}
    numHits = -1
    with open(filename, "r+") as blastfile:
        for line in blastfile:
            if line[0] == '#':
                numHits = getNumHits(line, numHits)
            else:
                newQuery = Query(line, genome, numHits)
                newQuery.__setValues__()
                newQuery.__findBYearFamilies__(bYearFamilies)
                if newQuery.query not in genomeDB:
                    genomeDB[newQuery.query] = []
                genomeDB[newQuery.query].append(newQuery)
    return genomeDB


# TODO: remove
def queriesToJSON(filepath, flankseq, queriesWorkingSet):
    '''
    Converts dictionary to JSON object and writes it to a file.
    Appears in:
        GetPrimers.py
        ParseSangerResults.py
        SequencesFromBlast.py
        VerifyPrimers.py
    '''
    pass


#    filename = f"{filepath}/JSONfiles/WorkingSet_{flankseq}.json"
#   jsonObject = {}
#  for q in queriesWorkingSet:
#     if queriesWorkingSet[q].query not in jsonObject:
#        jsonObject[queriesWorkingSet[q].query] = []
#   jsonObject[queriesWorkingSet[q].query] = dict(queriesWorkingSet[q])
# newJSONobject = json.dumps(jsonObject, indent=4)

#    with open(filename, "w") as outfile:
#       outfile.write(newJSONobject)

#  workingSetCSVFormat(f"{filepath}/CSVfiles/WorkingSet_{flankseq}.csv", queriesWorkingSet)
# workingSetGFFFormat(f"{filepath}/GFFfiles/WorkingSet_{flankseq}.gff", queriesWorkingSet)

# return


def readPrimer3Output(filename, task, queriesWorkingSet):
    '''
    Reads Primer3 Output files, creates Primer3 objects, and appends data to Query objects in queriesWorkingSet.
    Appears in:
        GetPrimers.py
        VerifyPrimers.py
    '''
    with open(filename, "r") as p3file:
        outputDict = {}
        for line in p3file:
            # A single = indicates the end of a P3 response. Next line will be a new response.
            if line[0] != "=":
                lineItems = line.split("=")  # key:value pairs for OUTPUTLABEL:data
                outputDict[lineItems[0].strip()] = lineItems[1].strip()
            elif line[0] == "=":
                queriesWorkingSet.__updateQueryWithPrimer3Data__(outputDict, task)
            else:
                print("ERROR")


# TODO: not yet in use
def readSangerBlastFiles():
    sangerDict = {}
    sangerDir = "DIGIToutput/SangerSequences"
    for subdir, dirs, files in os.walk(sangerDir):
        for oneFile in files:
            filepath = os.path.join(subdir, oneFile)
            if filepath.find(".tab") != -1:
                f = filepath.split("/")[-1]
                genome = f.split("_")[0]
                with open(filepath, "r") as tabFile:
                    for line in tabFile:
                        if line.find("# Query:") != -1:
                            metadata = line.split("/")[-1]
                        if line[0] != "#":
                            query = line.split("\t")[0]
                            if query not in sangerDict:
                                sangerDict[query] = SangerSeq(query, genome)
                            sangerDict[query].__addBlastLine__(genome, line)
    print(sangerDict.keys())


def runFilterfasta(queriesWorkingSet):
    '''
    Appears in:
        SequencesFromBlast.py
    '''
    coorA, coorB, coorW = queriesWorkingSet.__getCoordinatesForFilterFasta__()

    seqDataA = filterFasta("DIGITfiles/Genomes/A188v1/Zm-A188-REFERENCE-KSU-1.0.fa", coorA)
    seqDataB = filterFasta("DIGITfiles/Genomes/B73v5/Zm-B73-REFERENCE-NAM-5.0.fa", coorB)
    seqDataW = filterFasta("DIGITfiles/Genomes/W22v2/Zm-W22-REFERENCE-NRGENE-2.0.fa", coorW)

    queriesWorkingSet.__assignFilterFastaData__(seqDataA)
    queriesWorkingSet.__assignFilterFastaData__(seqDataB)
    queriesWorkingSet.__assignFilterFastaData__(seqDataW)


def runPrimer3(inputFile, flankseq, outputFilePath):
    '''
    Call Primer3 and submit an input file.
    Appears in:
        GetPrimers.py
        SequencesFromBlast.py
    '''
    with open(outputFilePath, "w") as outfile:
        subprocess.run(["primer3_core", inputFile], stdout=outfile)


def sameQuery(queryS, queryOG):
    '''
    Given a Query object from the sanger set and a query object from the original set, check if the metadata of the
    Query is the same for genome, database, or sequence start position. Print if anythin is mismatched and return False.
    Return True if the Query objects match.
    Appears in:
        Utils.py
    '''
    if (queryS.genome != queryOG.genome) or (queryS.chromosome != queryOG.chromosome) or (
            queryS.sStart != queryOG.sStart):
        return False
    return True


def setBestForGenome(allQueries):
    '''
    Create dictionary identifying the best BLAST hit for each query, within each genome.
    Appears in:
        ParseSangerResults.py
        SequencesFromBlast.py
    '''
    listQueries = []
    for genome in allQueries:
        for allele in allQueries[genome]:
            if allele not in listQueries:
                listQueries.append(allele)
            getBestQuery(genome, allele, allQueries)

    return listQueries


def setPrimerFoundNearFlankSeq(primerDict, queriesWorkingSet):
    '''
    Appears in:
        GetPrimerIncidenceRate.py
    '''
    for genome in primerDict:
        for query in primerDict[genome]:
            for i in range(0, len(primerDict[genome][query])):
                # wildtypeCoordinates = queriesWorkingSet.__getQueryWTCoordinates__(query)
                # queryObject = queriesWorkingSet[query]
                primerDict[genome][query][i].primerFoundNearFlankSeq = checkCoordinates(
                    primerDict[genome][query][i], queriesWorkingSet.workingSet[query])


def sortSangerQueries(originalData, bestSangerQueries):
    '''
    Appears in:
        ParseSangerResults.py
    '''
    queryData, badQueryData, queryDataInTwoSeqs = [], [], []

    for query in bestSangerQueries.workingSet:
        queryOG = query.split("_")[0]
        if queryOG in originalData.workingSet:
            if sameQuery(bestSangerQueries.workingSet[query], originalData.workingSet[queryOG]):
                queryData.append(query)
            else:
                badQueryData.append(query)

        # If the query does not appear in originalData, check if it appears as part of the twoSeq datasets
        else:
            queryA = queryOG + "a"
            queryB = queryOG + "b"

            if (queryA in originalData.workingSet) or (queryB in originalData.workingSet):
                queryDataInTwoSeqs.append(query)

            # If it still doesn't appear, its not there
            else:
                print(f"{query} does not exist in original dataset")

    return queryData, badQueryData, queryDataInTwoSeqs


def splitSangerQueriesAllIntoMatchingAndNonMatchingSets(sangerQueriesAll):
    '''
    Appears in:
        ParseSangerResults.py
    '''
    matchingSeqs, nonMatchingSeqs = {}, {}

    for genome in sangerQueriesAll:
        if genome not in matchingSeqs:
            matchingSeqs[genome] = {}
        if genome not in nonMatchingSeqs:
            nonMatchingSeqs[genome] = {}
        for sangerSeq in sangerQueriesAll[genome]:
            if sangerSeq.find("_ImperfectMatchDsGfp") == -1:
                matchingSeqs[genome][sangerSeq] = sangerQueriesAll[genome][sangerSeq]
            else:
                nonMatchingSeqs[genome][sangerSeq] = sangerQueriesAll[genome][sangerSeq]
    return matchingSeqs, nonMatchingSeqs


def workingQuerySelection(genome, allele, allQueries):
    '''
    Set ID for best_query --> indicates that this will belong to working set
    Appears in:
        Utils.py
    '''
    for i in range(0, len(allQueries[genome][allele])):
        if allQueries[genome][allele][i].bestAlleleForGenome == True:
            allQueries[genome][allele][i].bestHitForAllele = True
            return


# TODO: not in use
def workingSetCSVFormat(filename, qWorkSet):
    '''
    Print queriesWorkingSet in csv format.
    Appears in:
        Utils.py
    '''
    pass
    # csvName = filename.split(".")[0] + ".csv"
    # with open(csvName, "w") as csvFile:
    #    csvFile.write(f"Query,Genome,BYearFamily,Chromosome,PercentIdentity,AlignmentLength,Mismatches,GapOpens,qStar" +
    #                 f"t,qEnd,sStart,sEnd,eValue,BitScore,NumHitsForGenome,PercentDifferenceToNextHit,Strand,qStartS" +
    #                f"tatus,PrimerNameLeft,PrimerSequenceLeft,PrimerLeftProductSize,PrimerPenaltyLeft,PrimerNameRig" +
    #               f"ht,PrimerSequenceRight,PrimerRightProductSize,PrimerPenaltyRight,PrimerPairProductSize,Primer" +
    #              f"PairPenalty\n")

    # for q in qWorkSet:
    #   csvFile.write(f"{qWorkSet[q].query},{qWorkSet[q].genome},{qWorkSet[q].bYearFamilies}," +
    #                f"{qWorkSet[q].chromosome},{qWorkSet[q].perIdentity},{qWorkSet[q].alignmentLength}," +
    #               f"{qWorkSet[q].mismatches},{qWorkSet[q].gapOpens},{qWorkSet[q].qStart},{qWorkSet[q].qEnd}," +
    #              f"{qWorkSet[q].sStart},{qWorkSet[q].sEnd},{qWorkSet[q].eValue},{qWorkSet[q].bitScore}," +
    #             f"{qWorkSet[q].numHits},{qWorkSet[q].percentDiff},{qWorkSet[q].strand}," +
    #            f"{qWorkSet[q].qStartStatus},{qWorkSet[q].primerNameLeft},{qWorkSet[q].primerSequenceLeft}," +
    #           f"{qWorkSet[q].primerLeftProductSize},{qWorkSet[q].primerPenaltyLeft}," +
    #          f"{qWorkSet[q].primerNameRight},{qWorkSet[q].primerSequenceRight}," +
    #         f"{qWorkSet[q].primerRightProductSize},{qWorkSet[q].primerPenaltyRight}," +
    #        f"{qWorkSet[q].primerPairProductSize},{qWorkSet[q].primerPairPenalty},\n")


# TODO: not in use
def workingSetGFFFormat(filename, queriesWorkingSet):
    pass


#    gffName = filename.split(".")[0] + ".gff"
#   newGffFile = GffFile("")
#  newGffFile.__formatFromQueriesWorkingSet__(queriesWorkingSet)
# newGffFile.__writeToGffFile__(gffName)

# TODO: this is not used, do I need it at all?
def writePrimerNumHitsToCSV(dirname, flankseq, numHitsDir):
    '''
    Appears in:
    '''
    with open(dirname + flankseq + "_PrimerOccuranceInGenome.csv", "w") as csvfile:
        csvfile.write(f"PrimerName,NumHits,Genome,\n")
        for genome in numHitsDir:
            for primer in numHitsDir[genome]:
                csvfile.write(f"{primer},{numHitsDir[genome][primer]},{genome},\n")


def writeToBestQueriesFile(listQueries, flankseq, allQueries, filename):
    '''
    Write to file the best bit scores per genome per query, for manual review if needed.
    Appears in:
        ParseSangerResults.py
        SequencesFromBlast.py
    '''
    with open(filename, "w+") as newfile:
        newfile.write(
            f"Query,A188_BitScore,A188_NumHits,A188_QStartStatus,B73_BitScore,B73_NumHits,B73_QStartStatus," +
            f"W22_BitScore,W22_NumHits,W22_QStartStatus,BestGenomes\n")

        for allele in listQueries:
            bestBitScores = {"A188v1": 0, "B73v5": 0, "W22v2": 0}
            newfile.write(allele + "," +
                          makeStrBestInGenome("A188v1", allele, bestBitScores, allQueries) +
                          makeStrBestInGenome("B73v5", allele, bestBitScores, allQueries) +
                          makeStrBestInGenome("W22v2", allele, bestBitScores, allQueries) +
                          bestGenomesWrite(bestBitScores) + "\n")
            pickGenome(allele, bestBitScores, allQueries)
    return


def writePrimerIncidenceDataToFile(primerDict, flankseq):
    '''
    Appears in:
        GetPrimerIncidenceRate.py
    '''
    with open(
            f"DIGIToutput/FlankingSequences/{flankseq}/PrimerData/{flankseq}_PrimerIncidenceRate.csv",
            "w") as outfile:
        outfile.write("PrimerName,NumberHitsInGenome,Genome,PrimerFoundNearFlankingSequence")
        for genome in primerDict:
            for query in primerDict[genome]:
                for i in range(0, len(primerDict[genome][query])):
                    outfile.write(
                        f"{primerDict[genome][query][i].primerName},{primerDict[genome][query][i].numHits}," +
                        f"{primerDict[genome][query][i].genome}," +
                        f"{primerDict[genome][query][i].primerFoundNearFlankSeq}\n")

                ### TODO: rewrite this it looks ridiculous


def writeSangerComparisonFile(filename, dataList, originalData, sangerData):
    '''
    Appears in:
        Utils
    '''
    with open(filename, "w+") as f:
        f.write(
            f"Query,RefGenomeSanger,ChromosomeSanger,RefGenomeOriginal,ChromosomeOriginal,SeqStartSanger,SeqEndSa" +
            f"nger,SeqStartOriginal,SeqEndOriginal,EValueSanger,EValueOriginal,BitScoreSanger,BitScoreOriginal,Nu" +
            f"mHitsSanger,NumHitsOriginal,StrandSanger,StrandOriginal,QStartStatusSanger,QStartStatusOriginal\n")
        for query in dataList:
            queryOG = query.split("_")[0]
            f.write(
                f"{sangerData.workingSet[query].query},{sangerData.workingSet[query].genome},{sangerData.workingSet[query].chromosome}," +
                f"{originalData.workingSet[queryOG].genome},{originalData.workingSet[queryOG].chromosome},{sangerData.workingSet[query].sStart}," +
                f"{sangerData.workingSet[query].sEnd},{originalData.workingSet[queryOG].sStart},{originalData.workingSet[queryOG].sEnd}," +
                f"{sangerData.workingSet[query].eValue},{originalData.workingSet[queryOG].eValue},{sangerData.workingSet[query].bitScore}," +
                f"{originalData.workingSet[queryOG].bitScore},{sangerData.workingSet[query].numHits},{originalData.workingSet[queryOG].numHits}," +
                f"{sangerData.workingSet[query].strand},{originalData.workingSet[queryOG].strand},{sangerData.workingSet[query].qStartStatus}," +
                f"{originalData.workingSet[queryOG].qStartStatus}\n")


### TODO: should probably redo this as well
def writeSangerComparisonFileTwoSeqs(filename, dataList, originalData, sangerData):
    '''
    Appears in:
        ParseSangerResults.py
    '''
    with open(filename, "w+") as f:
        f.write(
            f"Query,RefGenomeSanger,ChromosomeSanger,RefGenomeOriginalA,ChromosomeOriginalA,RefGenomeOriginalB,Ch" +
            f"romosomeOriginalB,SeqStartSanger,SeqEndSanger,SeqStartOriginalA,SeqEndOriginalA,SeqStartOriginalB,S" +
            f"eqEndOriginalB,EValueSanger,EValueOriginalA,EValueOriginalB,BitScoreSanger,BitScoreOriginalA,BitSco" +
            f"reOriginalB,NumHitsSanger,NumHitsOriginalA,NumHitsOriginalB,StrandSanger,StrandOriginalA,StrandOrig" +
            f"inalB,QStartStatusSanger,QStartStatusOriginalA,QStartStatusOriginalB\n")

        for query in dataList:
            queryOGA = query.split("_")[0] + "a"
            queryOGB = query.split("_")[0] + "b"
            aData = thisIsDumb(queryOGA, originalData)
            bData = thisIsDumb(queryOGB, originalData)
            f.write(
                f"{sangerData.workingSet[query].query},{sangerData.workingSet[query].genome},{sangerData.workingSet[query].chromosome}," +
                f"{aData['genome']},{aData['chromosome']},{bData['genome']},{bData['chromosome']}," +
                f"{sangerData.workingSet[query].sStart},{sangerData.workingSet[query].sEnd},{aData['sStart']},{aData['sEnd']}," +
                f"{bData['sStart']},{bData['sEnd']},{sangerData.workingSet[query].eValue},{aData['eValue']}," +
                f"{bData['eValue']},{sangerData.workingSet[query].bitScore},{aData['bitScore']},{bData['bitScore']}," +
                f"{sangerData.workingSet[query].numHits},{aData['numHits']},{bData['numHits']},{sangerData.workingSet[query].strand}," +
                f"{aData['strand']},{bData['strand']},{sangerData.workingSet[query].qStartStatus},{aData['qStartStatus']}," +
                f"{bData['qStartStatus']}\n")


def writeToFile(dataList, ordernum, genomeDetail, seqType, originalData, bestSangers):
    '''
    Appears in:
        ParseSangerResults.py
    '''
    if len(dataList) > 0:
        writeSangerComparisonFile(
            f"DIGIToutput/SangerSequences/{ordernum}/CSVfiles/{genomeDetail}_{ordernum}_" +
            f"{seqType}.csv", dataList, originalData, bestSangers)


