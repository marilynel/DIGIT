
from difflib import SequenceMatcher
import json
import os
import re
import subprocess

from FilterFasta import *
from GFF import *
from Primer import Primer
from Query import Query
from QueriesWorkingSet import QueriesWorkingSet
from SangerSeq import SangerSeq


def getAorB(query, originalData):
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
    return {"genome": None, "chromosome": None, "sStart": None, "sEnd": None, "eValue": None, "bitScore": None, "numHits":
        None, "strand": None, "qStartStatus": None}


# TODO: not yet in use
def addToTestedPrimersList(testedPrimers):
    with open("PutOrderedPrimersHere/SucessfulPrimersFromSanger", "r") as testedFile:
        for primer in testedPrimers:
            testedFile.write(primer)


# TODO: is there a way to put these next two scripts together into one file?
def callPrimerBlastScript(pathToFastaInput, pathToDesitinationDir):
    '''
    Calls script to blast primers against the A188, B73, and W22 genomes. Specialized script for short sequences.
    Appears in:
        VerifyPrimers.py
    '''
    subprocess.run \
        (["sh", f"./DIGITfiles/RunBlastPrimers.sh", pathToFastaInput, pathToDesitinationDir])

    print \
        (f"SGE may take a while to run. To check the status of your jobs, enter 'qstat' in the command line. Running" +
          f" time may vary wildly.\n")


def callBlastScript(pathToFastaInput, pathToDesitinationDir, fastaContents):
    '''
    Calls script to blast the initial flanking sequences against the A188, B73, and W22 genomes.
    Appears in:
        BlastSangerOutput.py
    '''
    subprocess.run \
        (["sh", f"./DIGITfiles/RunBlastInitial.sh", pathToFastaInput, pathToDesitinationDir, fastaContents])

    print \
        (f"SGE may take a while to run. To check the status of your jobs, enter 'qstat' in the command line. Running" +
          f" time may vary wildly.\n")


# TODO: move this to Query class?
def checkCoordinates(primer, query):
    '''
    Appears in:
        Utils.py
    '''
    if matchChrAndGenome(primer, query):
        if query.strand == 1:
            if (primer.sStart >= query.wildtypeCoordinates[0]) and \
                    (primer.sEnd <= query.wildtypeCoordinates[1]):
                return True
        elif query.strand == -1:
            if (primer.sStart <= query.wildtypeCoordinates[0]) and (
                    primer.sEnd >= query.wildtypeCoordinates[1]):
                return True
    return False


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


def getQueryName(line, qName):
    blastData = line.split(' ')
    if blastData[1] == "Query:":
        return blastData[2].strip()
    return qName


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
                outputDict.clear()
            else:
                print("ERROR")


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

