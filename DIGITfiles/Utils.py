import json
import subprocess
import time
import os

from Query import Query
from Primer3Object import Primer3Object


def callBlastScript(pathToFastaInput, pathToDesitinationDir, fastaContents):
    '''
    Calls script to blast the initial flanking sequences against the A188, B73, and W22 genomes.
    Appears in:
        DIGIT.py
        BlastSangerOutput.py
    '''
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


def callPrimerBlastScript(pathToFastaInput, pathToDesitinationDir):
    '''
    Calls script to blast primers against the A188, B73, and W22 genomes. Specialized script for short sequences.
    Appears in:
        VerifyPrimers.py
    '''
    now = time.time()
    subprocess.run(
        [
            "sh",
            f"./DIGITfiles/RunBlastPrimers.sh",
            pathToFastaInput,
            pathToDesitinationDir  # ,
            # fastaContents
        ]
    )
    print(
        f"SGE may take a while to run. To check the status of your jobs, enter 'qstat' in the command line. Runni" +
        f"ng time may vary wildly.\n"
    )


def readPrimer3Output(filename, task, queriesWorkingSet):
    '''
    Reads Primer3 Output files, creates Primer3 objects, and appends data to Query objects in queriesWorkingSet.
    Appears in:
        GetPrimers.py
        VerifyPrimers.py
    '''
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
    '''
    Creates queriesWorkingSet from JSON file. Initializes Query objects and populates the dict queriesWorkingSet.
    Appears in:
        GetPrimers.py
        VerifyPrimers.py
    '''
    dataFromJSON = getDataFromJSON(filename)
    for allele in dataFromJSON:
        if allele not in queriesWorkingSet:
            queriesWorkingSet[allele] = Query(None, None, None)
            queriesWorkingSet[allele].__QueryFromJSON__(dataFromJSON[allele])


def getDataFromJSON(filename):
    '''
    Parse a JSON file and return a JSON object (key:value format)
    Appears in:
        Utils.py
    '''
    jsonFile = open(filename)
    dataFromJSON = json.load(jsonFile)
    jsonFile.close()

    return dataFromJSON


def writeP3InputFile(inputFilepath, queriesWorkingSet, flankseq, task):
    '''
    Initializes a Primer3Object for each query in queriesWorkingSet, and creates the input file that will be used with
    Primer3. Return value is the path to the Primer3 input file (str).
    Appears in:
        GetPrimers.py
    '''
    with open(inputFilepath, "w") as p3file:
        for allele in queriesWorkingSet:
            p3obj = Primer3Object(task)
            p3obj.__initValsFromQueryObj__(queriesWorkingSet[allele])

            if task == "generic":
                p3obj.__buildInputStrGeneric__()
            else:
                p3obj.__buildInputStrValidate__()
            p3file.write(p3obj.inputStr)

    if task == "generic":
        return f"DIGIToutput/FlankingSequences/{flankseq}/QueryData/Primer3/FindingPrimers/Primer3Input_{flankseq}.txt"
    else:
        return f"DIGIToutput/FlankingSequences/{flankseq}/QueryData/Primer3/VerifyingPrimers/Primer3VerificationInput_{flankseq}.txt"


def queriesToJSON(filepath, flankseq, queriesWorkingSet):
    '''
    Converts dictionary to JSON object and writes it to a file.
    Appears in:
        GetPrimers.py
        ParseSangerResults.py
        SequencesFromBlast.py
        VerifyPrimers.py
    '''
    # filename = f"DIGIToutput/FlankingSequences/{flankseq}/QueryData/JSONfiles/WorkingSet_{
    # flankseq}.json"
    filename = f"{filepath}/JSONfiles/WorkingSet_{flankseq}.json"
    jsonObject = {}
    for q in queriesWorkingSet:
        if queriesWorkingSet[q].query not in jsonObject:
            jsonObject[queriesWorkingSet[q].query] = []
        jsonObject[queriesWorkingSet[q].query] = dict(queriesWorkingSet[q])
    newJSONobject = json.dumps(jsonObject, indent=4)

    with open(filename, "w") as outfile:
        outfile.write(newJSONobject)

    workingSetCSVFormat(
        f"{filepath}/CSVfiles/WorkingSet_{flankseq}.csv",
        queriesWorkingSet
    )
    return


def workingSetCSVFormat(filename, queriesWorkingSet):
    '''
    Print queriesWorkingSet in csv format.
    Appears in:
        Utils.py
    '''
    csvName = filename.split(".")[0] + ".csv"
    with open(csvName, "w") as csvFile:
        csvFile.write(
            f"Query,Genome,Chromosome,PercentIdentity,AlignmentLength,Mismatches,GapOpens,qStart,qEnd,sStart,sEnd,eVa" +
            f"lue,BitScore,NumHitsForGenome,PercentDifferenceToNextHit,Strand,qStartStatus,PrimerNameLeft,PrimerSeque" +
            f"nceLeft,PrimerLeftProductSize,PrimerPenaltyLeft,PrimerNameRight,PrimerSequenceRight,PrimerRightProductS" +
            f"ize,PrimerPenaltyRight,PrimerPairProductSize,PrimerPairPenalty\n"
        )

        for q in queriesWorkingSet:
            csvFile.write(
                f"{queriesWorkingSet[q].query},{queriesWorkingSet[q].genome},{queriesWorkingSet[q].chromosome}," +
                f"{queriesWorkingSet[q].perIdentity},{queriesWorkingSet[q].alignmentLength}," +
                f"{queriesWorkingSet[q].mismatches},{queriesWorkingSet[q].gapOpens},{queriesWorkingSet[q].qStart}," +
                f"{queriesWorkingSet[q].qEnd},{queriesWorkingSet[q].sStart},{queriesWorkingSet[q].sEnd}," +
                f"{queriesWorkingSet[q].eValue},{queriesWorkingSet[q].bitScore},{queriesWorkingSet[q].numHits}," +
                f"{queriesWorkingSet[q].percentDiff},{queriesWorkingSet[q].strand}," +
                f"{queriesWorkingSet[q].qStartStatus},{queriesWorkingSet[q].primerNameLeft}," +
                f"{queriesWorkingSet[q].primerSequenceLeft},{queriesWorkingSet[q].primerLeftProductSize}," +
                f"{queriesWorkingSet[q].primerPenaltyLeft},{queriesWorkingSet[q].primerNameRight}," +
                f"{queriesWorkingSet[q].primerSequenceRight},{queriesWorkingSet[q].primerRightProductSize}," +
                f"{queriesWorkingSet[q].primerPenaltyRight},{queriesWorkingSet[q].primerPairProductSize}," +
                f"{queriesWorkingSet[q].primerPairPenalty},\n"
            )


def checkIfPrimerHasID():
    '''
    Create and return a list of primers that have already been developed.
    TODO: do I need this in here????

    '''
    listPreviousPrimers = []
    with open("DIGITfiles/compListPreviousPrimers", "r") as prevFile:
        for line in prevFile:
            listPreviousPrimers.append(line.strip())
    return listPreviousPrimers


def runPrimer3(inputFile, flankseq, outputFilePath):
    '''
    Call Primer3 and submit an input file.
    Appears in:
        GetPrimers.py
        SequencesFromBlast.py
    '''
    with open(outputFilePath, "w") as outfile:
        subprocess.run(
            [
                "primer3_core",
                inputFile
            ],
            stdout=outfile
        )


def composePrimerListOutput(filename, queriesWorkingSet):
    '''
    Composes a list of primers from queriesWorkingSet, organizes them, and prints them to csv.
    Appears in:
        VerifyPrimers.py
    '''
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


def processBlastOutput(filename, genome):
    '''
    Parse blast output files line by line, gathering data that will be used to create Query objects, which will be added
    to the allQueries dataset. Functions get_num_hits() and get_genome() are called to help parse hash (#) lines.
    Appears in:
        Utils.py
    '''
    genomeDB = {}
    num_hits = -1
    with open(filename, "r+") as blastfile:
        for line in blastfile:
            if line[0] == '#':
                num_hits = getNumHits(line, num_hits)
            else:
                newQuery = Query(line, genome, num_hits)
                newQuery.__setValues__()
                if newQuery.query not in genomeDB:
                    genomeDB[newQuery.query] = []
                genomeDB[newQuery.query].append(newQuery)
    return genomeDB


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


def getBestQuery(genome, query, allQueries):
    '''
    Sort hits to identify best hit for a query. Find the percent difference between the best and second best hits.
    Appears in:
        Utils.py
    '''
    bestQuery = allQueries[genome][query][0]  # Query
    secondBestQuery = None  # Query
    for i in range(1, len(allQueries[genome][
                              query])):  # iterate through Query objects in allQueries[genome][allele]
        ### NEW: changed > to >= make sure it works before deleting this comment!!!! 5/22 ###
        ### NOTE: this does change the output from before. yikes. still valid but the results are inconsistant...going to put the = on the "less than" comparison
        ### note again: that fixed it.
        if allQueries[genome][query][
            i].bitScore > bestQuery.bitScore:  # if the new bit score is better than bestQuery's bit score
            secondBestQuery = bestQuery  # redefine secondBestQuery as Query in bestQuery
            bestQuery = allQueries[genome][query][i]  # redefine bestQuery as the current Query

        elif allQueries[genome][query][i].bitScore <= bestQuery.bitScore:
            if not secondBestQuery or allQueries[genome][query][
                i].bitScore >= secondBestQuery.bitScore:  # if either (there is no secondBestQuery) OR (new query's bit score is better than the secondBestWuery bit score)
                secondBestQuery = allQueries[genome][query][
                    i]  # redefine second best query as the new query
    perDiff = -1
    if secondBestQuery:  # if there is a second-best query
        perDiff = abs(
            (
                    ((bestQuery.bitScore - secondBestQuery.bitScore) / (
                                bestQuery.bitScore + secondBestQuery.bitScore)) / 2
            ) * 100
        )

    for i in range(0, len(allQueries[genome][query])):
        if allQueries[genome][query][i] == bestQuery:
            allQueries[genome][query][i].bestAlleleForGenome = True
            allQueries[genome][query][i].percentDiff = perDiff

    return


def allQueriesToJSON(filename, allQueries):
    '''
    Converts dictionary to JSON object and writes it to a file.
    Appears in:
        GetPrimerIncidenceRate.py
        ParseSangerResults.py
        SequencesFromBlast.py
    '''
    jsonObject = {
        "A188v1": {},
        "B73v5": {},
        "W22v2": {}
    }

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


def buildWorkingSetFromAllQueries(allQueries, queriesWorkingSet):
    '''
    Iterates through entire dataset, and creates a new smaller working dataset with just the "best" hit for each query.
    Appears in:
        ParseSangerResults.py
        SequencesFromBlast.py
    '''
    for gen in allQueries:
        for q in allQueries[gen]:
            for i in range(0, len(allQueries[gen][q])):
                if allQueries[gen][q][i].bestHitForAllele:
                    queriesWorkingSet[allQueries[gen][q][i].query] = allQueries[gen][q][i]


def pickGenome(allele, bestBitScore, allQueries):
    '''
    Selects best genome.
    Appears in:
        Utils.py
    '''
    if bestBitScore["B73v5"] != 0:
        workingQuerySelection("B73v5", allele, allQueries)
        return
    elif bestBitScore["W22v2"] != 0:
        workingQuerySelection("W22v2", allele, allQueries)
        return
    elif bestBitScore["A188v1"] != 0:
        workingQuerySelection("A188v1", allele, allQueries)
    else:
        print("you really messed something up")


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


def writeToBestQueriesFile(listQueries, flankseq, allQueries, filename):
    '''
    Write to file the best bit scores per genome per query, for manual review if needed.
    Appears in:
        ParseSangerResults.py
        SequencesFromBlast.py
    '''
    with open(filename, "w+") as newfile:
        newfile.write(
            f"query,A188_bit_score,A188_num_hits,A188_qstart_status,B73_bit_score,B73_num_hits,B73_qstart_status,W22_" +
            f"bit_score,W22_num_hits,W22_qstart_status,best_genomes\n"
        )
        for allele in listQueries:
            bestBitScores = {
                "A188v1": 0,
                "B73v5": 0,
                "W22v2": 0
            }
            newfile.write(
                allele + "," +
                makeStrBestInGenome("A188v1", allele, bestBitScores, allQueries) +
                makeStrBestInGenome("B73v5", allele, bestBitScores, allQueries) +
                makeStrBestInGenome("W22v2", allele, bestBitScores, allQueries) +
                bestGenomesWrite(bestBitScores) +
                "\n"
            )
            pickGenome(allele, bestBitScores, allQueries)
    return


def bestGenomesWrite(bestBitScore):
    '''
    Find the best of the best hits, return strinigied list for writing to file.
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


def makeDirectories(necessaryDirs):
    '''
    Appears in:
        BlastSangerOutput.py
        ParseSangerResults.py
        SequencesFromBlast.py
    '''

    for dir in necessaryDirs:
        if not os.path.exists(dir):
            os.makedirs(dir)

