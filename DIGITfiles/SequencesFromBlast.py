import subprocess
from Query import Query
from FilterFasta import filterFasta
import os
import sys
import json

# all queries is a dict with structure:
#   {
#       databaseA:
#           {
#               query1:[QueryObj, QueryObj, QueryObj],
#               query2:[QueryObj]
#           }
#       databaseB:
#           {
#               query1:[QueryObj, QueryObj, QueryObj],
#               query2:[QueryObj]
#           }
#   }   etc.


allQueries = {
    "A188v1": {},
    "B73v5": {},
    "W22v2": {}
}

queriesWorkingSet = {}


def findBlastOutputFiles(dirname):
    '''
    This function goes to the specified directory, locates relevant (.tab)
    files, and sends the file names to function process_blast_output() for
    parsing.
    '''
    for subdir, dirs, files in os.walk(dirname):
        # TODO: really? i can make this prettier
        for oneFile in files:
            filename = os.path.join(subdir, oneFile)
            if filename.find(".tab") != -1:
                processBlastOutput(filename)

    return


def processBlastOutput(filename):
    '''
    Parse blast output files line by line, gathering data that will be used to
    create Query objects, which will be added to the all_queries dictionary.
    Functions get_num_hits() and get_genome() are called to help parse hash (#)
    lines.
    '''
    genome = ""
    num_hits = -1
    with open(filename, "r+") as blastfile:
        for line in blastfile:
            if line[0] == '#':
                # hash lines may hold data that will be necessary for Query
                # object initialization. genome and num_hits will remain
                # unchanged until after all of the hits for that query are read
                genome = getGenome(line, genome)
                num_hits = getNumHits(line, num_hits)
            else:
                newQuery = Query(line, genome, num_hits)
                newQuery.__setValues__()
                if newQuery.query not in allQueries[genome]:
                    # Add that query to all_queries[genome] as a key, value
                    # being an empty list
                    allQueries[genome][newQuery.query] = []
                # Place Query object in all_queries object
                allQueries[genome][newQuery.query].append(newQuery)
    return


def getGenome(line, genome):
    '''
    Parse a line for string indicating the genome name. If line does not have
    that string, return genome as original value.
    '''
    blastData = line.split(' ')
    if blastData[1] == "Database:":
        return blastData[2].split('/')[-1].strip()
    return genome


def getNumHits(line, numHits):
    '''
    Parse a line for string indicating the number of hits for that query in a
    genome. If line does not have that string, return num_hits as original
    value.
    '''
    blastData = line.split(' ')
    if blastData[2] == "hits":
        return int(blastData[1])
    return numHits


def setBestForGenome():
    '''
    Create dictionary identifying the best BLAST hit for each query, within
    each genome.
    '''
    listQueries = []

    for genome in allQueries:
        for allele in allQueries[genome]:
            if allele not in listQueries:
                listQueries.append(allele)
            getBestQuery(genome, allele)
    writeToBestQueriesFile(listQueries)

    return listQueries


def makeStrBestInGenome(genome, allele, best_bit_score):
    '''
    Find the specific hit that was specified as the best for that query ID and
    that genome (Query.best_for_genome == True) and return a string reporting
    data for the best_queries_by_genome.csv outfile.
    '''
    if allele in allQueries[genome]:
        for i in range(0, len(allQueries[genome][allele])):
            if allQueries[genome][allele][i].bestAlleleForGenome == True:
                best_bit_score[genome] = allQueries[genome][allele][i].bitScore
                return allQueries[genome][allele][i].__makeBestGenomeString__()
    return "none,none,none,"


def bestGenomesWrite(bestBitScore):
    '''
    Find the best of the best hits, return strinigied list for writing to file.
    '''
    genomeList = []
    bestScore = max(bestBitScore.values())
    for genome in bestBitScore:
        if bestBitScore[genome] == bestScore:
            genomeList.append(genome)
        else:
            bestBitScore[genome] = 0
    return str(genomeList)


def workingQuerySelection(genome, allele):
    '''
    set ID for best_query --> indicates that this will belong to working set
    '''
    for i in range(0, len(allQueries[genome][allele])):
        if allQueries[genome][allele][i].bestAlleleForGenome == True:
            allQueries[genome][allele][i].bestHitForAllele = True
            return


def pickGenome(allele, bestBitScore):
    if bestBitScore["B73v5"] != 0:
        workingQuerySelection("B73v5", allele)
        return
    elif bestBitScore["W22v2"] != 0:
        workingQuerySelection("W22v2", allele)
        return
    elif bestBitScore["A188v1"] != 0:
        workingQuerySelection("A188v1", allele)
    else:
        print("you really messed something up")
        exit()


def writeToBestQueriesFile(listQueries):
    with open("BestQueriesByGenome.csv", "w+") as newfile:
        newfile.write(
            "query,A188_bit_score,A188_num_hits,A188_qstart_status,B73_bit_sc" +
            "ore,B73_num_hits,B73_qstart_status,W22_bit_score,W22_num_hits,W2" +
            "2_qstart_status,best_genomes\n"
        )
        for allele in listQueries:
            bestBitScore = {}
            newfile.write(
                allele + "," +
                makeStrBestInGenome("A188v1", allele, bestBitScore) +
                makeStrBestInGenome("B73v5", allele, bestBitScore) +
                makeStrBestInGenome("W22v2", allele, bestBitScore) +
                bestGenomesWrite(bestBitScore) +
                "\n"
            )
            pickGenome(allele, bestBitScore)
    return


def getBestQuery(genome, query):
    '''
    Sort hits to identify best hit for a query. Find the percent difference
    between the best and second best hits.
    Returns:
    best_query              Query object
    '''
    bestQuery = allQueries[genome][query][0]
    secondBestQuery = None
    for i in range(1, len(allQueries[genome][query])):
        if allQueries[genome][query][i].bitScore > bestQuery.bitScore:
            secondBestQuery = bestQuery
            bestQuery = allQueries[genome][query][i]
        elif allQueries[genome][query][i].bitScore < bestQuery.bitScore:
            if not secondBestQuery or allQueries[genome][query][
                i].bitScore > secondBestQuery.bitScore:
                secondBestQuery = allQueries[genome][query][i]
    # TODO: come back and make this work
    # if second_best_query:
    #    best_query.diff = abs(
    #        ((best_query.bit_score - second_best_query.bit_score) / (best_query.bit_score + second_best_query.bit_score)) / 2) * 100

    for i in range(0, len(allQueries[genome][query])):
        if allQueries[genome][query][i] == bestQuery:
            allQueries[genome][query][i].bestAlleleForGenome = True

    return


def alleleInGenomeComp(genome, allele):
    if allele in allQueries[genome]:
        for i in range(0, len(allQueries[genome][allele])):
            if allQueries[genome][allele][i].bestAlleleForGenome == True:
                return allQueries[genome][allele][i].bitScore
    return -1


def setBestQuery(listQueries):
    for query in listQueries:
        # todo: the hell is this???????????
        bestBitScores = [
            alleleInGenomeComp("A188v1", query),
            alleleInGenomeComp("B73v5", query),
            alleleInGenomeComp("W22v2", query),
        ]


def queriesToJSON(filename):
    '''
    Converts dictionary to JSON object and writes it to a file.
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


def workingQueriesToJSON(filename):
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


def runFilterfasta(filelist):
    for gen in allQueries:
        for q in allQueries[gen]:
            for i in range(0, len(allQueries[gen][q])):
                if allQueries[gen][q][i].bestHitForAllele:
                    queriesWorkingSet[allQueries[gen][q][i].query] = allQueries[gen][q][i]

                # if allQueries[gen][q][i].query in filelist:
                #    print(f"{allQueries[gen][q][i].query} already exists in filterfasta folder")
                '''
                else:
                    if allQueries[gen][q][i].bestHitForAllele:
                        subprocess.run(
                            [
                                "SGE_Batch",
                                "-c",
                                f"./DIGITfiles/filterfasta.sh {allQueries[gen][q][i].query} {allQueries[gen][q][i].chromosome} {allQueries[gen][q][i].wildtypeCoordinates[0]} {allQueries[gen][q][i].wildtypeCoordinates[1]} {allQueries[gen][q][i].upperCoordinates[0]} {allQueries[gen][q][i].upperCoordinates[1]} {allQueries[gen][q][i].lowerCoordinates[0]} {allQueries[gen][q][i].lowerCoordinates[1]} {gen[:-2]}",
                                "-q",
                                "bpp",
                                "-P",
                                "8",
                                "-r",
                                f"sge.{allQueries[gen][q][i].query}"
                            ]
                        )
                '''
    # coordinates = {
    #    "test1" : ["chr1", [500, 600]],
    #    "test2" : ["chr1", [550, 650]]
    # }
    coorA, coorB, coorW = {}, {}, {}

    for q in queriesWorkingSet:
        # print(queriesWorkingSet[q].genome)
        if queriesWorkingSet[q].genome.strip() == "A188v1":
            coorA[queriesWorkingSet[q].query + "_wt"] = [queriesWorkingSet[q].chromosome,
                                                         queriesWorkingSet[q].wildtypeCoordinates]
            coorA[queriesWorkingSet[q].query + "_up"] = [queriesWorkingSet[q].chromosome,
                                                         queriesWorkingSet[q].upperCoordinates]
            coorA[queriesWorkingSet[q].query + "_lo"] = [queriesWorkingSet[q].chromosome,
                                                         queriesWorkingSet[q].lowerCoordinates]
        elif queriesWorkingSet[q].genome.strip() == "B73v5":
            coorB[queriesWorkingSet[q].query + "_wt"] = [queriesWorkingSet[q].chromosome,
                                                         queriesWorkingSet[q].wildtypeCoordinates]
            coorB[queriesWorkingSet[q].query + "_up"] = [queriesWorkingSet[q].chromosome,
                                                         queriesWorkingSet[q].upperCoordinates]
            coorB[queriesWorkingSet[q].query + "_lo"] = [queriesWorkingSet[q].chromosome,
                                                         queriesWorkingSet[q].lowerCoordinates]
        elif queriesWorkingSet[q].genome.strip() == "W22v2":
            coorW[queriesWorkingSet[q].query + "_wt"] = [queriesWorkingSet[q].chromosome,
                                                         queriesWorkingSet[q].wildtypeCoordinates]
            coorW[queriesWorkingSet[q].query + "_up"] = [queriesWorkingSet[q].chromosome,
                                                         queriesWorkingSet[q].upperCoordinates]
            coorW[queriesWorkingSet[q].query + "_lo"] = [queriesWorkingSet[q].chromosome,
                                                         queriesWorkingSet[q].lowerCoordinates]
        else:
            print(f"I don't even know what went wrong with {queriesWorkingSet[q].query}")

    seqDataA = filterFasta("DIGITfiles/Genomes/Zm-A188-REFERENCE-KSU-1.0.fa", coorA)
    seqDataB = filterFasta("DIGITfiles/Genomes/Zm-B73-REFERENCE-NAM-5.0.fa", coorB)
    seqDataW = filterFasta("DIGITfiles/Genomes/Zm-W22-REFERENCE-NRGENE-2.0.fa", coorW)
    print(seqDataA)


def main():
    findBlastOutputFiles(f"DIGITfiles/BlastOutput/{sys.argv[1]}")
    list_of_queries = setBestForGenome()
    setBestQuery(list_of_queries)
    queriesToJSON(f"DIGITfiles/AllBlastData_{sys.argv[1]}.json")

    filelist = os.listdir("DIGITfiles/FilterfastaFiles")

    filelist = [filelist[i][:-6] for i in range(len(filelist))]

    runFilterfasta(filelist)

    workingQueriesToJSON(f"DIGITfiles/WorkingQuerySet_{sys.argv[1]}.json")


if __name__ == '__main__':
    main()
