'''
1. make a fasta file out of the .seq files
2. blast the fasta against all 3 ref genomes
3. process blast output
4. not sure....compare to original results?
'''

import os
import sys
import re
from difflib import SequenceMatcher
import subprocess
import time

from Utils import *

dsgfpEnd = "TTTTACCGACCGTTCCCGACCGTTTTCATCCCTA"

def prepFasta(orderDir):
    order = orderDir.split("/")[-1]
    processedSeqFilesContents = ""
    for subdir, dirs, files in os.walk(sys.argv[1]):
        for oneFile in files:
            filename = os.path.join(subdir, oneFile)
            if filename.find(".seq") != -1:
                #print(f".seq file: {filename}")
                sequence, idx, notExact = getSequence(filename, dsgfpEnd)
                processedSeqFilesContents += writeToFasta(orderDir, filename[:-4], sequence[idx+len(dsgfpEnd):-22], order, notExact)
                #filename = writeToFasta(orderDir, filename[:-4], sequence[idx+len(dsgfpEnd):-22], order, notExact)

    #if not os.path.exists("DIGITfiles/BlastOutputSangerSeqs/" + order):
    #    os.makedirs("DIGITfiles/BlastOutputSangerSeqs/" + order)
    #filename = "DIGITfiles/BlastOutputSangerSeqs/" + order + "/" + order + ".fasta"
    filename = "PutSangerOutputFilesHere/" + order + "/" + order + ".fasta"
    with open(filename, "w") as outfile:
        outfile.write(processedSeqFilesContents)
    return filename

def checkConsecutiveNs(sequence):
    previousBasePairs = sequence[0]
    for i in range(0,len(sequence)):
        if previousBasePairs == "N" and sequence[i] == "N":
            return i-1
        else:
            previousBasePairs = sequence[i]
    return -2

# rename --> not writing to fasta just making a big str to write to fasta
def writeToFasta(folder, metadata, sequence, setstring, notExact):
    #with open(folder + "SangerSeq_" + setstring + ".fasta", "a+") as outfile:
    seqID = getSequenceID(metadata, notExact)
    #    outfile.write(">" + seqID + " " + metadata + "\n")
    fastaFormatStr = ">" + seqID + " " + metadata + "\n"
    endIdx = checkConsecutiveNs(sequence)
    if endIdx == -2:
    #        outfile.write(sequence + "\n")
        fastaFormatStr += sequence + "\n"
    elif endIdx == -1 or endIdx == 0:
    #        outfile.write("ERROR\n")
        fastaFormatStr += "ERROR\n"
    else:
    #        outfile.write(sequence[:endIdx] + "\n")
        fastaFormatStr += sequence[:endIdx] + "\n"
    #return(folder + "SangerSeq_" + setstring + ".fasta")
    return fastaFormatStr

def getSequence(filename, dsgfpEnd):
    sequence = ""
    with open(filename, "r+") as sangerFile:
        for line in sangerFile:
            sequence += line.strip()
    #idx = 0
    if sequence.find(dsgfpEnd) != -1:
        return sequence, sequence.find(dsgfpEnd), False
    else:
        seqMatchObj = SequenceMatcher(lambda x: x=="ACGT", sequence, dsgfpEnd)
        matchObj = seqMatchObj.find_longest_match(0, len(sequence), 0, len(dsgfpEnd))
        return sequence, matchObj.a, True


def getSequenceID(metadata, notExact):
    #print(metadata)
    seqID = re.findall (r"R\d{1,4}[A-Z]*\d{1,4}", metadata)
    #print(seqID)
    if notExact:
        return seqID[0] + "_ImperfectMatchDsGfp"
    return seqID[0]

#def runBlast(sangerFasta, order):
#    now = time.time()
    # Input:
    # $1 = Folder where fasta lives
    # $2 = name of fasta file s
    # $3 = output file location


#    print(f"arg 1: DIGITfiles/BlastOutputSangerSeqs/{order}")
##    print(f"arg 2: {sangerFasta}")
 #   print(f"arg 3: BlastOutputSangerSeqs")


    #subprocess.run(
    ##    [
      ##      "sh",
       ##     f"./DIGITfiles/RunBlastInitial.sh",
         #   f"DIGITfiles/BlastOutputSangerSeqs/{order}",
         #   f"{sangerFasta}",
         #   f"BlastOutputSangerSeqs"
        #]
   # )


def main():

    orderDir = sys.argv[1]
    #print(f"orderDir: {orderDir}")

    filename = prepFasta(orderDir)
    #print(f"filename: {filename}")
    #print(f"orderDir: {orderDir}")
    #print(f"order: {orderDir.split('/')[-1]}")
    order = orderDir.split("/")[-1]
    if not os.path.exists("DIGITfiles/BlastOutputSangerSeqs/" + order):
        os.makedirs("DIGITfiles/BlastOutputSangerSeqs/" + order)

    #runBlast(filename.split("/")[-1], orderDir.split("/")[-1])
    callBlastScript(
        filename,
        f"DIGITfiles/BlastOutputSangerSeqs/{order}",
        order
    )

if __name__ == '__main__':
    main()