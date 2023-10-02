'''
BlastSangerOutput.py creates a fasta file out of the .seq files returned from sanger sequencing and blasts those
sequences against the three reference genomes, A188, B73, and W22.

Written by: Marilyn Leary 2023

    Input:
    orderDir (sys.argv[1])  Path to directory containing sanger sequencing .seq files. Should be a subdirectory within
                            PutSangerOutputFilesHere/.

    Output:
    Three Blast output files with path and naming convention:
        DIGITfiles/SangerSequences/<orderNumber>/BlastOutputSangerSeqs/<referenceGenome_vs_orderNumber>.tab
'''

import os
import sys
import re
from difflib import SequenceMatcher

from Utils import *

# 34 BP tail at the end of the DsGFP insertion sequence
dsgfpEnd = "TTTTACCGACCGTTCCCGACCGTTTTCATCCCTA"


def prepFasta(orderDir):
    '''
    Create a fasta-formatted file from the .seq files contained in the given directory (orderDir). Returns path to
    resultant fasta file.
    '''
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
                    orderDir, filename[:-4], sequence[idx+len(dsgfpEnd):-22], order, notExact)

    filename = "PutSangerOutputFilesHere/" + order + "/" + order + ".fasta"
    with open(filename, "w") as outfile:
        outfile.write(processedSeqFilesContents)
    return filename


def checkConsecutiveNs(sequence):
    '''
    Look for 2 N's in a row in a sequence. N indicates an unknown base pair. Returns an index to the consectutive N's if
    found or nonsense value (-2) if not.
    '''
    previousBasePairs = sequence[0]
    for i in range(0, len(sequence)):
        if previousBasePairs == "N" and sequence[i] == "N":
            return i-1
        else:
            previousBasePairs = sequence[i]
    return -2


def makeFastaString(folder, metadata, sequence, setstring, notExact):
    '''
    Aggregates .seq data into a fasta-formatted string to be written to file. Returns that string.
    '''
    seqID = getSequenceID(metadata, notExact)
    fastaFormatStr = ">" + seqID + " " + metadata + "\n"
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


def getSequence(filename, dsgfpEnd):
    '''
    Read a .seq file and parse the data from it. Returns genomic sequence, the index where the last base pair of
    dsgfpEnd is found, and a boolean value where False indicates an exact match (not notExact) and True indicates a non-
    exact match (notExact).
    '''
    sequence = ""
    with open(filename, "r+") as sangerFile:
        for line in sangerFile:
            sequence += line.strip()
    if sequence.find(dsgfpEnd) != -1:
        return sequence, sequence.find(dsgfpEnd), False
    else:
        seqMatchObj = SequenceMatcher(lambda x: x=="ACGT", sequence, dsgfpEnd)
        matchObj = seqMatchObj.find_longest_match(0, len(sequence), 0, len(dsgfpEnd))
        return sequence, matchObj.a, True


def getSequenceID(metadata, notExact):
    '''
    Parse sequence ID (elsewhere known as query or allele) from the file name. Append string indicating its an imperfect
    match if relevant. Returns sequence ID.
    '''
    seqID = re.findall (r"R\d{1,4}[A-Z]*\d{1,4}", metadata)
    if notExact:
        return seqID[0] + "_ImperfectMatchDsGfp"
    return seqID[0]

#TODO: file str
def main():
    orderDir = sys.argv[1]
    filename = prepFasta(orderDir)
    order = orderDir.split("/")[-1]
    newDirs = [
        "DIGIToutput/SangerSequences",
        f"DIGIToutput/SangerSequences/{order}",
        f"DIGIToutput/SangerSequences/{order}/BlastOutput"
    ]
    makeDirectories(newDirs)
    # callBlastScript(pathToFastaInput, pathToDesitinationDir, fastaContents)
    callBlastScript(
        filename,
        f"DIGIToutput/SangerSequences/{order}/BlastOutput",
        order
    )


if __name__ == '__main__':
    main()