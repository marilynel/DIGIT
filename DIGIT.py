import os
import subprocess
import time
import logging
import datetime

from DIGITfiles.Utils import *

log = "log.log"


def okGo(dirname, fileSubstrings):
    '''
    Function okGo() ensures that the correct and necessary files exist in a directory in order to call another function
    or subprocess. Returns boolean value. No connection to the band.
    params:
    -   dirname         string, directory to be searched
    -   fileSubstrings  list of strings,files necessary to proceed
    '''
    ok = False
    listDirs = os.listdir(dirname)
    for fileSubstring in fileSubstrings:
        for dir in listDirs:
            if dir.find(fileSubstring) != -1:
                ok = True
                break
        else:
            print(
                f"\nYou do not have the correct files available in {dirname}/ to complete this task. Check to make su" +
                f"re you have completed the previous steps in order to proceed. See README for more details.\n"
            )
            exit()
    return ok


def dirSelect(dirname):
    '''
    User selects directory or folder within param <dirname> (string) to continue program with. Presented as numbered
    options in alphanumeric order. Exits if user selects an invalid int or noninteger value. Returns name of selected
    folder or directory as a string.
    '''
    listDirs = os.listdir(dirname)
    print(
        f"\nSelect a directory to work with. Make sure the folder you need is in the {dirname} folder.\n")

    listDirs.sort()
    for i in range(0, len(listDirs)):
        print(f"({i + 1}) {listDirs[i]}")

    print()
    sel = input()
    try:
        n = True
        while n:
            if listDirs[int(sel) - 1].find("BlastOutput") != -1:
                print("Invalid Directory. Pick another.")
                sel = input()
            elif int(sel) >= 1 and int(sel) <= len(listDirs):
                return listDirs[int(sel) - 1]

    except:
        print(f"\nPlease don't be weird, just select a real directory. Exiting.\n")
        exit()


def runBlast():
    listDirs = os.listdir("PutFlankingSequenceFilesHere/")
    if not listDirs:
        print(
            f"\nThere are no fasta files of flanking sequences available to Blast. Upload the flanking sequence file " +
            f"to the subdirectory PutFlankingSequenceFilesHere/ then restart this program.\n"
        )
        print(
            f"For assistance in using SFTP to upload files to CQLS please go to: https://www.hostinger.com/tutorial/h" +
            f"ow-to-use-sftp-to-safely-transfer-files/\n"
        )
        print(f"Exiting.\n")
        exit()

    fastaFile = dirSelect("PutFlankingSequenceFilesHere/")
    flankseq = fastaFile.split(".")[0]

    try:
        initialDirs = [
            f"DIGIToutput/FlankingSequences/{flankseq}",
            f"DIGIToutput/FlankingSequences/{flankseq}/BlastOutput",
            f"DIGIToutput/FlankingSequences/{flankseq}/BlastOutput/FlankingSequenceBlast"
            f"DIGIToutput/FlankingSequences/{flankseq}/QueryData"
        ]

        for dir in initialDirs:
            if not os.path.exists(dir):
                os.makedirs(dir)

        if len(os.listdir(f"DIGIToutput/FlankingSequences/{flankseq}/BlastOutput")) == 3:
            print("d")
            print(
                f"\nYou may already have blast results for this input file. Would you like to continue anyways? y or n")
            resp = input()
            if resp != 'y':
                print("e")
                exit()
        callBlastScript(
            f"PutFlankingSequenceFilesHere/{fastaFile}",
            f"DIGIToutput/FlankingSequences/{flankseq}/BlastOutput/FlankingSequenceBlast",
            flankseq
        )

        print(
            f"Return later to continue to continue the primer making process with building the predicted insertion se" +
            f"quences and running them against Primer3.\n"
        )
        print(
            f"Output files for Blast results will be located in the DIGITfiles/BlastOutput directory.\n")
        print(f"Goodbye.\n")
        logging.info(f"\tRunning callBlastScript() with {flankseq}")

    except:
        print("Process failed. Goodbye.")
        logging.info(f"\tFailed to run callBlastScript() {flankseq}")
        exit()


def runSequencesFromBlast():
    flankseq = dirSelect("DIGIToutput/FlankingSequences")

    # okGo(dirname, fileSubstrings)
    if okGo(f"DIGIToutput/FlankingSequences/{flankseq}", ["BlastOutput"]):
        print(f"\nFinding genomic sequences for alleles in {flankseq}.")
        ok = True
        while ok:
            print("Select from following:\n")
            print("(1) Make working set from genomes A188, B73, and W22")
            print("(2) Make working set from B73 only\n")
            b73only = input()
            if b73only == "1" or b73only == "2":
                ok = False

        now = time.time()
        subprocess.run(
            [
                "SGE_Batch",
                "-c",
                f"python3 DIGITfiles/SequencesFromBlast.py {flankseq} {b73only}",
                "-q",
                "bpp",
                "-P",
                "8",
                "-r",
                f"sge.runSequencesFromBlast_{flankseq}_{now}"
            ]
        )

        print(
            f"\nReturn later to continue to continue the primer making process with parsing the Primer3 output result" +
            f"s and submitting verification input to Primer3. SGE may take a while to run. To check the status of you" +
            f"r jobs, enter 'qstat' in the command line. Running time may vary wildly.\n"
        )

        print(
            f"\nResults will be saved in DIGIToutput/FlankingSequences/{flankseq}/QueryData"
        )
        logging.info(f"\tRunning SequencesFromBlast.py with {flankseq}")


def runGetPrimersAndVerify():
    # Select working directory
    flankseq = dirSelect("DIGIToutput/FlankingSequences")

    # make sure P3 outputfiles exiist in that dir
    # okGo(dirname, fileSubstrings)
    if okGo(f"DIGIToutput/FlankingSequences/{flankseq}/QueryData/", ["Primer3"]):
        print(
            f"Verifying DsGFP insertion sequences for alleles in {flankseq} using wildtype sequences and Primer3 Ouput."
        )

        # run subprocess as sge batch job
        subprocess.run(
            [
                "SGE_Batch",
                "-c",
                f"python3 DIGITfiles/GetPrimers.py {flankseq}",
                "-q",
                "bpp",
                "-P",
                "8",
                "-r",
                f"sge.runGetPrimersAndVerify_{flankseq}"
            ]
        )

        print(
            f"\nReturn later to continue to continue the primer making process with check verification output for iss" +
            f"ues, producing the primer dataset, and Blasting the primer dataset against all three genomes. SGE may t" +
            f"ake a while to run. To check the status of your jobs, enter 'qstat' in the command line. Running time m" +
            f"ay vary wildly.\n"
        )

        print(
            f"\nResults will be saved in DIGIToutput/FlankingSequences/{flankseq}/QueryData/Primer3/VerifyingPrimers/" +
            f"and in updated WorkingSet file in DIGIToutput/FlankingSequences/{flankseq}/QueryData/JSONfiles/"
        )
        logging.info(f"\tRunning GetPrimers.py with {flankseq}")


def runVerifyPrimers():
    flankseq = dirSelect("DIGIToutput/FlankingSequences")
    # okGo(dirname, fileSubstrings)
    if okGo(
            f"DIGIToutput/FlankingSequences/{flankseq}/QueryData/Primer3/VerifyingPrimers",
            ["Primer3VerificationOutput"]):
        print(
            f"Parsing Primer3 wildtype verification results."
        )
        subprocess.run(
            [
                "python3",
                "DIGITfiles/VerifyPrimers.py",
                f"{flankseq}"
            ]
        )
        print(
            f"Output will be in DIGIToutput/FlankingSequences/{flankseq}/ directory. For a consise list of the primer" +
            f"s with relevant information, go to DIGIToutput/FlankingSequences/{flankseq}/PrimerData/)."
        )
        logging.info(f"\tRunning VerifyPrimers.py with {flankseq}")


def findNumPrimersInDB():
    flankseq = dirSelect("DIGIToutput/FlankingSequences")
    # okGo(dirname, fileSubstrings)
    if okGo(f"DIGIToutput/FlankingSequences/{flankseq}/BlastOutput/PrimerBlast/", ["PrimerBlast"]):
        print(
            f"Primers and stuff."
        )
        subprocess.run(
            [
                "python3",
                "DIGITfiles/GetPrimerIncidenceRate.py",
                f"{flankseq}"
            ]
        )
        logging.info(f"\tRunning GetPrimerIncidenceRate.py with {flankseq}")


def showSangerMenu():
    print(
        f"Welcome to DIGIT's Sanger Sequence Menu, designed to guide you through processing and blasting the output o" +
        f"f sanger sequencing PCR results, as well as compare those results to your initial dataset.\n"
    )
    print(
        f"Please select from the following menu options:\n"
    )
    logging.info(f"\tSanger Menu requested")

    print(f"(1) Blast Sanger Sequencing results against maize genomes A188v1, B73v5, and W22v2")
    print(f"(2) Parse Sanger Sequence Blast results and compare with original sequence data\n")

    sel = input()
    print()
    logging.info(f"\tuser selected {sel}")

    if sel == "1":
        blastSangerSeqs()
    elif sel == "2":
        parseSangerResults()


def parseSangerResults():
    sangerseq = dirSelect(f"DIGIToutput/SangerSequences")
    # okGo(dirname, fileSubstrings)
    if okGo(f"DIGIToutput/SangerSequences/{sangerseq}/BlastOutput", ["A188", "B73", "W22"]):
        print(f"Parsing Sanger Blast results:")
        subprocess.run(
            [
                "python3",
                "DIGITfiles/ParseSangerResults.py",
                f"{sangerseq}"
            ]
        )
        logging.info(f"\tRunning ParseSangerResults.py with {sangerseq}")


def blastSangerSeqs():
    sangerseq = dirSelect("PutSangerOutputFilesHere/")
    # okGo(dirname, fileSubstrings)
    if okGo("PutSangerOutputFilesHere/" + sangerseq, [".seq"]):
        print("process and blast sanger data")
        now = time.time()
        subprocess.run(
            [
                "python3",
                f"DIGITfiles/BlastSangerOutput.py",
                f"PutSangerOutputFilesHere/{sangerseq}"
            ]
        )
        logging.info(f"\tRunning BlastSangerOutput.py with {sangerseq}")


def main():
    time = datetime.datetime.now()
    logging.basicConfig(filename=f"DIGIToutput/log", level=logging.INFO)
    logging.info(f"\t{time}")
    print(
        f"Welcome to DIGIT, the premiere tool for predicting DsGFP insertion sequences in maize and building primers " +
        f"for those sequences.\n"
    )
    print(f"Please select from the following menu options:\n")
    print(
        f"(1) Blast a collection of flanking sequences against maize genomes A188v1, B73v5, and W22v2")
    print(f"(2) Parse Blast results and find primers")
    print(f"(3) Parse primer data and run verification")
    print(
        f"(4) Check verification for any issues, produce primer dataset, and Blast primer dataset against all three" +
        f" genomes")
    print(f"(5) Check frequency of primers in genomes")
    print(f"(6) Go to Sanger Sequence Parsing menu")
    print(f"(7) Remove extraneous sge files from directory")
    print(f"(8) Use 'qstat' to check job status\n")

    sel = input()

    logging.info(f"\tuser selected {sel}")

    print()
    if sel == "1":
        runBlast()
    elif sel == "2":
        runSequencesFromBlast()
    elif sel == "3":
        runGetPrimersAndVerify()
    elif sel == "4":
        runVerifyPrimers()
    elif sel == "5":
        findNumPrimersInDB()
    elif sel == "6":
        showSangerMenu()
    elif sel == "7":
        subprocess.run(
            [
                "sh",
                f"./DIGITfiles/CleanUpDirectory.sh",
            ]
        )
    else:
        print("Bye")
    logging.info("")


if __name__ == "__main__":
    main()