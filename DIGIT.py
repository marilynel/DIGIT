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
    if os.path.exists(dirname):
        listDirs = os.listdir(dirname)
        for fileSubstring in fileSubstrings:
            for dir in listDirs:
                if dir.find(fileSubstring) != -1:
                    ok = True
                    break

    else:
        print("Files do not exist to support this action. Please try something else.")
        # main()
    return ok


def dirSelect(dirname):
    '''
    User selects directory or folder within param <dirname> (string) to continue program with. Presented as numbered
    options in alphanumeric order. Exits if user selects an invalid int or noninteger value. Returns name of selected
    folder or directory as a string.
    '''
    listDirs = os.listdir(dirname)
    print(
        f"\nSelect a directory or file to work with. Make sure the folder or file you need is in the {dirname} fold" +
        f"er. If you do not see the folder or file you need, type 'n' to exit.\n")

    listDirs.sort()
    for i in range(0, len(listDirs)):
        print(f"({i + 1}) {listDirs[i]}")

    print()
    sel = input()
    try:
        n = True
        while n:
            if int(sel) >= 1 and int(sel) <= len(listDirs):
                n = False
                return listDirs[int(sel) - 1]
            else:
                print("Please select a valid directory.")
                sel = input()
    except:
        print(f"\nPlease don't be weird, just select a real directory. Exiting.\n")
        logging.info(f"\tUser decided to be weird and tried to select {sel}")
        exit()


def filesExist(dirname, fileSubstrings):
    pass
    # if okGo(dirname, fileSubstrings):
    #    listDirs = os.listdir(dirname)
    #    print(f"\nThis step may already have been completed. The files in {dirname} include:\n")
    #    for item in listDirs:
    #        print(f" - {item}")
    #    print(f"\nWould you like to continue this process anyway or return to the main menu?")
    #    print(f"\n(0) Continue")
    #    print(f"(1) Return to main menu\n")

    #    sel = input()
    #    if int(sel) == 0:
    #        return
    #    elif int(sel) == 1:
    #        main()
    #    else:
    #        print("Why would you do this?")
    #    exit()


def runBlast():
    if okGo("PutFlankingSequenceFilesHere/", ["fasta"]):
        fastaFile = dirSelect("PutFlankingSequenceFilesHere/")
        flankseq = fastaFile.split(".")[0]

        initialDirs = [
            f"DIGIToutput/FlankingSequences/{flankseq}",
            f"DIGIToutput/FlankingSequences/{flankseq}/BlastOutput",
            f"DIGIToutput/FlankingSequences/{flankseq}/BlastOutput/FlankingSequenceBlast",
            f"DIGIToutput/FlankingSequences/{flankseq}/QueryData"
        ]

        makeDirectories(initialDirs)

        callBlastScript(
            f"PutFlankingSequenceFilesHere/{fastaFile}",
            f"DIGIToutput/FlankingSequences/{flankseq}/BlastOutput/FlankingSequenceBlast",
            flankseq
        )

        print(
            f"Return later to continue to continue the primer making process with building the predicted insertion " +
            f"sequences and running them against Primer3.\n")
        print(
            f"Output files for Blast results will be located in the DIGIToutput/FlankingSequences/{flankseq}/BlastO" +
            f"utput/FlankingSequenceBlast directory.\n")
        print(f"Goodbye.\n")
        logging.info(f"\tRunning callBlastScript() with {flankseq}")
    else:
        logging.info(f"\tFailed to run callBlastScript()")
        main()


def runSequencesFromBlast():
    flankseq = dirSelect("DIGIToutput/FlankingSequences")

    if okGo(f"DIGIToutput/FlankingSequences/{flankseq}/BlastOutput/FlankingSequenceBlast", ["tab"]):
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
        # filesExist(f"DIGIToutput/FlankingSequences/{flankseq}/QueryData/JSONfiles", ["json"])
        subprocess.run([
            "SGE_Batch",
            "-c",
            f"python3 DIGITfiles/SequencesFromBlast.py {flankseq} {b73only}",
            "-q",
            "bpp",
            "-P",
            "8",
            "-r",
            f"sge.runSequencesFromBlast_{flankseq}_{now}"
        ])

        print(
            f"\nReturn later to continue to continue the primer making process with parsing the Primer3 output resu" +
            f"lts and submitting verification input to Primer3. SGE may take a while to run. To check the status of" +
            f"your jobs, enter 'qstat' in the command line. Running time may vary wildly.\n")

        print(f"\nResults will be saved in DIGIToutput/FlankingSequences/{flankseq}/QueryData")
        print(f"\nParsed data for all blast results:\t/JSONfiles/AllBlastHits_{flankseq}.json")
        print(f"Working set of \"best\" queries:\t\t/JSONfiles/WorkingSet_{flankseq}.json")
        print(f"CSV version of working set:\t\t/CSVfiles/WorkingSet_{flankseq}.csv")
        print(f"GFF of working set:\t\t\t/GFFfiles/WorkingSet_{flankseq}.gff")
        print(f"Primer3 input file:\t\t\t/Primer3/FindingPrimers/Primer3Input_{flankseq}.txt")
        print(
            f"Primer3 output will appear as:\t\t/Primer3/FindingPrimers/Primer3Output_{flankseq}.txt")
        logging.info(f"\tRunning SequencesFromBlast.py with {flankseq}")


def runGetPrimersAndVerify():
    flankseq = dirSelect("DIGIToutput/FlankingSequences")

    if okGo(f"DIGIToutput/FlankingSequences/{flankseq}/QueryData/Primer3/FindingPrimers",
            ["Output"]):
        # filesExist(f"DIGIToutput/FlankingSequences/{
        # flankseq}/QueryData/Primer3/VerifyingPrimers", ["Output"])

        print(
            f"Verifying DsGFP insertion sequences for alleles in {flankseq} using wildtype sequences and Primer3 Ou" +
            f"put.")

        subprocess.run([
            "SGE_Batch",
            "-c",
            f"python3 DIGITfiles/GetPrimers.py {flankseq}",
            "-q",
            "bpp",
            "-P",
            "8",
            "-r",
            f"sge.runGetPrimersAndVerify_{flankseq}"
        ])

        print(
            f"\nReturn later to continue to continue the primer making process with check verification output for i" +
            f"ssues, producing the primer dataset, and Blasting the primer dataset against all three genomes. SGE m" +
            f"ay take a while to run. To check the status of your jobs, enter 'qstat' in the command line. Running " +
            f"time may vary wildly.\n")

        print(f"\nResults will be saved in DIGIToutput/FlankingSequences/{flankseq}/QueryData/")

        print(f"Updated working set of \"best\" queries:\t\t/JSONfiles/WorkingSet_{flankseq}.json")
        print(f"Updated CSV version of working set:\t\t/CSVfiles/WorkingSet_{flankseq}.csv")
        print(
            f"Primer3 verification input file:\t\t/Primer3/VerifyingPrimers/Primer3VerificationInput_{flankseq}.txt")
        print(
            f"Primer3 verification output will appear as:\t/Primer3/VerifyingPrimers/Primer3VerificationOutput_{flankseq}.txt")

        logging.info(f"\tRunning GetPrimers.py with {flankseq}")


def runVerifyPrimers():
    flankseq = dirSelect("DIGIToutput/FlankingSequences")
    if okGo(
            f"DIGIToutput/FlankingSequences/{flankseq}/QueryData/Primer3/VerifyingPrimers",
            ["Primer3VerificationOutput"]):
        filesExist(f"DIGIToutput/FlankingSequences/{flankseq}/QueryData/Primer3/VerifyingPrimers",
                   ["Output"])
        print(f"Parsing Primer3 wildtype verification results.")

        subprocess.run(["python3", "DIGITfiles/VerifyPrimers.py", f"{flankseq}"])
        print(f"Output will be in DIGIToutput/FlankingSequences/{flankseq}/ directory.")

        print(f"\nPrimer data:\t\t\t\t/PrimerData/PrimerResults_{flankseq}.csv")
        print(f"Failed primer data:\t\t\t/PrimerData/FailedPrimers_{flankseq}.csv")
        print(f"Blast output files:\t\t\t/BlastOutput/PrimerBlast/")
        print(
            f"Updated working set of \"best\" queries:\t/QueryData/JSONfiles/WorkingSet_{flankseq}.json")
        print(f"Updated CSV version of working set:\t/QueryData/CSVfiles/WorkingSet_{flankseq}.csv")

        logging.info(f"\tRunning VerifyPrimers.py with {flankseq}")


def findNumPrimersInDB():
    flankseq = dirSelect("DIGIToutput/FlankingSequences")
    # okGo(dirname, fileSubstrings)
    if okGo(f"DIGIToutput/FlankingSequences/{flankseq}/BlastOutput/PrimerBlast/", ["PrimerBlast"]):
        print(f"Primers and stuff.")
        filesExist(f"DIGIToutput/FlankingSequences/{flankseq}/PrimerData/", ["csv"])
        subprocess.run(["python3", "DIGITfiles/GetPrimerIncidenceRate.py", f"{flankseq}"])
        logging.info(f"\tRunning GetPrimerIncidenceRate.py with {flankseq}")

        print(f"Output will be in DIGIToutput/FlankingSequences/{flankseq}/ directory.")
        print(f"\nPrimer incidence rates per genome:\t\t/PrimerData/IncidenceRate.csv")
        print(
            f"Updated working set of \"best\" queries:\t\t/QueryData/JSONfiles/WorkingSet_{flankseq}.json")
        print(
            f"Updated CSV version of working set:\t\t/QueryData/CSVfiles/WorkingSet_{flankseq}.csv")
    else:
        print(
            f"\nFiles do not exist to support this action. Please go back and complete step 4 with Primer blast opt" +
            f"ion.")


def showSangerMenu():
    print(
        f"Welcome to DIGIT's Sanger Sequence Menu, designed to guide you through processing and blasting the output" +
        f"of sanger sequencing PCR results, as well as compare those results to your initial dataset.\n")
    print(f"Please select from the following menu options:\n")
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
    if okGo(f"DIGIToutput/SangerSequences/{sangerseq}/BlastOutput", ["A188", "B73", "W22"]):
        print(f"Parsing Sanger Blast results:")
        subprocess.run(["python3", "DIGITfiles/ParseSangerResults.py", f"{sangerseq}"])
        logging.info(f"\tRunning ParseSangerResults.py with {sangerseq}")


def blastSangerSeqs():
    sangerseq = dirSelect("PutSangerOutputFilesHere/")
    if okGo("PutSangerOutputFilesHere/" + sangerseq, [".seq"]):
        filesExist(f"DIGIToutput/SangerSequences/{sangerseq}/BlastOutput/", ["tab"])
        print("Process and blast sanger data")
        subprocess.run(["python3", f"DIGITfiles/BlastSangerOutput.py",
                        f"PutSangerOutputFilesHere/{sangerseq}"])
        logging.info(f"\tRunning BlastSangerOutput.py with {sangerseq}")


def main():
    time = datetime.datetime.now()
    logging.basicConfig(filename=f"DIGIToutput/log", level=logging.INFO)
    logging.info(f"\t{time}")
    print(
        f"Welcome to DIGIT, the premiere tool for predicting DsGFP insertion sequences in maize and building primer" +
        f"s for those sequences.\n")
    print(f"Please select from the following menu options:\n")
    print(
        f"(1) Blast a collection of flanking sequences against maize genomes A188v1, B73v5, and W22v2")
    print(f"(2) Parse Blast results and create and submit input files for Primer3")
    print(f"(3) Parse Primer3 output data and begin Primer3 verification")
    print(
        f"(4) Check Primer3 verification for any issues, produce primer dataset, and (optional) Blast primer datase" +
        f"t against all three genomes")
    print(
        f"(5) Check frequency of primers in genomes and make sure primer occurs within predicted sequence wildtype " +
        f"coordinates.")
    print(f"(6) Go to Sanger Sequence Parsing menu")
    print(f"(7) Look up data for specific allele, primer, or sanger sequence")
    print(f"(8) Remove extraneous sge files from directory")
    print(f"(9) Use 'qstat' to check job status\n")

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
        subprocess.run(["python3", "DIGITfiles/Lookup.py"])
    elif sel == "8":
        subprocess.run(["sh", f"./DIGITfiles/CleanUpDirectory.sh"])
    elif sel == "9":
        print(f"qstat\n")
        subprocess.run(["qstat"])
    else:
        print("Bye")
    logging.info("")


if __name__ == "__main__":
    main()