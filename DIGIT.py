import os
import subprocess
import time


def okGo(dirname, fileSubstrings):
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
    listDirs = os.listdir(dirname)
    print(
        f"\nSelect a directory to work with. Make sure the folder you need is in the {dirname} folder.\n")

    listDirs.sort()
    for i in range(0, len(listDirs)):
        print(f"({i + 1}) {listDirs[i]}")

    print()
    sel = input()
    if int(sel) < 1 or int(sel) > len(listDirs) + 1:
        print(f"\nPlease don't be weird, just select a real directory. Exiting.\n")
        exit()

    return listDirs[int(sel) - 1]


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

    flankseq = dirSelect("PutFlankingSequenceFilesHere/")

    try:
        now = time.time()
        print(flankseq)
        subprocess.run(
            [
                "sh",
                f"./DIGITfiles/RunBlastInitial.sh",
                f"{flankseq}"
            ]
        )
        print(
            f"SGE may take a while to run. To check the status of your jobs, enter 'qstat' in the command line. Runni" +
            f"ng time may vary wildly.\n"
        )
        print(
            f"Return later to continue to continue the primer making process with building the predicted insertion se" +
            f"quences and running them against Primer3.\n"
        )
        print(
            f"Output files for Blast results will be located in the DIGITfiles/BlastOutput directory.\n")
        print(f"Goodbye.\n")
    except:
        print("Well that didn't work. Exiting.")
        exit()


def runSequencesFromBlast():
    flankseq = dirSelect("DIGITfiles/BlastOutput")

    if okGo("DIGITfiles/BlastOutput/" + flankseq, ["A188", "B73", "W22"]):
        print(f"\nFinding genomic sequences for alleles in {flankseq}.")
        now = time.time()
        subprocess.run(
            [
                "SGE_Batch",
                "-c",
                f"python3 DIGITfiles/SequencesFromBlast.py {flankseq}",
                "-q",
                "bpp",
                "-P",
                "8",
                "-r",
                f"sge.runSequencesFromBlast_{flankseq}_{now}"
            ]
        )

        print(
            f"\nReturn later to continue to continue the primer making process with building the predicted insertion " +
            f"sequences. SGE may take a while to run. To check the status of your jobs, enter 'qstat' in the command " +
            f"line. Running time may vary wildly.\n"
        )
        exit()


def runGetPrimersAndVerify():
    # Select working directory
    flankseq = dirSelect("DIGIToutput")

    # make sure P3 outputfiles exiist in that dir
    if okGo("DIGIToutput/" + flankseq + "/Primer3Files/Output", ["Primer3Output"]):
        print(
            f"Verifying DsGFP insertion sequences for alleles in {flankseq} using wildtype sequences and Primer3 Ouput."
        )
        now = time.time()

        # run subprocess as sge batch job
        subprocess.run(  # ["python3", "DIGITfiles/GetPrimers.py", flankseq])
            [
                "SGE_Batch",
                "-c",
                f"python3 DIGITfiles/GetPrimers.py {flankseq}",
                "-q",
                "bpp",
                "-P",
                "8",
                "-r",
                f"sge.runGetPrimersAndVerify_{flankseq}_{now}"
            ]
        )

        print(
            f"Return later to continue to continue the primer making process with building the predicted insertion se" +
            f"quences. SGE may take a while to run. To check the status of your jobs, enter 'qstat' in the command li" +
            f"ne. Running time may vary wildly.\n"
        )


def main():
    print(
        f"Welcome to DIGIT, the premiere tool for predicting DsGFP insertion sequences in maize and building primers " +
        f"for those sequences. Please select from the following menu options:\n\n"
    )
    print(
        f"(1) Blast a collection of flanking sequences against maize genomes A188v1, B73v5, and W22v2")
    print(f"(2) Parse Blast results and find primers")
    print(f"(3) Parse primer data and run verification")
    print(f"(4) Check verification for any issues and produce primer dataset")
    print(f"(5) Remove extraneous sge files from directory")
    print(f"(6) Use 'qstat' to check job status\n")

    sel = input()
    if sel == "1":
        runBlast()
    elif sel == "2":
        runSequencesFromBlast()
    elif sel == "3":
        runGetPrimersAndVerify()
    elif sel == "4":
        print("TODO: check primer verification")
    elif sel == "5":
        subprocess.run(
            [
                "sh",
                f"./DIGITfiles/CleanUpDirectory.sh",
            ]
        )
    elif sel == "6":
        print("qstat")
        subprocess.run(
            ["qstat"]
        )
    else:
        print("Bye")


if __name__ == "__main__":
    main()