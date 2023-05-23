#!/bin/bash

# This script is called from within the DIGIT program to perform a Blast search of a selected fasta file against three
# maize reference genomes, A188, B73, and W22. Results will be saved to the specifed output folder. Sge directories will
# be created when these jobs are submitted to the queuing system. The contents of these directories are useful in the
# event of the Blast search not running as expected. Most of the time they will not be needed by the user and may be
# safely deleted with the menu option in the DIGIT main menu.

# Usage: ./RunBlastInistial.sh <arg1> <arg2> <arg3>
# Command line arguments:
#   arg1    $1          path to fasta input file
#   arg2    $2          path to output destination directory
#   arg3    $3          fasta file name/type (eg, FlankingSequences_SingleGenomicSequences; used to give unique
#                       identifying names to the outputfiles)

now="$(date +'%s')"

SGE_Batch -c 'blastn -db DIGITfiles/Genomes/A188v1/A188v1 -outfmt 7 -query '$1' -num_threads 8 -max_target_seqs 10 -out '$2'/A188v1_vs_'$3'.tab' -q bpp -P 8 -r sge.blastn_A188_$now
SGE_Batch -c 'blastn -db DIGITfiles/Genomes/B73v5/B73v5 -outfmt 7 -query '$1' -num_threads 8 -max_target_seqs 10 -out '$2'/B73v5_vs_'$3'.tab' -q bpp -P 8 -r sge.blastn_B73_$now
SGE_Batch -c 'blastn -db DIGITfiles/Genomes/W22v2/W22v2 -outfmt 7 -query '$1' -num_threads 8 -max_target_seqs 10 -out '$2'/W22v2_vs_'$3'.tab' -q bpp -P 8 -r sge.blastn_W22_$now
