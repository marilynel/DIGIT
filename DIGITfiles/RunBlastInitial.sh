#!/bin/bash

# Blast Script for CGRB
# blast a fasta file of sequences against a database
#
# Script accepts 2 command line arguments:
#     $1 is the database
#     $2 is the name of the flanking sequence fasta file
#
#     example usage:
#     ./run_blast.sh A188v1 flanking_sequence_fasta_files/FlankingSequences_SingleGenomicSequences.fasta


#DIR="$(dirname "${1}")" ; FILE="$(basename "${1}")"
#fafile=${FILE%.*}

fafile=${1%.*}
now="$(date +'%s')"

mkdir -p DIGITfiles/BlastOutput/$fafile

SGE_Batch -c 'blastn -db DIGITfiles/Genomes/A188v1/A188v1 -outfmt 7 -query PutFlankingSequenceFilesHere/'$1' -num_threads 8 -max_target_seqs 10 -out DIGITfiles/BlastOutput/'$fafile'/A188v1_vs_'$fafile'.tab' -q bpp -P 8 -r sge.blastn_A188_$now
SGE_Batch -c 'blastn -db DIGITfiles/Genomes/B73v5/B73v5 -outfmt 7 -query PutFlankingSequenceFilesHere/'$1' -num_threads 8 -max_target_seqs 10 -out DIGITfiles/BlastOutput/'$fafile'/B73v5_vs_'$fafile'.tab' -q bpp -P 8 -r sge.blastn_B73_$now
SGE_Batch -c 'blastn -db DIGITfiles/Genomes/W22v2/W22v2 -outfmt 7 -query PutFlankingSequenceFilesHere/'$1' -num_threads 8 -max_target_seqs 10 -out DIGITfiles/BlastOutput/'$fafile'/W22v2_vs_'$fafile'.tab' -q bpp -P 8 -r sge.blastn_W22_$now


#cd data/part1/input
#ln -s ~/nfs7/BPP/Fowler_Lab/marilyn_folder/insertion_location_id_program/data/part0/output/$fafile
#cd ../../..


#echo "Results will be in data/part0/output/$fafile"
#echo "If applicable, error information will be in sge.blastn_<database>_$now"
#fi

