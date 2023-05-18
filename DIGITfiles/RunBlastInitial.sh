#!/bin/bash

# TODO: redo this
# Blast Script for CGRB
# blast a fasta file of sequences against a database
#
# Script accepts 2 command line arguments:
#     $1 is the database
#     $2 is the name of the flanking sequence fasta file
#
#     example usage:
#     ./run_blast.sh A188v1 flanking_sequence_fasta_files/FlankingSequences_SingleGenomicSequences.fasta




now="$(date +'%s')"

SGE_Batch -c 'blastn -db DIGITfiles/Genomes/A188v1/A188v1 -outfmt 7 -query '$1' -num_threads 8 -max_target_seqs 10 -out '$2'/A188v1_vs_'$3'.tab' -q bpp -P 8 -r sge.blastn_A188_$now
SGE_Batch -c 'blastn -db DIGITfiles/Genomes/B73v5/B73v5 -outfmt 7 -query '$1' -num_threads 8 -max_target_seqs 10 -out '$2'/B73v5_vs_'$3'.tab' -q bpp -P 8 -r sge.blastn_B73_$now
SGE_Batch -c 'blastn -db DIGITfiles/Genomes/W22v2/W22v2 -outfmt 7 -query '$1' -num_threads 8 -max_target_seqs 10 -out '$2'/W22v2_vs_'$3'.tab' -q bpp -P 8 -r sge.blastn_W22_$now





