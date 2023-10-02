#!/bin/bash -l

now="$(date +'%s')"

#SGE_Batch -c "blastn -db ../../../DIGITfiles/Genomes/B73v5/B73v5 -query '$1' -outfmt 7 -max_target_seqs 50000 -evalue 1 -gapopen 10 -gapextend 2 -reward 1 -penalty -1 -word_size 5 -out '$2'" -q bpp -P 8 -r sge.blastn_W22_$now
#SAME AS FIRST SGE_Batch -c "blastn -db ../../../DIGITfiles/Genomes/W22v2/W22v2 -query '$1' -outfmt 7 -max_target_seqs 50000 -evalue 1 -gapopen 10 -gapextend 2 -penalty -1 -word_size 5 -out '$2'" -q bpp -P 8 -r sge.blastn_W22_$now
#SGE_Batch -c "blastn -db ../../../DIGITfiles/Genomes/W22v2/W22v2 -query '$1' -outfmt 7 -max_target_seqs 50000 -evalue 1 -gapopen 10 -gapextend 2 -word_size 5 -out '$2'" -q bpp -P 8 -r sge.blastn_W22_$now
# FIFTH IS FOURTH IN JSON< NOT THAT HELPFUL SGE_Batch -c "blastn -db ../../../DIGITfiles/Genomes/W22v2/W22v2 -query '$1' -outfmt 12 -max_target_seqs 50000 -gapopen 10 -gapextend 2 -word_size 5 -out '$2'" -q bpp -P 8 -r sge.blastn_W22_$now
#SGE_Batch -c "blastn -db ../../../DIGITfiles/Genomes/W22v2/W22v2 -query '$1' -outfmt 7 -max_target_seqs 50000 -word_size 5 -out '$2'" -q bpp -P 8 -r sge.blastn_W22_$now
SGE_Batch -c "blastn -db DIGITfiles/Genomes/A188v1/A188v1 -query '$1' -outfmt 7 -max_target_seqs 50000 -word_size 5 -out '$2'/A188v1_PrimerBlast.tab" -q bpp -P 8 -r sge.blastn_A188primers_$now
SGE_Batch -c "blastn -db DIGITfiles/Genomes/B73v5/B73v5 -query '$1' -outfmt 7 -max_target_seqs 50000 -word_size 5 -out '$2'/B73v5_PrimerBlast.tab" -q bpp -P 8 -r sge.blastn_B73primers_$now
SGE_Batch -c "blastn -db DIGITfiles/Genomes/W22v2/W22v2 -query '$1' -outfmt 7 -max_target_seqs 50000 -word_size 5 -out '$2'/W22v2_PrimerBlast.tab" -q bpp -P 8 -r sge.blastn_W22primers_$now

# to get only one mismatch
# q start < 3?
# also filter results by mismatch

