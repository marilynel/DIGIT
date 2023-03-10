#!/bin/bash

# $1 = query
# $2 = chromosome
# $3 = wildtype start
# $4 = wildtype end
# $5 = upper start
# $6 = upper end
# $7 = lower start
# $8 = lower end
# $9 = genome

mkdir filterfasta_files

echo '>'$1'_'$2'_'$9'_wildtype' >> filterfasta_files/$1.fasta
filterfasta.pl --match $2 --start $3 --end $4 genomes/Zm-$9-* >> filterfasta_files/$1.fasta

echo '>'$1'_'$2'_'$9'_upper' >> filterfasta_files/$1.fasta
filterfasta.pl --match $2 --start $5 --end $6 genomes/Zm-$9-* >> filterfasta_files/$1.fasta

echo '>'$1'_'$2'_'$9'_lower' >> filterfasta_files/$1.fasta
filterfasta.pl --match $2 --start $7 --end $8 genomes/Zm-$9-* >> filterfasta_files/$1.fasta

#    SGE_Batch -c 'part1/filterfasta.pl --match '$chromosome' --start '$hg_start' --end '$hg_end' resources/part0/reference_genomes/'$db'/Zm* >> data/filterfasta/filterfasta_output/'$temp'/'$query'/'$hg_name'.fasta' -q bpp -P 8 -r sge.${temp}_hg_${query}_${now}

