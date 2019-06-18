#!/bin/bash
#PBS -lselect=1:ncpus=4:mem=20gb
#PBS -lwalltime=10:00:00

cd $1

shopt -s nullglob
numfiles=(*)
numfiles=${#numfiles[@]}

for ((i=1; i<=${numfiles}; i++)); do for filename in res-$i.txt; do if [ ! -f $filename ]; then echo "$filename does not exist"; exit ; fi; done; done 

echo 1. Checked that all files exist.

cat res-*.txt >> GWAS_$1_imp.txt

echo 2. Concatenated the result files. Number of SNPs: 
wc -l GWAS_$1_imp.txt

#mkdir res_files
#mv res-*.txt res_file

#echo 3. Tidied up

awk '{ if ($9 >= 0.5) print $0 }' GWAS_$1_imp.txt > GWAS_$1_imp_qual05.txt

echo 4. Filtered by quality with 0.5 threshold. Number of SNPs remaining:
wc -l GWAS_$1_imp_qual05.txt

mv ../GWAS_$1.txt GWAS_$1.txt
mv ../$1_imprange.txt  $1_imprange.txt
mv ../ssimp_$1_array.sh ssimp_$1_array.sh

#module load anaconda3/personal
#Rscript ../insert_originalGWAS_into_SSimp.R GWAS_$1.txt GWAS_$1_imp_qual05.txt

cd ..
mv $1 $1_processed

echo GWAS_$1_imp_qual05.txt
