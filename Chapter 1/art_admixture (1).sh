#!/usr/bin/env bash
#PBS -N admixture_art
#PBS -q 
#PBS -l ncpus=28
#PBS -l mem=60gb
#PBS -l walltime=48:00:00
#PBS -o /98_log_files/admixture_art.log

#Global variables
DATADIRECTORY=./population_genetic
DATAINPUT=$DATADIRECTORY/02_data/puce_final
DATAOUTPUT=$SCRATCH/population_genetic/admixture
ADMIXTURENV=". /env.sh"
PLINK="/plink/1.9/plink"
NCPU=28

mkdir -p $DATAOUTPUT

# first we will transforme the vcf file into the bed format require by admixture, select the good populations and filtrate again
cd $DATAINPUT

$PLINK --allow-extra-chr --vcf $DATAINPUT/chip_pinctadapt_data_cleaned_outgroup.vcf --make-bed --out $DATAINPUT/admixture_file_art

rm $SCRATCH/population_genetic/tmp*

# ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just exchange the first column by 0
awk '{$1="0";print $0}' admixture_file_art.bim > admixture_file_art.bim.tmp
mv admixture_file_art.bim.tmp admixture_file_art.bim


# Now perform the analyses from k=1 to k=10 with statistiques implement in the software
$ADMIXTURENV

cd $DATAOUTPUT
admixture --cv $DATAINPUT/admixture_file_art.bed 1 -j28 | tee log_admixture_art_1.log
admixture --cv $DATAINPUT/admixture_file_art.bed 2 -j28 | tee log_admixture_art_2.log
admixture --cv $DATAINPUT/admixture_file_art.bed 3 -j28 | tee log_admixture_art_3.log
admixture --cv $DATAINPUT/admixture_file_art.bed 4 -j28 | tee log_admixture_art_4.log
admixture --cv $DATAINPUT/admixture_file_art.bed 5 -j28 | tee log_admixture_art_5.log
admixture --cv $DATAINPUT/admixture_file_art.bed 6 -j28 | tee log_admixture_art_6.log
admixture --cv --seed=50 $DATAINPUT/admixture_file_art.bed 7 -j28 | tee log_admixture_art_7.log
admixture --cv --seed=50 $DATAINPUT/admixture_file_art.bed 8 -j28 | tee log_admixture_art_8.log
admixture --cv --seed=50 $DATAINPUT/admixture_file_art.bed 9 -j28 | tee log_admixture_art_9.log
admixture --cv --seed=50 $DATAINPUT/admixture_file_art.bed 10 -j28 | tee log_admixture_art_10.log


