#!/usr/bin/bash

#declare variables used in the script
DATADIRECTORY=demographic_chapter
SCRIPT=$DATADIRECTORY/00_scripts/wall_genome/vcftools
HEADER=$DATADIRECTORY/00_scripts/headerP.txt
VCFTOOLS=". /vcftools/0.1.16/env.sh"
DATAOUTPUT=/vcftools

#create the output files
mkdir -p $SCRIPT
mkdir -p $DATAOUTPUT

#go to the folder containing the data in your project folder
cd $SCRATCH/06_freebayes/partial_filter
ls -d /partial_filter/*.vcf.vcf> /vcftools_files.txt

#list the file to filter
NAME='cat /vcftools_files.txt'
MAF=0.05
MISSING=0.9
MIN_MEAN_DEPTH=30
QUAL=30

#start the loop creating individual scripts for each VCF files located in the data folder
for FILE in $($NAME) #list all the files with the fna extension and store it in the variable FILE
do
        cp $HEADER $SCRIPT/vcftools_${FILE##*/}.sh ; # copy the header file in a new script file called vcftools_{genome_part}.sh
        echo "#PBS -N vcftools_${FILE##*/}" >> $SCRIPT/vcftools_${FILE##*/}.sh ;
        echo "#PBS -o $DATADIRECTORY/98_log_files/10_vcftools_${FILE##*/}.log" >> $SCRIPT/vcftools_${FILE##*/}.sh ;
        echo "cd $SCRATCH/06_freebayes/partial_filter" >> $SCRIPT/vcftools_${FILE##*/}.sh ; # append the echoed line in the script file (we go to the data folder of our project)
        echo "$VCFTOOLSENV" >> $SCRIPT/vcftools_${FILE##*/}.sh ; # preparing the vcftools
        echo "vcftools --temp $SCRATCH --max-missing $MISSING --min-meanDP $MIN_MEAN_DEPTH --minQ $QUAL --vcf ${FILE} --recode --out $DATAOUTPUT/filter_NA_${FILE##*/}" >> $SCRIPT/vcftools_${FILE##*/}.sh ;
        qsub $SCRIPT/vcftools_${FILE##*/}.sh ; # append the echoed line in the script file (here we ask to submit our script to the PBS claculation nodes for execution)
done ; # we finish the loop
