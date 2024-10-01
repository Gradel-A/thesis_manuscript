#!/usr/bin/bash

#declare variables used in the script
DATADIRECTORY=./demographic_chapter
SCRIPT=$DATADIRECTORY/00_scripts/wall_genome/ind_selection
HEADER=$DATADIRECTORY/00_scripts/headerS.txt
BCFTOOLSENV=". /bcftools/1.17/env.sh"

#create the output files
mkdir -p $SCRIPT
mkdir -p $DATAOUTPUT

#go to the folder containing the data in your project folder
cd $SCRATCH/06_freebayes/partial_vctools
ls -d /*.vcf.gz> /ind_selection.txt

#list the file to filter
NAME='cat /ind_selection.txt'

#start the loop creating individual scripts for each VCF files located in the data folder
for FILE in $($NAME) #list all the files with the fna extension and store it in the variable FILE
do
        cp $HEADER $SCRIPT/ind_selection_${FILE##*/}.sh ; # copy the header file in a new script file called ind_selection_{genome_part}.sh
        echo "#PBS -N ind_selection_${FILE##*/}" >> $SCRIPT/ind_selection_${FILE##*/}.sh ;
        echo "#PBS -o $DATADIRECTORY/98_log_files/0_1_ind_selection_${FILE##*/}.log" >> $SCRIPT/ind_selection_${FILE##*/}.sh ;
        echo "$BCFTOOLSENV" >> $SCRIPT/ind_selection_${FILE##*/}.sh;
        echo "cd $SCRATCH/wall_genome_vcf" >> $SCRIPT/ind_selection_${FILE##*/}.sh ;
        echo "bcftools view -o ind_selected_${FILE##*/} -O z -S /wall_genome_vcf/list_ind_selected_vcf.txt ${FILE}" >> $SCRIPT/ind_selection_${FILE##*/}.sh ;
      #  qsub $SCRIPT/ind_selection_${FILE##*/}.sh ; # append the echoed line in the script file (here we ask to submit our script to the PBS claculation nodes for execution)
done ; # we finish the loop 
