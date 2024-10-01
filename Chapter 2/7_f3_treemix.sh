#!/usr/bin/env bash

#PBS -q mpi_1
#PBS -l mem=60gb
#PBS -l walltime=03:00:00
#PBS -o /home1/datawork/agradel/demographic_chapter/98_log_files/treemix.log

TREEMIX=". /treemix/1.13/env.sh"
STACK="./stack/env.sh"
DIRECTORY=l/wall_genome_vcf
INPUT=/concatenation
OUTPUT=$SCRATCH/treemix
TREERESULTS=$OUTPUT/chap_demo_tree

mkdir -p $OUTPUT
$STACK
populations --in-vcf $INPUT/dadi_data_superclean_20240918.vcf --treemix -O ./ -M $INPUT/pop_map_dadi_data_superclean_20240918.txt
gzip $INPUT/pop_map_dadi_data_superclean_20240918.txt
#now we can do the analyses base on the files previously generated
$TREEMIX

##F3 stats first
threepop -i $INPUT/dadi_data_superclean_20240918.p.treemix.gz -k 500

mkdir -p $TREERESULTS/migration_0
treemix -i $INPUT/dadi_data_superclean_20240918.p.treemix.gz -k 1 -bootstrap -o $TREERESULTS/migration_0/treemix_results_replicate1
treemix -i $INPUT/dadi_data_superclean_20240918.p.treemix.gz -k 1 -bootstrap -o $TREERESULTS/migration_0/treemix_results_replicate2
treemix -i $INPUT/dadi_data_superclean_20240918.p.treemix.gz -k 1 -bootstrap -o $TREERESULTS/migration_0/treemix_results_replicate3
treemix -i $INPUT/dadi_data_superclean_20240918.p.treemix.gz -k 1 -bootstrap -o $TREERESULTS/migration_0/treemix_results_replicate4
treemix -i $INPUT/dadi_data_superclean_20240918.p.treemix.gz -k 1 -bootstrap -o $TREERESULTS/migration_0/treemix_results_replicate5
treemix -i $INPUT/dadi_data_superclean_20240918.p.treemix.gz -k 2 -bootstrap -o $TREERESULTS/migration_0/treemix_results_replicate2.6
treemix -i $INPUT/dadi_data_superclean_20240918.p.treemix.gz -k 2 -bootstrap -o $TREERESULTS/migration_0/treemix_results_replicate2.7
treemix -i $INPUT/dadi_data_superclean_20240918.p.treemix.gz -k 2 -bootstrap -o $TREERESULTS/migration_0/treemix_results_replicate2.8
treemix -i $INPUT/dadi_data_superclean_20240918.p.treemix.gz -k 2 -bootstrap -o $TREERESULTS/migration_0/treemix_results_replicate2.9
treemix -i $INPUT/dadi_data_superclean_20240918.p.treemix.gz -k 2 -bootstrap -o $TREERESULTS/migration_0/treemix_results_replicate2.10
treemix -i $INPUT/dadi_data_superclean_20240918.p.treemix.gz -k 3 -bootstrap -o $TREERESULTS/migration_0/treemix_results_replicate3.6
treemix -i $INPUT/dadi_data_superclean_20240918.p.treemix.gz -k 3 -bootstrap -o $TREERESULTS/migration_0/treemix_results_replicate3.7
treemix -i $INPUT/dadi_data_superclean_20240918.p.treemix.gz -k 3 -bootstrap -o $TREERESULTS/migration_0/treemix_results_replicate3.8
treemix -i $INPUT/dadi_data_superclean_20240918.p.treemix.gz -k 3 -bootstrap -o $TREERESULTS/migration_0/treemix_results_replicate3.9
treemix -i $INPUT/dadi_data_superclean_20240918.p.treemix.gz -k 3 -bootstrap -o $TREERESULTS/migration_0/treemix_results_replicate3.10
