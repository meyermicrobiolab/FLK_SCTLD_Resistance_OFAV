#!/bin/bash
#SBATCH --job-name=coverm_%j
#SBATCH --output=logs/coverm_%j.log
#SBATCH --error=logs/coverm_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=cynthiabecker@ufl.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=30gb
#SBATCH --time=4-00:00:00
#SBATCH --account=juliemeyer
#SBATCH --qos=juliemeyer-b

export OMP_NUM_THREADS=5

# GOAL: Get abundances for contigs that are likely producing antimicrobial-related compounds or toxins.

# Code for this script is from the blue/juliemeyer/cynthiabecker/ofav_mags/ folder
# Run the code in the location where the fasta file of putatively antimicrobial-encoding contigs are

# load packages
module load coverm

# coverm version - print to log file
coverm --version

# Save a file path to the cleaned fastqs needed for mapping to the genomes
# Base this on the file that contains a list of all the file paths. 
PATH_TO_R1=$(tr '\n' ' ' < FLK_1_noneuk_reads.txt)
PATH_TO_R2=$(tr '\n' ' ' < FLK_2_noneuk_reads.txt)

# Run coverM on all the reads (host-removed), to map them to the dereplicated genomes
# Generate separate output files for each method

coverm contig -1 ${PATH_TO_R1} -2 ${PATH_TO_R2} \
    --reference antimicrobial.fasta \
    --methods mean \
    --min-covered-fraction 10 \
    --output-file antimicro_mean.tsv --threads 5

coverm contig -1 ${PATH_TO_R1} -2 ${PATH_TO_R2} \
    --reference antimicrobial.fasta \
    --methods rpkm \
    --min-covered-fraction 10 \
    --output-file antimicro_rpkm.tsv --threads 5

####### PARAMETERS
# -1 is a list of all forward reads, formatted as full paths with a single space between each path
# -2 is a list of all reverse reads, formatted as full paths with a single space between each path
# --methods provides the methods to calculate. rpkm and mean are chosen.
# --min-covered-fraction 10 refers that the contig with less covered bases than this are reported as having zero coverage
# --output-file is the file to write the abundances to 
# --threads is the number of threads or CPUs to use, where the default is 1 
