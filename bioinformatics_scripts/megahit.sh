#!/bin/bash
#SBATCH --job-name=megahit_%j
#SBATCH --output=megahit_%j.log
#SBATCH --error=megahit_%j.err
#SBATCH --mem-per-cpu=10gb
#SBATCH --cpus-per-task=100
#SBATCH --partition=bigmem
#SBATCH --time=4-00:00:00

export OMP_NUM_THREADS=100

# Run this script from the folder with the reads to be assembled:

## Load packages which are already loaded on the HPC as modules
module load megahit/1.1.4
megahit --version

module load quast
quast --version

## Make variables from R1.txt and R2.txt files that contain the comma-separated list of files
R1=$(cat R1.txt)
R2=$(cat R2.txt)

## make sure the lists look correct
echo $R1 
echo $R2 

megahit -1 $R1 -2 $R2 --min-contig-len 1000 -o coassembly_FLK_OFAV -t 100 --presets meta-large

### PARAMETERS ###
# -1 This is a comma-separated list of the forward reads (gzipped) which I included in the R1.txt file
# -2 This is a comma-separated list of the reverse reads (gzipped) which I included in the R2.txt file
# --min-contig-len means all contigs will be over 1000 basepairs
# -o is the output folder name
# -t is the number of cpu threads
# --presets meta-large is for large and complex metagenomes such as soil microbiomes, whihc I apply here because coral microbiomes are complex

quast.py ./coassembly_FLK_OFAV/final.contigs.fa