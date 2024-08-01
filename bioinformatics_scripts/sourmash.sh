#!/bin/bash
#SBATCH --job-name=sour_%j
#SBATCH --output=logs/sour_%j.log
#SBATCH --error=logs/sour_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=100gb

# Run this code assuming sourmash is already loaded as a module on the HPC
module load sourmash/4.5.0
sourmash --version

sourmash sketch dna -p k=31 \
    /path/to/noneuk_reads/*1.fastq.gz \
    --output-dir /path/to/output/sourmash/

 sourmash compare /path/to/output/sourmash/*.sig \
    -o /path/to/output/sourmash/cmp.dist

sourmash plot /path/to/output/sourmash/cmp.dist \
    --labels