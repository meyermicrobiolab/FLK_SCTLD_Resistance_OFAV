#!/bin/bash
#SBATCH --job-name=bt2-build_%j
#SBATCH --output=logs/bt2-build_%j.log
#SBATCH --error=logs/bt2-build_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=20gb


# Use bowtie2, which is a module already loaded on the HPC:

module load bowtie2/2.4.2

bowtie2 --version # get version for good reporting

bowtie2

# First I need to save the metagenome assemblies to a folder where I can put the indexed files.
# Make a directory for each library
# Copy in the assembly into the new directory. This is where the bwa index will live as well

for assembly in ./output/eukrep/coral_set/*_eukrep/euk.final.contigs.fa
do
    CORAL=$(echo "$assembly" | awk -F'/' '{split($5, a, "_"); print a[1]}')
    bowtie2-build -f ${assembly} ./output/eukrep/coral_set/${CORAL}_eukrep/${CORAL}_build
done
