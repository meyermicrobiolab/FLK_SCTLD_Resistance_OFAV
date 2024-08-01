#!/bin/bash
#SBATCH --job-name=bt2_%j
#SBATCH --output=logs/bt2_%j.log
#SBATCH --error=logs/bt2_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --mem-per-cpu=4gb

export OMP_NUM_THREADS=50


# Use bowtie2 to remove eukaryote reads from symbiont host-removed reads
# bowtie2 is a module already loaded on the HPC

module load bowtie2/2.4.2

bowtie2 --version # get version for good reporting

bowtie2

for read1 in $(cat scripts/FLK_1_reads.txt)
do 
    CORAL=$(echo "$read1" | awk -F'/' '{split($7, a, "_"); print a[1]}')
    bowtie2 -q --phred33 -x ./output/bowtie2/indices/${CORAL}_build \
        -1 ${read1} \
        -2 /path/to/the/hostclean_reads/data/${CORAL}_*_2.fastq.gz \
        -S ${CORAL}_noneuk.sam \
        --un-conc-gz /path/to/noneuk_reads/${CORAL}_noneuk_%.fastq.gz \
        -p 50
    rm ${CORAL}_noneuk.sam
done

#### PARAMTER EXPLANATION ####
# -q reads are fastq
# -phred33 is the quality encoding. illumina uses phred+33 encoding
# --un-conc-gz makes the paired-end reads that DID NOT align write to the specified file. % will bereplaced with either 1 or 2 for the corresponding read

#### FOR LOOP EXPLANATION ####
# I am looping through each line in a file "FLK_1_reads.txt", which lists the full paths to read1 files used for metagenome assemblies for Florida Keys samples. 
# These files paths are to the host/symbiont cleaned reads
# For each file path, I extract a variable, "CORAL", that corresponds to the coral ID. This happens by printing the file path and parsing it by "/", taking the 7th part of it (the file name), parsing that by an underscore and taking the first bit
# I then start using bowtie2 on each sample, providing the path to the basename of the bowtie2-build index "CORAL_build"
# I tell bowtie to use both the forward and reverse reads. THe forward reads are specified with "read1"
# The reverse reads are specified as the file path and regular expressions (*) that fill in variable bits except the CORAL ID
# bowtie also outputs a sam file, that I will name based on the CORAL ID, but this isn't needed for downstream applications, so remove it
# Bowtie uses the --un-conc-gz flag to save paired end reads that did NOT align. I will include that flag and save all files to the same folder for ease of doing the co-assembly analyses. 
