#!/bin/bash
#SBATCH --job-name=logs/bwa_mem_%j
#SBATCH --output=logs/bwa_mem_%j.log
#SBATCH --error=logs/bwa_mem_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=cynthiabecker@ufl.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --mem-per-cpu=5gb
#SBATCH --time=4-00:00:00
#SBATCH --account=juliemeyer
#SBATCH --qos=juliemeyer-b

export OMP_NUM_THREADS=36

# execute script from the /blue/juliemeyer/cynthiabecker/ofav_mags/ folder using "sbatch scripts/bwa_mem.sh"

# load in bwa and samtools
module load bwa
module load samtools

# check versions of packages (useful for papers)
bwa
samtools --version

# Save a file path to the cleaned fastqs needed for mapping to the assemblies
# PATH_TO_FASTQ="/blue/juliemeyer/share/rrc_ofav/extra_fastqs"

# This script should only be run if you already did the bwa index command

# Loop through each assembly, found in sample-specific folder containing the assembly and index files from bwa index
# For each assembly, make a new directory for the mapping results and map the reads back to the assembly
# Pipe the results directly to samtools so that you can generate the sorted bam file

#for assembly in ./output/bwa_index/setExtra/*_assembly/final.contigs.fa
#do
#    CORAL=$(echo "$assembly" | awk -F'/' '{split($5, a, "_"); print a[1]}')
#    mkdir ./output/bwa_mem/setExtra/${CORAL}_map
#    bwa mem -t 5 ${assembly} \
#        ${PATH_TO_FASTQ}/${CORAL}_*1.fastq.gz ${PATH_TO_FASTQ}/${CORAL}_*2.fastq.gz | \
#        samtools sort -o ./output/bwa_mem/setExtra/${CORAL}_map/${CORAL}.bam -
#done

###### BWA OPTIONS ######
# -t is the number of threads to specify. I made sure it matches the --cpus-per-task and the OMP_NUM_THREADS that was set above

####### SAMTOOLS OPTIONS #####
# -o is the output file to specify
# - just indicates to use the current file, which is what is piped from bwa mem via |


########################################
###### BWA MEM FOR CO-ASSEMBLY ########
########################################

mkdir ./output/bwa_mem/coassembly_map # output folder
cd ./output/bwa_mem/coassembly_map/ # move to output folder

# set location of reads and coassembly
COASSEMBLY=/blue/juliemeyer/share/coassembly_rrc_ofav/coassembly_FLK_OFAV/final.contigs.fa
PATH_TO_READS="/blue/juliemeyer/share/coassembly_rrc_ofav"

# loop through all set of reads and generate mapping files for each sample
for file in $PATH_TO_READS/*1.fastq.gz
do
    fastq=$(basename $file)
    tail1=_1.fastq.gz
    tail2=_2.fastq.gz
    coral_id=${fastq%$tail1}
    echo "Processing sample ${coral_id}"
    bwa mem -t 36 ${COASSEMBLY} \
        ${PATH_TO_READS}/${coral_id}${tail1} \
        ${PATH_TO_READS}/${coral_id}${tail2} | \
        samtools sort -o ${coral_id}.bam -
done