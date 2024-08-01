#!/bin/bash

#SBATCH --job-name=salmon_quant_%j
#SBATCH --output=logs/salmon_quant_%j.log
#SBATCH --error=logs/salmon_quant_%j.err
#SBATCH --mem-per-cpu=5gb
#SBATCH --cpus-per-task=36
#SBATCH --time=4-00:00:00

export OMP_NUM_THREADS=36

## Load salmon and get version since it is a module already on the HPC
module load salmon/1.10.1
salmon --version

## move to the folder where you'd like the output
cd ./output/salmon/

## Start by indexing the nucleotide fasta file for mapping
salmon index -t /path/to/coassembly/coassembly_rrc_ofav/prodigal_output/FLK_OFAV_MG_pred.fna \
    -i MG_index \
    -k 31

## Next, run the main salmon command in a for loop to quantify all samples
PATH_TO_READS="/path/to/reads"

for file in $PATH_TO_READS/*1.fastq.gz
do
        fastq=$(basename $file)
        tail1=_1.fastq.gz
        tail2=_2.fastq.gz
        coral_id=${fastq%$tail1}
        echo "Processing sample ${coral_id}"
        salmon quant --meta -i MG_index --libType A \
                -1 ${PATH_TO_READS}/${coral_id}${tail1} \
                -2 ${PATH_TO_READS}/${coral_id}${tail2} \
                -o ${coral_id}.quant \
                -p 36
done

## INFO ON FLAGS
## -k is length of kmers
## -i is index name
## Salmon quant info:
## --meta changes the salmon quant command to be more specific for metagenomic data rather than transcriptomic data
##  --libType A tells the program to automatically decide the library type 
## -o is the output name
## -p is the nubmer of threads/cpus