#!/bin/bash
#SBATCH --job-name=prodigal_%j
#SBATCH --output=logs/prodigal_%j.log
#SBATCH --error=logs/prodigal_%j.err
#SBATCH --mem=100gb

# load prodigal which is a module on the HPC and print version
module load prodigal/2.6.3
prodigal -v

# make a directory to house results 
mkdir /path/to/prodigal_output
cd /path/to/prodigal_output

# run prodigal on the complete FLK co-assembly to get predicted genes
prodigal -i /path/to/coassembly/coassembly_FLK_OFAV/final.contigs.fa \
    -o FLK_OFAV_MG_coords.gff \
    -a FLK_OFAV_MG_pred.faa \
    -d FLK_OFAV_MG_pred.fna \
    -f gff \
    -p meta 

### PARAMETERS ###
# -i is input fasta file
# -o is output file, I wanted a gff file since I am more used to having those for reference
# -f gff is needed to make sure the output file is a gff not gbk
# -a is the output predicted protein translations
# -d is the nucleotide sequences of the predicted genes
# -p meta indicates that the input is a metagenome