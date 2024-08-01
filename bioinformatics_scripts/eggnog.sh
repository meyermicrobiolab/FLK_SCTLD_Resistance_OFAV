#!/bin/bash
#SBATCH --job-name=eggnog_%j
#SBATCH --output=logs/eggnog_%j.log
#SBATCH --error=logs/eggnog_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user cynthiabecker@ufl.edu
#SBATCH --mem-per-cpu=2gb
#SBATCH --cpus-per-task=100
#SBATCH --time=4-00:00:00
#SBATCH --account=juliemeyer
#SBATCH --qos=juliemeyer-b

export OMP_NUM_THREADS=100

# run this script from /blue/juliemeyer/cynthiabecker/ofav_mags/ with "sbatch scripts/eggnog.sh"

# load eggnog-mapper and print version
module load eggnog-mapper/2.1.6
emapper.py --version

# make a directory to house results 
mkdir /blue/juliemeyer/share/coassembly_rrc_ofav/eggnog_output/prodigal
cd /blue/juliemeyer/share/coassembly_rrc_ofav/eggnog_output/prodigal

# run eggnog-mapper on the complete FLK co-assembly
emapper.py -i /blue/juliemeyer/share/coassembly_rrc_ofav/prodigal_output/FLK_OFAV_MG_pred.faa \
    --itype proteins \
    --data_dir /data/reference/eggnog-mapper/2.1.6 \
    -o FLK_OFAV_MG \
    --decorate_gff yes \
    --excel --cpu 100

### PARAMETERS ###
# -i is input fasta file
# -itype specifies that i input a protein fasta file (output from prodigal)
# --data_dir is eggnog's database
# -o is the prefix for all output files (I chose FLK_OFAV_MG to represent florida keys, orbicella faveolata, and metagenome)
# --decorate_gff adds annotations to the fasta headers