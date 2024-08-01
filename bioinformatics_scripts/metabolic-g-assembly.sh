#!/bin/bash
#SBATCH --job-name=metabolic_%j
#SBATCH --output=logs/metabolic_%j.log
#SBATCH --error=logs/metabolic_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=25gb
#SBATCH --time=4-00:00:00

export OMP_NUM_THREADS=10

# load packages, which are modules already on the HPC
module load metabolic/4.0

# print the version
METABOLIC-G.pl -version

############## RUN ON THE CO-ASSEMBLY PREDICTED PROTEINS ###############

METABOLIC-G.pl -in /path/to/prodigal_output/ -o ./output/metabolic/metabolic-coassembly -m-cutoff 0.75 -t 10

#### NOTES ON FLAGS
# -in requires the folder to have a .faa file - the predicted proteins I generated with prodigal
# -o is a new folder