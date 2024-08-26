#!/bin/bash
#SBATCH --job-name=dRep_%j
#SBATCH --output=logs/dRep_%j.log
#SBATCH --error=logs/dRep_%j.err
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

# Run this script from the blue/juliemeyer/cynthiabecker/ofav_mags/ folder

# load packages
module load drep/2.3.2

# Use the dRep dereplicate workflow
dRep dereplicate output/dRep/dRep_5.14.24/ -g output/gtdb/coassembly_genomes/*.fa -p 5 -comp 50 

##### PARAMETERS
# -g this is the path to genomes
# -p this is the number of processors. I changed to 5, but the default is 6
# -comp this is the cutoff for completeness of the genomes. The default is 75, but I have lots of low-quality genomes, so I set the default to 50%. 