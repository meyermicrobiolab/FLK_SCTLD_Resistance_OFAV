#!/bin/bash
#SBATCH --job-name=rgi_%j
#SBATCH --output=logs/rgi_%j.log
#SBATCH --error=logs/rgi_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=cynthiabecker@ufl.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=5gb
#SBATCH --time=4-00:00:00
#SBATCH --account=juliemeyer
#SBATCH --qos=juliemeyer-b

export OMP_NUM_THREADS=20

# execute script from the /blue/juliemeyer/cynthiabecker/ofav_mags/ folder using "sbatch scripts/rgi.sh"

module load rgi/6.0.2
rgi main --version

mkdir /blue/juliemeyer/share/coassembly_rrc_ofav/FLK_coassembly/output/rgi
cd /blue/juliemeyer/share/coassembly_rrc_ofav/FLK_coassembly/output/rgi

# load the CARD database (antibiotic resistance genes) into the local working directory
# Note the database is located in my home folder at (~/databases/)

rgi load --card_json /home/cynthiabecker/databases/card.json --local

rgi main --input_sequence ../prodigal_prok/FLK_OFAV_MG_prok_pred.faa \
    --output_file ./FLK_OFAV_MG_prok_AMR \
    --local --clean --include_loose \
    -t protein \
    -n 20

## PARAMETER DETAILS ##
# --input_sequence - I used an amino acid fasta file, so that it will skip Prodigal ORF identification
# -t protein - I did this becuasee I included an amino acid fasta file
# --output_file - Makes a directory named this and this is the basename for output files'
# --local uses the database I loaded with the "rgi load" command. I loaded the db in my home folder (details on this are on the GITHUB documentation)
# --clean gets rid of temp files
# --include_loose Means I will include loose hits as well as the strict and perfect hist. Since this is from an environmental sample, I expect the AMR genes will differ. This setting is better for detecting new or distant homologs of AMR genes, which is what I expect for my data. 
# -n is number of threads