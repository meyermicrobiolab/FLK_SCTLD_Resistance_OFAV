#!/bin/bash
#SBATCH --job-name=eukrep_%j
#SBATCH --output=logs/eukrep_%j.log
#SBATCH --error=logs/eukrep_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50gb

# load in eukrep. Note, eukrep is not currently loaded as part of the modules on the HPC. 
# install via conda following the instructions at https://github.com/patrickwest/EukRep
# activate the conda environment PRIOR to executing the script

# Execute the following after you've installed the program with conda
# conda activate eukrep-env
# submit script with "sbatch /path/to/script.sh"

# Get version for papers
EukRep --version

# move to the output folder
cd /path/to/coassembly/FLK_coassembly/output/eukrep

# 7/29/24 - I need to run eukrep on the co-assembly
EukRep -i //path/to/coassembly/FLK_coassembly/coassembly_FLK_OFAV/final.contigs.fa \
    -o ./euk.final.contigs.fa \
    --prokarya ./prok.final.contigs.fa \
    --min 1000

echo "number of TOTAL contigs in FLK Co-assembly is: `grep -c "^>" /path/to/coassembly/FLK_coassembly/coassembly_FLK_OFAV/final.contigs.fa`" 
echo "number of EUKARYOTIC contigs in FLK Co-assembly is: `grep -c "^>" ./euk.final.contigs.fa`" 
echo "number of PROKARYOTIC contigs in FLK Co-assembly is: `grep -c "^>" ./prok.final.contigs.fa`"