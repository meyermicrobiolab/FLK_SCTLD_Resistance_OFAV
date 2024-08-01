#!/bin/bash
#SBATCH --job-name=eukrep_%j
#SBATCH --output=logs/eukrep_%j.log
#SBATCH --error=logs/eukrep_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50gb
#SBATCH --time=4-00:00:00

# load in eukrep. Note, eukrep is not currently loaded as part of the modules
# install via conda following the instructions at https://github.com/patrickwest/EukRep
# activate the conda environment PRIOR to executing the script

# Execute the following from the folder where the script is with:
# conda activate eukrep-env
# sbatch scripts/eukrep.sh 

# Get version for papers
EukRep --version

# Run eukrep to output a prokaryote assembly and a eukaryote assembly. Need to make sure conda environment is active. 
for assembly in ./output/assemblies/coral_set/*_assembly/final.contigs.fa
do
    CORAL=$(echo "$assembly" | awk -F'/' '{split($5, a, "_"); print a[1]}')
    mkdir ./output/eukrep/coral_set/${CORAL}_eukrep
    EukRep -i ${assembly} \
        -o ./output/eukrep/coral_set/${CORAL}_eukrep/euk.final.contigs.fa \
        --prokarya ./output/eukrep/coral_set/${CORAL}_eukrep/prok.final.contigs.fa \
        --min 1000
done

# Count the number of contigs in each eukaryote and prokaryote assembly for reporting in the MS
for assembly in ./output/assemblies/coral_set/*_assembly/final.contigs.fa
do
    CORAL=$(echo "$assembly" | awk -F'/' '{split($5, a, "_"); print a[1]}')
    echo "number of TOTAL contigs in ${CORAL} assembly is: `grep -c "^>" ${assembly}`" 
    echo "number of EUKARYOTIC contigs in ${CORAL} assembly is: `grep -c "^>" ./output/eukrep/coral_set/${CORAL}_eukrep/euk.final.contigs.fa`" 
    echo "number of PROKARYOTIC contigs in ${CORAL} assembly is: `grep -c "^>" ./output/eukrep/coral_set/${CORAL}_eukrep/prok.final.contigs.fa`"
done
