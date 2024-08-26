#!/bin/bash
#SBATCH --job-name=logs/index_%j
#SBATCH --output=logs/index_%j.log
#SBATCH --error=logs/index_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=cynthiabecker@ufl.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10gb
#SBATCH --time=4-00:00:00
#SBATCH --account=juliemeyer
#SBATCH --qos=juliemeyer-b

# execute script from the /blue/juliemeyer/cynthiabecker/ofav_mags/ folder using "sbatch scripts/bwa_index.sh"

# load in bwa and check version
module load bwa
bwa

# First I need to save the metagenome assemblies to a folder where I can put the indexed files.
# Make a directory for each library
# Copy in the assembly into the new directory. This is where the bwa index will live as well

#for path in /blue/juliemeyer/share/rrc_ofav/extra_fastqs/*megahit_out_metalarge*
#do
#    CORAL=$(echo "$path" | awk -F'/' '{split($NF, a, "_"); print a[1]}')
#    mkdir ./output/bwa_index/setExtra/${CORAL}_assembly
#    cp ${path}/final.contigs.fa ./output/bwa_index/setExtra/${CORAL}_assembly/
#done

# Next I need to make an index of each assembly. 
# I will loop through every directory containing an assembly, index it, 
#for assembly in ./output/bwa_index/setExtra/*_assembly/final.contigs.fa
#do
#    bwa index ${assembly}
#done

########################################
###### BWA MEM FOR CO-ASSEMBLY ########
########################################

bwa index /blue/juliemeyer/share/coassembly_rrc_ofav/coassembly_FLK_OFAV/final.contigs.fa