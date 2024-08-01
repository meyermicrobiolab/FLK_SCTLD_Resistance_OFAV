#!/bin/bash
#SBATCH --job-name=mmseqs_%j
#SBATCH --output=logs/mmseqs_%j.log
#SBATCH --error=logs/mmseqs_%j.err
#SBATCH --mem-per-cpu=10gb
#SBATCH --cpus-per-task=50
#SBATCH --time=4-00:00:00
#SBATCH --partition=bigmem


export OMP_NUM_THREADS=50

# Run the script from anywhere because I provide complete paths

# load mmseqs and info like the version
module load mmseqs2/14
mmseqs

mkdir /path/to/coassembly/coassembly_rrc_ofav/mmseqs_output/prodigal
cd /path/to/coassembly/coassembly_rrc_ofav/mmseqs_output/prodigal

# Run mmseqs on the coding domains (aa fasta) annotated by prodigal
mmseqs easy-taxonomy /path/to/coassembly/coassembly_rrc_ofav/prodigal_output/FLK_OFAV_MG_pred.faa \
    /data/reference/mmseqs/14/GTDB/GTDB \
    /path/to/coassembly/coassembly_rrc_ofav/mmseqs_output/prodigal/FLK_OFAV_MG_tax \
    tempdir \
    --search-type 1 --tax-lineage 1 --threads 50 --compressed 1

### PARAMETERS ####
# The first argument is the predicted proteins (amino acid fasta)
# The second argument is the reference database that should be loaded in the compute cluster somewhere
# The third argument is the output 
# --search-type 1 indicates to do an amino acid search which matches the input file (.faa)
# --tax-lineage 1 will add a column with the full lineage names, prefixed with the first letter of the rank 
# --compressed 1 compresses the output, which can take up SO MANY GB