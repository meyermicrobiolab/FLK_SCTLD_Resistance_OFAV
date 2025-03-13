#!/bin/bash
#SBATCH --job-name=antismash_bact_%j
#SBATCH --output=logs/antismash_bact_%j.log
#SBATCH --error=logs/antismash_bact_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=cynthiabecker@ufl.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=4-00:00:00
#SBATCH --account=juliemeyer
#SBATCH --qos=juliemeyer-b

export OMP_NUM_THREADS=25

# execute script from the /blue/juliemeyer/cynthiabecker/ofav_mags/ folder using "sbatch scripts/antismash_bact.sh"

module load antismash/7.0.0
antismash --version

cd /blue/juliemeyer/share/coassembly_rrc_ofav/FLK_coassembly/output/antismash

# I will test out antismash on both a prokaryotic and eukaryotic genome assemblies
# I expect to find little to nothing, but it cannot hurt to try!

# bacterial mode:
#for prok_assembly in ./output/eukrep/set4/*eukrep/prok.final.contigs.fa
#do
#    CORAL=$(echo "$prok_assembly" | awk -F'/' '{split($5, a, "_"); print a[1]}')
#    antismash ${prok_assembly} \
#        --taxon bacteria \
#        --tigrfam --pfam2go --cc-mibig --rre \
#        --output-dir ./output/antismash/bacteria/${CORAL}_antismash \
#        --output-basename ${CORAL}_prok \
#        --genefinding-tool prodigal-m \
#        --cpus 25
#done

# Run through the prokaryotic portion of the co-assembly
antismash ../eukrep/prok.final.contigs.fa \
        --taxon bacteria \
        --tigrfam --pfam2go --cc-mibig --rre \
        --output-dir ./ \
        --output-basename prok_coassembly \
        --genefinding-tool prodigal-m \
        --cpus 25