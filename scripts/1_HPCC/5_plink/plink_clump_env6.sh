#!/bin/sh

#################################################################
#Population Genomic Analysis pipeline for Mimulus cardinalis

#Author: Daniel Anstett
#################################################################


#################################################################
# Plink clump of baseline 
#SBATCH --job-name="clump_env6"
#SBATCH --account=st-angert-1
#SBATCH -t 0-01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=170G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniel.anstett@botany.ubc.ca
#SBATCH --output=clump.output.txt
#SBATCH --error=clump.error.txt


################################################################################

cd $SLURM_SUBMIT_DIR

/home/danstett/plink/plink --bfile /scratch/st-angert-1/24_plink/baseline --clump /scratch/st-angert-1/24_plink/plink_input_env6.tsv --clump-snp-field chr_snp --clump-field empirical_p --clump-kb 250 --clump-r2 0.4 --allow-extra-chr --out /scratch/st-angert-1/24_plink/baseline_env6 --clump-p1 0.9 --clump-p2 0.9

