#!/bin/sh

#################################################################
#Population Genomic Analysis pipeline for Mimulus cardinalis

#Author: Daniel Anstett
#################################################################


#################################################################
# Plink vcf to bed 
#SBATCH --job-name="vcf_bed"
#SBATCH --account=st-angert-1
#SBATCH -t 1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=170G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniel.anstett@botany.ubc.ca
#SBATCH --output=vcf_bed.output.txt
#SBATCH --error=vcf_bed.error.txt


################################################################################

cd $SLURM_SUBMIT_DIR

/home/danstett/plink/plink --vcf /scratch/st-angert-1/10_select_final_snps/baseline_filtered_variants.vcf.gz  \
    --keep-allele-order \
	--make-bed \
	--set-missing-var-ids @_# \
	--allow-extra-chr \
	--double-id \
	--out /scratch/st-angert-1/24_plink/baseline

