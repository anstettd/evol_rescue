#!/bin/sh
##################################################################
##Population Genomic Analysis pipeline for Mimulus cardinalis
#
##Author: Daniel Anstett
##Last updated Feb 23 2024
##################################################################
#
#
##################################################################
## Prep data for BayPass 


#!/bin/bash

#SBATCH --job-name="BP_Prep"
#SBATCH --account=st-angert-1
#SBATCH -t 1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=170G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniel.anstett@botany.ubc.ca
#SBATCH --output=BP_Prep.output.txt
#SBATCH --error=BP_Prep.error.txt

################################################################################

source ~/miniconda3/etc/profile.d/conda.sh
conda activate BP_PREP

cd $PBS_O_WORKDIR

Rscript /scratch/st-angert-1/scripts/11a_baypass_prep/PB_data_prep.R /scratch/st-angert-1/scripts/11a_baypass_prep/baseline_pop_id.txt /scratch/st-angert-1/10_select_final_snps/baseline_filtered_variants.vcf.gz /scratch/st-angert-1/11_baypass/ /scratch/st-angert-1/scripts/11a_baypass_prep/vcf2baypass.pl

conda deactivate


