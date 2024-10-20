#!/bin/sh
##################################################################
##Population Genomic Analysis pipeline for Mimulus cardinalis
#
##Author: Daniel Anstett
##Last updated Feb 22 2024
##################################################################
#
#
##################################################################
# Reformat timeseries vcf into abundance table per site/year 


#!/bin/bash

#SBATCH --job-name="BP_Prep_time_pop"
#SBATCH --account=st-angert-1
#SBATCH -t 1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=170G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniel.anstett@botany.ubc.ca
#SBATCH --output=BP_Prep_time_pop.output.txt
#SBATCH --error=BP_Prep_time_pop.error.txt
 
################################################################################

source ~/miniconda3/etc/profile.d/conda.sh
conda activate BP_PREP

cd $PBS_O_WORKDIR

Rscript /scratch/st-angert-1/scripts/11a_baypass_prep/PB_data_prep_time.R /scratch/st-angert-1/scripts/11a_baypass_prep/timeseries_pop_id.txt /scratch/st-angert-1/10_select_final_snps/timeseries_filtered_variants.vcf.gz /scratch/st-angert-1/11_baypass/timeseries_pop/ /scratch/st-angert-1/scripts/11a_baypass_prep/vcf2baypass.pl

conda deactivate


