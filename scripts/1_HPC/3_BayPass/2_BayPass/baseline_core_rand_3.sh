#!/bin/sh
#################################################################
#Population Genomic Analysis pipeline for Mimulus cardinalis

#Author: Daniel Anstett
#Last updated March 7 2024
#################################################################


#################################################################
# BayPass Core Model 



#!/bin/bash

#SBATCH --job-name="core_r3"
#SBATCH --account=st-angert-1
#SBATCH -t 2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=170G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniel.anstett@botany.ubc.ca
#SBATCH --output=baseline_core_rand_3.output.txt
#SBATCH --error=baseline_core_rand_3.error.error.txt

################################################################################ 
################################################################################

cd $PBS_O_WORKDIR

/home/danstett/baypass_2.2/sources/g_baypass -npop 55 -gfile /scratch/st-angert-1/11_baypass/10k_random_genotype/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table10K_random3  -nthreads 8 -outprefix /scratch/st-angert-1/12_baypass_results/core_model_rand_3




