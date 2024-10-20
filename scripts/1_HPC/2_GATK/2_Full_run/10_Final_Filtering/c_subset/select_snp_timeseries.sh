#!/bin/sh
#################################################################
#Population Genomic Analysis pipeline for Mimulus cardinalis

#Author: Daniel Anstett
#Last updated Feb 22 2024
#################################################################


#################################################################
# Generate VCF for only timeseries populations 



#!/bin/bash

#SBATCH --job-name="select_snp_timeseries"
#SBATCH --account=st-angert-1
#SBATCH -t 1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniel.anstett@botany.ubc.ca
#SBATCH --output=select_snp_timeseries.output.txt
#SBATCH --error=select_snp_timeseries.error.txt

################################################################################

source ~/miniconda3/etc/profile.d/conda.sh
conda activate gatk4

cd $PBS_O_WORKDIR

#Timeseries
gatk SelectVariants -V /scratch/st-angert-1/10_select_final_snps/all_vqsr_filtered_variants.vcf.gz  -R /scratch/st-angert-1/mc_ref/CE10g_v2.0.fa -select-type SNP -O /scratch/st-angert-1/10_select_final_snps/timeseries_filtered_variants.vcf.gz --arguments_file /project/st-angert-1/list/timeseries_names.csv

gatk IndexFeatureFile -I /scratch/st-angert-1/10_select_final_snps/timeseries_filtered_variants.vcf.gz

conda deactivate


