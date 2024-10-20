#!/bin/bash
#################################################################
#Population Genomic Analysis pipeline for Mimulus cardinalis

#Author: Daniel Anstett
#Last updated March 7 2024
#################################################################


#################################################################
# Make Pop VCF


#!/bin/bash

#SBATCH --job-name="pi_pop_vcf_batch"
#SBATCH --account=st-angert-1
#SBATCH -t 00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniel.anstett@botany.ubc.ca
#SBATCH --output=pi_pop_vcf_batch.output.txt
#SBATCH --error=pi_pop_vcf_batch.error.txt

################################################################################

source ~/miniconda3/etc/profile.d/conda.sh
conda activate gatk4

cd $PBS_O_WORKDIR
for f1 in /scratch/st-angert-1/16_popVCF/space_time/basetime/*csv; 

do
base_vcf=`basename $f1 .csv`


vcftools --gzvcf /scratch/st-angert-1/20_parse_beagle_basetime/$base_vcf".vcf.gz" --site-pi --out /scratch/st-angert-1/21_beagle_pi/$base_vcf

sleep 10


done

conda deactivate

