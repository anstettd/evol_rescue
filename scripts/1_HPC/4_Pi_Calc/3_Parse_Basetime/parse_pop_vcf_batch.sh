#!/bin/sh
#################################################################
#Population Genomic Analysis pipeline for Mimulus cardinalis

#Author: Daniel Anstett
#Last updated March 7 2024
#################################################################


#################################################################
# Make Pop VCF


#!/bin/bash

#SBATCH --job-name="make_pop_vcf_batch"
#SBATCH --account=st-angert-1
#SBATCH -t 2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=50G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniel.anstett@botany.ubc.ca
#SBATCH --output=make_pop_vcf_batch.output.txt
#SBATCH --error=make_pop_vcf_batch.error.txt

################################################################################

source ~/miniconda3/etc/profile.d/conda.sh
conda activate samtools

cd $PBS_O_WORKDIR
for f1 in /scratch/st-angert-1/16_popVCF/space_time/basetime/*csv; 

do
base_vcf=`basename $f1 .csv`

bcftools view -S $f1 -Oz -o /scratch/st-angert-1/20_parse_beagle_basetime/$base_vcf".vcf.gz" /scratch/st-angert-1/18_beagle/beagle_output.vcf.gz --threads 16
sleep 20

done

conda deactivate

