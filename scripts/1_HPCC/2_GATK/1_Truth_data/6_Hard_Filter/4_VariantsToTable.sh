# Variants to Table script
# Daniel Anstett
# May 12 2021 


#Run this ins in directory with gold vcf's before hard filtering to get variant stats
# see R script for getting stats

#indels
gatk VariantsToTable -V gold_chr1_indels.vcf.gz -F CHROM -F POS -F DP -O cov_indel2_chr1.tsv
gatk VariantsToTable -V gold_chr2_indels.vcf.gz -F CHROM -F POS -F DP -O cov_indel2_chr2.tsv
gatk VariantsToTable -V gold_chr3_indels.vcf.gz -F CHROM -F POS -F DP -O cov_indel2_chr3.tsv
gatk VariantsToTable -V gold_chr4_indels.vcf.gz -F CHROM -F POS -F DP -O cov_indel2_chr4.tsv
gatk VariantsToTable -V gold_chr5_indels.vcf.gz -F CHROM -F POS -F DP -O cov_indel2_chr5.tsv
gatk VariantsToTable -V gold_chr6_indels.vcf.gz -F CHROM -F POS -F DP -O cov_indel2_chr6.tsv
gatk VariantsToTable -V gold_chr7_indels.vcf.gz -F CHROM -F POS -F DP -O cov_indel2_chr7.tsv
gatk VariantsToTable -V gold_chr8_indels.vcf.gz -F CHROM -F POS -F DP -O cov_indel2_chr8.tsv




#snps
gatk VariantsToTable -V gold_chr1_snp.vcf.gz -F CHROM -F POS -F DP -O cov_snp2_chr1.tsv
gatk VariantsToTable -V gold_chr2_snp.vcf.gz -F CHROM -F POS -F DP -O cov_snp2_chr2.tsv
gatk VariantsToTable -V gold_chr3_snp.vcf.gz -F CHROM -F POS -F DP -O cov_snp2_chr3.tsv
gatk VariantsToTable -V gold_chr4_snp.vcf.gz -F CHROM -F POS -F DP -O cov_snp2_chr4.tsv
gatk VariantsToTable -V gold_chr5_snp.vcf.gz -F CHROM -F POS -F DP -O cov_snp2_chr5.tsv
gatk VariantsToTable -V gold_chr6_snp.vcf.gz -F CHROM -F POS -F DP -O cov_snp2_chr6.tsv
gatk VariantsToTable -V gold_chr7_snp.vcf.gz -F CHROM -F POS -F DP -O cov_snp2_chr7.tsv
gatk VariantsToTable -V gold_chr8_snp.vcf.gz -F CHROM -F POS -F DP -O cov_snp2_chr8.tsv




