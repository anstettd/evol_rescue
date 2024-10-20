#!/bin/bash

FILE="chr_size.csv"

while IFS=, read -r chr max_in;

do
for a in $(seq -f  "%07.0f" 1 5000000 $max_in)
do
        echo "A ->" $a

        b=$(expr $a + 4999999 )
        echo "B ->" $b

        OUT="-I /scratch/st-angert-1/3_genotypeGVCF/gvcf_all_chr"$chr"_"$a"_"$b".vcf.gz"
        echo $OUT

        echo; echo

if [ $b -lt $max_in ]
then
	echo $OUT >> "gather_VCF_In_chr"$chr".csv"
else
	OUT="-I /scratch/st-angert-1/3_genotypeGVCF/gvcf_all_chr"$chr"_"$a"_"$max_in".vcf.gz"
	echo $OUT >> "gather_VCF_In_chr"$chr".csv"
fi
done

done < $FILE
