#!/bin/bash

FILE="chr_size.csv"

while IFS=, read -r chr max_in;

do
for a in $(seq -f  "%07.0f" 1 5000000 $max_in)
do
        echo "A ->" $a

        b=$(expr $a + 4999999 )
        echo "B ->" $b

        OUT="genotypeGVCF_all_chr"$chr"_"$a"_"$b
        echo $OUT

        echo; echo

if [ $b -lt $max_in ]
then

cat genotypeGVCF_template_all.pbs | \
       sed "s|XX|$chr|g" | \
        sed "s|YY|$a|g" | \
        sed "s|ZZ|$b|g" >$OUT.pbs
else
OUT="genotypeGVCF_all_chr"$chr"_"$a"_"$max_in
cat genotypeGVCF_template_all.pbs | \
       sed "s|XX|$chr|g" | \
        sed "s|YY|$a|g" | \
        sed "s|ZZ|$max_in|g" >$OUT.pbs

fi
done

done < $FILE
