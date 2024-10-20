#!/bin/bash

FILE="chr_size.csv"

while IFS=, read -r chr max_in;

do

OUT="gather_VCF_chr"$chr

cat gather_VCF_template.pbs | sed "s|XX|$chr|g" >$OUT.pbs 

done < $FILE
