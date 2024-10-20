#!/bin/bash

FILE="/project/st-angert-1/list/chr_size.csv"

while IFS=, read -r chr max_in;

do

OUT="make_sites_snp_chr"$chr

cat make_sites_snp_template.pbs| sed "s|XX|$chr|g" >$OUT.pbs 

done < $FILE 

rm *chr0*

