#!/bin/bash

FILE="/project/st-angert-1/list/snp_cov2_summary.csv"

while IFS=, read -r chr d1 d2 ;

do

OUT="hard_filter_snp_chr"$chr

cat hard_filter_snp_template.pbs | sed "s|XX|$chr|g"|\
	sed "s|YY|$d1|g" |\
	sed "s|ZZ|$d2|g" >$OUT.pbs 

done < $FILE
