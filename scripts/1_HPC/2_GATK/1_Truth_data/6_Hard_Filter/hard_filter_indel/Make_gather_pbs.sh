#!/bin/bash

FILE="/project/st-angert-1/list/indel_cov_summary.csv"

while IFS=, read -r chr d1 d2 ;

do

OUT="hard_filter_indel_chr"$chr

cat hard_filter_indel_template_2.pbs | sed "s|XX|$chr|g"|\
	sed "s|YY|$d1|g" |\
	sed "s|ZZ|$d2|g" >$OUT.pbs 

done < $FILE
