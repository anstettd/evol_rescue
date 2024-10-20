#!/bin/bash

FILE="/project/st-angert-1/list/numbered_marked_rg_bam.csv"

while IFS=, read -r sample_num sample ;

do

OUT="apply_bqsr_num_"$sample_num

cat apply_bqsr_template.pbs | sed "s|XX|$sample|g"|\
	sed "s|YY|$sample_num|g" >$OUT.pbs

done < $FILE
