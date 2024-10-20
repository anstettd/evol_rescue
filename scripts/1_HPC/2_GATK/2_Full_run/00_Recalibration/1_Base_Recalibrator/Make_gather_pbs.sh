#!/bin/bash

FILE="/project/st-angert-1/list/numbered_marked_rg_bam.csv"

while IFS=, read -r sample_num sample ;

do

OUT="base_recal_num"$sample_num

cat base_recal_template.pbs | sed "s|XX|$sample|g"|\
	sed "s|YY|$sample_num|g" >$OUT.pbs

done < $FILE
