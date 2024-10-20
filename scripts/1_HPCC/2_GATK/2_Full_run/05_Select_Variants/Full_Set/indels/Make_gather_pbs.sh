#!/bin/bash

FILE="/project/st-angert-1/list/chr_size.csv"

while IFS=, read -r chr max_in;

do

OUT="select_indel_chr"$chr

cat select_indel_template.pbs | sed "s|XX|$chr|g" >$OUT.pbs 

done < $FILE
