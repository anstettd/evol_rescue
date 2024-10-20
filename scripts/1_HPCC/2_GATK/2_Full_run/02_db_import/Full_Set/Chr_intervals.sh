#!/bin/bash

FILE="/project/st-angert-1/list/chr_size.csv"

while IFS=, read -r chr max_in;

do
for a in $(seq -f  "%07.0f" 1 5000000 $max_in)
do
        echo "A ->" $a

        b=$(expr $a + 4999999 )
        echo "B ->" $b

        OUT="db_import_chr"$chr"_"$a"_"$b
        echo $OUT

        echo; echo

if [ $b -lt $max_in ]
then

cat db_import_template.pbs | \
       sed "s|XX|$chr|g" | \
        sed "s|YY|$a|g" | \
        sed "s|ZZ|$b|g" >$OUT.pbs
else
OUT="db_import_chr"$chr"_"$a"_"$max_in
cat db_import_template.pbs | \
       sed "s|XX|$chr|g" | \
        sed "s|YY|$a|g" | \
        sed "s|ZZ|$max_in|g" >$OUT.pbs

fi
done

done < $FILE
