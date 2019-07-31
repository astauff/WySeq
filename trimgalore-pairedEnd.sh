#!/bin/bash
for files in *R1*q.gz
do
file2="$(echo $files | sed 's/R1/R2/g')"
echo "Trimming read pair including $files and $file2."
trim_galore --paired -j 2 -o trimGalored/ $files $file2
done
