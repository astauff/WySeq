#!/bin/bash
for files in *.fastq.gz
do 
trim_galore -j 2 -o trimGalored/ $files
done
