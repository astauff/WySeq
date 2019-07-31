#!/bin/bash
cd subreadCounts

for files in *_counts
do
outprefix="$(echo $files | awk -F 'Aligned' '{print $1}')"
tail -n +2 $files > ${outprefix}_processedCounts.tab


done
