mkdir -p alignmentOutputs

for files in *R1*q.gz
do
file2="$(echo $files | sed 's/R1/R2/g')"
file2="$(echo $file2 | sed 's/val_1/val_2/g')"

echo "Beginning alignment for read pair including $files and $file2."

outprefix="$(echo $files | awk -F 'R1' '{print $1}')"

STAR --runThreadN 4 --genomeDir $1 --readFilesIn  $files $file2 --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 100000 --alignMatesGapMax 100000 --outFileNamePrefix "./alignmentOutputs/$outprefix" --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts

done
