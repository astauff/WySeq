mkdir subreadCounts


for files in *.bam
do
featureCounts -T 8 -M -a $1 -o ./subreadCounts/${files%_Aligned.sortedByCoord.out.bam}_counts $files 
done

