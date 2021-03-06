#!/bin/bash

#keep this in mind
#export PATH='/home/wyattlab/Desktop/AndrewRNAseq_Pipeline':$PATH

#source_dire="/home/wyattlab/Desktop/AndrewRNAseq_Pipeline/trimGalored"
#anne_source_dir='/media/wyattlab/Data/DataDrive/MSU_miRNA_Proj/mRNA_Sternberger_MSU/trimGalored'

#the name of the file used in the process
meta_data_filename=$3
addtodir="/trimGalored/alignmentOutputs/subreadCounts"
r_dir=$2
subread_dir="${r_dir}${addtodir}" 
echo $subread_dir
vpub='/media/wyattlab/Data/DataDrive/MSU_miRNA_Proj/VpubSTAR' 
vpubGTF='/media/wyattlab/Data/DataDrive/MSU_miRNA_Proj/vpubescens_v1.1_maker_standard_no_FASTA.gtf'
arab='/media/wyattlab/Data/DataDrive/pilotPostTrim/AraportTair/Araport11'
arabGTF='/home/wyattlab/Desktop/AndrewRNAseq_Pipeline/Athaliana_447_Araport11.gene_exons.gtf'
mouse='/home/wyattlab/Desktop/AndrewRNAseq_Pipeline/MouseGenome/STARMouse_M20_GRCm38'
mouseGTF='/home/wyattlab/Desktop/AndrewRNAseq_Pipeline/MouseGenome/gencode.vM20.primary_assembly.annotation.gtf'

if [ $1 == "--help" ]
then
	echo "Help command executed."
	echo "==============================================="
	echo "Command Line Arguments:"
	echo "First Arg: desired species (viola, arabidopsis, mouse)"
	printf "\n"
	echo "Second Arg: desired working directory for WyattSeq"
	printf "\n"	
	echo "Third Arg: metadata file (.csv)"
	printf "\n"	
	echo "Additional Info:"	
	echo "The input for the R output directory is just the second argument followed by /trimGalored/alignmentOutputs/subreadCounts"
	printf "\n"
	echo "When entering metadata file name, use global address. Ex. /media/wyattlab/Data/DataDrive/mccall-kelly/190701-transfer/McCallMetadata.csv"
	printf "\n"
	echo "Do not capitalize species name in first argument."
	printf "\n\n"
	printf "Compatible Species:\n- Arabidopsis\n- Viola\n- Mouse\n"
	echo "==============================================="
fi

#loop to set genome and gtf files

if [ $1 == 'viola' ]
then
	echo "Genome set to violet."	
	genome=$vpub
	gtf=$vpubGTF

elif [ $1 == 'arabidopsis' ]
then
	echo "Genome set to arabidopsis."	
	genome=$arab
	gtf=$arabGTF

elif [ $1 == 'mouse' ]
then
	echo "Genome set to mouse."
	genome=$mouse
	gtf=$mouseGTF

elif [ $1 == '--help' ]
then
	exit

else
	echo Please specify the desired genome.
	echo "User provided: $1"
	exit
fi



echo $genome
echo $gtf


#loop for the source directory this script will run in

if [ $2 == "" ]
then
	source_dire='./'
else
	source_dire=$2
	echo "Directory set to $source_dire"
fi

echo "cd $source_dire"
cd $source_dire

#finding out if data is paired

echo -n "Is your read data paired (y/n)?"
read paired

if [ $paired == "y" ]
then
	echo "Running trimgalore-pairedEnd.sh"
#	sh /home/wyattlab/Desktop/AndrewRNAseq_Pipeline/trimgalore-pairedEnd.sh 
	cd ${source_dire}/trimGalored/

elif [ $paired == "n" ]
then
	echo "Running trimgalore-singleEnd.sh"
	sh /home/wyattlab/Desktop/AndrewRNAseq_Pipeline/trimgalore-singleEnd.sh
	cd ${source_dire}/trimGalored/

fi
echo -n "Current directory is: "
pwd


#moving into directory where trimmed files were moved to start running Star Alignment
#This will use the argument for source directory


#where genome goes
#running Star Alignment, may want to use this in case files are too big: --limitBAMsortRAM 6091545403
echo "Running StarAlign"
#sh /media/wyattlab/Data/DataDrive/MSU_miRNA_Proj/mRNA_Sternberger_MSU/StarAlign-General.sh $genome

#moving to directory alignmentOutputs
echo "cd alignmentOutputs"
cd alignmentOutputs

#where gtf is used
#To make the BamtoCounts script work correctly and for multi-species, we need to use the argument for
#species to set both a GTF file for use and the star formatted genome, the GTF will have to passed to 
#this next script
#running BamtoCount.sh
echo "Running BamtoCount.sh"
#sh /home/wyattlab/Desktop/AndrewRNAseq_Pipeline/BamtoCounts.sh $gtf

#running rmvFirstLine.sh
echo "Running rmvFirstLine.sh"
#sh /home/wyattlab/Desktop/AndrewRNAseq_Pipeline/rmvFirstLine.sh

echo "cd $subread_dir"
cd $subread_dir

#Running the R script
echo "Running PilotDE.R"
Rscript /home/wyattlab/Desktop/AndrewRNAseq_Pipeline/PilotDE_AutomatedRNAseq.R $meta_data_filename $subread_dir $1 --no-save #> RLog.txt

#Running second R script stored within new bash script (might not want to run second script every time)
#sh /home/wyattlab/Desktop/AndrewRNAseq_Pipeline/secondaryR.sh


