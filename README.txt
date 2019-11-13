=== WySeq ===

Developed in The Wyatt Lab. Ohio University, Athens, OH.
July, 2019.

Contributors: Colin Kruse, Andrew Stauffer

Version: 1.2

Usage Notes
===========
- WySeq is an automated RNAseq pipeline that utilizes the StarAlign and TrimGalore packages.

- The program uses a mixture of bash scripts, Python and R.

- To function properly, WySeq takes in three command line arguments. Those arguments are, in order: species, working directory and metadata file.

- Users must change pathways to .gtf and genome files at the top of the main script.

- The help arguments for WySeq is --help as the first argument in the command line.

Required packages
=================
- StarAlign
- TrimGalore

Example Call of WyattSeq
========================

./wyseq.sh mouse '/media/wyattlab/Data/DataDrive/mccall-kelly/190701-transfer'  /media/wyattlab/Data/DataDrive/mccall-kelly/190701-transfer/McCallMetadata.csv


*** WySeq only works for Arabidopsis, Viola and Mouse ***
