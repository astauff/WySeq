=== WyattSeq ===

Developed in The Wyatt Lab. Ohio University, Athens, OH.
July, 2019.

Contributors: Colin Kruse, Andrew Stauffer

Version: 1.0

Usage Notes
===========
- WyattSeq is an automated RNAseq pipeline that utilizes the StarAlign and TrimGalore packages.

- The program uses a mixture of bash scripts, Python and R.

- To function properly, WyattSeq takes in four command line arguments. Those arguments are, in order: species, working directory, metadata file and subread count output directory.

- Users must change pathways to .gtf and genome files at the top of the main script.

- The help arguments for WyattSeq is --help as the first argument in the command line.

Required packages
=================
- StarAlign
- TrimGalore

Example Call of WyattSeq
========================

./wyattseq.sh mouse '/media/wyattlab/Data/DataDrive/mccall-kelly/190701-transfer' '/media/wyattlab/Data/DataDrive/mccall-kelly/190701-transfer/McCallMetadata.csv' '/media/wyattlab/Data/DataDrive/mccall-kelly/190701-transfer/trimGalored/alignmentOutputs/subreadCounts'

*** WyattSeq only works for Arabidopsis, Viola and Mouse ***
