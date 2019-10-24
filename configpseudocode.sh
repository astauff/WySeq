#Check that trim_galore, cutadapt, R, STAR, and subread can be found in the path

#Get the current directory that WyattSeq is in

#add command to export the directory to path at the end of the .bashrc profile 
#(might check that it exists first, then give the user the option to add to it)
export PATH="DIR:$PATH"
#cat a source command of the GO_BP functions then append the base GO_BP functions without the previous source command in it
echo source($DIR/GO_BP_Functions.r) > PilotDE_automated.r
cat ASourceScriptofPilotDEMinusSourceCommand.r >> PilotDE_automated.r
#NOTE: this is so that users don't have to deal with changing anything in the Rscript

#Run Rscript to install all R packages required
#I've added the basics of that to github


