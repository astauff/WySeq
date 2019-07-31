library(edgeR)
library(ggplot2)
library(grid)
library(limma)
library(NMF)
library(plyr)
library(RColorBrewer)
library(reshape2)
library(scales)
library(gridExtra)
library(dplyr)
library(pheatmap)
library(tidyr)
library(cowplot)
library(gplots)
library(WGCNA)
options(stringsAsFactors = F)
allowWGCNAThreads()
library(genefilter)

#***WHERE WE USE THE ARGUMENTS TAKEN IN THE MAIN SCRIPT FROM COMMAND LINE***
args <- commandArgs()

#setting the working directory from command line argument
dir<-args[7]
setwd(dir)
#Printing the name of our working directory
print(paste("Working Directory: ",dir,sep = ""))

#Taking the command line argument for the name of the libsizes file (.csv)
meta_data_filename <- args[6]
#Printing the name of libsizes file to confirm this works
print(paste("Metadata File: ",meta_data_filename,sep = ""))

#Taking plant species from command line argument
plant_species <- args[8]
#Printing plant species
print(paste("Species: ",plant_species,sep = ""))
### Create CPM Matrix:

#Feed in an argument for the directory containing the counts files (merged or not, we can handle it either way)
#setwd("/media/wyattlab/Data/DataDrive/pilotPostTrim/BamFiles/MergedBamFiles")
#tweak

#working directory would be set to the directory containing the counts'
#dir <- "/home/wyattlab/Desktop/AndrewRNAseq_Pipeline/example_dir/alignmentOutputs/subreadCounts/"
#dir <- '/media/wyattlab/Data/DataDrive/MSU_miRNA_Proj/mRNA_Sternberger_MSU/example_dir/alignmentOutputs/subreadCounts'
#setwd("/media/wyattlab/Data/DataDrive/pilotPostTrim/BamFiles/MergedBamFiles/OutputCounts/")

#This will pull all files in the directory and merge them together by the column titled "ensembl_geneID"
#We may have to refine this function slightly based on the output you get from your test run.
#dataset<-NULL
if(exists("dataset")){
  rm(dataset)
}
file_list <- list.files(pattern = "\\.tab$")
i=1
for(file in file_list){
  if (exists("dataset")){
    temp_dataset <-read.table(file_list[i], header=TRUE, sep="\t")
    #This should need to be done but may have already been done for this test dataset outside of R.  We'll reinstate this upon the test run with the sternberger data
    temp_dataset<-temp_dataset[,-(2:6)]
    colnames(temp_dataset)<-c("ensembl_geneID",file_list[i])
    dataset<-merge(dataset,temp_dataset,by.x="ensembl_geneID",by.y="ensembl_geneID")
    rm(temp_dataset)
  }
  if(!exists("dataset")){
    dataset <- read.table(file_list[i], header=TRUE, sep="\t")
    #This should need to be done but may have already been done for this test dataset outside of R.  We'll reinstate this upon the test run with the sternberger data
    dataset<-dataset[,-(2:6)]
    colnames(dataset)<-c("ensembl_geneID",file_list[i])
  }

  # if the merged dataset does exist, append to it

  i=i+1
}

#Removed in the final version but can be handy if you don't want to keep importing and merging all the files in a directory.
postimport<-dataset
postimport->dataset

#RT = _1, MR = _2, HY = _3, CO = _4

rownames(dataset) <- dataset[,1]
dataset<-dataset[,-1]
#dataset<-dataset[,grep(pattern = "_4",x = colnames(dataset))]
colnames(dataset)<-sapply(strsplit(colnames(dataset),split = "_S"),'[', 1)

#excludes COL_T_RT4 (largest outlier of COL_T_RT) Might actually be COL_T_RT1 that is biggest outlier
#dataset$'14_1'<-NULL
#Excludes PGM_T_RT1 (largest outlier of PGM_T_RT)
#dataset$'1_1'<-NULL
#excludes COL_V_RT1 (largest outlier of COL_V_RT)
#dataset$'11_1'<-NULL
#Excludes PGM_V_RT1 (largest outlier of PGM_V_RT)
#dataset$'7_1'<-NULL
#Excludes PGM_V_RT4 (largest outlier of PGM_V_RT)
#dataset$'8_1'<-NULL

#makes final_data from data in dataset, excluding rows 1-5



#It doesn't seem like the alignment log stats are in the files, this means we probably don't need to run the line below.
#final_data<-dataset[-(1:5),]
final_data<-dataset

#rownames(final_data)<-strtrim(rownames(final_data),9) #Could be generalized
#creates count_data from data in final_data


count_data <- final_data

### MDS Plot:
#libsize was renamed to meta_data, lookout for libsize in following lines



#this whole metadata part will likely be read in as a file. Where we use file Ann gave us?
#not really sure that the two lines below are needed since we are reading meta_data in as a file
meta_data <- data.frame(colSums(count_data))
#meta_data$Tissue <- strtrim(colnames(count_data),4) #could be generalized




#think of new name for lib_sizes, "mata_datas" sounds goofy
#reading in the file for metadata, THIS WILL NEED TO BE PASSED AS AN ARGUMENT

#here is the path to the meta_data_filename: '/media/wyattlab/Data/DataDrive/MSU_miRNA_Proj/mRNA_Sternberger_MSU/TemplateForRNAseqPipeline_Sternberger.csv'
#here is the path until we figure out datadrive issue: '/home/wyattlab/Desktop/AndrewRNAseq_Pipeline/TemplateForRNAseqPipeline_Sternberger.csv'


lib_sizes<-read.csv(meta_data_filename, stringsAsFactors = FALSE,header = TRUE)
row.names(lib_sizes)<-sapply(strsplit(lib_sizes$forwardfile,split = "_S"),'[', 1)
lib_sizes["sample"] <- row.names(lib_sizes)
lib_sizes<-merge(lib_sizes,meta_data,by.y="row.names",by.x="sample")

#lib_sizes$sample <- NULL
row.names(lib_sizes)<-lib_sizes$sample
#lib_sizes$Row.names <- NULL
#lib_sizes$Tissue<-NULL

# 
#lib_sizes$block<-sapply(strsplit(lib_sizes$forwardfile,split = "_"),'[', 1) #could be generalized?
# lib_sizes$Tissue<-strtrim(lib_sizes$ExpID,8) #could be generalized

 
# 
# LOOK HERE
test<-as.data.frame(t(count_data))
test<-merge(lib_sizes,test,by.x="sample",by.y = "row.names")
row.names(test)<-test$sampleid
test<-test[,-(1:9)]
count_data<-as.data.frame(t(test))
# END LOOK HERE 


# #AT this point we can remove outliers as needed
# #count_data$COL_T_CO4<-NULL
# 
# #tweak
# 
# #creating model matrix
expDesign <- model.matrix(~0 + lib_sizes$group)

# #tweak, check for correct order
# #colnames(expDesign) <- unique(sort(strtrim(colnames(count_data),8)))
colnames(expDesign)<-sapply(strsplit(colnames(expDesign),split = "group"),'[', 2) #make sure Tissue is still in use

# 
# #tweak (check for correct order)

 
y <- DGEList(counts= count_data[rowSums(cpm(as.matrix(count_data)) > 1) >= 3,])
y<- calcNormFactors(y, method = "TMM")
# 
v<-voomWithQualityWeights(y,design = expDesign, method = "genebygene",plot=TRUE)


if(!all(isUnique(lib_sizes$block)) & length(unique(lib_sizes$block)) > 1){
b=(lib_sizes$block)
corfit<-duplicateCorrelation(v,expDesign,block = b)
vv<-voomWithQualityWeights(y,design=expDesign,plot = TRUE, block=b, correlation = corfit$consensus.correlation)
v<-vv
}

# 
mds <- plotMDS(v, ndim=3)
# 
mds.out <- as.data.frame(mds$cmdscale.out)
# 
mds.out$Tissue <- strtrim(colnames(count_data),9)
mds.out$Sample <- factor(rownames(lib_sizes), levels = rownames(lib_sizes))
# 
cols <- colorRampPalette(brewer.pal(7, "Paired"))
# 
ggplot(mds.out, aes(x=V1, y=V2, color=Tissue)) +
   geom_point(size=5) +
   ylim(-6,6) +
   xlim(-6,6) +
   ylab(bquote('Leading' ~Log[2]~ 'Fold Change Dim 2')) +
   xlab(bquote('Leading' ~Log[2]~ 'Fold Change Dim 1')) +
   theme_bw() +
   theme(
     text = element_text(size = 16)
   )
# 
 fit<-lmFit(v,expDesign)
# #fit<-lmFit(v,expDesign)

 
# #tweak, create contrast matrix
design.pairs<-function(levels){
  n <- length(levels)
  design <- matrix(0,n,choose(n,2))
  rownames(design) <- levels
  colnames(design) <- 1:choose(n,2)
  k <- 0
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      k <- k+1
      design[i,k] <- 1
      design[j,k] <- -1
      colnames(design)[k] <- paste(levels[i],"_v_",levels[j],sep="")
    }
  }
  design
}

#we can start again from here
levels=colnames(expDesign)
contMatrix<-design.pairs(levels = levels)
DE<-contrasts.fit(fit,contMatrix)
DE<-eBayes(DE)

#This part will need to be done separately for each species

if(plant_species == "viola"){
  Annotations <- read.csv(file = '/home/wyattlab/Desktop/AndrewRNAseq_Pipeline/VPannotations_w_TAIR_2.csv', header = TRUE, stringsAsFactors = FALSE) 
  row.names(Annotations)<-Annotations$Vpub_asc
}

if(plant_species == "arabidopsis"){
require(biomaRt)
require(xml2)
dbthale<-useEnsembl(biomart="plants_mart",dataset = "athaliana_eg_gene",host="plants.ensembl.org")
AthalianaAnnotations<-getBM(attributes = c('ensembl_gene_id','description','external_gene_name','entrezgene_id'),mart=dbthale)
Annotations<-AthalianaAnnotations
}

if(plant_species == "mouse"){
require(biomaRt)
require(xml2)
dbmus<-useMart("ensembl", dataset="mmusculus_gene_ensembl")
Atributes<-listAttributes(dbmus)
MouseAnnotations<-getBM(attributes=c('ensembl_gene_id','description','mgi_symbol','entrezgene_id'),mart=dbmus)
Annotations<-MouseAnnotations
}

#Add if exists comments before these rm commands to prevent warnings
if(exists("deframe")){
rm(deframe)
  }
if(exists("dframe")){
  rm(dframe)
  }
for(contrasts in colnames(contMatrix)){
  assign("temp",topTable(DE,coef =paste(contrasts),adjust.method = "BH",p.value = 1,number=length(DE$coefficients)))
  #print(paste(contrasts,"_DE",sep = ""))
  #print(c("The number of significantly DE genes for ",paste(contrasts),sum(topTable(DE,coef =paste(contrasts),adjust.method = "BH",p.value = 1,number=length(DE$coefficients))[5]<0.05)))
  #Add if exists command here to check if Annotations exists, and only run if that is true
  if(exists("Annotations")){
    temp<-merge(Annotations,temp,by.x="row.names",by.y="row.names")
    temp$Row.names <- NULL
    assign(paste(contrasts,"_DE",sep = ""),temp)
  }else{
    assign(paste(contrasts,"_DE",sep = ""),temp)
  }
  write.csv(temp,file=paste(dir,contrasts,".csv",sep = ""))
  if(exists("dframe")){
    deframe=data.frame(contrasts,sum(topTable(DE,coef =paste(contrasts),adjust.method = "BH",p.value = 1,number=length(DE$coefficients))[5]<0.05),sum(topTable(DE,coef =paste(contrasts),adjust.method = "BH",p.value = 1,number=length(DE$coefficients))[4]<0.05))
    names(deframe)=c("contrast","SigDEGenes_adj","SigDEGenes")
    dframe<-rbind(dframe,deframe)
  }
  if(!exists("dframe")){
    dframe=data.frame(contrasts,sum(topTable(DE,coef =paste(contrasts),adjust.method = "BH",p.value = 1,number=length(DE$coefficients))[5]<0.05),sum(topTable(DE,coef =paste(contrasts),adjust.method = "BH",p.value = 1,number=length(DE$coefficients))[4]<0.05))
    names(dframe)<-c("contrast","SigDEGenes_adj","SigDEGenes")
  }
}

save.image(paste(dir,"/","DE_workspace.Rdata",sep = ""))

# 
# 
# 
# 
# #contMatrixPvsC<-makeContrasts(dRT=(COL_T_RT-COL_V_RT)-(PGM_T_RT-PGM_V_RT), dHY=(COL_T_HY-COL_V_HY)-(PGM_T_HY-PGM_V_HY),dMR=(COL_T_MR-COL_V_MR)-(PGM_T_MR-PGM_V_MR),dCO=(COL_T_CO-COL_V_CO)-(PGM_T_CO-PGM_V_CO),levels=expDesign)
# #contMatrixPvsC<-makeContrasts(doubleDECOLvPGM=((COL_T_RT-COL_V_RT)-(PGM_T_RT-PGM_V_RT)),levels=expDesign)
# #DE<-contrasts.fit(fit,contMatrixPvsC)
# DE<-contrasts.fit(fit,contMatrix)
# 
# DE<-eBayes(DE)
# 

# 
# 
# #Tweak outputs

# 
# 
# COL_RT_TvV<-merge(AthalianaAnnotations,COL_T_RT_v_COL_V_RT_DE,by.x="ensembl_gene_id",by.y="row.names")
# COL_RT_TvV<-merge(COL_RT_DE, geneInfo, by.x="ensembl_gene_id",by.y="Gene")
# 
# PGM_RT_TvV<-merge(AthalianaAnnotations,PGM_T_RT_v_PGM_V_RT_DE,by.x="ensembl_gene_id",by.y="row.names")
# PGM_RT_TvV<-merge(PGM_RT_DE, geneInfo, by.x="ensembl_gene_id",by.y="Gene")
# 
# COL_v_PGM_DE_RT_V<-merge(AthalianaAnnotations,COL_V_RT_v_PGM_V_RT_DE,by.x="ensembl_gene_id",by.y="row.names")
# COL_v_PGM_DE_RT_V<-merge(PGM_v_COL_DE_RT_V, geneInfo, by.x="ensembl_gene_id",by.y="Gene")
# 
# 
# COL_v_PGM_DE_RT_T<-merge(AthalianaAnnotations,COL_T_RT_v_PGM_T_RT_DE,by.x="ensembl_gene_id",by.y="row.names")
# COL_v_PGM_DE_RT_T<-merge(PGM_v_COL_DE_RT_T, geneInfo, by.x="ensembl_gene_id",by.y="Gene")
# 
# doubleDECOLvPGM<-topTable(DE,coef = "doubleDECOLvPGM",adjust.method = "BH",p.value = 1,number=length(DE$coefficients))
# doubleDECOLvPGM<-merge(AthalianaAnnotations,doubleDECOLvPGM,by.x="ensembl_gene_id",by.y="row.names")
# write.csv(doubleDECOLvPGM,file = "/home/wyattlab/Dropbox/Colins_Stuff/Wyatt Lab/PilotDataOutputs/doubleDECOLvPGM_RT.csv",row.names = FALSE)
# 
# write.csv(COL_RT_TvV,file="/home/wyattlab/Dropbox/Colins_Stuff/Wyatt Lab/PilotDataOutputs/COL_10minTvVert_RT.csv",row.names = FALSE)
# 
# write.csv(PGM_RT_TvV,file="/home/wyattlab/Dropbox/Colins_Stuff/Wyatt Lab/PilotDataOutputs/PGM_10minTvVert_RT.csv",row.names = FALSE)
# 
# write.csv(COL_v_PGM_DE_RT_V,file="/home/wyattlab/Dropbox/Colins_Stuff/Wyatt Lab/PilotDataOutputs/Comp_COLvPGM_Vert_RT.csv",row.names = FALSE)
# 
# write.csv(COL_v_PGM_DE_RT_T,file="/home/wyattlab/Dropbox/Colins_Stuff/Wyatt Lab/PilotDataOutputs/Comp_COLvPGM_treat_RT.csv",row.names = FALSE)
# 
# #MDS plotting, at some point I wanna try this out
# mds <- plotMDS(v, ndim=3)
# 
# mds.out <- as.data.frame(mds$cmdscale.out)
# 
# mds.out$Tissue <- strtrim(colnames(count_data),9)
# mds.out$Sample <- factor(rownames(lib_sizes), levels = rownames(lib_sizes))
# 
# cols <- colorRampPalette(brewer.pal(7, "Paired"))
# 
# ggplot(mds.out, aes(x=V1, y=V2, color=Tissue)) +
#   geom_point(size=5) +
#   ylim(-6, 6) +
#   xlim(-6, 6) +
#   ylab(bquote('Leading' ~Log[2]~ 'Fold Change Dim 2')) +
#   xlab(bquote('Leading' ~Log[2]~ 'Fold Change Dim 1')) +
#   theme_bw() +
#   theme(
#     text = element_text(size = 16)
#   )
# 
# 
# WGCNA_data<-as.data.frame(t(v$E))
# rownames(WGCNA_data)
# 
# goodgenescheck = goodSamplesGenes(WGCNA_data, verbose = 3)
# goodgenescheck$allOK
# hclust_tree = hclust(dist(WGCNA_data), method = "average")
# 
# par(cex = 0.6)
# par(mar = c(0,4,2,0))
# plot(hclust_tree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
#      cex.axis = 1.5, cex.main = 2)
# 
# traits <- model.matrix(~0 + factor(lib_sizes$Tissue))
# colnames(traits) <- sapply(strsplit(colnames(traits),split = "Tissue)"),'[', 2)
# 
# hclust_tree2 = hclust(dist(WGCNA_data), method = "average")
# 
# traitColors = numbers2colors(traits, signed = FALSE, colors=blueWhiteRed(100))
# # Plot the sample dendrogram and the colors underneath.
# plotDendroAndColors(hclust_tree2, traitColors,
#                     groupLabels = names(traits),
#                     main = "Sample dendrogram and trait heatmap")
# 
# 
# powers = c(c(1:10), seq(from = 12, to=80, by=2))
# # Call the network topology analysis function
# sft = pickSoftThreshold(WGCNA_data, powerVector = powers, verbose = 5,blockSize = 20000)
# # Plot the results:
# par(mfrow = c(1,2))
# cex1 = 0.9
# 
# # Scale-free topology fit index as a function of the soft-thresholding power
# plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
#      main = paste("Scale independence"))
# text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      labels=powers,cex=cex1,col="red");
# # this line corresponds to using an R^2 cut-off of h
# abline(h=c(0.90,0.80),col=c("red","blue"))
# # Mean connectivity as a function of the soft-thresholding power
# plot(sft$fitIndices[,1], sft$fitIndices[,5],
#      xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
#      main = paste("Mean connectivity"))
# text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# 
# net = blockwiseModules(WGCNA_data, randomSeed=5678, power = 7,
#                        TOMType = "signed", minModuleSize = 30,
#                        reassignThreshold = 1e-6, mergeCutHeight = 0.15,
#                        numericLabels = TRUE, pamRespectsDendro = FALSE,
#                        saveTOMs = F,verbose = 3,deepSplit = 3,maxBlockSize = 20000)
# 
# 
# table(net$colors)
# 
# # Convert labels to colors for plotting
# colors = labels2colors(net$colors)
# # Plot the dendrogram and the module colors underneath
# plotDendroAndColors(net$dendrograms[[1]], colors[net$blockGenes[[1]]],
#                     "Module colors",
#                     dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05)
# table(colors)
# 
# genes = ncol(WGCNA_data)
# samples = nrow(WGCNA_data)
# # Recalculate MEs with color labels
# MEs0 = moduleEigengenes(WGCNA_data, colors)$eigengenes
# MEs = orderMEs(MEs0)
# moduleTraitCor = cor(MEs, traits, use = "p");
# moduleTraitPvalue = corPvalueStudent(moduleTraitCor, samples)
# 
# textMatrix = paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "");
# dim(textMatrix) = dim(moduleTraitCor)
# par(mar = c(6, 8.5, 3, 3));
# # Display the correlation values within a heatmap plot
# labeledHeatmap(Matrix = moduleTraitCor,
#                xLabels = colnames(traits),
#                yLabels = substring(names(MEs),3),
#                ySymbols = names(MEs),
#                colorLabels = TRUE,
#                colors = blueWhiteRed(50),
#                textMatrix = textMatrix,
#                setStdMargins = FALSE,
#                cex.text = 0.75,
#                zlim = c(-.7,.7),
#                main = paste("Module-trait relationships"))
# 
# geneInfo=data.frame(Gene=colnames(WGCNA_data),Module=colors)
# 
# 
# write.csv(geneInfo,file="/home/wyattlab/Dropbox/Colins_Stuff/Wyatt Lab/PilotDataOutputs/GenesInModules_2ndWGCNA_RT.csv")
