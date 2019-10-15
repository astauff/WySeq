#This file contains all three GO functions needed for WyattSeq and any supplementary functions 
suppressMessages(require(gage))
suppressMessages(require(biomaRt))
suppressMessages(require(xml2))

suppressMessages(require(gage))
suppressMessages(require(gdata))
suppressMessages(require(ggplot2))
suppressMessages(require(ParamHelpers))
suppressMessages(require(biomaRt))
suppressMessages(require(xml2))
suppressMessages(require(gage))
suppressMessages(require(pathview))

trim<-function (x) 
{
  gsub("^\\s+|\\s+$", "", x)
}
#biological process enrichment function will use the next two functions
wrap.it <- function(x, len)
{ 
  sapply(x, function(y) paste(strwrap(y, len), 
                              collapse = "\n"), 
         USE.NAMES = FALSE)
}


# Call this function with a list or vector
wrap.labels <- function(x, len)
{
  if (is.list(x))
  {
    lapply(x, wrap.it, len)
  } else {
    wrap.it(x, len)
  }
}

Arabidopsis_GOBP<-function(x, out.prefix = "GO_Biological_Process",wraplen=35){
  goDB<-go.gsets(species = "arabidopsis", id.type = "entrez")
  print(paste("Beginning ontological analysis for ", out.prefix,".",sep=""))
  prev_dir <- getwd()
  exp.fc <- x$logFC
  if(!exists("exp.fc")){
    print("exp.fc does not exist")
    print(x$logFC)
  }
  names(exp.fc) <- x$ensembl_gene_id
  fc.go.p <- gage(exp.fc, gsets = goDB$go.sets[goDB$go.subs$BP], 
                  ref = NULL, samp = NULL)
  greater <- as.data.frame(fc.go.p$greater)
  greater$Term_ID <- trim(rownames(greater))
  less <- as.data.frame(fc.go.p$less)
  less$Term_ID <- trim(rownames(less))
  GO_out <- merge(greater, less, by = "Term_ID")
  names(GO_out) <- c("Term_ID", "p.geomean.greater", "stat.mean.greater", 
                     "p.val.greater", "q.val.greater", "set.size.greater", 
                     "exp1.greater", "p.geomean.less", "stat.mean.less", 
                     "p.val.less", "q.val.less", "set.size.less", "exp1.less")
  GO_out$p.value <- ifelse(GO_out$stat.mean.greater > 0, GO_out$p.val.greater, 
                           GO_out$p.val.less)
  GO_out$FDR <- ifelse(GO_out$stat.mean.greater > 0, GO_out$q.val.greater, 
                       GO_out$q.val.less)
  GO_out <- GO_out[, c(1, 6, 3, 14, 15)]
  names(GO_out) <- c("Term_ID", "set.size", "mean_logFC", 
                     "P.value", "FDR")
  #**********this next line causes all the data to be deleted from GO_out**********#
  GO_out <- GO_out[!is.na(GO_out$P.value), ]
  row.names(GO_out) <- GO_out$Term_ID
  GO_out <- GO_out[order(GO_out$P.value), ]
  #GO_out$Genes <- "NULL"
  x<-x[order(abs(x$logFC),decreasing = TRUE),]
  for (i in GO_out$Term_ID) {
    GO_entrezgenes <- goDB$go.sets[[i]]
    GO_out[i, "Genes"] <- paste(x[x$ensembl_gene_id %in% GO_entrezgenes,"external_gene_name"], collapse = ",")
  }
  if(length(GO_out$P.value)>0){
    GO_out$plotval<- -log(GO_out$P.value,base = 10)
    #return(GO_out)
    print(paste(out.prefix,"BP.tiff",sep = ""))
    tiff(filename = paste(out.prefix,"BP.tiff",sep = ""), height = 4, width = 5.5, units = "in", res = 400, compression = "lzw")
    
    par(mai = c(0.45, 2.2, 0.08, 0.08))
    plot.new()
    lowx = 0
    highx = max(GO_out$plotval)
    #More specific GOs only
    #GO_out<-GO_out[GO_out$set.size <= 300,]
    number=min(c(10,nrow(GO_out[GO_out$P.value<=1.00,])))
    #number=min(nrow(GO_out))
    print(lowx)
    print(highx)
    plot.window(xlim = c(lowx, highx*1.04), ylim = c(1 - 0.20, number + 0.20), xaxs = "i")
    
    # This part adds little horizontal lines
    abline(v = 0, col = 1, lty = 1, lwd = 2)
    
    for (i in 1:number) {
      abline(h = i, lty = 2, lwd = 0.80, col = "grey70") }
    # This part generates the bars. It iterates through each row in object "GO_out"
    vsize = 0.37
    lt = par()$usr[1]; rt = par()$usr[2]
    for (i in 1:number) {
      theP = GO_out[i, ]$plotval
      thelab = GO_out[i, ]$Term_ID
      thelab = as.character(unlist(strsplit(thelab, split = "\ ", fixed = T)))
      thelab = thelab[2:length(thelab)]
      thelab = paste(thelab, collapse = " ")
      thegenes = GO_out[i, ]$Genes
      thegenes = as.character(unlist(strsplit(thegenes, split = ",", fixed = T)))
      take = min(length(thegenes), 4)
      thegenes = thegenes[1:take]
      thegenes = paste(thegenes, collapse = ", ")
      #Male Colors
      #if (i %% 2 == 0) thecol = "#e8a9e5" # choose bar color here
      #if (i %% 2 != 0) thecol = "#dd1c77" # choose bar color here
      #Female colors
      if (i %% 2 == 0) thecol = "#fff5a3" # choose bar color here
      if (i %% 2 != 0) thecol = "#e2cc00" # choose bar color here
      if (i %% 2 == 0) acol = "black" # choose bar color here
      if (i %% 2 != 0) acol = "midnightblue" # choose bar color here
      rect(0, i - vsize, theP, i + vsize, col = thecol, border = 1, lwd = 1.25)
      thelab<-wrap.labels(thelab,wraplen)
      axis(2, at = i, labels = thelab, cex.axis = 0.7, las = 1, tick=F, line = -0.80, font=2,col.axis = acol)
      text(rt, i, thegenes, cex = 0.6, font = 4, pos = 2, col = "black", offset = 0.10) 
    }
    # This part labels the bottom axis and margin
    box(lwd = 3)
    axis(1, line = -0.90, lwd = 0, tick = F, col = "white", cex.axis = 0.80, font.axis = 2)
    axis(1, labels = FALSE, tick = T, tcl = -0.2, lwd.ticks = 3)
    mtext("-log(P-value)", side = 1, line = 1.10, cex = .6, font = 2)
    
    # This adds the bottom-left label
    par(mai = c(0, 0, 0, 0), new = T)
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))
    text(0.14, 0.02, paste("(",wrap.labels(out.prefix,20),")",sep = ""), font = 2, cex = 1,col = "red")
    # May also need to add your setwd statement here.
    graphics.off()
  }
  return(GO_out)
}

Vpub_GOBP<-function (x, out.prefix = "GO_Biological_Process",wraplen=35){
  goDB<-go.gsets(species = "arabidopsis", id.type = "entrez")
  prev_dir <- getwd()
  exp.fc <- x$logFC
  x$TAIR_ID<-sapply(strsplit(x$Arab_asc,split = "\\."),'[', 1)
  names(exp.fc) <- x$TAIR_ID
  fc.go.p <- gage(exp.fc, gsets = goDB$go.sets[goDB$go.subs$BP], 
                  ref = NULL, samp = NULL)
  greater <- as.data.frame(fc.go.p$greater)
  greater$Term_ID <- trim(rownames(greater))
  less <- as.data.frame(fc.go.p$less)
  less$Term_ID <- trim(rownames(less))
  GO_out <- merge(greater, less, by = "Term_ID")
  names(GO_out) <- c("Term_ID", "p.geomean.greater", "stat.mean.greater", 
                     "p.val.greater", "q.val.greater", "set.size.greater", 
                     "exp1.greater", "p.geomean.less", "stat.mean.less", 
                     "p.val.less", "q.val.less", "set.size.less", "exp1.less")
  GO_out$p.value <- ifelse(GO_out$stat.mean.greater > 0, GO_out$p.val.greater, 
                           GO_out$p.val.less)
  GO_out$FDR <- ifelse(GO_out$stat.mean.greater > 0, GO_out$q.val.greater, 
                       GO_out$q.val.less)
  GO_out <- GO_out[, c(1, 6, 3, 14, 15)]
  names(GO_out) <- c("Term_ID", "set.size", "mean_logFC", 
                     "P.value", "FDR")
  #**********this next line causes all the data to be deleted from GO_out**********#
  GO_out <- GO_out[!is.na(GO_out$P.value), ]
  row.names(GO_out) <- GO_out$Term_ID
  GO_out <- GO_out[order(GO_out$P.value), ]
  GO_out$Genes <- "NULL"
  x<-x[order(abs(x$logFC),decreasing = TRUE),]
  for (i in GO_out$Term_ID) {
    GO_entrezgenes <- goDB$go.sets[[i]]
    GO_out[i, "Genes"] <- paste(row.names(x)[row.names(x) %in% GO_entrezgenes], collapse = ",")
  }
  GO_out$plotval<- -log(GO_out$P.value,base = 10)
  #return(GO_out)
  print(paste(out.prefix,"BP.tiff",sep = ""))
  tiff(filename = paste(out.prefix,"BP.tiff",sep = ""), height = 4, width = 5.5, units = "in", res = 400, compression = "lzw")
  
  par(mai = c(0.45, 2.2, 0.08, 0.08))
  plot.new()
  lowx = 0
  highx = max(GO_out$plotval)
  #More specific GOs only
  #GO_out<-GO_out[GO_out$set.size <= 300,]
  number=min(c(10,nrow(GO_out[GO_out$P.value<=1.00,])))
  #number=min(nrow(GO_out))
  print(lowx)
  print(highx)
  plot.window(xlim = c(lowx, highx*1.04), ylim = c(1 - 0.20, number + 0.20), xaxs = "i")
  
  # This part adds little horizontal lines
  abline(v = 0, col = 1, lty = 1, lwd = 2)
  
  for (i in 1:number) {
    abline(h = i, lty = 2, lwd = 0.80, col = "grey70") }
  # This part generates the bars. It iterates through each row in object "GO_out"
  vsize = 0.37
  lt = par()$usr[1]; rt = par()$usr[2]
  for (i in 1:number) {
    theP = GO_out[i, ]$plotval
    thelab = GO_out[i, ]$Term_ID
    thelab = as.character(unlist(strsplit(thelab, split = "\ ", fixed = T)))
    thelab = thelab[2:length(thelab)]
    thelab = paste(thelab, collapse = " ")
    thegenes = GO_out[i, ]$Genes
    thegenes = as.character(unlist(strsplit(thegenes, split = ",", fixed = T)))
    take = min(length(thegenes), 4)
    thegenes = thegenes[1:take]
    thegenes = paste(thegenes, collapse = ", ")
    #Male Colors
    #if (i %% 2 == 0) thecol = "#e8a9e5" # choose bar color here
    #if (i %% 2 != 0) thecol = "#dd1c77" # choose bar color here
    #Female colors
    if (i %% 2 == 0) thecol = "#fff5a3" # choose bar color here
    if (i %% 2 != 0) thecol = "#e2cc00" # choose bar color here
    if (i %% 2 == 0) acol = "black" # choose bar color here
    if (i %% 2 != 0) acol = "midnightblue" # choose bar color here
    rect(0, i - vsize, theP, i + vsize, col = thecol, border = 1, lwd = 1.25)
    thelab<-wrap.labels(thelab,wraplen)
    axis(2, at = i, labels = thelab, cex.axis = 0.7, las = 1, tick=F, line = -0.80, font=2,col.axis = acol)
    text(rt, i, thegenes, cex = 0.6, font = 4, pos = 2, col = "black", offset = 0.10) 
  }
  # This part labels the bottom axis and margin
  box(lwd = 3)
  axis(1, line = -0.90, lwd = 0, tick = F, col = "white", cex.axis = 0.80, font.axis = 2)
  axis(1, labels = FALSE, tick = T, tcl = -0.2, lwd.ticks = 3)
  mtext("-log(P-value)", side = 1, line = 1.10, cex = .6, font = 2)
  
  # This adds the bottom-left label
  par(mai = c(0, 0, 0, 0), new = T)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  text(0.14, 0.02, paste("(",wrap.labels(out.prefix,20),")",sep = ""), font = 2, cex = 1,col = "red")
  # May also need to add your setwd statement here.
  graphics.off()
  return(GO_out)
}

mouse_GOBP<-function (x, out.prefix = "GO_Biological_Process",wraplen=35) 
{
  goDB<-go.gsets(species = "mouse", id.type = "entrez")
  prev_dir <- getwd()
  exp.fc <- x$logFC
  names(exp.fc) <- x$entrezgene
  fc.go.p <- gage(exp.fc, gsets = goDB$go.sets[goDB$go.subs$BP], 
                  ref = NULL, samp = NULL)
  greater <- as.data.frame(fc.go.p$greater)
  greater$Term_ID <- trim(rownames(greater))
  less <- as.data.frame(fc.go.p$less)
  less$Term_ID <- trim(rownames(less))
  GO_out <- merge(greater, less, by = "Term_ID")
  names(GO_out) <- c("Term_ID", "p.geomean.greater", "stat.mean.greater", 
                     "p.val.greater", "q.val.greater", "set.size.greater", 
                     "exp1.greater", "p.geomean.less", "stat.mean.less", 
                     "p.val.less", "q.val.less", "set.size.less", "exp1.less")
  GO_out$p.value <- ifelse(GO_out$stat.mean.greater > 0, GO_out$p.val.greater, 
                           GO_out$p.val.less)
  GO_out$FDR <- ifelse(GO_out$stat.mean.greater > 0, GO_out$q.val.greater, 
                       GO_out$q.val.less)
  GO_out <- GO_out[, c(1, 6, 3, 14, 15)]
  names(GO_out) <- c("Term_ID", "set.size", "mean_logFC", 
                     "P.value", "FDR")
  GO_out <- GO_out[!is.na(GO_out$P.value), ]
  row.names(GO_out) <- GO_out$Term_ID
  GO_out <- GO_out[order(GO_out$P.value), ]
  GO_out$Genes <- "NULL"
  x<-x[order(abs(x$logFC),decreasing = TRUE),]
  for (i in GO_out$Term_ID) {
    GO_entrezgenes <- goDB$go.sets[[i]]
    GO_out[i, "Genes"] <- paste(x[x$entrezgene_id %in% GO_entrezgenes, 3], collapse = ",")
  }
  GO_out$plotval<- -log(GO_out$P.value,base = 10)
  #return(GO_out)
  print(paste(out.prefix,"BP.tiff",sep = ""))
  tiff(filename = paste(out.prefix,"BP.tiff",sep = ""), height = 4, width = 5.5, units = "in", res = 400, compression = "lzw")
  
  par(mai = c(0.45, 2.2, 0.08, 0.08))
  plot.new()
  lowx = 0
  highx = max(GO_out$plotval)
  #More specific GOs only
  #GO_out<-GO_out[GO_out$set.size <= 300,]
  number=min(c(10,nrow(GO_out[GO_out$P.value<=1.00,])))
  #number=min(nrow(GO_out))
  print(lowx)
  print(highx)
  plot.window(xlim = c(lowx, highx*1.04), ylim = c(1 - 0.20, number + 0.20), xaxs = "i")
  
  # This part adds little horizontal lines
  abline(v = 0, col = 1, lty = 1, lwd = 2)
  
  for (i in 1:number) {
    abline(h = i, lty = 2, lwd = 0.80, col = "grey70") }
  # This part generates the bars. It iterates through each row in object "GO_out"
  vsize = 0.37
  lt = par()$usr[1]; rt = par()$usr[2]
  for (i in 1:number) {
    theP = GO_out[i, ]$plotval
    thelab = GO_out[i, ]$Term_ID
    thelab = as.character(unlist(strsplit(thelab, split = "\ ", fixed = T)))
    thelab = thelab[2:length(thelab)]
    thelab = paste(thelab, collapse = " ")
    thegenes = GO_out[i, ]$Genes
    thegenes = as.character(unlist(strsplit(thegenes, split = ",", fixed = T)))
    take = min(length(thegenes), 4)
    thegenes = thegenes[1:take]
    thegenes = paste(thegenes, collapse = ", ")
    #Male Colors
    #if (i %% 2 == 0) thecol = "#e8a9e5" # choose bar color here
    #if (i %% 2 != 0) thecol = "#dd1c77" # choose bar color here
    #Female colors
    if (i %% 2 == 0) thecol = "#fff5a3" # choose bar color here
    if (i %% 2 != 0) thecol = "#e2cc00" # choose bar color here
    if (i %% 2 == 0) acol = "black" # choose bar color here
    if (i %% 2 != 0) acol = "midnightblue" # choose bar color here
    rect(0, i - vsize, theP, i + vsize, col = thecol, border = 1, lwd = 1.25)
    thelab<-wrap.labels(thelab,wraplen)
    axis(2, at = i, labels = thelab, cex.axis = 0.7, las = 1, tick=F, line = -0.80, font=2,col.axis = acol)
    text(rt, i, thegenes, cex = 0.6, font = 4, pos = 2, col = "black", offset = 0.10) 
  }
  # This part labels the bottom axis and margin
  box(lwd = 3)
  axis(1, line = -0.90, lwd = 0, tick = F, col = "white", cex.axis = 0.80, font.axis = 2)
  axis(1, labels = FALSE, tick = T, tcl = -0.2, lwd.ticks = 3)
  mtext("-log(P-value)", side = 1, line = 1.10, cex = .6, font = 2)
  
  # This adds the bottom-left label
  par(mai = c(0, 0, 0, 0), new = T)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  text(0.14, 0.02, paste("(",wrap.labels(out.prefix,20),")",sep = ""), font = 2, cex = 1,col = "red")
  # May also need to add your setwd statement here.
  graphics.off()
  return(GO_out)
}


GO_Wrapper<-function(x, out.prefix = "GO_Biological_Process",wraplen=35,SpeciesInput=plant_species){
  #Add an if statement for each of the three organisms
  if(SpeciesInput == "mouse"){
    GO_out<-mouse_GOBP(x, out.prefix, wraplen)
  }
  
  if(SpeciesInput == "arabidopsis"){
    GO_out<-Arabidopsis_GOBP(x, out.prefix, wraplen)
  }
    
  if(SpeciesInput == "viola"){
    GO_out<-Vpub_GOBP(x, out.prefix, wraplen)
  }
    #store each output to GO_out
    return(GO_out)
}


PathwayD<-function (x, SpeciesInput=plant_species, out.suffix="Kegg") {
  print(paste("Beginning pathways analysis for ", out.suffix,".",sep=""))
  prev_dir <- getwd()
  out.suffix <- gsub(" ", "", out.suffix)
  new_dir <- paste(out.suffix, "KEGGPathways", sep = "_")
  if(!exists(new_dir)){
    dir.create(new_dir)
  }
  setwd(new_dir)
  if(SpeciesInput == "mouse"){
    keggDB<-kegg.gsets(species = "mouse", id.type = "entrez")
    species="mmu"
  }
  
  if(SpeciesInput == "arabidopsis"){
    keggDB<-kegg.gsets(species = "ath", id.type = "entrez")
    species="ath"
  }
  
  if(SpeciesInput == "viola"){
    keggDB<-kegg.gsets(species = "ath", id.type = "entrez")
    species="ath"
  }  
  exp.fc <- x$logFC
  names(exp.fc) <- x$entrezgene_id
  #the two lines from above
  #keggDB<-kegg.gsets(species = "mmu", id.type = "entrez")
  #keggDB<-keggDB$kg.sets
  #end of lines
  fc.kegg.p <- gage(exp.fc, gsets = keggDB$kg.sets, ref = NULL, samp = NULL)
  greater <- as.data.frame(fc.kegg.p$greater)
  greater$Term_ID <- trim(rownames(greater))
  less <- as.data.frame(fc.kegg.p$less)
  less$Term_ID <- trim(rownames(less))
  path_out <- merge(greater, less, by = "Term_ID")
  names(path_out) <- c("Term_ID", "p.geomean.greater", "stat.mean.greater", 
                       "p.val.greater", "q.val.greater", "set.size.greater", 
                       "exp1.greater", "p.geomean.less", "stat.mean.less", 
                       "p.val.less", "q.val.less", "set.size.less", "exp1.less")
  path_out$p.value <- ifelse(path_out$stat.mean.greater > 
                               0, path_out$p.val.greater, path_out$p.val.less)
  path_out$FDR <- ifelse(path_out$stat.mean.greater > 0, path_out$q.val.greater, 
                         path_out$q.val.less)
  path_out <- path_out[, c(1, 6, 3, 14, 15)]
  names(path_out) <- c("Term_ID", "set.size", "mean_logFC", 
                       "P.value", "FDR")
  path_out <- path_out[!is.na(path_out$P.value), ]
  row.names(path_out) <- path_out$Term_ID
  path_out <- path_out[order(path_out$P.value), ]
  path_out[path_out$Term_ID==""]
  #path_out$Genes <- "NULL"
  for (i in path_out$Term_ID) {
    path_entrezgenes <- keggDB$kg.sets[[i]]
    path_out[i, "Genes"] <- paste(x[x$entrezgene %in% path_entrezgenes, 3], collapse = ",")
  }
  write.table(path_out, file = paste(out.suffix, "_results.xls", 
                                     sep = "_"), sep = "\t", quote = F, row.names = F)
  # try(plotSigTerms_barplot(results = path_out, order = "top", 
  #                          num = 25, name = "Top25_significant_singleDirection_KEGG_barplot"))
  # try(plotSigTerms_barplot(results = path_out, order = "upregulated", 
  #                          num = 25, name = "Top25_upregulated_significant_singleDirection_KEGG_barplot"))
  # try(plotSigTerms_barplot(results = path_out, order = "downregulated", 
  #                          num = 25, name = "Top25_downregulated_significant_singleDirection_KEGG_barplot"))
  
  if(length(path_out$Term_ID[path_out$FDR<=0.05 & path_out$Term_ID != "ath01100 Metabolic pathways"])>0){
    for (pid in path_out$Term_ID[path_out$FDR<=0.05 & path_out$Term_ID != "ath01100 Metabolic pathways" & path_out$Term_ID != "mmu01100 Metabolic pathways"]) {
      pathway.id <- substr(pid, 1, 8)
      term_name <- gsub("\\-","",gsub("/", ":", gsub(" ", "_", substr(pid,10, nchar(pid)))))
      print(paste("Trying pathview command for loop on ", pid))
      try(pathview(gene.data = exp.fc, pathway.id = pathway.id, 
                   species = species, out.suffix = paste(term_name,out.suffix, sep = "."), limit = list(gene = 2, cpd = 2), low = list(gene = "blue", cpd = "blue"), 
                   mid = list(gene = "gray", cpd = "gray"), high = list(gene = "orange",cpd = "orange")),silent=FALSE)
      print(paste("Done pathview command for loop on .", pid))
      pid_entrezgenes <- keggDB$kg.sets[[pid]]
      pid_table <- x[x$entrezgene %in% pid_entrezgenes, ]
      write.table(pid_table, file = paste(pathway.id, term_name,out.suffix, "Differential_Expression_subset.xls", sep = "."), sep = "\t", quote = F, row.names = F)
      file.remove(paste(pathway.id, ".png", sep = ""))
      file.remove(paste(pathway.id, ".xml", sep = ""))
    }
  }
  setwd(prev_dir)
  print(paste("Completed pathway analysis for ", out.suffix,".",sep=""))
  if(length(path_out$Term_ID[path_out$FDR<=0.05])>0){
    return(path_out)
  }else{
    return(NULL)
  }
  
}
