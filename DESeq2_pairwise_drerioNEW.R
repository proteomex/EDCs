########################################################################
### DESeq2 RNA Seq Analysis Script for multiple factors / treatments ###
########################################################################

# Author: Hannes Reinwald
# Contact: hannes.reinwald@ime.fraunhofer.de

### README ### -----------------------------------------------------------
# This script is desgined to run DESeq2 normalization and statistical testing on RNAseq experiments
# in Ecotox testing for multiple concentrations. This script will automatically adopt to the k numbers
# of samples and N numbers of conditions (High,Mid,Low or Treat1, Treat2, Treat3, ... , Control) and 
# will generate Log2FC values with respect to control groups.

# Requiered Input: 
# - CountMatrix (txt or csv is fine, if csv read.csv2 is used)
# - coldata.csv

# Main Analysis Steps are: 
# 1) RLE normalization (DESeq2) across all samples element of CountMatrix
# 2) Calc LFC (with respect to control) and LFcutOff (upper 90% quantile of abs(LFC values))
# 3) apeglm shrinking on LFC

# 4) Multiple t-testing with Benjamin-Hochberg correction (padj < 0.05) 
#    and independet hypothesis weighing (IHW) to identify DEGs
#    for LFC and apeglm(LFC) values for H0: LFC = 0; apeglm(LFC) = 0
# 4.1) Annotation of Genelists via org.Dr.eg.db [https://www.bioconductor.org/packages/release/data/annotation/html/org.Dr.eg.db.html]
# 4.2) DEGs with LFC / apeglm(LFC) > < abs(LFcutOff) are kept as potential molecular marker features

# 5) Plotting
# 5.1 MA- & Vulcano plots for N-1 conditions for LFC [LFcut & padj cut]
# 5.2 MA- & Vulcano plots for N-1 conditions for apeglm(LFC) [LFcut & padj cut]
# 5.3 DataQC
#     - Correlation bio. replicates (rlog(meanCounts))
#     - MeanSD plots (rlog(meanCounts))
#     - RLE_normalisation
#     - normCounts_rlogTransformation

# 6) Output: 
# 6.1 DESeq2 Results (xcel; still working on html output though ... )
#     - N-1 DESeq2 result files for LFC
#     - N-1 DESeq2 result files for apeglm(LFC)
# 6.2 DEGs (padj < 0.05 and LFcut for: 
#     - N-1 files for LFC
#     - N-1 files for apeglm(LFC)
##############

### LOAD PACKAGES ### -------------------------------------------------
require(DESeq2) #installed
require(IHW)    #installed
require(apeglm) #installed
require(qvalue)
#require(locfdr)
#require(RColorBrewer)
require(ggplot2)
require(ggpubr)
require(hexbin)
#####################

### DATA IMPORT & Color settings ### ---------------------------------------------------
message("\nImporting CountMatrix and coldata files ...")

## coldata import ##
tmp = list.files(full.names = T, pattern = "[Cc]oldata")
if(grepl(".csv",tmp)==T){
  coldata = read.csv2(file = tmp, row.names = 1, header = T)
  # In case csv file was not differently exported run read.csv
  if (ncol(coldata) <= 1) {
    coldata = read.csv(file = tmp, row.names = 1, header = T)
  }
} else {
  # if not ending with csv import with read.delim2
  coldata = read.delim2(file = tmp, row.names = 1, header = T)
  if (ncol(coldata) <= 1) {
    # if file was differently formated use read.delim
    coldata = read.delim(file = tmp, row.names = 1, header = T)
  }
}
coldata <- droplevels(coldata[which(coldata$Condition != ""), # filtering for non empty rows
                              which(!(grepl("X.",colnames(coldata))) == T)]) # filtering cols
coldata$Condition <- gsub(" ","",coldata$Condition) #removing any spaces in Condition levels
coldata$Tank <- gsub("-","_",coldata$Tank)

# set factor levels
coldata$Tank      <- factor(coldata$Tank)
coldata$Condition <- factor(coldata$Condition)
coldata$Substance <- factor(coldata$Substance)
substance         <- levels(coldata$Substance) #Extract the name(s) of the tested substance(s)

# relevel factors in Condition in right order:
condition = levels(coldata$Condition)
tmp = as.character(coldata$Condition[grep("[Cc]ontrol",coldata$Condition)][1]) #extract control factor level
if(all(grepl("[Ee]xposure",condition[2:4]))==T & length(condition)==4) {# Control, LE, ME, HE
  coldata$Condition <- factor(coldata$Condition, levels = condition[c(1,3,4,2)])
} else if (all(grepl("[Ee]xposure",condition[2:3]))==T & length(condition)==3) { # Control, LE, HE
  coldata$Condition <- factor(coldata$Condition, levels = condition[c(1,3,2)])
} else { # just specifying the reference level --> Control
  coldata$Condition <- relevel(coldata$Condition, ref = tmp)
}
message(paste0("Reference level:\t",levels(coldata$Condition)[1]),"\nSorting ",length(condition)-1," Treatments:\t",paste(levels(coldata$Condition)[-1], collapse = " "))
condition <- levels(coldata$Condition)

## CountMatrix import ##
# Import all provided count mtx files and merge them together;
# then filter according to the provided sample ID in coldata
importCountMtx = function(file){
  if(any(grepl(".txt",file))==T) { # import txt
    tmp1 = file[grepl(".txt",file)]
    count.matrix = read.table(file = tmp1, row.names = 1, header = T)
    if(ncol(count.matrix) <= 1) {
      count.matrix = read.delim2(file = tmp1, row.names = 1, header = T)
      if(ncol(count.matrix) <= 1) {
        count.matrix = read.delim(file = tmp1, row.names = 1, header = T)
      }
    }
  } else { # import csv
    tmp1 = file[grepl(".csv",file)]
    count.matrix = read.csv2(tmp1, header = T, row.names = 1, fill = F)
    if(ncol(count.matrix) <= 1) {
      count.matrix = read.csv(tmp1, header = T, row.names = 1, fill = F)
    }
  }
  return(count.matrix)
}
countMtxMerge = function(cmtx.ls){
  mergeMtx = function(df1,df2){merge(df1,df2, by=0, all = T)}
  mtx = Reduce(mergeMtx, cmtx.ls)
  row.names(mtx) <- mtx$Row.names
  return(mtx[,-1])
}

cmtx.ls = list()
tmp = list.files(full.names = F, pattern = "[Cc]ount[Mm]atrix")
for(i in tmp){ cmtx.ls[[sub("[Cc]ount[Mm]atrix","",i)]] <- importCountMtx(i) }
if(length(cmtx.ls) == 1){
  count.matrix = cmtx.ls[[1]]
} else {
  count.matrix <- countMtxMerge(cmtx.ls)
}
# subset cols in count.matrix after rows in coldata. So there is only the need to trim the coldata file manually
count.matrix <- subset(count.matrix, select = row.names(coldata))

# Sort & check correct data import:
count.matrix <- count.matrix[,rownames(coldata)] # reduce the to analyse samples in mtx by the samples listed in coldata
stopifnot(all(rownames(coldata) == colnames(count.matrix))) # Must be TRUE!
message("Done!\n")
rm(i,tmp, cmtx.ls)

### Set annotation colors for downstream plotting ###
col.Cond = colorRampPalette(RColorBrewer::brewer.pal(n=8, name="YlOrRd"))(length(levels(coldata$Condition)))
col.Cond[1] <- "#56B4E9" #Defines color for control
col.Tank = colorRampPalette(c("gray95","gray50"))(length(levels(coldata$Tank)))
ann_colors = list(Tank = setNames(col.Tank, levels(coldata$Tank)),
                  Condition = setNames(col.Cond, levels(coldata$Condition)))
####################################

### Setup output environment ###
dir.create("DESeq2_Pairwise", showWarnings = F)
setwd("DESeq2_Pairwise")
home = getwd()

##################
###   DESeq2   ###
##################
## Create DESeq2 object -------------------------------------------------------
p  = 0.05 # padj cutoff
Qt = 0.75  # %-Quantile for effect size cut off (i.e. 0.9 => 90% = top 10% abs(LFC))
message(paste0("\n Starting DESeq2 Analysis - Pairwise Wald's t-test and IHW \n padj < ",p,"\n Tested Substance: \t",substance,"\n"))

id.ls = list()
for(i in condition) {id.ls[[i]] = rownames(coldata[coldata$Condition %in% i,])}
stopifnot(length(id.ls) > 1) # Stop script here if there is only one condition provided
# We will need this object later ... 

# Check if we have a balanced test design.
t1 = all(table(coldata$Condition) == table(coldata$Condition)[1])
t2 = all(table(coldata$Tank) == table(coldata$Tank)[1])
if(all(c(t1,t2) == T)){
  # if TRUE we have a balanced test design between Conditions and Tanks.
  # Hence we can apply a multifactor model as: ~ Tank + Condition
  message("Test design of experiment is BALANCED.\nApplying multi-factor model:\t~ Tank + Condition")
  dds <- DESeqDataSetFromMatrix(countData = count.matrix, colData = coldata,
                                design = ~ Tank + Condition) 
  # Multifactor model controlling for variability in Tanks
} else {
  # if FALSE test design is not balanced and we can not directly apply the batch effect correction here.
  # Hence the model design in this case is only ~ Tank
  message("Test design of experiment is NOT BALANCED.\nApplying uni-factor model:\t~ Condition")
  dds <- DESeqDataSetFromMatrix(countData = count.matrix, colData = coldata,
                                design = ~ Condition)
}
rm(t1,t2)

## Count Matrix filtering -----------------------------------------------------
message(paste0("\nRemoving low abundant counts from count matrix.",
               "\nMinimum number of gene counts per row: \t\t",ncol(count.matrix)))

# Filter count matrix - removing low count genes
dds  = dds[row.names(counts(dds)[rowSums(count.matrix) >= ncol(count.matrix), ]),]

# Get a list of different disp estimates for later QC plotting
fitType = c("parametric", "local") #additonally "mean" possible
estDisp.ls <- list()
for(i in fitType){
  message(paste("\nEstimate Dispersions for fitType: ",i))
  x <- estimateSizeFactors(dds)
  x <- estimateDispersions(x, fitType = i)
  estDisp.ls[[i]] <- x
}
rm(x,i,fitType)

## RUN DESeq() ----------------------------------------------------------------
message("\nRunning DESeq2 pairwise comparison ...")
dds <- DESeq(dds, test = "Wald")                   # Pairwise comparison
#dds <- DESeq(dds, test = "LRT", reduced = ~ Tank) # ANOVA-like approach

## Extract norm. & transformed count.matrix ------------------------------------
normMtx =list()
normMtx[["norm"]]   = round(counts(dds, normalized = T),3) # [[1]]  normalized read counts
normMtx[["ntd"]]    = round(assay(normTransform(dds)),3)   # [[2]] (n+1)log2 transformed norm counts
normMtx[["nrl"]]    = round(assay(rlog(dds, blind = F)),3) # [[3]] rlog transformed mean read counts
normMtx[["nrl.bl"]] = round(assay(rlog(dds, blind = T)),3) # [[4]] rlog transformed mean read counts; BLIND!

message("\nExporting DESeq2 normalized count matrix. This might take a while...\nSaving to txt table ...")
write.table(normMtx[[1]], paste0(substance,"_DESeqNormCounts.txt")) # Extract DESEq2 normalized count matrix
#write.table(df.rld, paste0(substance,"_rlogDESeqNormCounts.txt")) # Extract rlog(norm counts) matrix
message("Done!\n")

## Extract DESeq2 result tables to list objects --------------------------------
resNames = resultsNames(dds)[grepl("^[Cc]ondition_",resultsNames(dds))]
## Non-shrunk results
res.ls <- list() #create an empty list to store objects in
for(i in resNames) {
  x = gsub("[Cc]ondition_","", i)
  message(paste0("Extracting DESeq2 results for: ",x," [IHW, non-shrunk Log2FC] ..."))
  res.ls[[gsub("_.+","",x)]] <- results(dds, name = i, filterFun = ihw, alpha = p)
  message("Done!\n")
}
## apeglm shunk results
resLfs.ls <- list() # save log2FC shrunk results in different list
for(i in resNames) {
  x = gsub("[cC]ondition_","", i)
  message(paste0("Extracting DESeq2 results for: ",x," [IHW, apeglm shrunk Log2FC] ..."))
  resLfs.ls[[gsub("_.+","",x)]] <- lfcShrink(dds, coef = i, type = "apeglm", 
                                             res = res.ls[[gsub("_.+","",x)]])
  message("Done!\n")
}
rm(x,i)

## Append qvalue to res.ls / resLfs.ls -----------------------------------------
message("\nCalculating and appending qvalues to DESeq2 results ...")
# res.ls
res.ls = lapply(res.ls, function(x){
  x$qval = qvalue(x$pvalue)$qvalue
  x })
# resLfs.ls
resLfs.ls = lapply(resLfs.ls, function(x){
  x$qval = qvalue(x$pvalue)$qvalue
  x })
message("Done!\n")
###################


###################
###   biomaRt   ###
###################
## Create annotation object for DESeq2 result tables ---------------------------
message("\nAnnotating ENSEMBL Gene IDs in result tables using biomaRt ...")
# if rerio mart object found locally load it, create new mart object from host 
if(file.exists("~/biomaRt/drerio_mart.Robj")) {
  message("\nDanio rerio mart object found locally. \nLoading 'drerio_mart.Robj' into R session ...")
  load("~/biomaRt/drerio_mart.Robj")
} else if (file.exists("S:/data/biomaRt/drerio_mart.Robj")){
  message("\nDanio rerio mart object found in S:/data/biomaRt/ \nLoading 'drerio_mart.Robj' into R session ...")
  load("S:/data/biomaRt/drerio_mart.Robj")
} else {
  message("\nCould not find 'drerio_mart.Robj'. Creting new mart object from 'www.ensembl.org'",
          "\nMake sure you have a working interent connection. Otherwise this will fail!\nConnecting to server ...")
  rerio = biomaRt::useMart("ENSEMBL_MART_ENSEMBL",dataset="drerio_gene_ensembl",host="https://www.ensembl.org")
}

## Build gene annotation object & merge ----------------------------------------
id     = rownames(res.ls[[1]])
idType = "ensembl_gene_id"
attr   = c("ensembl_gene_id","external_gene_name","description","gene_biotype",
           "entrezgene_id")
message("\nAnnotating Ensembl gene IDs ...")
GeneAnno = biomaRt::getBM(attributes = attr, mart = rerio, uniqueRows = F,
                          filters = idType, values = id, useCache = F)
stopifnot(length(intersect(GeneAnno$ensembl_gene_id,id)) == length(id)) #Check Sum

# !!! In a few cases there can be more than one entrez id for a single ensembl gene id !!!
# Here we remove duplicated entries (~ 360 out of ~ 25 000, so it's really a minor fraction)
# As it turns out higher integer values in the entrez ID correspond to a deeper level of organization
# i.e. ID 999 -> Protein A, ID 1000348 -> Protein A subunit a1. For our purpose it is good enough 
# to retrieve the overall gene / protein information.
# => We sort GeneAnno by entrezgene_id and then remove the duplicates with !duplicate().
# That way, in case multiple entrez ids are assigned to a single ensembl id the lower entrez id
# will be kept for downstream analysis. 
GeneAnno = dplyr::arrange(GeneAnno, entrezgene_id)
GeneAnno = GeneAnno[!duplicated(GeneAnno$ensembl_gene_id),] #keep only non duplicated
stopifnot(length(id) == nrow(GeneAnno)) #Check Sum
row.names(GeneAnno) = GeneAnno$ensembl_gene_id
# To ensure compatability with downstream GSEA scripts rename external_gene_name & entrezgene_id
colnames(GeneAnno) = c("ensembl_gene_id","SYMBOL","description","biotype","ENTREZID")

# Now merge GeneAnno with DESeq2 result tables
AnnoFun = function(x){
  tmp = merge(as.data.frame(x), GeneAnno, by=0)[,-1]
  row.names(tmp) = tmp$ensembl_gene_id
  tmp[,-which(colnames(tmp) %in% "ensembl_gene_id")]
}
res.ls    = lapply(res.ls, FUN = AnnoFun)
resLfs.ls = lapply(resLfs.ls, FUN = AnnoFun)
# Finally sort res.ls / resLfs.ls after the padj value
res.ls    = lapply(res.ls, function(x){dplyr::arrange(x, padj)})
resLfs.ls = lapply(resLfs.ls, function(x){dplyr::arrange(x, padj)})
rm(id, rerio)
message("Done!\n")
###################


###  EXPORT DESeq2 RESULTS  ### -------------------------------------------------------
message("Exporting DESeq2 result tables. This might take a while :)")
dir.create("Results", showWarnings = F)
for(k in names(res.ls)){
  n = gsub("[Ee]xposure","",k)
  message(paste("Saving DESeq2 result table",k,"[IHW, non-shrunk Log2FC] to csv ..."))
  write.csv2(res.ls[[k]], file = paste0("Results/",substance,"_res_",n,".csv"))
  message(paste("Saving DESeq2 result table",k,"[IHW, apeglm shrunk Log2FC] to csv ..."))
  write.csv2(resLfs.ls[[k]], file = paste0("Results/",substance,"_reslfs_",n,".csv"))
}
message("Done!\n")

## Select only DEGs for export and downstream plotting ##
deg.ls    <- lapply(res.ls, function(x){subset(x, padj <= p)}) # pcut
degCut.ls <- lapply(res.ls, function(X){                       # pcut & log2FC cutoff 
  lfcut = quantile(abs(X$log2FoldChange), Qt) #(90% quantile of abs(log2FC))
  if(lfcut > 1){lfcut = 1} # in case lfcut greater 1, set to 1
  x = subset(X, padj <= p)
  subset(x, abs(log2FoldChange) >= lfcut)
})
# Selecting only DEGs with padj < p & abs(apeglm shrunk Log2FC) > lfcut (90% quantile from abs(non shrunk log2FC))
degLfs.ls = list()
for(k in names(resLfs.ls)){
  # determine lfc cut off based on non-shrunk values
  lfcut = quantile(abs(res.ls[[k]]$log2FoldChange), Qt)
  if(lfcut > 1){lfcut = 1} # in case lfcut greater 1, set to 1
  message(paste("Log2FC cut off for:\t",k,"\t----->\t",round(lfcut,2),"(Top",(1-Qt)*100,"% Quantile)"))
  # Select padj < p from resLfs.ls and subset by log2FC
  degLfs.ls[[k]] <- subset(subset(resLfs.ls[[k]], padj <= p), abs(log2FoldChange) >= lfcut)
}

## Exporting DEGs in separate result csv tables ##
dir.create("DEGs", showWarnings = F)
out = c("unshrink","unshrinkFCcut","lfsFCcut") # output dir. for different levels of filtering
for(k in out){dir.create(paste0("DEGs/",k), showWarnings = F)}
# Export to respective dir
for(k in names(deg.ls)) {
  message(paste("Saving selection of DEGs for",k,"to csv ..."))
  n = gsub("[Ee]xposure","",k)
  write.csv2(deg.ls[[k]],paste0("DEGs/",out[1],"/",substance,"_pcut_",n,".csv"))   # unshrink
  write.csv2(degCut.ls[[k]],paste0("DEGs/",out[2],"/",substance,"_pcut_LFcut_",n,".csv"))# unshrink + FCcut
  write.csv2(degLfs.ls[[k]],paste0("DEGs/",out[3],"/",substance,"_lfs_pcut_LFcut_",n,".csv"))# shrink + FCcut
}
message("Done!\n")
rm(n,k,lfcut,out)
###############################


######################################   PLOTTING   ##########################################
dir.create(paste0(home,"/Plots"), showWarnings = F)
dir.create(paste0(home,"/Plots/DataQC"), showWarnings = F)
message("\nStarting Data Composition and QC plotting ...")

### QC - DESeq2 norm. & count distr. & countMtx filtering  ### -----------
message("DESeq2 normalization & gene count distribution ...")
while (!is.null(dev.list()))  dev.off()
pdf(file = paste0("Plots/DataQC/",substance,"_CountNorm_mtxFiltering.pdf"),
    width = 9, height = 11.7)
par(mfrow=c(4,2))

## Normalization bar plots ##
barplot(colSums(counts(dds)),
        ylab= "Total gene counts",
        main= "Raw counts",
        las= 3, #rotating sample labels 90
        ylim= c(0,1.2*max(colSums(counts(dds)))))
barplot(colSums(counts(dds, normalized=T)),  #plot mean normalized counts
        ylab= "Total gene counts",
        main= "DESeq2 norm. counts",
        las= 3, #rotating sample labels 90
        ylim= c(0,1.2*max(colSums(counts(dds)))))
boxplot(log10(counts(dds)+1),
        ylab= "Log10(Gene counts + 1)",
        las = 3, #rotating sample labels 90
        main= "Raw counts")
boxplot(log10(counts(dds, normalized=T)+1),
        ylab= "Log10(Gene counts + 1)",
        las = 3, #rotating sample labels 90
        main= "DESeq2 norm. counts")

## Gene count distr. & countMtx filtering ##
cut = ncol(count.matrix) # min counts per row
# quantile distribution of gene count matrix row sums
tmp = apply(count.matrix, 1, sum)
tmp = quantile(log10(tmp+1), probs = seq(0, 1, .05))
tmp = data.frame(x=seq(0,1,.05)*100, y=tmp)
plot(tmp, xlab = "Quantile-%", ylab = "log10(gene row sum + 1)", main = "Gene row sum quantile distribution")
abline(h = log10(cut+1), lwd=1.5, lty=2, col = "firebrick")
text(x = 0, y = log10(cut+1)*1.3, pos = 4, offset = 0.5,labels = paste("Cutoff:",cut))

# rowsum cutoff vs remaining number of genes in count matrix
cutOffPlot <- function(countMtx, cut) {
  n = ncol(countMtx)
  if(missing(cut)){cut = n}
  if(n < 6) {stop("Number of count matrix columns < 6")}
  X = seq(n-6,n+50,1) 
  X[X == 0] <- 1
  X = sort(union(X, c(1:5)))
  
  rNbr = c() # empty vector to store Nbr of genes per cutoff in
  for(x in X) {
    tmp  = nrow(countMtx[which(rowSums(countMtx) >= x), ])
    rNbr = append(rNbr, tmp)
  }
  df <- data.frame(x = X, y = rNbr)
  plot(df, xlab = "Row sum cut off", ylab = "Nbr of genes",
       main = "Genes in count matrix")
  abline(v=cut, lwd=1.5, lty=2, col = "firebrick")
  #abline(h=df[cut,2], lwd=1.5, lty=2, col = "blue")
  text(x = cut, y = df[cut,2], pos = 4, offset = 1.5,
       labels = paste("Cutoff:",cut,"=",df[cut,2],"genes"))
}
cutOffPlot(count.matrix, cut)

# gene count distribution before AND after row sum cutoff (based on on row means)
tmp  = apply(count.matrix,1,mean)
tmp1 = apply(count.matrix[which(rowSums(count.matrix) >= cut), ], 1, mean)
hist(log2(tmp+1), breaks = 40, main = "countMtx", 
     ylim = c(0,2500), xlim = c(0,20),
     xlab = "log2(Row mean count +1)")
hist(log2(tmp1 +1), breaks = 40, main = paste("filtered countMtx - min. rowSum:",cut),
     ylim = c(0,2500), xlim = c(0,20),
     xlab ="log2(Row mean count +1)") #histogram of filtered counts

dev.off()
rm(cut,tmp, tmp1)
##############################################################

### QC - DESeq2's dispersion estimate models ### ---------------------------------------------
message("DESeq2's dispersion estimate models (estimateDispersions) ...")
while (!is.null(dev.list()))  dev.off()
pdf(file = paste0("Plots/DataQC/",substance,"_DispEstimates.pdf"),
    width = 7, height = 5.9, onefile = T, bg = "transparent")
for(i in names(estDisp.ls)){
  plotDispEsts(estDisp.ls[[i]], main = paste(substance,'- DESeq2 fit type:',i))
}
dev.off()
rm(estDisp.ls, i)
gc() #free some memory 
################################################

### QC - pvalue and LFC distribution ### -----------------------------------------------------
message("pvalue and LFC distribution ...")
while (!is.null(dev.list()))  dev.off()
n = length(condition)-1
pdf(file = paste0("Plots/DataQC/",substance,"_pvalue_LFC_distr.pdf"),
    width = 3*n, height = 13.33, onefile = T, bg = "transparent")
#png(filename = paste0("Plots/DataQC/",substance,"_pvalue_LFC_distr.png"),
#    width = 6*n, height = 18, units = "cm", bg = "white", pointsize = 7, res = 450)
par(mfrow=c(4,n))
## pval distr
for(i in names(res.ls)){
  hist(res.ls[[i]]$pvalue, main= paste(substance,i),col = "gray50", 
       border = "gray50", xlab = "p-value", breaks = 40)
}
## pval vs padj & qval
pVSpadjPlot <- function(res, pval=.05, q=.05, p=.05, title = ""){
  plot(res$pvalue, res$padj, xlab="pvalue", ylab="conversion of pvalue",
       xlim= c(0,0.25), ylim = c(0,1), main= paste(substance,title),
       col = "dodgerblue1",pch = 20)
  points(res$pvalue, res$qval, col="springgreen3", pch = 20)
  abline(h=p, lwd=1, lty=2)
  abline(v=pval, lwd=1, lty=2)
  xpos <- 0.1
  text(x = .25, y = .2, pos = 2, offset = 0, col = "black",
       labels = paste("Genes with pvalue <",pval,":\t",sum(res$pvalue <= pval, na.rm = T)))
  text(x = .25, y = .15, pos = 2, offset = 0, col = "springgreen3",
       labels = paste("Genes with Qvalue <",q,":\t",sum(res$qval <= q, na.rm = T)))
  text(x = .25, y = .1, pos = 2, offset = 0, col = "dodgerblue",
       labels = paste("Genes with padj (BH) <",p,":\t",sum(res$padj <= p, na.rm = T)))
}
for(i in names(res.ls)) { pVSpadjPlot(res.ls[[i]], title = i) }
## Log2FC distr.
lfcDistPlot <- function(res, q=Qt, LFcut=quantile(abs(res$log2FoldChange),q), 
                        title="", Ylim = c(0,2000), xlab="Log2(fc)") {
  m = max(abs(res$log2FoldChange))
  if(m > 5){m = 5}
  hist(res$log2FoldChange, breaks= 500, main= paste(substance,title), col= "gray50",
       border= "gray50", xlab=xlab, xlim= c(-m, m), ylim = Ylim)
  abline(v= c(LFcut,-LFcut), col= "dodgerblue", lty= 2, lwd= 1)
  legend(x="topright", text.col = "dodgerblue", bty = "n",
         legend = paste0("LFcut: +/-",round(LFcut, digits = 2),
                         "\n",(1-q)*100,"% Quantile"))
}
for(i in names(res.ls)) { lfcDistPlot(res.ls[[i]], title = i) }
for(i in names(res.ls)) {
  lfcut = quantile(abs(res.ls[[i]]$log2FoldChange),Qt)
  lfcDistPlot(resLfs.ls[[i]], title = i, LFcut = lfcut, Ylim = c(0,200), xlab = "apgelm(lfc)")
}
dev.off()
rm(n,i,lfcut)
########################################

### QC - Norm. gene count transformation ### -------------------------------------------------
message("Norm. count matrix transformation ...")
while (!is.null(dev.list()))  dev.off()
pdf(file = paste0("Plots/DataQC/",substance,"_NormCount_Transformation.pdf"),
    width = nrow(coldata)-2, height = 4, onefile = T, bg = "transparent")
par(mfrow=c(1,3))
boxplot(normMtx$norm , notch = TRUE,las = 3, #rotating sample labels 90
        main = "Normalized read counts", ylab = "norm. read counts")
boxplot(normMtx$ntd, notch = TRUE, las = 3, #rotating sample labels 90
        main = "log2 Transformation", ylab = "log2 (norm. read counts + 1)")
boxplot(normMtx$nrl.bl, notch = TRUE, las = 3, #rotating sample labels 90
        main = "rlog Transformation", ylab = "rlog (norm. read counts)")
dev.off()
############################################

### QC - PearsonCor of biol. replicates ### ----------------------------------------------------
message("Biological replicates correlation ...")
# get all possible combinations of samples for each condition
comb = lapply(id.ls, function(x){gtools::combinations(length(x),2,x)})
# plot dim.
w = nrow(comb[[1]])   #width
h = length(condition) #height
# plotting function
corplot.1 = function(df, tmp, Lab="", title="") {
  d  = densCols(df[,tmp[1]], df[,tmp[2]], colramp = colorRampPalette(rev(rainbow(10, end = 4/6))))
  ggplot(df) + geom_point(aes(df[,tmp[1]], df[,tmp[2]], col = d), size = .8, alpha = .4) +
    labs(x=paste(Lab,tmp[1]),y=paste(Lab,tmp[2]),
         subtitle = (paste0(title," [",coldata[tmp[1],"Tank"]," vs ",coldata[tmp[2],"Tank"],"]"))) +
    annotate("text", label = paste0("Pearson: ",round(cor(df[,tmp[1]], df[,tmp[2]]),2)),
             x = (max(df[,tmp[1]])), y = (min(df[,tmp[2]])), hjust=1, vjust=0) + 
    annotate("text", label = paste0("R2: ",round((cor(df[,tmp[1]], df[,tmp[2]]))^2,2)),
             x = (min(df[,tmp[1]])), y = (max(df[,tmp[2]])), hjust=0, vjust=1) +
    scale_color_identity() + coord_equal(ratio=1) + theme_bw()
}
multiCorPlot = function(mtx, label="") {
  gg.ls = list()
  for(k in names(comb)){
    x = comb[[k]]
    for(i in 1:nrow(x)){
      tmp = x[i,] #IDs to plot with
      df = as.data.frame(mtx[ ,tmp])
      gg.ls[[paste0(k,i)]] <- corplot.1(df, tmp, Lab = label,
                                        title = gsub("[Ee]xposure","",k))
    }
  }
  gg.ls
}
## log2 ##
gg = multiCorPlot(normMtx$ntd, label = "log2(norm.counts)")
while (!is.null(dev.list()))  dev.off()
png(filename = paste0("Plots/DataQC/",substance,"_Correlation_log2.png"),
    width = 7.33*w, height = 7.6*h, units = "cm", bg = "white", 
    pointsize = .5, res = 500)
print(ggpubr::ggarrange(plotlist = gg, ncol = w, nrow = h, 
                        labels = LETTERS[1:length(gg)]))
dev.off()
## rlog ##
gg = multiCorPlot(normMtx$nrl.bl, label = "rlog(norm.counts)")
while (!is.null(dev.list()))  dev.off()
png(filename = paste0("Plots/DataQC/",substance,"_Correlation_rlog.png"),
    width = 7.33*w, height = 7.6*h, units = "cm", bg = "white", 
    pointsize = .5, res = 500)
print(ggpubr::ggarrange(plotlist = gg, ncol = w, nrow = h, 
                        labels = LETTERS[1:length(gg)]))
dev.off()
# clear env
rm(w,h,gg,comb)
gc()
###########################################

### QC - SD varriance plot ### -------------------------------------------------
my_meanSdPlot <- function(normMtx){
  sdp1 <- vsn::meanSdPlot(normMtx$norm, ranks = T, plot = F)
  sdp1 <- sdp1$gg + ggtitle("Mean counts - Ranked")+scale_y_continuous(trans='log10')+ylab("sd (log scale)") +
    theme_bw()
  
  sdp2 <- vsn::meanSdPlot(normMtx$ntd, ranks = T, plot = F)
  sdp2 <- sdp2$gg + ggtitle("log2 (mean counts) - Ranked") +
    theme_bw() #+ ylim(0,3)
  
  sdp3 <- vsn::meanSdPlot(normMtx$nrl.bl, ranks = T, plot = F)
  sdp3 <- sdp3$gg + ggtitle("rlog (mean counts) - Ranked") +
    theme_bw() #+ ylim(0,3)
  
  sdp1b <- vsn::meanSdPlot(normMtx$norm, ranks = F, plot = F) #original scale with ranks=F
  sdp1b <- sdp1b$gg + ggtitle("Mean counts") +
    theme_bw()
  
  sdp2b <- vsn::meanSdPlot(normMtx$ntd, ranks = F, plot = F) #original scale with ranks=F
  sdp2b <- sdp2b$gg + ggtitle("log2 (mean counts)") +
    theme_bw() #+ ylim(0,3)
  
  sdp3b <- vsn::meanSdPlot(normMtx$nrl.bl, ranks = F, plot = F) #original scale with ranks=F
  sdp3b <- sdp3b$gg + ggtitle("rlog (mean counts)") + 
    theme_bw()#+ ylim(0,3)
  
  print(ggarrange(sdp1, sdp2, sdp3, sdp1b, sdp2b, sdp3b,
                  labels = c("A1", "B1", "C1","A2", "B2", "C2"),
                  ncol = 3, nrow =2))
}

message("Plotting standard deviation composition ...")
while (!is.null(dev.list()))  dev.off()
pdf(file = paste0("Plots/DataQC/",substance,"_meanSD.pdf"), #print directly to pdf
    width = 17, height = 9.6)
my_meanSdPlot(normMtx)
dev.off()
message("Done!\n")
##############################

### Sample Distance Matrix ### --------------------------------------------------
message("Sample distance matrix plotting ...")
myDismtx <- function(mtx, method="euclidean", top=500, title="") {
  if(nrow(mtx) < top){top <- nrow(mtx)}
  message(paste0("Plotting Sample Dist for Var.Top:\t\t",top))
  
  # Subsetting the top genes with max variance
  Xvar <- apply(mtx,1,var)
  X <- t(mtx[names(sort(Xvar, decreasing = T)[1:top]),])
  
  # transpose input, calculate sample distance and create a distance matrix
  sampleDist <- dist(X, method = method) #must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
  distMtx <- as.matrix(sampleDist)
  #rownames(distMtx) <- paste(coldata$Condition, coldata$Tank, sep="-")
  
  # create annotation object for heatmap
  ann_col <- subset(coldata, select = c("Condition"))#,"Substance"))
  ann_row <- subset(coldata, select = c("Tank"))
  #rownames(ann_row) <- paste(coldata$Mixture, coldata$BioReplicate, sep="-")
  
  colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(15)
  pheatmap::pheatmap(distMtx, angle_col = "90", display_numbers = T,
                     #treeheight_col = 40, #default 50
                     #fontsize = 12, #default 10
                     #fontsize_number = 0.7*12,
                     #cellwidth = 26,
                     #cellheight = 26,
                     drop_levels = T,
                     number_format = if(method == "manhattan"){"%1.0f"}else{"%.1f"},
                     main = paste(title,"- Dist:",method,"- topVar:",top),
                     clustering_distance_rows = sampleDist,
                     clustering_distance_cols = sampleDist,
                     annotation_col = ann_col, # Condition!
                     annotation_row = ann_row, # Tank!
                     annotation_colors = ann_colors,
                     col = colors)
}

mtx    = normMtx$nrl.bl #Input for disMtx
topVar = c(2500,5000,nrow(mtx)) # Number top N varying genes/proteins

while(!is.null(dev.list()))  dev.off()
pdf(file = paste0(home,"/Plots/",substance,"_SampleDistMtx.pdf"), width = 8, height = 6.5)
for(dis in c("euclidean","canberra","manhattan")){
  message(paste0("\nComputing Sample Dist with dist measure:\t",dis))
  for(i in topVar){ myDismtx(mtx, method=dis, top=i, title=substance) }
}
dev.off()
rm(i,dis,topVar,mtx)
message("Done!\n")

##############################

### PCA & t-SNE ### -----------------------------------------------------------------
message("PCA and t-SNE clustering ...")
# Add other k-mean cluster related type => See ENSEMBL Course

myPCA  <- function(mtx, pcaM = "svd", top = 500, title = "", IDlabs = F) {
  if(nrow(mtx) < top){top <- nrow(mtx)}
  message(paste0("Plotting PCA for ",title," Var.Top:\t\t\t",top))
  
  # Subsetting the top genes with max variance
  Xvar <- apply(mtx,1,var)
  X <- t(mtx[names(sort(Xvar, decreasing = T)[1:top]),])
  pc <- pcaMethods::pca(X, method = pcaM, center = T, nPcs=2)
  pcaDf <- merge(coldata, pcaMethods::scores(pc), by=0)
  
  g = ggplot(pcaDf, aes(PC1, PC2, colour = Condition, shape=Tank)) +
    geom_point(size = 3, alpha = .65) + scale_colour_manual(values = ann_colors$Condition) +
    ggtitle(paste0(title," - exVar: ",round(pc@R2cum[2]*100,1),"%; ",pcaM,"; Top:",pc@nVar)) +
    xlab(paste0("PC1: ",round((pc@R2)[1]*100,1),"% variance")) +
    ylab(paste0("PC2: ",round((pc@R2)[2]*100,1),"% variance")) + #stat_ellipse() +
    theme_bw() + theme(aspect.ratio = 1, legend.position = "none") +
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12))
  if(IDlabs == T){
    return(g + ggrepel::geom_text_repel(aes(x=PC1, y=PC2, label=Row.names), nudge_x = .5, nudge_y = .5))
  } else {
    return(g)
  }
}
mytSNE <- function(mtx, top = 500, title = "", IDlabs = F) {
  if(nrow(mtx) < top){top <- nrow(mtx)}
  message(paste0("Plotting t-SNE for ",title," Var.Top:\t\t\t",top))
  
  # Subsetting the top genes with max variance
  Xvar <- apply(mtx,1,var)
  X <- t(mtx[names(sort(Xvar, decreasing = T)[1:top]),])
  
  # perform Rtsne calc.
  perplex = (nrow(X)-1)/3
  set.seed(42) #to make results reproducable! 
  tsne <- Rtsne::Rtsne(X, dims=2, initial_dims=nrow(X), check_duplicates=F, num_threads=2,
                       pca = TRUE, perplexity = perplex, #(should not be bigger than 3 * perplexity < nrow(X) - 1, see details for interpretation)
                       theta = 0.0, #Speed/accuracy trade-off (increase for less accuracy), set to 0.0 for exact TSNE (default: 0.5)
                       partial_pca = F, #(requires the irlba package). This is faster for large input matrices (default: FALSE)
                       max_iter = 10000, #number of iterations
                       is_distance = FALSE, pca_center = TRUE, pca_scale = FALSE, verbose = F,
                       normalize = F) #Default True; Set to F as DESeq's RLE norm was performed prior!
  row.names(tsne$Y) <- row.names(X)
  colnames(tsne$Y) <- c("tSNE_1","tSNE_2")
  
  # ggplot
  ggtsne <- merge(coldata, tsne$Y, by=0)
  g = ggplot(ggtsne) + geom_point(aes(x=tSNE_1, y=tSNE_2, color = Condition, shape = Tank), size=3, alpha = .7) +
    ggtitle(paste(title,"; Top:",top)) + xlab('t-SNE 1') + ylab('t-SNE 2') +
    scale_colour_manual(values = ann_colors$Condition) + theme_bw() + theme(aspect.ratio = 1)
  if(IDlabs == T){
    return(g + ggrepel::geom_text_repel(aes(x=tSNE_1, y=tSNE_2, label=Row.names), nudge_x = .5, nudge_y = .5))
  } else {
    return(g)
  }
}

topVar = c(500,2000,10000,nrow(normMtx$nrl)) # Number top N varying genes/proteins
gg = list() #to store plots in
for(i in topVar){ # ref.norm - ref channel & global mean normalized
  gg[[paste0("pca",i)]]  <- myPCA(normMtx$nrl.bl, top=i, title=paste(substance,"rlog[norm].bl"))
  gg[[paste0("tsne",i)]] <- mytSNE(normMtx$nrl.bl, top=i, title=paste(substance,"rlog[norm].bl")) }
gg1 = list() #to store plots in
for(i in topVar){ # norm - only global mean normalized
  gg1[[paste0("pca",i)]]  <- myPCA(normMtx$nrl, top=i, title=paste0(substance,"rlog[norm]"))
  gg1[[paste0("tsne",i)]] <- mytSNE(normMtx$nrl, top=i, title=paste0(substance,"rlog[norm]")) }
gg2 = list()
for(i in topVar){
  gg2[[paste0("pca",i)]]  <- myPCA(normMtx$nrl.bl, top=i, title=paste(substance,"rlog[norm].bl"), IDlabs=T)
  gg2[[paste0("tsne",i)]] <- mytSNE(normMtx$nrl.bl, top=i, title=paste(substance,"rlog[norm].bl"), IDlabs=T) }
gg3 = list()
for(i in topVar){
  gg3[[paste0("pca",i)]]  <- myPCA(normMtx$nrl, top=i, title=paste(substance,"rlog[norm]"), IDlabs=T)
  gg3[[paste0("tsne",i)]] <- mytSNE(normMtx$nrl, top=i, title=paste(substance,"rlog[norm]"), IDlabs=T) }
# Export to pdf
while(!is.null(dev.list())) dev.off()
pdf(file = paste0(home,"/Plots/",substance,"_PCA.tSNE.rlogNormCounts.pdf"),
    width = 12, height = length(topVar)*4.2 )
print(ggpubr::ggarrange(plotlist=gg, ncol=2, nrow=length(topVar)))
print(ggpubr::ggarrange(plotlist=gg2, ncol=2, nrow=length(topVar)))
print(ggpubr::ggarrange(plotlist=gg1, ncol=2, nrow=length(topVar)))
print(ggpubr::ggarrange(plotlist=gg3, ncol=2, nrow=length(topVar)))
dev.off()
message("Done!\n")
rm(gg,gg1,gg2,gg3,topVar,i)

###################

### MA plots & Vulcano plots ### ------------------------------------------------------
message("Plotting MA and Vulcano plots ...")
# MA FUN
MAfun <- function(res, title="", topN=5, Symbol=T, LFcut= quantile(abs(res$log2FoldChange),Qt)){
  
  x = res
  if(!any(colnames(x)=="baseMean" | colnames(x)=="log2FoldChange")){
    x$mean.ProtAbund <- 2^(x$mean.ProtAbund)
    colnames(x)[grep("mean.ProtAbund", colnames(x))] <- "baseMean"
    colnames(x)[grep("log2FC", colnames(x))] <- "log2FoldChange"
    datType = "prot"
  } else {datType = "rna"}
  
  my_theme <- theme(plot.title = element_text(size = 10, face = 'plain'), #line = element_line(size = .1),
                    axis.text = element_text(size = 10),
                    axis.title = element_text(size = 10),
                    legend.text = element_text(size = 10),
                    text = element_text(size = 7, face = 'plain'))
  
  ma <- ggpubr::ggmaplot(x, fdr = .05, fc = 2^(LFcut), size = 1, top = topN, legend = "bottom",
                         select.top.method = 'fc', # fc or padj
                         palette = c("#B31B21", "#1465AC", "darkgray"),
                         genenames = if(Symbol == F){NULL}else{as.vector(res$SYMBOL)},
                         xlab = if(datType == "prot"){"Mean intensity"}else{bquote(~log[2]~ "(mean expression)")},
                         #font.legend = "bold", font.main = "bold", font.label = c("bold", 11), label.rectangle = F
                         ylab = bquote(~log[2]~ "(fold change)")) +
    ggtitle(paste0(title)) + theme_bw() + theme(legend.position = "bottom")
  
  if(LFcut != 0){
    ma <- ma +
      geom_text(aes(x = max(log2(x$baseMean)*.99), y = LFcut*1.5, label = round( LFcut,2)), size = 3.5) + 
      geom_text(aes(x = max(log2(x$baseMean)*.99), y =-LFcut*1.5, label = round(-LFcut,2)), size = 3.5)
  }
  return(ma + my_theme)
}
# Vulcano FUN
VulcFun <- function(res, title="", topN=5, Symbol=T, LFcut= quantile(abs(res$log2FoldChange),Qt)){
  if(LFcut == 0){LFcut = 0.0000001}
  
  if(Symbol == T){
    select <- res[order(res$padj),"SYMBOL"]
  }else{
    select <- rownames(res[order(res$padj),])
  }
  
  my_theme <- theme(plot.title = element_text(size = 10, face = 'plain'), line = element_line(size = .6),
                    axis.text = element_text(size = 10),
                    axis.title = element_text(size = 10),
                    legend.text = element_text(size = 8),
                    text = element_text(size = 8, face = 'plain'))
  
  vu <- EnhancedVolcano::EnhancedVolcano(res, x = "log2FoldChange", y = "padj",
                                         title = paste0(title," [LFcut:",round(LFcut,2),"]"), subtitle = NULL,
                                         lab = if(Symbol == F){rownames(res)}else{res$SYMBOL},
                                         selectLab = if(topN == 0){select[NULL]}else{select[1:topN]}, # select topgenes based on padj value
                                         legendLabels = c('NS','Log2 FC','padj','padj & Log2 FC'),
                                         xlab = bquote(~log[2]~ "(fold change)"), ylab = bquote(~-log[10]~"("~italic(padj)~")"),
                                         FCcutoff = LFcut,   #default 1
                                         pCutoff = .05,      #default 10e-6
                                         labSize = 3.0,
                                         pointSize = 1, #default 0.8
                                         col = c("grey30", "grey30", "royalblue", "red2"),
                                         #shape = c(1, 0, 17, 19),   #default 19, http://sape.inf.usi.ch/quick-reference/ggplot2/shape for details
                                         colAlpha = 0.4, #default 0.5
                                         hline = c(0.01, 0.001), hlineCol = c('grey40','grey55'), hlineType = 'dotted', hlineWidth = 0.6,
                                         gridlines.major = T, gridlines.minor = T, drawConnectors = F,
                                         #widthConnectors = 0.2, colConnectors = 'grey15'
  ) + xlim(-max(abs(res$log2FoldChange))*1.05, max(abs(res$log2FoldChange)*1.05)) +
    ylim(0, max(-log10(res$padj))*1.05) + theme_bw() + theme(legend.position = "bottom")
  
  return(vu + my_theme)
  
}

# plot
ggMA <- list()
ggVU <- list()
while (!is.null(dev.list()))  dev.off()
#pdf(file = paste0("Plots/",substance,"_MAplots_Vulcano.pdf"), width = 5*length(res.ls), height = 10)

# res - non-shrunk
for(k in names(res.ls)){
  ggMA[[k]] <- MAfun(res.ls[[k]], title = paste(substance,k))
  ggVU[[k]] <- VulcFun(res.ls[[k]], title = paste(substance,k))
}
message("Printing non-shrunk results to png ...")
png(file = paste0("Plots/",substance,"_MAplots_Vulcano.png"),
    width = 8.5*length(res.ls), height = 20, units = "cm", res = 600)
print(ggpubr::ggarrange(plotlist = c(ggMA,ggVU),ncol = length(res.ls), nrow =2))
dev.off()

# resLfs - lfc shrunk
stopifnot(all(names(res.ls) == names(resLfs.ls)))
for(k in names(resLfs.ls)){
  ggMA[[k]] <- MAfun(resLfs.ls[[k]], title = paste(substance,k),
                     LFcut = quantile(abs(res.ls[[k]]$log2FoldChange),Qt))
  ggVU[[k]] <- VulcFun(resLfs.ls[[k]], title = paste(substance,k),
                       LFcut = quantile(abs(res.ls[[k]]$log2FoldChange),Qt))
}
message("Printing apeglm shrunk results to png ...")
png(file = paste0("Plots/",substance,"_MAplots_Vulcano.lfs.png"),
    width = 8.5*length(resLfs.ls), height = 20, units = "cm", res = 600)
print(ggpubr::ggarrange(plotlist = c(ggMA,ggVU),ncol = length(res.ls), nrow =2))
dev.off()

rm(ggMA, ggVU, k)
message("Done!\n")
################################

### Diff.Expr. Protein Correlation ### -----------------------------------------
DE.Cor.Venn = function(df.x, df.y, pcut = .05, typ.x = "DEG in HE", typ.y = "DEG in LE", title = "",
                       LFcut.x = quantile(abs(df.x$log2FC),Qt),
                       LFcut.y = quantile(abs(df.y$log2FC),Qt),
                       COL = c("red3","deepskyblue1","blue4"),
                       lfsSelect = F, lfsSet.x = NULL, lfsSet.y = NULL){
  # subset df.x/.y for plotting based on pcut and LFcut settings
  x = droplevels(subset(df.x[df.x$padj <= pcut, ], abs(log2FC) >= LFcut.x))[,"Gene"]
  y = droplevels(subset(df.y[df.y$padj <= pcut, ], abs(log2FC) >= LFcut.y))[,"Gene"]
  
  # Provide custom gene set selection based on lfs cut
  if(lfsSelect == T){
    message("LF shrunk DEG selection provided")
    x = lfsSet.x
    y = lfsSet.y
  }
  
  # get common DEGs
  ls = list(DEG.x = unique(x), DEG.y = unique(y))
  names(ls) <- c(typ.x, typ.y)
  int = intersect(ls[[1]], ls[[2]])
  uni = union(ls[[1]], ls[[2]]) # union set of x & y
  
  # write df for gg corrplot
  X = droplevels(df.x[df.x$Gene %in% uni,c("Gene","SYMBOL","log2FC")])
  Y = droplevels(df.y[df.y$Gene %in% uni,c("Gene","SYMBOL","log2FC")])
  stopifnot(nrow(X) == nrow(Y))
  df = merge(X,Y, by = "Gene")[,-4]
  
  # append Type info to df
  df[df$Gene %in% ls[[1]], "Type"] <- typ.x
  df[df$Gene %in% ls[[2]], "Type"] <- typ.y
  df[df$Gene %in% int, "Type"] <- "DEG Overlap"
  df$Type <- factor(df$Type, levels = c("DEG Overlap", typ.y, typ.x))
  
  col = setNames(COL,c("DEG Overlap", typ.y, typ.x))
  gg = list()
  
  ## Venn plot ##
  venn = eulerr::euler(ls[c(2,1)], shape = "ellipse")
  s = round(venn$stress,3)
  e = round(venn$diagError,3)
  tmp = sort(unlist(lapply(ls, length)))
  pct = round(length(int)/tmp[1]*100,1)
  gg[["VE"]] = plot(venn,
                    fills = list(fill = col[c(2,3,1)], alpha = .6), #labels = list(col = "black", font = 4),
                    legend = list(col = "black", font = 4),
                    main = paste0("Stress: ",s,"\nDiag.Er: ",e,"\nShared ",pct,"%"), #main = paste0(fn,": sign. terms [padj < ",pcut,"]"),
                    quantities = TRUE, shape = "ellipse", lty = 0)
  
  ## Corr plot ##
  # cor test
  cdf <- df[which(df$Type %in% 'DEG Overlap'),]
  if(nrow(cdf) > 2){
    res <- cor.test(cdf$log2FC.x, cdf$log2FC.y, method = "pearson", alternative = "greater")
    p <- signif(res$p.value, digits = 2)
    t <- round(res$statistic,1)
    d <- round(res$parameter,1)
    c <- round(res$estimate,2)
    if(p <= .001) {i = "***"
    }else if(p <= .01) {i = "**"
    }else if(p <= .05) {i = "*"
    }else {i = "NS"}
  } else {
    warning("Not enough common observations between HE & LE to run cor.test!\n")
    c <- round(cor(cdf$log2FC.x, cdf$log2FC.y, method = "pearson"),2)
    p <- t <- d <- i <- "NA"
  }
  cA <- round(cor(df$log2FC.x,df$log2FC.y),2)
  
  my_theme <- theme(axis.text = element_text(size = 12),
                    axis.title = element_text(size = 13))
  # cor plot
  gg[["COR"]] = ggplot(df, aes(x=log2FC.x, y=log2FC.y, color = Type)) +
    geom_point(size = 3, alpha = .4, shape = 16) +
    geom_point(data = df[which(df$Type %in% 'DEG Overlap'),], #2nd layer
               aes(x=log2FC.x, y=log2FC.y), shape = 16, size = 3, alpha = .5) +
    scale_color_manual(values = col) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = .4) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = .4) +
    theme_bw() + theme(legend.position = 'bottom') + geom_rug(alpha = .5) +
    labs(title = paste(substance,"- DEG cor:",title), 
         subtitle = paste0(typ.x," vs ", typ.y,": (t=",t," df=",d," p=",p,")",i), 
         x = paste0("Log2FC (",typ.x,")"), 
         y = paste0("Log2FC (",typ.y,")")) +
    annotate("text", label = paste0("Pearson Cor = ",c,"\nR2 = ",round(c^2,2)),
             color = COL[1],
             x = (min(df$log2FC.x)), 
             y = (max(df$log2FC.y)),
             hjust=0, vjust=.52) + 
    annotate("text", label = paste0("Pearson Cor = ",cA,"\nR2 = ",round(cA^2,2)),
             x = (max(df$log2FC.x)), 
             y = (min(df$log2FC.y)),
             hjust=.9, vjust=0) + my_theme
  
  return(gg)
}

# Add a gene column & replace 'log2FoldChange' with 'log2FC'; to fitt upper function
mod_res <- function(res.ls){
  RES <- lapply(res.ls, function(res){
    x = res
    x$Gene <- row.names(x)
    colnames(x)[grep("log2FoldChange", colnames(x))] <- "log2FC"
    return(x)
  })
  return(RES)
}

message("\nStart shiny DEG Correlation & Venn plotting ...")

gg <- list() #to store plots in
RES = mod_res(res.ls)
RESlfs = mod_res(resLfs.ls)
# No LFcut
for(k in names(RES)[-length(RES)]){
  n.x = gsub("[.].+$","", tail(names(RES),1))
  n.y = gsub("[.].+$","", k)
  message(paste("Plotting\t", n.y, "vs", n.x))
  tmp = DE.Cor.Venn(RES[[length(RES)]], RES[[k]], LFcut.x=0, LFcut.y=0, title = paste0("padj<",p),
                    typ.x = n.x, typ.y = n.y)
  gg[[paste0(n.y,".VE")]]  <- tmp$VE
  gg[[paste0(n.y,".COR")]] <- tmp$COR
}
# LFcut
for(k in names(RES)[-length(RES)]){
  n.x = gsub("[.].+$","", tail(names(RES),1))
  n.y = gsub("[.].+$","", k)
  message(paste("Plotting\t", n.y, "vs", n.x))
  tmp = DE.Cor.Venn(RES[[length(RES)]], RES[[k]], title = paste0("padj<",p," & LFcut (",(1-Qt)*100,"%Q)"),
                    typ.x = n.x, typ.y = n.y)
  gg[[paste0(n.y,".VEcut")]]  <- tmp$VE
  gg[[paste0(n.y,".CORcut")]] <- tmp$COR
}
# LFcut on apeglm shrunk lfc values
for(k in names(RES)[-length(RES)]){
  n.x = gsub("[.].+$","", tail(names(RES),1))
  n.y = gsub("[.].+$","", k)
  message(paste("Plotting\t", n.y, "vs", n.x))
  LFcut.x = quantile(abs(RES[[length(RES)]]$log2FC),Qt)
  LFcut.y = quantile(abs(RES[[k]]$log2FC),Qt)
  X = RESlfs[[length(RESlfs)]]
  Y = RESlfs[[k]]
  s.x = droplevels(subset(X[X$padj <= p, ], abs(log2FC) >= LFcut.x))[,"Gene"]
  s.y = droplevels(subset(Y[Y$padj <= p, ], abs(log2FC) >= LFcut.y))[,"Gene"]
  tmp = DE.Cor.Venn(RES[[length(RES)]], RES[[k]], title = paste0("padj<",p," & LFcut (",(1-Qt)*100,"%Q) on lfs"),
                    LFcut.x = LFcut.x, LFcut.y = LFcut.y,
                    typ.x = n.y, typ.y = n.x, lfsSelect = T,
                    lfsSet.x = s.x, lfsSet.y = s.y)
  gg[[paste0(n.y,".VEcut.LFS")]]  <- tmp$VE
  gg[[paste0(n.y,".CORcut.LFS")]] <- tmp$COR
}
# LFcut on apeglm shrunk lfc values with shrunk lfc
for(k in names(RES)[-length(RES)]){
  n.x = gsub("[.].+$","", tail(names(RES),1))
  n.y = gsub("[.].+$","", k)
  message(paste("Plotting\t", n.y, "vs", n.x))
  tmp = DE.Cor.Venn(RESlfs[[length(RESlfs)]], RESlfs[[k]], title = paste("padj <",p,"with lfs values"),
                    LFcut.x = 0, LFcut.y = 0,
                    typ.x = n.y, typ.y = n.x)
  gg[[paste0(n.y,".VEcut.LFS2")]]  <- tmp$VE
  gg[[paste0(n.y,".CORcut.LFS2")]] <- tmp$COR
}

while (!is.null(dev.list()))  dev.off()
pdf(file = paste0("Plots/",substance,"_DEG_CorrVenn.pdf"),
    width = 12, height = 6*(length(res.ls)-1) )
print(ggpubr::ggarrange(plotlist = gg, ncol = 2, nrow = length(res.ls)-1))
dev.off()
message("Done!\n")
rm(RES,RESlfs,tmp,gg,k,n.x,n.y,LFcut.x,LFcut.y,X,Y,s.x,s.y)
gc()
######################################

### Multi Venn Plot ### --------------------------------------------------------
myVenn <- function(deg, title = "", shape = "ellipse", ...) {
  # Venn fun
  set.seed(42)
  venn <- eulerr::euler(deg, shape = shape, ...)
  s <- round(venn$stress,3)
  e <- round(venn$diagError,3)
  
  # plot
  return(
    plot(venn,
         fills = list(fill = ann_colors$Condition[-1], alpha = .6),
         legend = list(col = "black", font = 4),
         main = paste0(title,"\nStress: ",s," - Diag.Er: ",e),
         quantities = TRUE, shape = shape, lty = 0)
  )
}

if(length(condition)-1 > 6){
  warning("Dataset contains ",length(condition)-1," treatment conditions.\n",
          "Multi venn diagrams can only be drawn for a maximum of 6 treatments.\nSkipping multi venn plotting.")
} else {
  message("Multi Venn plotting ...")
  # use deg.ls ; degCut.ls & degLfs.ls as input
  while (!is.null(dev.list()))  dev.off()
  pdf(file = paste0("Plots/",substance,"_VennPlots.pdf"),
      width = 13, height = 11)
  par(mfrow = c(2,2))
  gg <- list()
  for(i in c("deg.ls","degCut.ls","degLfs.ls")){
    deg <- lapply(get(i), row.names)
    if(i == "deg.ls"){n = paste0("padj<",p)}
    if(i == "degCut.ls"){n = paste0("padj<",p," & LFcut(",(1-Qt)*100,"% Qt)")}
    if(i == "degLfs.ls"){n = paste0("padj<",p," & LFcut(",(1-Qt)*100,"% Qt) on lfs")}
    # Venn1
    venn::venn(deg, ilab=TRUE, zcolor = ann_colors$Condition[-1], lty = 0, 
               ilcs = 1, sncs = 1, box = F)
    text(0,1000, labels = paste(substance,"DEGs [",n,"]"), pos = 4)
    # Venn 2
    gg[[i]] = myVenn(deg, paste(substance,"DEGs\n",n))
  }
  print(ggpubr::ggarrange(plotlist = gg, ncol = 2, nrow = 2))
  dev.off()
  rm(i,n, gg, deg)
}
message("Done!\n")
#######################

### Session Information & save Rdata  ### --------------------------------------
sink(paste0(home,"/DESeq2_SessionInfo.txt"))
print(date())
print(devtools::session_info())
sink()
message("Saving RData object. This might take a while ...")
save.image(paste0(home,"/DESeq2_pairwise.RData"))
message(paste0("\n\nFinished Wald's pairwise testing DESeq2's DEG analysis & data plotting.",
               "\nAll outputs were saved in:\n",home,
               "\n\nWhat a ride! Finally END of SCRIPT! :)\nJ.A.R.V.I.S. over and out\n"))
setwd("../")
gc()
###    END OF SCRIPT   ####
