# Title: "Methylation Array Analysis"
# Author: "YC"
# Date: "7/17/2020"

args<-commandArgs(T)
baseDir <- args[1]
outDir <- args[2]

baseDir <- "/Users/caiy02/Documents/8_SCLC_Kayla/idats"
outDir <- "/Users/caiy02/Documents/8_SCLC_Kayla/output"

## Load library
library("data.table")
library("minfi")
library("Biobase")
library("limma")
library("ggfortify")
library("ggplot2")
library("qdapTools")
library('plyr')
library('tidyverse')
library("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
library("IlluminaHumanMethylationEPICmanifest")
library("sva")
library("plotly")
library("factoextra")
library("IlluminaHumanMethylationEPICmanifest")
library("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
library("IlluminaHumanMethylation450kmanifest")
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library("RColorBrewer")
library("Rtsne")
library("plotly")
library("qdapTools")
library("gplots")
library("grid")
library("ComplexHeatmap")
library("circlize")
library("gplots")

##### Read in sheet and idat files
targets <- read.metharray.sheet(baseDir, pattern = "samplesheet.csv")

targets_850k <- targets[targets$Array_Type=="850K",]
targets_450k <- targets[targets$Array_Type=="450K",]

RGSet_850K <- read.metharray.exp(targets = targets_850k, force = T) 
RGSet_450K <- read.metharray.exp(targets = targets_450k, force = T)

RGSet <- combineArrays(RGSet_850K, RGSet_450K, outType =c("IlluminaHumanMethylation450k"), verbose = T)

##### QC and plots
qcReportpdf = qcReport(RGSet, sampNames = targets$Sample_ID, sampGroups = targets$Array_Type, pdf = "qcReport.pdf")
qcReportpdf

## Pre-process the samples to normalize using preprocessIllumina {minfi}
Mset.illumina <- preprocessIllumina(RGSet, bg.correct = T, normalize = "controls")
qc <- getQC(Mset.illumina)

## Generate QC plot as PNG
png(
  file = paste(outDir, "/QC.png", sep = ""), 
  width = 1280, 
  height = 1280, 
  pointsize = 50)
plotQC(qc)
dev.off()

##### Remove poor quality samples 
detP <- detectionP(RGSet)
colnames(detP) <- RGSet$Sample_ID

## examine mean detection p-values across all samples to identify any failed samples
#pal <- brewer.pal(10,"Paired") #adjust according to sample number
pal <- c('blue', "red", "cyan4", "orange", "green")
par(mfrow = c(1, 2))
png(
  file = paste(outDir, "/pvalue_sample_filter.png", sep = ""),
  width = 2048,
  height = 2048,
  pointsize = 50
)
barplot(
  colMeans(detP),
  col = pal[factor(targets$Staging)],
  las = 2,
  cex.names = 0.5,
  cex.axis = 0.5,
  names.arg = targets$Sample_ID,
  ylab = "Mean Detection p-values"
)
abline(h = 0.05, col = "red")
legend("topleft",
       legend = levels(factor(targets$Staging)),
       fill = pal,
       bg = "white")
dev.off()


## Filter out low QC samples and probes from all objects
keep <- colMeans(detP) < 0.05
RGSet <- RGSet[, keep]
RGSet
## Remove poor quality samples from targets data
targets <- targets[keep,]
detP <- detP[, keep]
dim(detP)
mSetSq <- preprocessQuantile(RGSet)

## Ensure probes are in the same order in the mSetSq and detP objects
detP <- detP[match(featureNames(mSetSq), rownames(detP)),]

## Remove any probes that have failed in one or more samples
keep_ <- rowSums(detP < 0.01) == ncol(mSetSq)
table(keep_)
gset.funnorm <- mSetSq[keep_,]
dim(gset.funnorm)

## Remove SNP
gset.funnorm <- addSnpInfo(gset.funnorm) 
gset.funnorm <- dropLociWithSnps(gset.funnorm, snps = c("SBE", "CpG"), maf = 0) 

## Get annotation and remove sex probes
annot = getAnnotation(gset.funnorm)
sex_probes = annot$Name[annot$chr %in% c("chrX", "chrY")]
gset.funnorm = gset.funnorm[!(rownames(gset.funnorm) %in% sex_probes),]

## Get beta values to obtain the data matrix
gset.funnorm.beta <- getBeta(gset.funnorm)
gset.funnorm.beta_df <- as.data.frame(gset.funnorm.beta)
colnames(gset.funnorm.beta_df) <- gset.funnorm$Sample_ID
colnames(gset.funnorm.beta) <- gset.funnorm$Sample_ID
#write.csv(gset.funnorm.beta_df,file = "beta_withSampleID.csv",quote = F)

## Sex prediction###
predictedSex <- getSex(mSetSq, cutoff = -2)
mSetSq_withsex <- addSex(mSetSq, sex = predictedSex)
head(predictedSex)
plotSex(mSetSq_withsex)
plotSex(mSetSq)

##### Remove batch effect 
batch = targets$Batch
modcombat = model.matrix(~ 1, data = targets)
combat_beta = ComBat(
  dat = gset.funnorm.beta,
  batch = batch,
  mod = modcombat,
  par.prior = T,
  prior.plots = F
)

combat_beta_df <- as.data.frame(combat_beta)
#write.csv(combat_beta_df,file = "beta_withSampleID_batch_correction.csv",quote = F)

##### heatmaps unsupervised heatmap mpute 0 and 1
var_probes <- apply(combat_beta, 1, var)
select_var <- names(sort(var_probes, decreasing=T)) #[1:10000]
top_variable_beta <- combat_beta[select_var,]

#after batch correction some beta >1 or < 0, so impute >1 is 1 and <0 is 0
top_variable_beta_imputed <- ifelse(top_variable_beta<0, 0, top_variable_beta)
top_variable_beta_imputed <- ifelse(top_variable_beta_imputed>1, 1, top_variable_beta_imputed)

#write.csv(top_variable_beta_imputed,file = "top_variable_beta_imputed.csv",quote = F)

#### Built function for variables color assignment
annotationlist_builder <- function(metatable, customcolorlist = NULL){
  annotlist <- list()
  ## Create the color list from manual distinct colors, expanded appropriately
  colorlist <- rep(c("dodgerblue2", "#E31A1C", "green4","#6A3D9A", "#FF7F00",
                     "black", "gold1", "skyblue2", "#FB9A99", "palegreen2",
                     "#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon", "orchid1",
                     "deeppink1", "blue1", "steelblue4", "darkturquoise", "green1",
                     "yellow4", "yellow3", "darkorange4", "brown"),10)
  colorcount <- 0
  for (colnum in seq_len(ncol(metatable))) {
    # in this case is FALSE
    if (!is.null(customcolorlist) && colnames(metatable[,colnum,drop=FALSE]) %in% names(customcolorlist)) {
      annotlist[colnames(metatable[,colnum,drop=FALSE])] <- customcolorlist[colnames(metatable[,colnum,drop=FALSE])]}
    
    if (!is.null(customcolorlist) && !colnames(metatable[,colnum,drop=FALSE]) %in% names(customcolorlist) | is.null(customcolorlist)) {
      if (!is.numeric(metatable[,colnum])) {
        
        annotlist[[colnum]] <- colorlist[(1+colorcount): (colorcount + length(na.omit(unique(metatable[,colnum]))))]
        colorcount <- colorcount + length(na.omit(unique(metatable[,colnum])))
        
        names(annotlist[[colnum]]) <- na.omit(unique(metatable[,colnum]))
        names(annotlist)[[colnum]] <- colnames(metatable)[colnum]
      }
      else {collist <- colorRamp2(seq( min(metatable[,colnum], na.rm = TRUE), max(metatable[,colnum], na.rm = TRUE), length = 3),
        brewer.pal(3, "Purples"))
      annotlist[[colnum]] <- collist
      names(annotlist)[colnum] <- colnames(metatable)[colnum]
      }
    }
  }
  
  return(annotlist)
}

## to select how many probes with most top variance
top_10000 <- top_variable_beta_imputed[1:10000,]

maptab <- as.matrix(top_10000)
rowdistance <- dist(na.omit(top_10000), method = "euclidean")
rowcluster <- hclust(rowdistance, method = "ward.D2")
rowclusterparam <- rowcluster

## Build Annotations from metatable
hatop <- NULL
metatable <- targets
colannotationlist <- annotationlist_builder(metatable)
colmetatable=metatable

temp1 <- vector("list", length(colannotationlist))
names(temp1) <- names(colannotationlist)
annotlegendlist <- lapply(temp1, function(x) 
  x[[1]] = list(title_gp = gpar(fontsize = 5, fontface = "bold"), labels_gp = gpar(fontsize=4)))

columnannotation_height = unit(min(34/length(colannotationlist), 5),"mm")

hatop <- HeatmapAnnotation(df = data.frame(colmetatable),
                           col = colannotationlist, #color followed by pre-defined annotation1
                           na_col = "white", # na color is white
                           show_annotation_name = TRUE,
                           annotation_name_gp = gpar(fontsize = 6, fontface="bold"),
                           annotation_name_side = "left",
                           simple_anno_size = columnannotation_height,
                           #height = unit(4,"cm"), #Height of the whole column annotations.
                           #annotation_height = unit(4,"mm"), #Height of each annotation
                           show_legend = T, #showlegendparam,
                           annotation_legend_param = c(annotlegendlist, 
                                                       list(title_gp = gpar(fontsize = 6, fontface = "bold"), labels_gp = gpar(fontsize = 6),
                                                            grid_height = unit(5, "mm"), grid_width = unit(5, "mm"))
                           )
)

## Define the Heatmap
col_fun <- colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))

ht <- Heatmap(maptab,
                col = col_fun,                ## Define the color scale
                cluster_columns = T,          ## Cluster the columns
                column_km = 6,   
                cluster_rows = rowclusterparam,             ## Cluster the rows
                
                show_column_names = T,        ## Show the Column Names (which is sample #)
                column_names_gp = gpar(fontsize = 5),       ## Column Name Size
                show_row_names = F,           ## Show Row names (which is probes)
                row_names_side = "left",      ## Place the row names
                row_names_gp = gpar(fontsize = 8),
                
                show_row_dend = T, #nrow(maptab) <=100,         ## Show row dendrogram
                show_column_dend = T,         ## Show col dendrogram
                
                show_heatmap_legend = T,
                top_annotation = hatop[,2:7],
                bottom_annotation = hatop[,8:14],
                heatmap_legend_param = list(title = "Value", #legend_height = unit(2, "cm"), 
                                            labels_gp = gpar(fontsize=6), title_gp = gpar(fontsize = 8, fontface = "bold")
                ),
                
                height = unit(18,"cm"),       ## Heatmap body height
                width = unit(14,"cm"),        ## Heatmap body weight
                #column_split = 3,
                #heatmap_width = unit(min(ncol(top_50), 20),"cm")
)

draw(ht, padding = unit(c(1, 1, 1, 1), "mm"), merge_legend =F)

