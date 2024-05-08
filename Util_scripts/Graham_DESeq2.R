#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(RUVSeq))
suppressPackageStartupMessages(library(EDASeq))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library("argparse"))

# Rewrite DESeq2 plotPCA function to specify which PC to plot
plotPCA_ = function(object, intgroup="condition", ntop=500, returnData=FALSE, pcX=1, pcY=2, ignoreReps=FALSE)
{
  # calculate the variance for each gene
  rv <- rowVars(assay(object))

  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]

  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))

  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )

  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }

  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])
  
  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=":"))
  } else {
    colData(object)[[intgroup]]
  }
    
  if (returnData) {
    attr(d, "percentVar") <- c(percentVar[pcX], percentVar[pcY])
    return(d)
  }
  rep <- colData(object)[['rep']]
#   rep <- row.names(colData(dds))
  if (!ignoreReps){
      # Add replicate info
#       rep <- colData(object)[['rep']]

        # assembly the data for the plot
      d <- data.frame(PCX=pca$x[,pcX], PCY=pca$x[,pcY], group=group, rep=rep, intgroup.df, name=colnames(object))

      ggplot(data=d, aes_string(x="PCX", y="PCY", color="group", shape="rep")) + geom_point(size=3) + 
        xlab(paste0("PC", pcX, ": ",round(percentVar[pcX] * 100),"% variance")) +
          ylab(paste0("PC", pcY, ": ",round(percentVar[pcY] * 100),"% variance")) +
      scale_shape_manual(values=seq(0,15))
  } else {
#       batch <- as.factor(colData(object)[['batch']])
        # assembly the data for the plot
      d <- data.frame(PCX=pca$x[,pcX], PCY=pca$x[,pcY], group=group, intgroup.df, name=colnames(object))

      ggplot(data=d, aes_string(x="PCX", y="PCY", color="group")) + 
          geom_point(size=3) +
          geom_text(aes(label=rep),hjust=0, vjust=0) +
          xlab(paste0("PC", pcX, ": ",round(percentVar[pcX] * 100),"% variance")) +
          ylab(paste0("PC", pcY, ": ",round(percentVar[pcY] * 100),"% variance")) +
      scale_shape_manual(values=seq(0,15))
  }
#     + coord_fixed()
}

# create parser object

parser <- ArgumentParser()

parser$add_argument("--counts", nargs="+", help="Countdata for all the samples")
parser$add_argument("--contrast", nargs="+", help="Contrast matrix indicating all the contrasts between the given samples")
parser$add_argument("--outdir", nargs="+", help="Path to output directory where all results will be saved")

args <- parser$parse_args()

count <- args$counts
contrast <- args$contrast
outdir <- args$outdir
model_name <- gsub("^.*/", "", outdir)

### countData

countData <- read.csv(count, sep='\t')
rownames(countData) <- countData[,1]
countData[,1] <- NULL


### colData

col1_colData <- as.vector(colnames(countData))
cond_colData <- gsub('(.*)\\..*', '\\1', col1_colData)
rep_colData <- gsub('.*\\.(.*)','\\1',col1_colData)
colData <- data.frame(condition = cond_colData,
                      rep = rep_colData, stringsAsFactors = TRUE)
rownames(colData) <- col1_colData


### contrast_matrix

contrast_matrix <- read.csv(contrast, sep='\t')
control <- as.character(contrast_matrix$Treatment2)
base_control <- unique(control)
base_control

## DESeq

dds <- DESeqDataSetFromMatrix(countData = countData, 
                                 colData = colData, 
                                 design = ~ condition + rep)
dds <- dds[ rowSums(fpm(dds, robust = FALSE)>=2) > 10, ]


for (i in 1:length(base_control)){
    dds$condition <- relevel(dds$condition, ref = base_control[i])
    dds <- estimateSizeFactors(dds)
    dds <- estimateDispersions(dds, fitType='local')
    dds <- nbinomWaldTest(dds)

    coeff <- coef(dds)
    
    treatment <- as.vector(contrast_matrix$Treatment1[contrast_matrix$Treatment2==base_control[i]])
#     print(treatment)
    
    for (j in 1:length(treatment)){
        result <- results(dds, contrast = c('condition',treatment[j],base_control[i]), alpha = 0.05, independentFiltering=F)
        result <- result[order(result$padj),]
        result_logFC <- result[(abs(result$log2FoldChange))>1,]
        
        
        pdf(paste0(outdir, '/plots/', treatment[j], '_vs_', base_control[i], '.MA_plot.filtered_log2FC.pdf'), width=6, height=6)
        DESeq2::plotMA(result_logFC, main=paste0('MA plot for ', treatment[j], '_vs_', base_control[i]), ylim = c(-8,8))
        dev.off()
    
        pdf(paste0(outdir, '/plots/', treatment[j], '_vs_', base_control[i], '.MA_plot.pdf'), width=6, height=6)
        DESeq2::plotMA(result, main=paste0('MA plot for ', treatment[j], '_vs_', base_control[i]), ylim = c(-8,8))
        dev.off()
    
        write.table(result, file=paste0(outdir, '/', treatment[j], '_vs_', base_control[i],'.txt'),
                quote=FALSE, row.names=TRUE, sep='\t')
        write.table(result_logFC, file=paste0(outdir, '/', treatment[j], '_vs_', base_control[i],'.filtered_log2FC.txt'),
                quote=FALSE, row.names=TRUE, sep='\t')
    
        }
    
    write.table(coeff, file=paste0(outdir, '/',base_control[i],'.beta_values.txt'),
                quote=FALSE, row.names=TRUE, sep='\t')
    
}
save.image(paste0(outdir, '/', model_name, '.RData'))
