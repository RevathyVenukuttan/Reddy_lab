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

parser <- ArgumentParser()

parser$add_argument("--counts", nargs="+", help="Countdata for all the samples")
parser$add_argument("--outdir", nargs="+", help="Path to output directory where all results will be saved")

args <- parser$parse_args()

count <- args$counts
outdir <- args$outdir


countData <- read.csv(count, sep='\t')
rownames(countData) <- countData[,1]
countData[,1] <- NULL

col1_colData <- as.vector(colnames(countData))
cond_colData <- gsub('(.*)\\..*', '\\1', col1_colData)
comp_colData <- gsub('.*_(.*)', '\\1', cond_colData)
treat_colData <- c(rep('no',31), rep('yes',33))

colData <- data.frame(compound = comp_colData, treatment = treat_colData,
                      stringsAsFactors = TRUE)

rownames(colData) <- col1_colData

dds <- DESeqDataSetFromMatrix(countData = countData, 
                                 colData = colData, 
                                 design = ~ treatment+compound+treatment:compound)

dds$compound <- relevel(dds$compound, ref = 'dex')
dds <- dds[ rowSums(fpm(dds, robust = FALSE)>=2) > 10, ]
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds, fitType='local')
dds <- nbinomWaldTest(dds)

# vsd <- tryCatch({
#     vst(dds, blind=FALSE)
# }, error = function(e){
#     varianceStabilizingTransformation(dds, blind=FALSE)
# })

# pdf(paste0(outdir, '/plots/model1_interaction.pca.1_vs_2.pdf'), width=5, height=5)
# plotPCA_(vsd, intgroup=c('compound'), ntop=500, pcX=1, pcY=2)
# dev.off()
# # Create PCA plot for the second and third PCs
# pdf(paste0(outdir, '/plots/model1_interaction.pca.2_vs_3.pdf'), width=5, height=5)
# plotPCA_(vsd, intgroup=c('compound'), ntop=500, pcX=2, pcY=3)
# dev.off()
# # Create PCA plot for the third and fourth PCs
# pdf(paste0(outdir, '/plots/model1_interaction.pca.3_vs_4.pdf'), width=5, height=5)
# plotPCA_(vsd, intgroup=c('compound'), ntop=500, pcX=3, pcY=4)
# dev.off()

coeff <- coef(dds)

write.table(coeff, file=paste0(outdir, '/model1_interaction.beta_values.txt'),
            quote=FALSE, row.names=TRUE, sep='\t')

contrast_names <- grep('treatmentyes.', resultsNames(dds), value=TRUE)

for (i in 1:length(contrast_names)){
    result <- results(dds, name=(contrast_names[i]), alpha = 0.05)
    result <- result[order(result$padj),]
    result_logFC <- result[(abs(result$log2FoldChange))>1,]
        
        
    pdf(paste0(outdir, '/plots/', contrast_names[i], '.MA_plot.filtered_log2FC.pdf'), width=6, height=6)
    DESeq2::plotMA(result_logFC, main=paste0('MA plot for ', contrast_names[i]), ylim = c(-8,8))
    dev.off()
    
    pdf(paste0(outdir, '/plots/', contrast_names[i], '.MA_plot.pdf'), width=6, height=6)
    DESeq2::plotMA(result, main=paste0('MA plot for ', contrast_names[i]), ylim = c(-8,8))
    dev.off()
    
    write.table(result, file=paste0(outdir, '/', contrast_names[i],'.txt'),
                quote=FALSE, row.names=TRUE, sep='\t')
    write.table(result_logFC, file=paste0(outdir, '/', contrast_names[i],'.filtered_log2FC.txt'),
                quote=FALSE, row.names=TRUE, sep='\t')
}

save.image(paste0(outdir, '/model1_interaction.RData'))