#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(RUVSeq))
suppressPackageStartupMessages(library(EDASeq))
suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-count", "--countData", nargs="+", help="Count Data matrix for the RUVSeq")
parser$add_argument("--outfile", help="Output file where result of RUVSeq will be saved")
args <- parser$parse_args()

### countData

countData_time_00 <- read.csv(args$count, sep='\t')
rownames(countData_time_00) <- countData_time_00[,1]
countData_time_00[,1] <- NULL
countData_time_00 <- as.matrix(subset(countData_time_00))

### colData

col1_colData_00 <- as.vector(colnames(countData_time_00))
cond_colData_00 <- gsub('(.*)\\..*', '\\1', col1_colData_00)
rep_colData_00 <- gsub('.*\\.(.*)','\\1',col1_colData_00)

colData_time_00 <- data.frame(condition = cond_colData_00, rep = rep_colData_00, stringsAsFactors = TRUE)
rownames(colData_time_00) <- col1_colData_00


## create expression set for RUVSeq

set_time_00 <- newSeqExpressionSet(counts=countData_time_00, phenoData = colData_time_00)
idx_time_00  <- rowSums(counts(set_time_00) > 5) >= 2
set_time_00 <- set_time_00[idx_time_00, ]

differences <- makeGroups(colData_time_00$condition)

set_00 <- RUVs(set_time_00, unique(rownames(set_time_00)), k=2, differences)
ruv_colData_time_00 <- pData(set_00)
write.table(ruv_colData_time_00, file=paste0(args$outfile), quote = FALSE, row.names=TRUE, sep='\t')
