# Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0
# Install required packages
install.packages("data.table")
install.packages("DESeq2")
install.packages("GEOquery")
install.packages("umap")
install.packages("limma")
install.packages("car")

#load library 
library(GEOquery)

#========= Load data from GEO============-
# load the GEO series matrix file
geo_data_path <- "/Users/alinaelahie/Desktop/GSE270472_series_matrix.txt.gz"
gse <- getGEO(filename = geo_data_path)

# extract the expression data from the GEO object
exprs_data <- exprs(gse[[1]])

# check the dimensions of the expression data
dim(exprs_data)

#======================= Load the Metadata (conditon)=================
# load the metadata (condition data)
library(readxl)

metadata_path <- "/Users/alinaelahie/Desktop/GSE270472_HD_KO_NSC.xlsx"
metadata <- read_excel(metadata_path)

# view the first few rows of the metadata
head(metadata)

#===================== Preprocess the data================
#   Data plots for selected GEO samples
# load counts table from GEO
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE270472", "file=GSE270472_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)

# pre-filter low count genes
# keep genes with at least 2 counts > 10
keep <- rowSums( tbl >= 10 ) >= 2
tbl <- tbl[keep, ]

# log transform raw counts
# instead of raw counts can display vst(as.matrix(tbl)) i.e. variance stabilized counts
dat <- log10(tbl + 1)

# box-and-whisker plot
par(mar=c(7,4,2,1))
boxplot(dat, boxwex=0.7, notch=T, main="GSE270472", ylab="lg(cnt + 1)", outline=F, las=2)

# UMAP plot (dimensionality reduction)
library(umap)
dat <- dat[!duplicated(dat), ] # first remove duplicates
ump <- umap(t(dat), n_neighbors = 5, random_state = 123)
plot(ump$layout, main="GSE270472 UMAP plot, nbrs =5", xlab="", ylab="", pch=20, cex=1.5)
library(car)
pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)

#=================== Differential Expression Analysis ==================
# load DESeq2 package
library(DESeq2)

col_data <- data.frame(condition = factor(metadata$condition))

# create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = exprs_matrix_filtered,
                              colData = col_data,
                              design = ~ condition)

# pre-filter low count genes
dds <- dds[rowSums(counts(dds)) > 1, ]

# run DESeq2 analysis
dds <- DESeq(dds)

# results for differential expression analysis
res <- results(dds)

# view results
head(res)
#===========Create Visulizations================
#volcano plot?

#==========Differential Abundance Analysis ===============
# install and load gprofiler2
install.packages("gprofiler2")
library(gprofiler2)

# Perform gProfiler analysis
gprofiler_results <- gprofiler(query = rownames(res)[which(res$pvalue < 0.05)], 
                               organism = "hsapiens")

# View the results
head(gprofiler_results)
################################################################


