#!/usr/bin/env Rscript

# load libs
library(DESeq2)
library(deconvSeq)
library(BisqueRNA)
library(Seurat)
library(SingleR)
library(SingleCellExperiment)
library(DWLS)
library(optparse)

# specify script inputs
option_list = list(
    make_option(c("-r", "--rdata_file"), type="character", default=NULL, 
              help="path to Rdata file from NextFlow", metavar="character"),
    make_option(c("-s", "--sc_data_file"), type="character", default=NULL, 
              help="path to scRNA-seq data to infer cell types", metavar="character"),
    make_option(c("-n", "--n_threads"), type="character", default='4', 
              help="number of threads to use", metavar="character"),
	make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# normalize the bulk data from the path
# given to the script
dds_path = opt$rdata_file
load(dds_path)

bulk <- DESeq(dds)
bulk <- counts(dds, normalized=TRUE)
#head(bulk)

# make a scRNA-seq eset
scrna = opt$sc_data_file
scrna <- SingleCellExperiment(scrna)


# get bisque results
bis_res <- BisqueRNA::ReferenceBasedDecomposition(bulk, scrna, scrna$markers)

# get deconvSeq results
dec_res <- 
