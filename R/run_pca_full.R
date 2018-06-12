missing.dataset <- as.matrix(read.table("~/Projects/PPCA/PCAMV/data4matlab.txt", sep="\t"))
inds            <- as.matrix(read.table("~/Projects/PPCA/PCAMV/inds4matlab.txt", sep="\t"))
source("~/Projects/PPCA/ppcaNet/R/pca_full_pdwk.R")
Rcpp::sourceCpp('~/Projects/PPCA/ppcaNet/src/pca_updates.cpp')
X <- missing.dataset
missing.dataset[inds] <- NaN
ncomp <- 2
tmp <- pca_full(missing.dataset, NA, "vb")
