#!/usr/bin/env Rscript

library(CompQuadForm)
setwd('~/Desktop/protein_clustering/Spring/Simulations/for_davies_alg')
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 
q_eigs <- as.vector(read.csv(args[1], header=FALSE)[['V1']])
q <- q_eigs[1]

eigs <- q_eigs[2:length(q_eigs)]
p_val_file<-file(args[2])
writeLines(as.character(c(davies(q, eigs, h=rep(1, length(eigs)))$Qq)), 
           p_val_file)
close(p_val_file)


