args <- commandArgs(TRUE)
inputfilename <- as.character(args[1])
chr <- as.character(args[2])
poolroot <- as.character(args[3])
founderfile <- as.character(args[4])
outdir <- as.character(args[5])

library(limSolve)
source("scripts/haplotyper.limSolve.code.R")
runscan(inputfilename, chr, poolroot, founderfile, outdir)

