library(fastGenoGAM)

## ================
## General workflow
## ================

## data
config <- "config_h3k4me3.txt"

## specify number of cores
BiocParallel::register(BiocParallel::SnowParam(worker = 1))

## Read in data, the chunk and the overhang give the tile,
## the units into which the entire genome will be cut for parallel computation
## While chunks are disjoint, tiles are overlapping.
chunkSize <- 20000 
overhangSize <- 1000 

## set custom HDF5 folder if needed, otherwise it will use the
## default dump folder of HDF5Array
settings <- GenoGAMSettings(hdf5Control = list(dir = "sacCer3_h5"))#,
                            #chromosomeList = c("any chromosome")) ## restrict to only some chromosomes if needed, otherwise leave argument out.

## with HDF5
ggd <- GenoGAMDataSet(config, design = ~ s(x) + s(x, by = genotype), directory = folder, hdf5 = TRUE, chunkSize = chunkSize, settings = settings, overhangSize = overhangSize)

## size factor
ggd <- computeSizeFactors(ggd)

## Run 
lambda <- NULL ## set to NULL to use Cross Validation
theta <- NULL ## set to NULL to use Cross Validation
result <- genogam(ggd, intervalSize = 200, lambda = lambda, theta = theta)



