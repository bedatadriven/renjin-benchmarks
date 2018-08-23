
source("../../harness/runners.R")


#' Times a benchmark on a single interpreter and variant
#' 
timeBenchmark <- function(benchmark, interpreter, variant = NA, inProcessIterations = 10) {
  
  dir.create("_results", showWarnings = FALSE)
  
  interpreterMetadata <- collectMetadata(interpreter)
  
  benchmarkVariant <- if(is.na(variant)) {
    benchmark
  } else {
    sprintf("%s_%s", benchmark, variant)
  }
  
  benchmarkScript <- sprintf("%s/%s.R", benchmark, benchmarkVariant)
  
  for(i in 1:3) {
    timingsFile <- file.path("_results", sprintf("%s_%s_%d.timings", benchmarkVariant, interpreter$id, i))
    metadataFile <- file.path("_results", sprintf("%s_%s_%d.run", benchmarkVariant, interpreter$id, i))
    
    startTime <- Sys.time()
    
    message(sprintf("Starting execution #%d of %s on %s...", i, benchmarkVariant, interpreter$id))
    
    interpreter$runScript("../../harness/time_microbenchmark.R", c(benchmarkScript, timingsFile, inProcessIterations))
    
    executionMetadata <- list(
      benchmark = benchmark,
      benchmarkVariant = variant,
      execution = i,
      startTime = as.character(startTime))
    
    write.dcf(file = metadataFile, c(interpreterMetadata, executionMetadata))
  }
}

collectMetadata <- function(interpreter) {
  metadataFile <- tempfile()
  interpreter$runScript("../../harness/detect_metadata.R", metadataFile)
  c(interpreter$metadata, as.list(read.dcf(metadataFile)[1,]))
}


compareBenchmark <- function(benchmark) {

  timeBenchmark("compute_pi", gnur())
  timeBenchmark("compute_pi", gnur(), "rcpp")
  timeBenchmark("compute_pi", renjinDev())
  
}