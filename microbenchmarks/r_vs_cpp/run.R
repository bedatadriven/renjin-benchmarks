
source("../../harness/runners.R")


#' Times a benchmark on a single interpreter and variant
#' 
timeBenchmark <- function(benchmark, interpreter, variant = NA) {
  
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
    
    interpreter$runScript("../../harness/time_microbenchmark.R", c(benchmarkScript, timingsFile))
    
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
  as.list(read.dcf(metadataFile)[1,])
}

writeMetadata <- function(metadataFile, interpreter, benchmark, variant, executionIndex) {
  m <- list(
    interpreter = intepreter$id,
    benchmark = benchmark,
    variant = variant,
    execution = executionIndex
  )
}

compareBenchmark <- function(benchmark) {

  # Ensure results folder exists
  dir.create("_results")
  
  
}