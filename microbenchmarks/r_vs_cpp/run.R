
source("../../harness/runners.R")


#' Times a benchmark on a single interpreter and variant
#' 
timeBenchmark <- function(benchmark, interpreter, variant = NA, executions = 3, inProcessIterations = 10) {
  
  dir.create("_results", showWarnings = FALSE)
  
  interpreterMetadata <- collectMetadata(interpreter)
  
  benchmarkVariant <- if(is.na(variant)) {
    benchmark
  } else {
    sprintf("%s_%s", benchmark, variant)
  }
  
  benchmarkScript <- sprintf("%s/%s.R", benchmark, benchmarkVariant)
  
  for(i in 1:executions) {
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

verify_benchmark <- function(benchmark, interpreter, variant = NA) {
  dir.create("_results", showWarnings = FALSE)
  
  benchmarkVariant <- if(is.na(variant)) {
    benchmark
  } else {
    sprintf("%s_%s", benchmark, variant)
  }
  
  benchmarkScript <- sprintf("%s/%s.R", benchmark, benchmarkVariant)
  outputFile <- file.path("_results", sprintf("%s_%s.rds", benchmarkVariant, interpreter$id))
  
  cat(sprintf("%20s %10s: ", benchmarkVariant, interpreter$id))
  
  success <- tryCatch({
    interpreter$runScript("../../harness/verify_microbenchmark.R", c(benchmarkScript, outputFile))
    cat("OK\n")
    TRUE
  },  error = function(e) {
    cat("ERROR\n")
    FALSE  
  })

  success
  
}

collectMetadata <- function(interpreter) {
  metadataFile <- tempfile()
  interpreter$runScript("../../harness/detect_metadata.R", metadataFile)
  c(interpreter$metadata, as.list(read.dcf(metadataFile)[1,]))
}


compareBenchmark <- function(benchmark, executions = 3, iterations = 10) {

  timeBenchmark(benchmark, gnur(), inProcessIterations = iterations)
  timeBenchmark(benchmark, gnur(), "rcpp", inProcessIterations = iterations)
  timeBenchmark(benchmark, renjinDev(), inProcessIterations = iterations)
  
}

suite <- c("curly", "dmvnorm", "fibonacci_seq", "fuzzycluster", "gibbs", "leibniz", )

verify_all <- function() {
  for(benchmark in suite) {
    verify_benchmark(benchmark, gnur())
  }
}

compareAll <- function(executions = 3, iterations = 10) {
  for(benchmark in suite) {}
}

