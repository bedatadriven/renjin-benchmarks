library(ggplot2)

read_microbenchmark_results <- function(dir = "_results") {
  
  runs <- list.files(path = dir, pattern = ".+\\.run")
  runTimings <- sapply(runs, function(runFile) paste0(gsub(runFile, pattern = "\\.run$", replacement=""), ".timings"))
  completedRuns <- file.exists(file.path(dir, runTimings))
  
  # Find all metadata names
  metadataNames <- unique(unlist(recursive = TRUE, lapply(runs, function(runFile) {
    dcf <- read.dcf(file.path(dir, runFile))
    colnames(dcf)
  })))
  
  # Combine runs into a single table
  runTables <- lapply(seq_along(runs)[completedRuns], function(runIndex) {
    runFile <- runs[runIndex]
    timingFile <- file.path(dir, runTimings[runIndex])
    
    # Ensure we have complete metadata
    
    dcf <- read.dcf(file.path(dir, runFile))[1, ]
    metadata <- matrix(dcf[metadataNames], nrow=1)
    colnames(metadata) <- metadataNames
    metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
    
    # Correct metadata types
    metadata[, 'execution'] <- as.integer(metadata[, 'execution'])
    
    # Read timings data if it exists
    timings <- if(file.exists(timingFile)) {
      nanoseconds <- as.double(readLines(con = timingFile))
      data.frame(iteration = seq_along(nanoseconds), nanoseconds = nanoseconds)
    } else {
      data.frame(iteration = integer(0), nanoseconds = double(0))
    }
    cbind(metadata, timings)
  })
  
  # Combine into a single table
  
  do.call(rbind, runTables)
}

plot_runs <- function(results) {
  
  gg <- ggplot(results, aes(x=iteration, y=nanoseconds)) + 
    geom_point(aes(col=InterpreterVersion)) + 
    facet_grid(rows = vars(execution)) +
    labs(subtitle="Area Vs Population", 
         y="Nanoseconds", 
         x="In-process iteration")
  

  gg
}

