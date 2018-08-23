
# Script to collect metadata about the environment in which 
# we are running 

# this script is run by the interpreter-under-test

# Command line args
# [ results file ] 

resultsFile <- commandArgs(TRUE)[1]

# Writes metadata to a DCF file


m <- list()


# Find the name of the interpereter in which this script is running

m$interpreterName <- if(is.null(version$engine)) {
  "GNU R"
} else {
  version$engine
}

# Find the version of the interpreter in which this script is running
m$interpreterVersion <- if(identical(version$engine, "Renjin")) {
  gsub(version$version.string, pattern = "Renjin version ", replacement = "", fixed = TRUE)
} else {
  paste(version$major, version$minor, sep=".")
}

m$jdkVersion <- if(identical(version$engine, "Renjin")) {
  sys <- import(java.lang.System)
  if(grepl(sys$getProperty("java.vm.name"), pattern="OpenJDK")) {
    vendor <- "OpenJDK"
  } else {
    vendor <- sys$getProperty("java.vendor")
    vendor <- gsub(vendor, pattern = "Oracle Corporation)", replacement = "Oracle")
  }
  sprintf("%s %s", vendor, sys$getProperty("java.version"))
} else {
  NA_character_
}

if(!is.na(resultsFile)) {
  write.dcf(m, resultsFile)
}
