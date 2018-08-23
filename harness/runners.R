
# Defines runners for GNU R and Renjin

gnur <- function() {
  
  bin <- file.path(R.home(), "bin", "Rscript")
  runScript <- function(scriptFile, scriptArgs) {
    system2(bin, c(scriptFile, scriptArgs))  
  }
  
  list(id = "gnu", runScript = runScript)
  
}

renjin <- function(bin = "renjin", loop.compiler = TRUE, vm.options) {
  
  env <- character(0)
  if(!missing(vm.options)) {
    env <- sprintf('RENJIN_OPTS="%s"', paste(vm.options, collapse = " "))
  }
  
  optimFlags <- if(loop.compiler) "--compile-loops" else NULL
  
  runScript <- function(scriptFile, scriptArgs) {
    system2(command = bin, args = c("-f", scriptFile, optimFlags, "--args", scriptArgs), env = env)
  }
  list(
    id = "renjin", 
    runScript = runScript, 
    metadata = c(LoopCompilerEnabled = loop.compiler))
}

renjinDev <- function(repo.path = "/home/alex/dev/renjin", ...) {
  renjin(bin = file.path(repo.path, "dist", "generic", "target", "verify", "renjin-1.0.0-SNAPSHOT", "bin", "renjin"), ...)
}

