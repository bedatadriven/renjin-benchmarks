

# Dumps version properties of the given interpreter

rv <- R.Version()

props = list()
props$major <- rv$major
props$minor <- rv$minor
props$svn <- R.Version()[['svn rev']]
props$arch <- rv$arch

# Identify the version of blas being used
openFileOutput <- system(sprintf('lsof -p %d', Sys.getpid()), intern = TRUE)
openFiles <- unlist(strsplit(openFileOutput,split="\\s+"))
props$blas <- grep(openFiles, pattern="blas.*\\.so", value=T, ignore.case = T)

writeLines(paste(names(props), props, sep="="))
