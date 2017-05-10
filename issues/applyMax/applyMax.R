
load("maxslow.RData")

rowMax <- apply(fprob[1:nrow(fprob),],1,max)

actualMean <- mean(rowMax) 
expectedMean <- 0.8735726
diff <- abs(actualMean - expectedMean)

stopifnot(diff < 0.0001)
