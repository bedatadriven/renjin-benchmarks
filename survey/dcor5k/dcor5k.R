
library(energy)

# Load a dataset from General Social Survey 2010
# Based on script from
# https://github.com/ajdamico/asdfree/tree/master/General%20Social%20Survey

load("GSS.2010.CS.rda")

# GNU R can't handle more than 5k observations, so we need to 
# sample from the complete file for the purposes of comparison
GSS.2010.CS.df <- GSS.2010.CS.df[ sample(1:nrow(GSS.2010.CS.df), size = 5000), ]

nr <- nrow(GSS.2010.CS.df)
vars <- names(GSS.2010.CS.df)

# Find correlated pairs of variables
# Skip the first 10 variables which are year, id, etc
for(i in 10:nr)  {
  for(j in (i+1):nr) {
    x <- GSS.2010.CS.df[, i]
    y <- GSS.2010.CS.df[, j]
    valid <- !(is.na(x) | is.na(y))
    result <- DCOR(as.double(x[valid]), as.double(y[valid]))
    if(result$dCor > 0.5) {
      cat(sprintf("%.3f %15s x %15s\n", result$dCor, vars[i], vars[j]))
    }
  }
}
