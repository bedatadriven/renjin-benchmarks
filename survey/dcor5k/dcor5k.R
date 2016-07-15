
library(energy)

# Load a dataset from General Social Survey 2010
# Based on script from
# https://github.com/ajdamico/asdfree/tree/master/General%20Social%20Survey

load("GSS.2010.CS.rda")

# GNU R can't handle more than 5k observations, so we need to 
# sample from the complete file for the purposes of comparison
GSS.2010.CS.df <- GSS.2010.CS.df[ sample(1:nrow(GSS.2010.CS.df), size = 5000), ]

nr <- nrow(GSS.2010.CS.df)
vars <- attr(GSS.2010.CS.df, "var.labels")

# clean up 

# Find correlated pairs of variables
# Skip the first 2 variables which are year, id, etc
# and limit to 400 variables otherwise the benchmark will
# run too long.
for(i in 3:400)  {
  for(j in (i+1):400) {
    x <- GSS.2010.CS.df[, i]
    y <- GSS.2010.CS.df[, j]
    valid <- !(is.na(x) | is.na(y))
    if(any(valid)) {
        result <- DCOR(as.double(x[valid]), as.double(y[valid]))
        if(!is.na(result$dCor)) {
            if(result$dCor > 0.6) {
              cat(sprintf("%.3f '%s' x '%s'\n", result$dCor, vars[i], vars[j]))
            }
        }
    }
  }
}
