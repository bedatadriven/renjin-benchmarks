
library(energy)

# Load a dataset from General Social Survey 2010
# Based on script from
# https://github.com/ajdamico/asdfree/tree/master/General%20Social%20Survey

load("GSS.2010.CS.rda")

# Few demographic variables such as sex, age, etc
#demo <- c(16, 4, 54, 55, 57, 67, 77, 78, 114)
demo <- c(44, 92)
# Attitude questions, not asked of all participants
attitudes <- c(146, 4459:4493, 4510:4656, 4692:4720, 4729:4789, 4805:4814, 4817:4823)


nr <- nrow(GSS.2010.CS.df)
vars <- attr(GSS.2010.CS.df, "var.labels")

# clean up 

# Find correlated pairs of variables
for(i in demo)  {
  x <- GSS.2010.CS.df[, i]
  cat(sprintf("finding associations between '%s' and attitudes:\n", vars[i]))
  for(j in attitudes) {
    y <- GSS.2010.CS.df[, j]
    valid <- !(is.na(x) | is.na(y))
    if(any(valid)) {
        result <- DCOR(as.double(x[valid]), as.double(y[valid]))
        if(!is.na(result$dCor)) {
            if(result$dCor > 0.25) {
              cat(sprintf("%.3f '%s' x '%s' (n=%d)\n", 
                result$dCor, vars[i], vars[j], sum(valid)))
            }
        }
    }
  }
}
