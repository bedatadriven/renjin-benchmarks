
library(survey)

svydata <- readRDS("california.rds")
names(svydata) <- tolower(names(svydata))
 
svydsgn <- svrepdesign(
  weight = ~pwgtp ,
  repweights = 'pwgtp[1-9]' ,
  scale = 4 / 80 ,
  rscales = rep( 1 , 80 ) ,
  mse = TRUE ,
  data = svydata) 
  
agep <- svymean(~agep, svydsgn, se=TRUE)
relp <- svymean(~adjinc, svydsgn, se=TRUE)

print(agep)
print(relp)

