# human liver cohort
# synapse account
# user=synapsedata
# pwd=synapse
# ieuan.clay@gmail.com
# June 2015

### set up session
#rm(list=ls())

# reproducibility
set.seed(8008)
## packages
# CRAN
library(e1071)
## Blocks for timing


data <- lapply(dir(".", pattern = ".txt$", full.names = TRUE), read.delim, stringsAsFactors=FALSE)
names(data) <- dir(".", pattern = ".txt$", full.names = FALSE)

#file.remove(dir(".", pattern = ".txt$", full.names = TRUE))

# check loaded data
cat(str(data))

## reformat loaded data
# bind columns together (each row in file is one column)
tmp <- data.frame(t(rbind(data$phenotype.txt)), stringsAsFactors = FALSE)
tmp$indvidual_id <- rownames(tmp)
names(tmp)<-tmp[1,] # first line in file is column names
tmp <- tmp[2:nrow(tmp),] # drop header

# convert numeric columns from character
numeric_cols <- c(1, 5:18)
for (numeric_col in numeric_cols){
  tmp[,numeric_col] <- as.numeric(tmp[,numeric_col])
  tmp[is.na(tmp[,numeric_col]), numeric_col] <- 0 # replace all NAs with 0s, i know this is horrible.
}

liverdata <- tmp


do.svm <- function(liverdata){
  cat("> START: do.svm()\n")
  ### run some simple predictive modelling on liver cohort clinical data
  # see integration/humanLiverCohort.R for more advanced calculations

  results <- list() # placeholder

  ## SVM on gender and liver stats using liver enzyme data
  # set test/train, sampling 1/3 rows as test
  traintest <- rep(TRUE, nrow(liverdata))
  traintest[sample(1:length(traintest), size = round(nrow(liverdata)/3), replace = FALSE)] <- FALSE

  train <- liverdata[traintest,]
  test <- liverdata[!traintest,]
  livercols <- 9:18
  model <- svm(x = as.matrix(train[,livercols]), # all liver enzme activiy stats
               y=factor(train$GENDER, levels = c("Male", "Female")),
               scale = TRUE, type = "C")
  cat(">>> DONE: svm()\n")
  # test model on training set
  pred <- predict(model, as.matrix(train[,livercols]))
  cat(">>> DONE: predict()\n")
  res <- classAgreement(table(pred, factor(train$GENDER, levels = c("Male", "Female"))))
  results <- append(results,
                    list(data.frame(
                      dat="gendertrain",
                      var=names(res),
                      coeff=unlist(res)
                    ))
  )

  # and on the testset
  pred <- predict(model,as.matrix(test[,livercols]))
  cat(">>> DONE: predict()\n")
  res <- classAgreement(table(pred, test$GENDER))

  results <- append(results,
                    list(data.frame(
                      dat="gendertest",
                      var=names(res),
                      coeff=unlist(res)
                    ))
  )

  ## liver triglycerides
  model <- svm(x = as.matrix(train[,livercols]), # all liver enzme activiy stats
               y=train$`Liver_Triglyceride_(mg_per_dL)`,
               scale = TRUE, type = "eps-regression")
  cat(">>> DONE: svm()\n")
  # test model on training set
  train$pred <- predict(model, as.matrix(train[,livercols]))
  # and on the testset
  test$pred <- predict(model,as.matrix(test[,livercols]))

  # return results
  results <- append(results,
                    list(data.frame(
                      dat=c("triglyceridesTrain", "triglyceridesTest"),
                      var="spearman",
                      coeff=c(
                        cor(x=train$`Liver_Triglyceride_(mg_per_dL)`, y=train$pred, method = "spearman"),
                        cor(x=test$`Liver_Triglyceride_(mg_per_dL)`, y=test$pred, method = "spearman")
                        )
                    ))
  )

  cat("> END: do.svm()\n")
  return(do.call("rbind",results))


}

### reporting
res <- do.svm(liverdata)
tail(res)

# final clean up
#rm(list=ls())
#gc()
