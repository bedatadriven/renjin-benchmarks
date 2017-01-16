# Copyright (c) 2015 Andre Verissimo (andre.verissimo@tecnico.ulisboa.pt)
# Copyright (c) 2016 BeDataDriven B.V.
# License: ...
#
## Survival analysis

##### set up session #####
rm(list = ls())
START_WORKFLOW <- as.numeric(Sys.time())
# reproducibility
set.seed(8008)

##### loading packages #####
library(glmnet)
library(survival)
##### Set global vars #####
INPUT <- "survival_tcga_v1.csv"

# parameters for survival analysis
params = list();
params$nlambda <- 10000 # number of different lambdas to be tested
#
# commented out as they are not being used
#params$lam_max <- 1e-7
#params$lam_rat <- 0.001
#params$lambda  <- seq(params$lam_max, params$lam_rat * params$lam_max, length.out=params$nlambda)
#
params$alpha   <- seq(0, 1, 0.1)
params$nalpha  <- length(params$alpha)

##### Blocks for timing #####

  dat <- read.table(INPUT, sep = ",", row.names = 1, header = TRUE)
  #
  var_arr = 3:ncol(dat)
  # get the variables from data
  #  ydata is a data.frame keeps the status of the patient and time of last follow-up
  ydata     <- cbind(time = dat$time, status = dat$status)
  #  xdata keeps the gene expression for each patient
  xdata     <- dat[ , var_arr]
  #
  surv_data <- list(ydata = ydata, xdata = xdata)

  cat('Starting calculation...\n')
  #
  alpha_vec  <- array(0, params$nalpha)
  #
  xdata <- as.matrix(surv_data$xdata)
  ydata <- surv_data$ydata
  #
  # create results structure, as it was not possible to determine before the first loop
  my_results = list()
  # for each alpha and lambda, determine the glmnet
  for (mm in 1:length(params$alpha)) {
    # set the alpha value
    alpha_v = params$alpha[mm]
    # get local results
    temp_results = glmnet(xdata, ydata, family = 'cox', alpha = alpha_v, nlambda = params$nlambda, standardize = FALSE )
    # save results
    item = list()
    item$lambda = temp_results$lambda
    item$beta   = temp_results$beta
    my_results[[mm]] <- item
  }
  
  alpha_vec  <- params$alpha
  cat("End calculation.\n")
  res <- list( my_results = my_results, alphas = alpha_vec)

  print(str(res))

END_WORKFLOW <- as.numeric(Sys.time())
TOTAL_TIME <- END_WORKFLOW - START_WORKFLOW
print(TOTAL_TIME)
write(TOTAL_TIME, file = "TIMINGS", append = TRUE)

# final clean up
rm(list = ls())
gc()
