

library(randomForest)

## Train the model

data_apprenti=read.csv2("data_apprenti_TUFA.csv",sep=";",dec=",")
var_eject=c(2,3,4,5,6,7,8)
var_apprenti=data_apprenti[,-var_eject]

TUFA <- randomForest(i_TYPO~.,data=var_apprenti,ntree=100,mtry=7, replace=TRUE)

## Predict 

data_commune=read.csv2("tableau_commune.csv",sep=";",dec=",")
var_pred_com=data_commune[,-c(1,2,3,4,5,6,7)]
set.seed(7) 
typologie=predict(TUFA,var_pred_com,type="class")

# Write to disk to ensure the operation is not optimized away
write.csv(typologie, "typologie.csv")


