
library(randomForest)

load("trajectoryModel.RData")

features <- read.csv("training.csv", row.names = 1)

results <- predict(model, features)

# Make sure we output the results so they
# aren't optimized away
print(table(results))

