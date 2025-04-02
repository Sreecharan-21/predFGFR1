library(randomForest)
library(caret)
load("data/random_forest_model.rda")
load("data/process1.rda")
command <- "java -jar data/padel/bin.jar -2d -fingerprints -removesalt -retainorder -detectaromaticity -standardizenitro -dir inputmols/inputmols.smi -file data/temp/out.csv"
system(command, ignore.stdout = TRUE)
descriptors <- read.csv("data/temp/out.csv", header = TRUE)
descriptors <- descriptors[, -1]
descriptors <- as.data.frame(sapply(descriptors, as.numeric))
descriptors <- descriptors[colnames(process1$ranges)[-1]]
descriptors <-  cbind(IC50 = sample(2,size = nrow(descriptors),replace = T),descriptors)
descriptors_normalized <- predict(process1, descriptors)
descriptors_normalized <- descriptors_normalized[colnames(rfmer$trainingData)[-18]]
predictions <- predict(rfmer, descriptors_normalized)
predictions <- as.character(predictions)
predictions[predictions=="1"] <- 'Inhibitor'
predictions[predictions=="0"] <- 'Non-Inhibitor'
write.csv(cbind(Molecules = paste('mol',1:nrow(descriptors)),predictions), 'Results/predictions.csv', row.names = F)
