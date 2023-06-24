# Import of descriptor dataset dataset
des = read.csv("C:\\Users\\Sreecharan\\Desktop\\Onlydes.csv")
View(des)
nod1 = des[,-1] # removal of character column
View(nod1)
anyNA(nod1) # check wether any null values are present or not

#splitting the data-set into training and test
set.seed(21)
sp1 = sample(nrow(nod1), nrow(nod1)*0.7)

train1 = nod1[sp1,] # training set data
test1 = nod1[-sp1,] # test set data
anyNA(train1)

#normalization
library(caret)
process1 = preProcess(train1, method = c("range"))

norm_train1 = predict(process1, train1) #Train normalization
norm_train1

norm_test1 = predict(process1, test1) #Test normalization
norm_test1

#check for Zero-variance
library(caret)
library(lattice)

v_train1 = preProcess(norm_train1,
                      method =  "zv",
                      thresh = 0.0)
#removing the zero-variance from training and test
rmv_train1 = predict(v_train1, norm_train1)
rmv_test1 = predict(v_train1, norm_train1)

#checking of correlation

cor_matrix1 = cor(rmv_train1)
cor_matrix1
cutoff1 = 0.85
highlycorrelated1 = findCorrelation(cor_matrix1, cutoff = cutoff1)
highlycorrelated1

#remove highly correlated variables for train data
new_cor_train1 = rmv_train1[,-highlycorrelated1]
new_cor_train1
#remove highly correlated variables for test data
new_cor_test1 = rmv_test1[,-highlycorrelated1]
new_cor_test1

#feature selection
library(Boruta)
library(dplyr)
set.seed(213)
Boruta_train13 = Boruta(as.factor(new_cor_train1$IC50) ~ ., data = new_cor_train1,
                        doTrace = 2, maxRuns = 500)
Boruta_train13
print(Boruta_train13)
sp33 = attStats(Boruta_train13)
sele_boruta1 = getSelectedAttributes(Boruta_train13,withTentative = FALSE)
Boruta_imps3 = select(new_cor_train1,MPC8,	nAromBond,	MDEN.23,	SpMax1_Bhp,
                      SCH.7,	naaaC,	AATS1i,	R_TpiPCTPC,	nAtomP,	ETA_BetaP_s,
                      SM1_Dzi,	nN,	AATS0i,	SpMin1_Bhi,	SpMax2_Bhv,	AATS5i,
                      VCH.5,	SpMin4_Bhi,	SpMin3_Bhi,	ATSC3i,	AATSC4e,
                      SpMin2_Bhe,	SpMax1_Bhm,	SpMax8_Bhp,	minsNH2,	SpMin2_Bhs,
                      GGI5,	maxHother,	ATSC0c,	MDEN.22,	IC2,	GATS1c,	AATS6i,
                      GGI3,	minsssN,	MDEN.33,	MDEC.13,	GGI7,	JGI8,	SRW9,
                      ETA_dPsi_A,	SpMax2_Bhm,	TopoPSA,	MDEC.23,	ETA_EtaP_F,
                      MDEC.33,	GATS4e,	SpMin2_Bhm,	ATSC2c,	VPC.6,	ATSC0i,
                      minHBa,	minaaN,	GATS3i,	MLFER_BH)

plot(Boruta_train13, main = "Boruta Important Variables",bg = "gray",
     ylim = c(-5,15),lwd = 0.5)

library('xlsx')
write.xlsx2(sp33,"C:\\Users\\Sreecharan\\Desktop\\sp33.xlsx")

# feature  selection by random forest 

library(randomForest)
rf_train1s<-randomForest(as.factor(new_cor_train1$IC50) ~ ., 
                         data = new_cor_train1,importance=TRUE, ntree = 500)
rf_train1s
rf.imp1s<-importance(rf_train1s)

rf_imp113 = select(new_cor_train1, MPC8,	nAromBond,	MDEN.23,	nN,	SpMax1_Bhp,
                   naaaC,	MDEN.33,	SCH.7,	R_TpiPCTPC,	SM1_Dzi,	VCH.5,
                   minsssN,	IC2,	nRing,	SaasN,	SpMax2_Bhv,	SpMin2_Bhe,	AATS1i,
                   GGI7,	AATS0i,	SpMin3_Bhi,	ATSC2c,	ETA_BetaP_s,	nAtomP,
                   ATS7s,	SpMin4_Bhi,	SpMax8_Bhp,	MDEC.23,	SpMin1_Bhi,	ATSC3i,
                   SRW9,	GGI5,	ATSC0i,	GGI9,	maxHother,	AATSC1m,	GATS4e,
                   hmin,	nHBAcc,	AATSC4e,	AATSC1c,	nHeteroRing,	SpMax1_Bhm,
                   GATS3i,	SpMin2_Bhs,	minssO,	ATS4m,	MDEC.13,	GGI3,	VPC.6,
                   MLFER_A,	ETA_BetaP_ns_d,	minHBa,	SpMin2_Bhm,	MDEN.22)

# Get importance values as a data frame
imp1s = as.data.frame(importance(rf_train1s))
imp1s = cbind(vars=rownames(imp1s), imp1s)
imp1s = imp1s[order(imp1s$MeanDecreaseGini),]
imp1s$vars = factor(imp1s$vars, levels=unique(imp1s$vars))

barplot(imp1s$MeanDecreaseGini, names.arg=imp1s$vars, 
        xlim = c(10,450),ylim = c(0,1.5),col = "skyblue",
        xlab = "Attributes",ylab = "MeanDecreaseGini")

# Writing data_frame to excel sheet
write.xlsx2(rf.imp1s,"C:\\Users\\Sreecharan\\Desktop\\rf.imp1s3.xlsx")

#Recursive feature elimination
library(caret)
set.seed(212)
# Performing RFE
rfe_train1s <- rfe(new_cor_train1,as.factor(new_cor_train1$IC50),
                   sizes = c(56),
                   rfeControl = rfeControl(functions = treebagFuncs))
result_rfe1s = print(rfe_train1s)
ref_imp1s= predictors(rfe_train1s)
rfes1s = subset(new_cor_train1,select = ref_imp1s[-1])
# Writing data_frame to excel sheet
write.xlsx2(rfes1s,"C:\\Users\\Sreecharan\\Desktop\\rfe113.xlsx")

common3 = select(new_cor_train1,MPC8,nAromBond,nN,MDEN.23,MDEN.33,minsssN,
                 SpMax1_Bhp,SM1_Dzi,naaaC,MDEC.23,SCH.7,IC2,VCH.5,SpMax2_Bhv,
                 SpMin4_Bhi,GGI7,SpMin3_Bhi,GGI5,SpMin2_Bhe,SpMax8_Bhp,
                 R_TpiPCTPC,AATS0i,AATS1i,AATS5i, AATS6i)

library(Boruta)
set.seed(213)
Boruta_train13 = Boruta(as.factor(new_cor_train1$IC50) ~ ., data = common3,
                        doTrace = 2)
Boruta_train13
kk1 = attStats(Boruta_train13) 
sele_boruta1 = getSelectedAttributes(Boruta_train13, withTentative = FALSE) # Only Imp
Boruta_imp13 = select(new_cor_train1,MPC8,nAromBond,nN,MDEN.23,MDEN.33,minsssN,
                      SpMax1_Bhp,SM1_Dzi,naaaC,MDEC.23,SCH.7,IC2,VCH.5,SpMax2_Bhv,
                      SpMin4_Bhi,GGI7,SpMin3_Bhi,GGI5,SpMin2_Bhe,SpMax8_Bhp,
                      R_TpiPCTPC,AATS0i,AATS1i,AATS5i, AATS6i)
Boruta_imp13

# model buildiing code
library(e1071)
library(caret)
B_training_data1 = select(new_cor_train1, MPC8,nAromBond,nN,MDEN.23,MDEN.33,minsssN,
                          SpMax1_Bhp,SM1_Dzi,naaaC,MDEC.23,SCH.7,IC2,VCH.5,SpMax2_Bhv,
                          SpMin4_Bhi,GGI7,SpMin3_Bhi,GGI5,SpMin2_Bhe,SpMax8_Bhp,
                          R_TpiPCTPC,AATS0i,AATS1i,AATS5i, AATS6i)

B_testing_data1 = select(new_cor_test1,MPC8,nAromBond,nN,MDEN.23,MDEN.33,minsssN,
                         SpMax1_Bhp,SM1_Dzi,naaaC,MDEC.23,SCH.7,IC2,VCH.5,SpMax2_Bhv,
                         SpMin4_Bhi,GGI7,SpMin3_Bhi,GGI5,SpMin2_Bhe,SpMax8_Bhp,
                         R_TpiPCTPC,AATS0i,AATS1i,AATS5i, AATS6i)

# Randomforest
library(randomForest)
controlmr <- trainControl(method="repeatedcv", number=10)
tunegridr <- expand.grid(.mtry= c (1:25), .ntree=c(100 : 500))
set.seed = 123
rfmer <- train(B_testing_data1,as.factor(new_cor_train1$IC50),method = "rf",
               tuneGRid=tunegridr, trControl = controlmr)
rfmer
plot(rfmer)
# Predicting tuned train data
rfmr = as.factor(new_cor_train1$IC50)
rfmr_train_pred = predict(rfmer, B_training_data1)
caret::confusionMatrix(rfmr_train_pred,rfmr)
F1_Score(rfmr_train_pred,as.factor(new_cor_train1$IC50))
mltools::mcc(rfmr_train_pred,as.factor(new_cor_train1$IC50))
# Predicting tuned test data
rfmtr = as.factor(new_cor_test1$IC50)
rfmtr_test_pred = predict(rfmer,B_testing_data1)
confusionMatrix(rfmtr_test_pred,rfmtr)
F1_Score(rfmtr_test_pred,as.factor(new_cor_test1$IC50))
mltools::mcc(rfmtr_test_pred,as.factor(new_cor_test1$IC50))


