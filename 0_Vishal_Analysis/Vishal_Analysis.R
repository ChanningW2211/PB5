library(ggplot2)
library(smotefamily)
library(dplyr)
library(scales)
library(e1071)
library(caret)
library(class)
library(xgboost)
library(zoo)
library(lubridate)
library(tseries)

accuracy <- function(x){sum(diag(x)/(sum(rowSums(x)))) * 100}
tpr <- function(x){x[2,2]/(x[2,1] + x[2,2])}
tnr <- function(x){x[1,1]/(x[1,1] + x[1,2])}
fnr <- function(x){x[2,1]/(x[2,1] + x[2,2])}
MinMax <- function(x,y) (x-min(y))/(max(y) - min(y))

water <- read.csv("./main/beach_main.csv", fileEncoding="UTF-8-BOM")

# water_raw <- read.csv("./clark.csv", fileEncoding="UTF-8-BOM")
# acf(waterq$Entero, na.action = na.pass)
# length(na.omit(waterq$Entero))
# adf.test(na.omit(waterq$Entero))

# plot(x = as.Date(waterq$Date, format="%d/%m/%Y"), y = waterq$Rain.24., pch=19,
#      type="l")
# 
# water1 <- water1[waterq$Class == 1,]
# 
# ggplot(waterq) + 
#   geom_point(aes(x = as.Date(Date, format="%d/%m/%Y"), y = Rain.24.,
#                  colour = factor(Class), alpha=factor(Class))) +
#   scale_alpha_discrete(range = c(0.8,1)) +
#   labs(x = "Date", y="Rainfall_24hrs")
# 
# 
# water1$Date = as.yearmon(water1$Date, "%d/%m/%Y")
# 
# ggplot(water1, 
#        aes(month(Date, label=TRUE, abbr=TRUE), 
#            Entero, group=factor(year(Date)), 
#            colour=factor(year(Date)))) +
#   geom_point() +
#   labs(x="Month", colour="Year")


nrow(water)
head(water)

#making test set
test <- water[water$Class == 1,]
test <- test[sample(nrow(test)),]
water_0 <- water[water$Class == 0,]
water_0 <- water_0[sample(nrow(water_0), nrow(test)), ]
test <-rbind(test, water_0)
test <- test[sample(nrow(test)),]

#create copy to remove from synthesised dataset
org_test <- data.frame(test)


# ggplot(test) +
#   geom_point(aes(x = Wspeed, y = Wdirection,
#                  colour = factor(Class), alpha=factor(Class))) +
#   scale_alpha_discrete(range = c(0.3,1))

#adasyn
adasyn = ADAS(water, water$Class,K=2)$data
adasyn$class <- NULL

adasyn = adasyn %>% anti_join(org_test)

#shuffle rows
set.seed(42)
rows <- sample(nrow(adasyn))
adasyn <- adasyn[rows, ]

#normalise

test[c(1,2,3,4,6)] <- MinMax(test[c(1,2,3,4,6)], adasyn[c(1,2,3,4,6)])
test$Wdirection <-  test$Wdirection / 360
test$SolarHours <- test$SolarHours / 24

#adasyn <- round(adasyn, digits = 1)

adasyn[c(1,2,3,4,6)] <- MinMax(adasyn[c(1,2,3,4,6)], adasyn[c(1,2,3,4,6)])
adasyn$Wdirection <-  adasyn$Wdirection / 360
adasyn$SolarHours <- adasyn$SolarHours / 24

labels <- adasyn[,8]

nrow(adasyn)
nrow(test)

#knn


pr <- knn(adasyn[,-8],test[,-8],cl=adasyn$Class,k=5)

#confusion matrix
tab <- table(pr,test$Class)
tab
  
tpr(tab)
tnr(tab)
fnr(tab)
accuracy(tab)

#svm method1

adasyn$Class <- as.factor(adasyn$Class)
actuals <- test[, 8]
test$Class <- as.factor(test$Class)

svmfit <- svm(Class ~ ., data = adasyn, 
              scale = FALSE, kernel = "radial", 
              cost = 1, gamma = 0.2622849, cross = 10)
pred <- predict(svmfit, test[1:7])
confusionMatrix(pred,test$Class)


#svm method 2

fitControl <- trainControl(method = "cv", number = 5)
C_range = sapply(seq(-1,3,0.0125), function(x) 10^x)
sigma_range = sapply(seq(-3,1,0.0125), function(x) 10^x)

fitGrid <- expand.grid(C = C_range, sigma = sigma_range)

svm_classifier = train(form = Class ~ ., data = adasyn, method = 'svmRadial')
                   #tuneGrid = fitGrid,
                   #trControl = fitControl,
                   #tuneLength = 10)

svm_classifier$bestTune

y_pred = predict(svm_classifier, newdata = test[-8])

cm <- table(test[, 8], y_pred)

tpr(cm)
tnr(cm)
fnr(cm)
accuracy(cm)

# cross-validation for svm

# folds = createFolds(adasyn$Class, k = 10)
# 
# cv = lapply(folds, function(x){
# 
#   training_fold = adasyn[-x, ]
#   test_fold = adasyn[x, ] 
# 
#   classifier = svm(formula = Class ~ .,
#                    data = training_fold,
#                    type = 'C-classification',
#                    kernel = 'radial')
#   
#   y_pred = predict(classifier, newdata = test_fold[-8])
#   cm = table(test_fold[, 8], y_pred)
#   accuracy = (cm[1,1] + cm[2,2]) / (cm[1,1] + cm[2,2] + cm[1,2] + cm[2,1])
#   return(accuracy)
# })  
# 
# accuracy = mean(as.numeric(cv))
# accuracy

# tune_out <- tune.svm(x = adasyn[1:7], y = adasyn$Class, 
#                      type = "C-classification", 
#                      kernel = "polynomial", degree = 2, cost = 10^(-1:2), 
#                      gamma = c(0.1, 2, 5), coef0 = c(0.1, 1, 5))

#xgboost - BDT

train_matrix <- as.matrix(adasyn[,-8])

bst <- xgboost(data = train_matrix, label = labels, 
               max.depth = 10, eta = 1, 
               nthread = 2, nrounds = 100,
               subsample = .8,
               objective = "binary:logistic",
               verbose = 0)

pred <- predict(bst, as.matrix(test[,-8]))

prediction <- as.numeric(pred > 0.5)

cm <- table(test[, 8], prediction)

tpr(cm)
tnr(cm)
fnr(cm)
accuracy(cm)
accuracy(tab)

#adaboost

adaboost <- ada(x=train_matrix, y=labels,
                loss = "logistic",
                type="gentle",
                iter=50,
                nu=1)

adaboost <- boosting(Class~., data=adasyn, boos=TRUE, mfinal=100)

pred <- predict(adaboost, test[,-8])

cm <- table(test[, 8], as.numeric(pred$class))

tpr(cm)
tnr(cm)
fnr(cm)
accuracy(cm)


# ANN

library(keras)
library(tensorflow)
use_condaenv("keras-tf", required = T)

train_x <- as.matrix(adasyn[,-8])
dimnames(train_x) <- NULL

test_x <- as.matrix(test[,-8])
dimnames(test_x) <- NULL

train_y <- adasyn[,8]
test_y <- test[,8]

train_y_oneh = to_categorical(train_y, num_classes = 2)
test_y_oneh = to_categorical(test_y, num_classes = 2)


model <- keras_model_sequential() 

model %>% 
  layer_dense(units = 8, activation = 'relu', 
              input_shape = c(7)) %>% 
  layer_dense(units = 2, activation = 'softmax')

model %>% compile(
  loss = 'categorical_crossentropy',
  optimizer = 'adam',
  metrics = 'accuracy'
)

history <- model %>% fit(
  train_x, 
  train_y_oneh, 
  epochs = 200, 
  batch_size = 5, 
  validation_split = 0.2
)

plot(history)

plot(history$metrics$loss, main="Model Loss", xlab = "epoch", ylab="loss", col="blue", type="l")
lines(history$metrics$val_loss, col="red")
legend("topright", c("train","test"), col=c("blue", "red"), lty=c(1,1))

plot(history$metrics$acc, main="Model Accuracy", xlab = "epoch", ylab="accuracy", col="blue", type="l")
lines(history$metrics$val_acc, col="red")
legend("bottomright", c("train","test"), col=c("blue", "red"), lty=c(1,1))


classes <- model %>% predict_classes(test_x, batch_size = 128)
table(test_y, classes)

score <- model %>% evaluate(test_x, test_y_oneh, batch_size = 128)
print(score)

save_model_hdf5(model, "ann_model.h5")

# ostsc
# lstm

# train.label <- 
# train.sample <- 
# test.label <- 
# test.sample <- 
# 
# 
# ostc_train <- OSTSC(train.sample, train.label, parallel = FALSE)
# ostc_test <- OSTSC(test.sample, test.label, parallel = FALSE)
# over_train.sample <- ostc_train$sample
# over_train.label <- ostc_train$label
# 
# over_test.sample <- ostc_test$sample
# over_test.label <- ostc_test$label
# 
# 
# model.over <- keras_model_sequential()
# model.over %>%
#     layer_lstm(10, input_shape = c(dim(over.x)[2], dim(over.x)[3])) %>%
#     layer_dropout(rate = 0.1) %>%
#     layer_dense(dim(over.y)[2]) %>%
#     layer_dropout(rate = 0.1) %>%
#     layer_activation("softmax")
#     layer_dropout(rate = 0.4) %>% 
#     layer_dense(units = 128, activation = "relu") %>%
#     layer_dropout(rate = 0.3) %>%
#     layer_dense(units = 10, activation = "softmax")
#     history.over <- LossHistory$new()
#     model.over %>% compile(
#     loss = "categorical_crossentropy",
#     optimizer = "adam",
#     metrics = c("accuracy",f1_score_0 = metric_f1_0, f1_score_1 = metric_f1_1)
# )
# 
# lstm.after <- model.over %>% fit(
#     x = over.x,
#     y = over.y,
#     validation_data=list(vali.x,vali.y),
#     callbacks = list(history.over),
#     epochs = 50
# )

#bls-smote

# blsmote = BLSMOTE(water,water$Class,K=5)$data
# blsmote$class <- NULL
# 
# blsmote = blsmote %>% anti_join(test)
# 
# set.seed(42)
# rows <- sample(nrow(blsmote))
# blsmote <- blsmote[rows, ]
# 
# table(blsmote$Class)
# 
# blsmote <- round(blsmote, digits = 2)
