library(pdist)
library(MASS)
library(fields)

water <- read.csv("./PB5_trainset.csv", fileEncoding="UTF-8-BOM")

colnames(water)

waterC = water[,-c(1,2,6,10, 12,13,14,15)]

colnames(waterC)

waterC$Class <- water$Entero_level


MinMax <- function(x,y) (x-min(y))/(max(y) - min(y))

waterC[c(1,2,3,5)] <- MinMax(waterC[c(1,2,3,5)], waterC[c(1,2,3,5)])

waterC$Wdirection <-  waterC$Wdirection / 360
waterC$Solarhours <- waterC$Solarhours / 24

#waterC <- abs(waterC)

# generate the new training samples


#OversampledTrain = JillOSTSC(as.matrix(waterC[,-9]), waterC$Class, parallel = FALSE, per=0.50)
OversampledTrain = OSTSC(as.matrix(waterC[,-9]), waterC$Class, parallel = FALSE, per=0.50)

osd <- OversampledTrain$sample

osd <-  as.data.frame(osd)

colnames(osd) <- colnames(waterC[,-9])

osd$Class <- OversampledTrain$label


# bin the categorical values into 0,1,2 - on_offshore, 0,1 - beachtype
osd[osd < 0] <- 0

osd$BeachType[osd$BeachType < 0.5] <- 0
osd$BeachType[0.5 < osd$BeachType] <- 1

osd$on_offshore[osd$on_offshore < 0.5] <- 0
osd$on_offshore[0.5 < osd$on_offshore & osd$on_offshore < 1.5] <- 1
osd$on_offshore[1.5 < osd$on_offshore] <- 2

head(osd)
nrow(osd)

write.csv(osd,"overSampledDataset.csv", row.names = FALSE)

a_ = do.call("paste", waterC[,1:9])
b_ = do.call("paste", osd[,1:9])


waterC$InB = a_ %in% b_

table(waterC$InB)

head(osd[b_ %in% a_,c(1,2,3,4,5,6,7)],20)

index <- as.numeric(rownames(osd[b_ %in% a_,c(1,2)]))

plot(1:length(index), index)

plot(1:nrow(osd[b_ %in% a_,]),osd[b_ %in% a_,]$Class)
plot(1:nrow(osd),osd$Class)
plot(1:nrow(waterC), waterC$Class)

plot(1:nrow(osd), osd$BeachType)