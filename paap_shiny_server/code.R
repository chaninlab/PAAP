#######set directory
setwd("~/Dropbox/Antihypertensive PAAP")
#######Load package
library(protr)
library(seqinr)
library(randomForest)


#######Building RF model
x <- read.fasta('train.fasta', seqtype="AA", as.string = TRUE)
D = read.csv("label.csv", header = TRUE) 
m = length(x)
paac <- matrix(nrow = m, ncol = 23)

for(i in 1:m){ 
paac[i, ]= extractPAAC(x[[i]][1],lambda = 3, w = 0.1 , props = c("Hydrophobicity", "Hydrophilicity", "SideChainMass"))
}

training = data.frame(paac,Class = D[,ncol(D)])

Model = randomForest(Class ~ ., training, ntree= 400,mtry =3)

#######Feature extraction for Testing set
xtest <- read.fasta('example.fasta', seqtype="AA", as.string = TRUE)###read data
xtest2 <- xtest[(sapply(xtest, protcheck))]###check special symbol
m2 = length(xtest2)
paactest <- matrix(nrow = m2, ncol = 23)
for(i in 1:m2){ 
paactest[i, ]= extractPAAC(xtest2[[i]][1],lambda = 3, w = 0.1 , props = c("Hydrophobicity", "Hydrophilicity", "SideChainMass"))
}
data <- data.frame(paactest)

#######Predicting unknown sequences
data.frame(Prediction= predict(Model,data))



####
saveRDS(Model, "Model.rds")

mod <- readRDS("Model.rds")

