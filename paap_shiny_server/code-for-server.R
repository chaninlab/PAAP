# Loads the Model to memory
mod <- readRDS("Model.rds")

# Feature extraction for Testing set
xtest <- read.fasta('example.fasta', seqtype="AA", as.string = TRUE)###read data
xtest2 <- xtest[(sapply(xtest, protcheck))]###check special symbol
m2 = length(xtest2)
paactest <- matrix(nrow = m2, ncol = 23)
for(i in 1:m2){ 
  paactest[i, ]= extractPAAC(xtest2[[i]][1],lambda = 3, w = 0.1 , props = c("Hydrophobicity", "Hydrophilicity", "SideChainMass"))
}
data <- data.frame(paactest)

# Predicting unknown sequences
data.frame(Prediction= predict(mod,data))