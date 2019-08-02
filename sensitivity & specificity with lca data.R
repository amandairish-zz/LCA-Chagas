# try to calc sens/spec

library(epiR)
library(tidyverse)

# load data
data <- read.csv("/Users/amandairish/Desktop/LCA/dxcomp_041019_csv.csv")

# keep relevant columns
data.epir <- data %>%
  select(INDEX_ID, O_result, INB1_CUT0, INB2_CUT0, H_7POS, W_7POS, ARC_STAT,
         NPOS, SERO2POS, OHWI, Birth)

# calculate terms for use with epiR package

# create indicator var for true/false pos/neg with ARC_STATUS as gold standard
data.epir$O.true.pos <- ifelse(data.epir$O_result==1 & data.epir$ARC_STAT==1, 1, 0) # 1 if cond met, 0 if not
data.epir$O.false.pos <- ifelse(data.epir$O_result==1 & data.epir$ARC_STAT==0, 1, 0)
data.epir$O.false.neg <- ifelse(data.epir$O_result==0 & data.epir$ARC_STAT==1, 1, 0)
data.epir$O.true.neg <- ifelse(data.epir$O_result==0 & data.epir$ARC_STAT==0, 1, 0) 

data.epir$I1.true.pos <- ifelse(data.epir$INB1_CUT0==1 & data.epir$ARC_STAT==1, 1, 0)
data.epir$I1.false.pos <- ifelse(data.epir$INB1_CUT0==1 & data.epir$ARC_STAT==0, 1, 0)
data.epir$I1.false.neg <- ifelse(data.epir$INB1_CUT0==0 & data.epir$ARC_STAT==1, 1, 0)
data.epir$I1.true.neg <- ifelse(data.epir$INB1_CUT0==0 & data.epir$ARC_STAT==0, 1, 0) 

data.epir$I2.true.pos <- ifelse(data.epir$INB2_CUT0==1 & data.epir$ARC_STAT==1, 1, 0)
data.epir$I2.false.pos <- ifelse(data.epir$INB2_CUT0==1 & data.epir$ARC_STAT==0, 1, 0)
data.epir$I2.false.neg <- ifelse(data.epir$INB2_CUT0==0 & data.epir$ARC_STAT==1, 1, 0)
data.epir$I2.true.neg <- ifelse(data.epir$INB2_CUT0==0 & data.epir$ARC_STAT==0, 1, 0) 

data.epir$H.true.pos <- ifelse(data.epir$H_7POS==1 & data.epir$ARC_STAT==1, 1, 0)
data.epir$H.false.pos <- ifelse(data.epir$H_7POS==1 & data.epir$ARC_STAT==0, 1, 0)
data.epir$H.false.neg <- ifelse(data.epir$H_7POS==0 & data.epir$ARC_STAT==1, 1, 0)
data.epir$H.true.neg <- ifelse(data.epir$H_7POS==0 & data.epir$ARC_STAT==0, 1, 0) 

data.epir$W.true.pos <- ifelse(data.epir$W_7POS==1 & data.epir$ARC_STAT==1, 1, 0)
data.epir$W.false.pos <- ifelse(data.epir$W_7POS==1 & data.epir$ARC_STAT==0, 1, 0)
data.epir$W.false.neg <- ifelse(data.epir$W_7POS==0 & data.epir$ARC_STAT==1, 1, 0)
data.epir$W.true.neg <- ifelse(data.epir$W_7POS==0 & data.epir$ARC_STAT==0, 1, 0) 


# sum up # true pos/neg
nO.true.pos <- sum(data.epir$O.true.pos)
nO.false.pos <- sum(data.epir$O.false.pos)
nO.false.neg <- sum(data.epir$O.false.neg)
nO.true.neg <- sum(data.epir$O.true.neg)

nI1.true.pos <- sum(data.epir$I1.true.pos)
nI1.false.pos <- sum(data.epir$I1.false.pos)
nI1.false.neg <- sum(data.epir$I1.false.neg)
nI1.true.neg <- sum(data.epir$I1.true.neg)

nI2.true.pos <- sum(data.epir$I2.true.pos)
nI2.false.pos <- sum(data.epir$I2.false.pos)
nI2.false.neg <- sum(data.epir$I2.false.neg)
nI2.true.neg <- sum(data.epir$I2.true.neg)

nH.true.pos <- sum(data.epir$H.true.pos)
nH.false.pos <- sum(data.epir$H.false.pos)
nH.false.neg <- sum(data.epir$H.false.neg)
nH.true.neg <- sum(data.epir$H.true.neg)

nW.true.pos <- sum(data.epir$W.true.pos)
nW.false.pos <- sum(data.epir$W.false.pos)
nW.false.neg <- sum(data.epir$W.false.neg)
nW.true.neg <- sum(data.epir$W.true.neg)


# create tables with ARC status as "gold standard"

# O test
dat.O <- as.table(matrix(c(nO.true.pos, nO.false.pos, nO.false.neg, nO.true.neg),
                         nrow = 2, byrow = TRUE))
colnames(dat.O) <- c("Dis+","Dis-")
rownames(dat.O) <- c("Test+","Test-")
sens.spec.0 <- epi.tests(dat.O, conf.level = 0.95)
print(sens.spec.0); summary(sens.spec.0)


# I1 test
dat.I1 <- as.table(matrix(c(nI1.true.pos, nI1.false.pos, nI1.false.neg, nI1.true.neg),
                         nrow = 2, byrow = TRUE))
colnames(dat.I1) <- c("Dis+","Dis-")
rownames(dat.I1) <- c("Test+","Test-")
sens.spec.I1 <- epi.tests(dat.I1, conf.level = 0.95)
print(sens.spec.I1); summary(sens.spec.I1)


# I2 test
dat.I2 <- as.table(matrix(c(nI2.true.pos, nI2.false.pos, nI2.false.neg, nI2.true.neg),
                          nrow = 2, byrow = TRUE))
colnames(dat.I2) <- c("Dis+","Dis-")
rownames(dat.I2) <- c("Test+","Test-")
sens.spec.I2 <- epi.tests(dat.I2, conf.level = 0.95)
print(sens.spec.I2); summary(sens.spec.I2)


# H test
dat.H <- as.table(matrix(c(nH.true.pos, nH.false.pos, nH.false.neg, nH.true.neg),
                          nrow = 2, byrow = TRUE))
colnames(dat.H) <- c("Dis+","Dis-")
rownames(dat.H) <- c("Test+","Test-")
sens.spec.H <- epi.tests(dat.H, conf.level = 0.95)
print(sens.spec.H); summary(sens.spec.H)


# H test
dat.W <- as.table(matrix(c(nW.true.pos, nW.false.pos, nW.false.neg, nW.true.neg),
                         nrow = 2, byrow = TRUE))
colnames(dat.W) <- c("Dis+","Dis-")
rownames(dat.W) <- c("Test+","Test-")
sens.spec.W <- epi.tests(dat.W, conf.level = 0.95)
print(sens.spec.W); summary(sens.spec.W)
