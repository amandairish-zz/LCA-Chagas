# Chagas LCA with CB/JW
# Last updated by AI 5/20/19

library(BayesLCA)
library(tidyverse)

# load & select data to use
data <- read.csv("/Users/amandairish/Desktop/LCA/dxcomp_041019_csv.csv")

# get rid of extraneous columns - for this, just use INB1_CUT0
data.lca <- data %>%
  select(O_result, INB1_CUT0, H_7POS, W_7POS)

# EM implementation

# This is the same method that SAS & poLCA use
# alpha, beta, delta are set to 1 for uninformative priors
# note that BayesLCA recommends using SD rather than SE; 
# see documentation for details (at least this is true for EM?)
set.seed(7)
lca.em <- blca.em(data.lca, 2, iter = 10000, restarts = 5, 
                  verbose = TRUE, sd = TRUE, conv = 1e-06, 
                  small = 1e-100)

lca.em


# bootstrap implementation

# blca.boot repeatedly samples from the data with replacement then utilises an 
# EM algorithm to find maximum posterior (MAP) and standard error estimates 
# of the parameters.
# alpha, beta, delta are set to 1 for uninformative priors
set.seed(13)
lca.boot <- blca.boot(data.lca, 2, iter = 10000, B = 10000, relabel = TRUE,
                      verbose = TRUE, verbose.update = 1000, small = 1e-100)
lca.boot


# See if data are normally distributed for calculation of CIs
dimnames(lca.boot$samples$itemprob[,,1:4])
# order O, I, H, W

# sensitivities
hist.O <- hist(lca.boot$samples$itemprob[,1,1])
hist.I <- hist(lca.boot$samples$itemprob[,1,2])
hist.H <- hist(lca.boot$samples$itemprob[,1,3])
hist.W <- hist(lca.boot$samples$itemprob[,1,4])

# fairly normally distributed except for I which is left-skewed

O.sens <- lca.boot$itemprob[1,1]
O.spec <- 1 - lca.boot$itemprob[2,1]

I.sens <- lca.boot$itemprob[1,2]
I.spec <- 1 - lca.boot$itemprob[2,2]

H.sens <- lca.boot$itemprob[1,3]
H.spec <- 1 - lca.boot$itemprob[2,3]

W.sens <- lca.boot$itemprob[1,4]
W.spec <- 1 - lca.boot$itemprob[2,4]

# quantile method of CIs

# O.sens
quantile((lca.boot$samples$itemprob[,1,1]),0.025)
quantile((lca.boot$samples$itemprob[,1,1]),0.975)

# O.spec
1 - as.numeric(quantile((lca.boot$samples$itemprob[,2,1]),0.025))
1 - as.numeric(quantile((lca.boot$samples$itemprob[,2,1]),0.975))

# O sens: 96.5 (94.2, 98.6)
# O spec: 98.8 (96.2, 1.0)


# +/- 1.96*SE method
# O.sens
O.sens + 1.96*lca.boot$itemprob.se[1,1]
O.sens - 1.96*lca.boot$itemprob.se[1,1]

# O.spec
O.spec + 1.96*lca.boot$itemprob.se[2,1]
O.spec - 1.96*lca.boot$itemprob.se[2,1]
