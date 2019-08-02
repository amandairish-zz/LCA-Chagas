library(BayesLCA)
library(tidyverse)

set.seed(123)

tau <- c(0.4, 0.6)
theta <- rbind(rep(c(0.8, 0.2), each = 2), rep(c(0.2, 0.8), each = 2))
X <- rlca(500, itemprob = theta, classprob = tau)
x <- data.blca(X)

data("Alzheimer")
alz <- data.blca(Alzheimer)

fit1 <- blca(X, 2, method = "em")
fit1

fit2 <- blca.em(x, 2, alpha = theta * 5, beta = (1 - theta) * 5,
                delta = tau * 5)
fit2

summary(fit2)

plot(fit1, which = 1:2)


####
# Gibbs sampling

sj3.gibbs <- blca.gibbs(alz, 3, relabel = FALSE)

par(mfrow = c(4, 2))
plot(sj3.gibbs, which = 5)


sj30 <- blca.gibbs(alz, 3, relabel = TRUE)
plot(sj30, which = 5)

raftery.diag(as.mcmc(sj30))

# try with Chagas data

data <- read.csv("/Users/amandairish/Desktop/LCA/dxcomp_041019_csv.csv")

# get rid of extraneous columns - for this, just use INB1_CUT0
data.lca <- data %>%
  select(O_result, INB1_CUT0, H_7POS, W_7POS)

gibbs <- blca.gibbs(data.lca, 2, thin = 1/10, iter = 50000, burn.in = 1000, relabel = TRUE)
plot(gibbs, which = 3)
plot(gibbs, which = 1:2)

raftery.diag(as.mcmc(gibbs))


