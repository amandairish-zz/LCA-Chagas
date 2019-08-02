# Chagas LCA with CB/JW
# Last updated by AI 6/17/19

library(BayesLCA)
library(tidyverse)

# load & select data to use
data <- read.csv("/Users/amandairish/Desktop/LCA/dxcomp_041019_csv.csv")

# get rid of extraneous columns - for this, just use INB1_CUT0
data_lca <- data %>%
  select(O_result, INB1_CUT0, H_7POS, W_7POS)

# EM implementation

# This is the same method that SAS & poLCA use
# alpha, beta, delta are set to 1 for uninformative priors
# note that BayesLCA recommends using SD rather than SE; 
# see documentation for details (at least this is true for EM?)
set.seed(7)
lca_em <- blca.em(data_lca, 2, iter = 10000, restarts = 5, 
                  verbose = TRUE, sd = TRUE, conv = 1e-06, 
                  small = 1e-100)

lca_em


# bootstrap implementation

# blca.boot repeatedly samples from the data with replacement then utilises an 
# EM algorithm to find maximum posterior (MAP) and standard error estimates 
# of the parameters.
# alpha, beta, delta are set to 1 for uninformative priors
set.seed(13)
lca_boot <- blca.boot(data_lca, 2, iter = 10000, B = 10000, relabel = TRUE,
                      verbose = TRUE, verbose.update = 1000, small = 1e-100)
lca_boot
summary(lca_boot)
plot.blca.boot(lca_boot, which = 1:4)

lca_boot$counts # number of times each unique datapoint occurred (e.g. how many 0000, 1111, etc.)

# Save sensitivity/specificity
O_sens <- as.numeric(lca_boot$itemprob[1,1])
O_spec <- 1 - as.numeric(lca_boot$itemprob[2,1])

I_sens <- as.numeric(lca_boot$itemprob[1,2])
I_spec <- 1 - as.numeric(lca_boot$itemprob[2,2])

H_sens <- as.numeric(lca_boot$itemprob[1,3])
H_spec <- 1 - as.numeric(lca_boot$itemprob[2,3])

W_sens <- as.numeric(lca_boot$itemprob[1,4])
W_spec <- 1 - as.numeric(lca_boot$itemprob[2,4])

# See if data are normally distributed for calculation of CIs
dimnames(lca_boot$samples$itemprob[,,1:4])
# order O, I, H, W

# sensitivities
hist_O_sens <- hist(lca_boot$samples$itemprob[,1,1])
hist_I_sens <- hist(lca_boot$samples$itemprob[,1,2])
hist_H.sens <- hist(lca_boot$samples$itemprob[,1,3])
hist_W_sens <- hist(lca_boot$samples$itemprob[,1,4])
# fairly normally distributed except for I which is left-skewed

# 1 - specificities
hist_O_spec <- hist(lca_boot$samples$itemprob[,2,1])
hist_I_spec <- hist(lca_boot$samples$itemprob[,2,2])
hist_H_spec <- hist(lca_boot$samples$itemprob[,2,3])
hist_W_spec <- hist(lca_boot$samples$itemprob[,2,4])
# O & H right-skewed (so specificity would be left-skewed)

# Use quantile method of CIs since some skewedness

# O_sens
O_sens_lower_CI <- as.numeric(quantile((lca_boot$samples$itemprob[,1,1]),0.025))
O_sens_upper_CI <- as.numeric(quantile((lca_boot$samples$itemprob[,1,1]),0.975))

# O_spec
O_spec_lower_CI <- 1 - as.numeric(quantile((lca_boot$samples$itemprob[,2,1]),0.975))
O_spec_upper_CI <- 1 - as.numeric(quantile((lca_boot$samples$itemprob[,2,1]),0.025))


# I_sens
I_sens_lower_CI <- as.numeric(quantile((lca_boot$samples$itemprob[,1,2]),0.025))
I_sens_upper_CI <- as.numeric(quantile((lca_boot$samples$itemprob[,1,2]),0.975))

# I_spec
I_spec_lower_CI <- 1 - as.numeric(quantile((lca_boot$samples$itemprob[,2,2]),0.975))
I_spec_upper_CI <- 1 - as.numeric(quantile((lca_boot$samples$itemprob[,2,2]),0.025))


# H_sens
H_sens_lower_CI <- as.numeric(quantile((lca_boot$samples$itemprob[,1,3]),0.025))
H_sens_upper_CI <- as.numeric(quantile((lca_boot$samples$itemprob[,1,3]),0.975))

# H_spec
H_spec_lower_CI <- 1 - as.numeric(quantile((lca_boot$samples$itemprob[,2,3]),0.975))
H_spec_upper_CI <- 1 - as.numeric(quantile((lca_boot$samples$itemprob[,2,3]),0.025))


# W_sens
W_sens_lower_CI <- as.numeric(quantile((lca_boot$samples$itemprob[,1,4]),0.025))
W_sens_upper_CI <- as.numeric(quantile((lca_boot$samples$itemprob[,1,4]),0.975))

# W_spec
W_spec_lower_CI <- 1 - as.numeric(quantile((lca_boot$samples$itemprob[,2,4]),0.975))
W_spec_upper_CI <- 1 - as.numeric(quantile((lca_boot$samples$itemprob[,2,4]),0.025))


# Make a dataframe to store and export results
Sensitivity <- c(O_sens, O_sens_lower_CI, O_sens_upper_CI, I_sens, I_sens_lower_CI, 
                 I_sens_upper_CI, H_sens, H_sens_lower_CI, H_sens_upper_CI, W_sens,
                 W_sens_lower_CI, W_sens_upper_CI)
Specificty <- c(O_spec, O_spec_lower_CI, O_spec_upper_CI, I_spec, I_spec_lower_CI,
                I_spec_upper_CI, H_spec, H_spec_lower_CI, H_spec_upper_CI, W_spec,
                W_spec_lower_CI, W_spec_upper_CI)

df <- as.data.frame(rbind(Sensitivity, Specificty))

names(df) <- c("O estimate", "O lower CI", "O upper CI", "I estimate", "I lower CI",
               "I upper CI", "H estimate", "H lower CI", "H upper CI", "W estimate",
               "W lower CI", "W upper CI")

df

write.csv(df, "/Users/amandairish/Desktop/LCA/LCA_bootstrap_with_CIs_6-4-19.csv")



# Robustness check - test for 3 latent classes
# bootstrap implementation
# alpha, beta, delta are set to 1 for uninformative priors
set.seed(16)
lca_boot_3 <- blca.boot(data_lca, 3, iter = 10000, B = 10000, relabel = TRUE,
                      verbose = TRUE, verbose.update = 1000, small = 1e-100)
lca_boot_3
summary(lca_boot_3)
plot.blca.boot(lca_boot, which = 1:5)



# +/- 1.96*SE method
# O.sens
#O.sens + 1.96*lca.boot$itemprob.se[1,1]
#O.sens - 1.96*lca.boot$itemprob.se[1,1]

# O.spec
#O.spec + 1.96*lca.boot$itemprob.se[2,1]
#O.spec - 1.96*lca.boot$itemprob.se[2,1]
