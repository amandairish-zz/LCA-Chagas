# LCA first attempt

# Using poLCA
library(poLCA)
library(tidyverse)

# load data
data <- read.csv("/Users/amandairish/Desktop/LCA/dxcomp_041019_csv.csv")

# get rid of extraneous columns - for this, just use INB1_CUT0
data.lca <- data %>%
  dplyr::select(O_result, INB1_CUT0, H_7POS, W_7POS)

# change coding from 0/1 to 1/2 for 0_result, INB1/INB2, H_7POS, and W_7POS
data.lca$O_result <- data.lca$O_result + 1
data.lca$INB1_CUT0 <- data.lca$INB1_CUT0 + 1
#data.lca$INB2_CUT0 <- data.lca$INB2_CUT0 + 1
data.lca$H_7POS <- data.lca$H_7POS + 1
data.lca$W_7POS <- data.lca$W_7POS + 1

# define function 
f <- cbind(O_result, INB1_CUT0, H_7POS, W_7POS)~1

# run lca
set.seed(1234)
lca.test <- poLCA(f, data.lca, nclass=2, maxiter = 10000, graphs = TRUE, nrep=5)

#######
# Bootstrap CIs

# first write a for-loop to repeat the LCA 1000 times to generate a sample of sens & spec

result <- data.frame(matrix(nrow = 1000, ncol = 8)) # create a df to store results of for-loop
colnames(result) <- c("Osens", "Ospec", "Isens", "Ispec", "Hsens", "Hspec", "Wsens", "Wspec")

for(i in 1:1000) {
  lca <- poLCA(f, data.lca, nclass=2, maxiter = 10000, graphs = FALSE)
  result$Osens[i] <- lca$probs$O_result[2,2]
  result$Ospec[i] <- lca$probs$O_result[1,1]
  result$Isens[i] <- lca$probs$INB1_CUT0[2,2]
  result$Ispec[i] <- lca$probs$INB1_CUT0[1,1]
  result$Hsens[i] <- lca$probs$H_7POS[2,2]
  result$Hspec[i] <- lca$probs$H_7POS[1,1]
  #result$Wsens[i] <- lca$probs$W7_POS[2,2]
  #result$Wspec[i] <- lca$probs$W7_POS[1,1]
}


#######
# now try using INB2
data.lca2 <- data %>%
  dplyr::select(O_result, INB2_CUT0, H_7POS, W_7POS)

# change coding from 0/1 to 1/2 for 0_result, INB1/INB2, H_7POS, and W_7POS
data.lca2$O_result <- data.lca2$O_result + 1
data.lca2$INB2_CUT0 <- data.lca2$INB2_CUT0 + 1
data.lca2$H_7POS <- data.lca2$H_7POS + 1
data.lca2$W_7POS <- data.lca2$W_7POS + 1

# define function 
f2 <- cbind(O_result, INB2_CUT0, H_7POS, W_7POS)~1

# run lca
lca2.2 <- poLCA(f2, data.lca2, nclass=2, graphs = TRUE)

