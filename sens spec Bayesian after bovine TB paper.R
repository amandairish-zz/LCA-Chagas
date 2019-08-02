# Attempt to calculate sensitivity and specificity using Bayesian analysis

# 4 tests, assume 1 population for now (not sure if this is valid?)
# Also assuming conditional independence of tests, which also may not be valid!

# test 1 === O
# test 2 === H
# test 3 === W
# test 4 === I

# OHWI Frequencies
# NNNN === 278
# NNNP === 30
# NNPN === 6
# NNPP === 9
# NPNN === 1
# NPNP === 4
# NPPP === 10
# PNNN === 1
# PNNP === 7 
# PNPN === 1
# PNPP === 28
# PPNN === 1
# PPNP === 6
# PPPN === 2
# PPPP === 416

library(runjags)
library(lattice)

pop2 <- matrix(c(416, 2, 6, 1, 28, 1, 7, 1, 10, 4, 1, 9, 6, 30, 278), 1, 15)

ns2 <- apply(pop2, 1, sum)

dd2 <- list(n = ns2, pop = pop2)

modelData2 <- dump.format(dd2)
