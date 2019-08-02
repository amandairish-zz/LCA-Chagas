# ================================
# = three tests and 1 population-South West Region animals
# ================================

# test 1 === Meat Insp (VI)
# test 2 === SCITT standard (fat1)
# test 3 === Bov05 (fat2)


# +++/++-/-++/-+-/+-+/+--/--+/---

pop <- matrix(c(166, 20, 96, 89, 134, 271, 1166, 13895), 1, 8)

ns <- apply(pop, 1, sum)

dd<-list(n = ns, pop = pop)

modelData <- dump.format(dd)
modelInit1 <- dump.format(list(Sefat1=0.5, Sefat2=0.6, SeVI=0.7, Spfat1=0.95, Spfat2=0.97, SpVI=0.99, pi=runif(length(ns), 0.1, 0.3)))
modelInit2 <- dump.format(list(Sefat1=0.4, Sefat2=0.5, SeVI=0.6, Spfat1=0.93, Spfat2=0.96, SpVI=0.94, pi=runif(length(ns), 0.1, 0.3)))
modelInit3 <- dump.format(list(Sefat1=0.7, Sefat2=0.7, SeVI=0.5, Spfat1=0.90, Spfat2=0.98, SpVI=0.98, pi=runif(length(ns), 0.1, 0.3)))
modelInits <- c(modelInit1, modelInit2,  modelInit3)

# ====================================
# =  define model string for runjags =
# ====================================

modelString <- "
model{
pop[1, 1:8] ~ dmulti(par[1, 1:8], n)

par[1,1] <- pi*SeVI*(Sefat1*Sefat2+covDp) + (1-pi)*(1-SpVI)*((1-Spfat1)*(1-Spfat2)+covDn)
par[1,2] <- pi*SeVI*(Sefat1*(1-Sefat2)-covDp) + (1-pi)*(1-SpVI)*((1-Spfat1)*Spfat2-covDn)
par[1,3] <- pi*(1-SeVI)*(Sefat1*Sefat2+covDp) + (1-pi)*SpVI*((1-Spfat1)*(1-Spfat2)+covDn)
par[1,4] <- pi*(1-SeVI)*(Sefat1*(1-Sefat2)-covDp) + (1-pi)*SpVI*((1-Spfat1)*Spfat2-covDn)
par[1,5] <- pi*SeVI*((1-Sefat1)*Sefat2-covDp) + (1-pi)*(1-SpVI)*(Spfat1*(1-Spfat2)-covDn)
par[1,6] <- pi*SeVI*((1-Sefat1)*(1-Sefat2)+covDp) + (1-pi)*(1-SpVI)*(Spfat1*Spfat2+covDn)
par[1,7] <- pi*(1-SeVI)*((1-Sefat1)*Sefat2-covDp) + (1-pi)*SpVI*(Spfat1*(1-Spfat2)-covDn)
par[1,8] <- pi*(1-SeVI)*((1-Sefat1)*(1-Sefat2)+covDp) + (1-pi)*SpVI*(Spfat1*Spfat2+covDn)

ls <- (Sefat1-1)*(1-Sefat2)
us <- min(Sefat1,Sefat2) - Sefat1*Sefat2
lc <- (Spfat1-1)*(1-Spfat2)
uc <- min(Spfat1,Spfat2) - Spfat1*Spfat2
rhoD <- covDp / sqrt(Sefat1*(1-Sefat1)*Sefat2*(1-Sefat2))
rhoDc <- covDn / sqrt(Spfat1*(1-Spfat1)*Spfat2*(1-Spfat2))

pi ~ dunif(0, 0.3) 
Sefat1 ~ dbeta(1,1) 
Spfat1 ~ dbeta(1545, 17) 
Sefat2 ~ dbeta(1, 1) 
Spfat2 ~ dbeta(1496, 66) 
SeVI ~ dbeta(1, 1) 
SpVI ~ dbeta(1, 1) 
# covDn ~ dunif(lc, uc)
# covDp ~ dunif(ls, us)
covDn <- 0
covDp <- 0
}
"

print(system.time(
  x <- run.jags(	model = modelString,
                 monitor = c("Sefat1", "Sefat2", "SeVI", "Spfat1", "Spfat2",  "SpVI", "rhoD", "rhoDc", "pi"),
                 data = modelData,
                 burnin = 20000,
                 sample = 30000,
                 thin = 10, 
                 n.chains=3, 
                 inits = modelInits
  )
))


x <- list(	result = as.mcmc.list(x[[1]]),
           model = modelString,
           data = modelData)


print(summary(x[[1]]))

# store results
# dput(x, file=paste("result",gsub(":","",date())))				


# # =====================================
# # = some exploratory plots and output =
# # =====================================
pstring<-rev(c("Sefat1", "Sefat2", "SeVI","Spfat1", "Spfat2",  "SpVI", "rhoD", "rhoDc", "pi"))
print(densityplot(	x[[1]][,pstring],
                   bw=0.005,
                   aspect=1,
                   layout=c(3,3)))

print(summary(x[[1]]))
quartz()
par(mfrow=c(9,1))
par(mai=c(0,0,0,0))
traceplot(x[[1]])

# crossplot
quartz()
par(las=2)  
crosscorr.plot(x[[1]],las=2)

gelman.diag(x[[1]])
gelman.plot(x[[1]])
