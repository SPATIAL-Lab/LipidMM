library(coda)
library(lattice)
library(rjags)
library(R2jags)
library(msm)
library(ggmcmc)
library(mcmcplots)
library(MASS)
library(viridisLite)
library(EnvStats)
library(bayestestR)

setwd("C:/Users/ydmag/Google Drive/U of U/Proxy project/LipidMM")

####make sure to run all the helper functions in "code/LipidMM-helper functions.R"
African <- read.csv("data/EA-3 data afr.csv")

#subset the data
GR <- African[which(African$Source == "GR"),]
SV <- African[which(African$Source == "SV"),]
RF <- African[which(African$Source == "RF"),]

###compile prior parameters######

##d13C
###means and vcovs
Afr.GR.d13C.mean.est <- c(mean(GR$d.n.C29, na.rm = T),
                        mean(GR$d.n.C31, na.rm = T), mean(GR$d.n.C33, na.rm = T))
Afr.GR.d13C.vcov <- var(data.frame(GR$d.n.C29, GR$d.n.C31, GR$d.n.C33), use = "complete.obs")

Afr.SV.d13C.mean.est <- c(mean(SV$d.n.C29, na.rm = T),
                        mean(SV$d.n.C31, na.rm = T), mean(SV$d.n.C33, na.rm = T))
Afr.SV.d13C.vcov <- var(data.frame(SV$d.n.C29, SV$d.n.C31, SV$d.n.C33), use = "complete.obs")

Afr.RF.d13C.mean.est <- c(mean(RF$d.n.C29, na.rm = T),
                        mean(RF$d.n.C31, na.rm = T), mean(RF$d.n.C33, na.rm = T))
Afr.RF.d13C.vcov <- var(data.frame(RF$d.n.C29, RF$d.n.C31, RF$d.n.C33), use = "complete.obs")

###rows are sources i
###columns are chain n
Afr.d13C.mu <- rbind(Afr.GR.d13C.mean.est, Afr.SV.d13C.mean.est, Afr.RF.d13C.mean.est)

Afr.d13C.vcov <- rbind(Afr.GR.d13C.vcov, Afr.SV.d13C.vcov, Afr.RF.d13C.vcov)

#the compiled prior parameters
Afr.d13C.mu
Afr.d13C.vcov

###concentrations are natrual log transformed
###means and vcovs
Afr.GR.prod.mean.est <- c(mean(log(GR$c.n.C29), na.rm = T),
                        mean(log(GR$c.n.C31), na.rm = T), mean(log(GR$c.n.C33), na.rm = T))
Afr.GR.prod.vcov <- var(data.frame(log(GR$c.n.C29), log(GR$c.n.C31), log(GR$c.n.C33)), use = "complete.obs")

Afr.SV.prod.mean.est <- c(mean(log(SV$c.n.C29), na.rm = T),
                        mean(log(SV$c.n.C31), na.rm = T), mean(log(SV$c.n.C33), na.rm = T))
Afr.SV.prod.vcov <- var(data.frame(log(SV$c.n.C29), log(SV$c.n.C31), log(SV$c.n.C33)), use = "complete.obs")

Afr.RF.prod.mean.est <- c(mean(log(RF$c.n.C29), na.rm = T),
                        mean(log(RF$c.n.C31), na.rm = T), mean(log(RF$c.n.C33), na.rm = T))
Afr.RF.prod.vcov <- var(data.frame(log(RF$c.n.C29), log(RF$c.n.C31), log(RF$c.n.C33)), use = "complete.obs")

###rows are sources i
###columns are chain n
Afr.prod.mu <- rbind(Afr.GR.prod.mean.est, Afr.SV.prod.mean.est, Afr.RF.prod.mean.est)

Afr.prod.vcov <- rbind(Afr.GR.prod.vcov, Afr.SV.prod.vcov, Afr.RF.prod.vcov)

#the compiled prior parameters
Afr.prod.mu
Afr.prod.vcov

#######model parameters#####
#initialize parameters
I <- 3  #number of sources
N <- 3  #number of chains
K <- 50 #number of grams of leaves to integrate per source

#common MCMC parameters
n.iter = 8e5
n.burnin = 2e5
n.thin = floor(n.iter-n.burnin)/2500
#average runtime is ~2 hours/sample

#model parameters to save
parameters <- c("d13C.mix.m","RA.mix.m", "f","f.sum.prod_n_i","exp.prod_k","d13C.k")

####first data point: Garcin et al. 2014 Rhum lake

#relative abundance among the tree chains from raw concentration values
#in the order of n-C27, n-C29, n-C31
RA.rhum.l <- c(0.317, 0.384, 0.299)

#d13C
d13C.rhum.l <- c(-28.6, -29.7, -24.8)

#analytical precision
d13C.sd.rhum.l <- c(0.3, 0.3, 0.3)

##Data to pass to the model
#prior parameters in the first two lines
#model parameters in the third
#data in the fourth
dat = list(d13C.mu.est = Afr.d13C.mu, d13C.sigma.est = Afr.d13C.vcov,
           prod.mu.est = Afr.prod.mu, prod.sigma.est = Afr.prod.vcov,  
           I = I, N = N, K = K, 
           RA.mix = RA.rhum.l, d13C.mix = d13C.rhum.l, d13C.mea.sd = d13C.sd.rhum.l)

#Start time
t1 = proc.time()

set.seed(t1[3])
#Run it
rhum.l.mix = do.call(jags.parallel,list(model.file = "code/LipidMM-JAGS-Multinorm.R", 
                                          parameters.to.save = parameters, 
                                          data = dat, n.chains = 3, n.iter = n.iter, 
                                          n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1

#use rhat to check convergence
rhum.l.mix$BUGSoutput$summary

save(rhum.l.mix, file = "out/rhum_l_results.RData")

####second data point: Garcin et al. 2014 Asso lake

#relative abundance among the tree chains from raw concentration values
#in the order of n-C27, n-C29, n-C31
RA.asso.l <- c(0.502, 0.289, 0.209)

#d13C
d13C.asso.l <- c(-32.9, -31.1, -23.8)

#analytical precision
d13C.sd.asso.l <- c(0.2, 0.8, 0.8)

##Data to pass to the model
#prior parameters in the first two lines
#model parameters in the third
#data in the fourth
dat = list(d13C.mu.est = Afr.d13C.mu, d13C.sigma.est = Afr.d13C.vcov,
           prod.mu.est = Afr.prod.mu, prod.sigma.est = Afr.prod.vcov,  
           I = I, N = N, K = K, 
           RA.mix = RA.asso.l, d13C.mix = d13C.asso.l, d13C.mea.sd = d13C.sd.asso.l)

#Start time
t1 = proc.time()

set.seed(t1[3])
#Run it
asso.l.mix = do.call(jags.parallel,list(model.file = "code/LipidMM-JAGS-Multinorm.R", 
                                        parameters.to.save = parameters, 
                                        data = dat, n.chains = 3, n.iter = n.iter, 
                                        n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1

#use rhat to check convergence
asso.l.mix$BUGSoutput$summary

save(asso.l.mix, file = "out/asso_l_results.RData")

####third data point: Garcin et al. 2014 Baro lake

#relative abundance among the tree chains from raw concentration values
#in the order of n-C27, n-C29, n-C31
RA.baro.l <- c(0.317,0.462,0.167)

#d13C
d13C.baro.l <- c(-35.3, -36.2, -32.7)

#analytical precision
#original data for n-c33 is 0, but this is not realistic
#so it is set at 0.04 as it rounds down to 0
d13C.sd.baro.l <- c(0.7, 0.4, 0.04)

##Data to pass to the model
#prior parameters in the first two lines
#model parameters in the third
#data in the fourth
dat = list(d13C.mu.est = Afr.d13C.mu, d13C.sigma.est = Afr.d13C.vcov,
           prod.mu.est = Afr.prod.mu, prod.sigma.est = Afr.prod.vcov,  
           I = I, N = N, K = K, 
           RA.mix = RA.baro.l, d13C.mix = d13C.baro.l, d13C.mea.sd = d13C.sd.baro.l)

#Start time
t1 = proc.time()

set.seed(t1[3])
#Run it
baro.l.mix = do.call(jags.parallel,list(model.file = "code/LipidMM-JAGS-Multinorm.R", 
                                        parameters.to.save = parameters, 
                                        data = dat, n.chains = 3, n.iter = n.iter, 
                                        n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1

#use rhat to check convergence
baro.l.mix$BUGSoutput$summary

save(baro.l.mix, file = "out/baro_l_results.RData")

#### MAPs, medians, and 89% HDIs ####

rhum.l.GR.map <- map_estimate(rhum.l.mix$BUGSoutput$sims.list$f[,1], method = "KernSmooth")
rhum.l.GR.median <- median(rhum.l.mix$BUGSoutput$sims.list$f[,1])
rhum.l.GR.hdi <- hdi(rhum.l.mix$BUGSoutput$sims.list$f[,1], ci = .89)

rhum.l.SV.map <- map_estimate(rhum.l.mix$BUGSoutput$sims.list$f[,2], method = "KernSmooth")
rhum.l.SV.median <- median(rhum.l.mix$BUGSoutput$sims.list$f[,2])
rhum.l.SV.hdi <- hdi(rhum.l.mix$BUGSoutput$sims.list$f[,2], ci = .89)

rhum.l.RF.map <- map_estimate(rhum.l.mix$BUGSoutput$sims.list$f[,3], method = "KernSmooth")
rhum.l.RF.median <- median(rhum.l.mix$BUGSoutput$sims.list$f[,3])
rhum.l.RF.hdi <- hdi(rhum.l.mix$BUGSoutput$sims.list$f[,3], ci = .89)

asso.l.GR.map <- map_estimate(asso.l.mix$BUGSoutput$sims.list$f[,1], method = "KernSmooth")
asso.l.GR.median <- median(asso.l.mix$BUGSoutput$sims.list$f[,1])
asso.l.GR.hdi <- hdi(asso.l.mix$BUGSoutput$sims.list$f[,1], ci = .89)

asso.l.SV.map <- map_estimate(asso.l.mix$BUGSoutput$sims.list$f[,2], method = "KernSmooth")
asso.l.SV.median <- median(asso.l.mix$BUGSoutput$sims.list$f[,2])
asso.l.SV.hdi <- hdi(asso.l.mix$BUGSoutput$sims.list$f[,2], ci = .89)

asso.l.RF.map <- map_estimate(asso.l.mix$BUGSoutput$sims.list$f[,3], method = "KernSmooth")
asso.l.RF.media <- median(asso.l.mix$BUGSoutput$sims.list$f[,3])
asso.l.RF.hdi <- hdi(asso.l.mix$BUGSoutput$sims.list$f[,3], ci = .89)

baro.l.GR.map <- map_estimate(baro.l.mix$BUGSoutput$sims.list$f[,1], method = "KernSmooth")
baro.l.GR.median <- median(baro.l.mix$BUGSoutput$sims.list$f[,1])
baro.l.GR.hdi <- hdi(baro.l.mix$BUGSoutput$sims.list$f[,1], ci = .89)

baro.l.SV.map <- map_estimate(baro.l.mix$BUGSoutput$sims.list$f[,2], method = "KernSmooth")
baro.l.SV.median <- median(baro.l.mix$BUGSoutput$sims.list$f[,2])
baro.l.SV.hdi <- hdi(baro.l.mix$BUGSoutput$sims.list$f[,2], ci = .89)

baro.l.RF.map <- map_estimate(baro.l.mix$BUGSoutput$sims.list$f[,3], method = "KernSmooth")
baro.l.RF.median <- median(baro.l.mix$BUGSoutput$sims.list$f[,3])
baro.l.RF.hdi <- hdi(baro.l.mix$BUGSoutput$sims.list$f[,3], ci = .89)

#####sensitivity tests#####
####test 1: using different priors####
#using a subset of CS2: western africa
W.African <- read.csv("data/EA-4 data w-afr.csv")

#subset data
W.GR <- W.African[which(W.African$Source == "GR"),]
W.SV <- W.African[which(W.African$Source == "SV"),]
W.RF <- W.African[which(W.African$Source == "RF"),]

###compile parameters
#concentration
W.GR.prod.mean.est <- c(mean(log(W.GR$c.n.C29), na.rm = T),
                      mean(log(W.GR$c.n.C31), na.rm = T), mean(log(W.GR$c.n.C33), na.rm = T))
W.GR.prod.vcov <- var(data.frame(log(W.GR$c.n.C29), log(W.GR$c.n.C31),
                               log(W.GR$c.n.C33)), use = "complete.obs")

W.SV.prod.mean.est <- c(mean(log(W.SV$c.n.C29),na.rm = T),
                      mean(log(W.SV$c.n.C31),na.rm = T), mean(log(W.SV$c.n.C33), na.rm = T))
W.SV.prod.vcov <- var(data.frame(log(W.SV$c.n.C29), log(W.SV$c.n.C31),
                               log(W.SV$c.n.C33)), use = "complete.obs")

W.RF.prod.mean.est <- c(mean(log(W.RF$c.n.C29), na.rm = T),
                      mean(log(W.RF$c.n.C31), na.rm = T), mean(log(W.RF$c.n.C33), na.rm = T))
W.RF.prod.vcov <- var(data.frame(log(W.RF$c.n.C29), log(W.RF$c.n.C31),
                               log(W.RF$c.n.C33)), use = "complete.obs")

W.prod.vcov <- rbind(W.GR.prod.vcov, W.SV.prod.vcov, W.RF.prod.vcov)

W.prod.mu <- rbind(W.GR.prod.mean.est, W.SV.prod.mean.est, W.RF.prod.mean.est)

#d13C
W.GR.d13C.mean.est <- c(mean(W.GR$d.n.C29, na.rm = T),
                      mean(W.GR$d.n.C31, na.rm = T), mean(W.GR$d.n.C33, na.rm = T))
W.GR.d13C.vcov <- var(data.frame(W.GR$d.n.C29, W.GR$d.n.C31, W.GR$d.n.C33), use = "complete.obs")

W.SV.d13C.mean.est <- c(mean(W.SV$d.n.C29, na.rm = T),
                      mean(W.SV$d.n.C31, na.rm = T), mean(W.SV$d.n.C33, na.rm = T))
W.SV.d13C.vcov <- var(data.frame(W.SV$d.n.C29, W.SV$d.n.C31, W.SV$d.n.C33), use = "complete.obs")

W.RF.d13C.mean.est <- c(mean(W.RF$d.n.C29, na.rm = T),
                      mean(W.RF$d.n.C31, na.rm = T), mean(W.RF$d.n.C33, na.rm = T))
W.RF.d13C.vcov <- var(data.frame(W.RF$d.n.C29, W.RF$d.n.C31, W.RF$d.n.C33), use = "complete.obs")

W.d13C.vcov <- rbind(W.GR.d13C.vcov, W.SV.d13C.vcov, W.RF.d13C.vcov)

W.d13C.mu <- rbind(W.GR.d13C.mean.est, W.SV.d13C.mean.est, W.RF.d13C.mean.est)

##data point 1: rhum
##Data to pass to the model
#prior parameters in the first two lines
#model parameters in the third
#data in the fourth
#only the priors have been changed to W.xxx, other parameters are the same
dat = list(d13C.mu.est = W.d13C.mu, d13C.sigma.est = W.d13C.vcov,
           prod.mu.est = W.prod.mu, prod.sigma.est = W.prod.vcov,  
           I = I, N = N, K = K, 
           RA.mix = RA.rhum.l, d13C.mix = d13C.rhum.l, d13C.mea.sd = d13C.sd.rhum.l)

#Start time
t1 = proc.time()

set.seed(t1[3])
#Run it
W.rhum.l.mix = do.call(jags.parallel,list(model.file = "code/LipidMM-JAGS-Multinorm.R", 
                                        parameters.to.save = parameters, 
                                        data = dat, n.chains = 3, n.iter = n.iter, 
                                        n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1

#use rhat to check convergence
W.rhum.l.mix$BUGSoutput$summary

save(W.rhum.l.mix, file = "out/W_rhum_l_results.RData")

##data point 2: asso
##Data to pass to the model
#prior parameters in the first two lines
#model parameters in the third
#data in the fourth
#only the priors have been changed to W.xxx, other parameters are the same
dat = list(d13C.mu.est = W.d13C.mu, d13C.sigma.est = W.d13C.vcov,
           prod.mu.est = W.prod.mu, prod.sigma.est = W.prod.vcov,  
           I = I, N = N, K = K, 
           RA.mix = RA.asso.l, d13C.mix = d13C.asso.l, d13C.mea.sd = d13C.sd.asso.l)

#Start time
t1 = proc.time()

set.seed(t1[3])
#Run it
W.asso.l.mix = do.call(jags.parallel,list(model.file = "code/LipidMM-JAGS-Multinorm.R", 
                                          parameters.to.save = parameters, 
                                          data = dat, n.chains = 3, n.iter = n.iter, 
                                          n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1

#use rhat to check convergence
W.asso.l.mix$BUGSoutput$summary

save(W.asso.l.mix, file = "out/W_asso_l_results.RData")

##data point 3: baro
##Data to pass to the model
#prior parameters in the first two lines
#model parameters in the third
#data in the fourth
#only the priors have been changed to W.xxx, other parameters are the same
dat = list(d13C.mu.est = W.d13C.mu, d13C.sigma.est = W.d13C.vcov,
           prod.mu.est = W.prod.mu, prod.sigma.est = W.prod.vcov,  
           I = I, N = N, K = K, 
           RA.mix = RA.baro.l, d13C.mix = d13C.baro.l, d13C.mea.sd = d13C.sd.baro.l)

#Start time
t1 = proc.time()

set.seed(t1[3])
#Run it
W.baro.l.mix = do.call(jags.parallel,list(model.file = "code/LipidMM-JAGS-Multinorm.R", 
                                          parameters.to.save = parameters, 
                                          data = dat, n.chains = 3, n.iter = n.iter, 
                                          n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1

#use rhat to check convergence
W.baro.l.mix$BUGSoutput$summary

save(W.baro.l.mix, file = "out/W_baro_l_results.RData")

####test 2: sensitivity to likelihood functions####
### a) RA likelihood functions completely ignored
##data point 1: rhum
##Data to pass to the model
dat = list(d13C.mu.est = Afr.d13C.mu, d13C.sigma.est = Afr.d13C.vcov,
           prod.mu.est = Afr.prod.mu, prod.sigma.est = Afr.prod.vcov, 
           I = I, N = N, K = K, 
           RA.mix = RA.rhum.l, d13C.mix = d13C.rhum.l, d13C.mea.sd = d13C.sd.rhum.l)

#Start time
t1 = proc.time()

set.seed(t1[3])
#Run it with the model that has the RA likelihood function masked
rhum.l.test.a = do.call(jags.parallel,list(model.file = "code/LipidMM-JAGS-Multinorm no RA.R", 
                                          parameters.to.save = parameters, 
                                          data = dat, n.chains = 3, n.iter = n.iter, 
                                          n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1

#use rhat to check convergence
rhum.l.test.a$BUGSoutput$summary

save(rhum.l.test.a, file = "out/rhum_l_test_a.RData")

##data point 2: asso
##Data to pass to the model
dat = list(d13C.mu.est = Afr.d13C.mu, d13C.sigma.est = Afr.d13C.vcov,
           prod.mu.est = Afr.prod.mu, prod.sigma.est = Afr.prod.vcov, 
           I = I, N = N, K = K, 
           RA.mix = RA.asso.l, d13C.mix = d13C.asso.l, d13C.mea.sd = d13C.sd.asso.l)

#Start time
t1 = proc.time()

set.seed(t1[3])
#Run it with the model that has the RA likelihood function masked
asso.l.test.a = do.call(jags.parallel,list(model.file = "code/LipidMM-JAGS-Multinorm no RA.R", 
                                          parameters.to.save = parameters, 
                                          data = dat, n.chains = 3, n.iter = n.iter, 
                                          n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1

#use rhat to check convergence
asso.l.test.a$BUGSoutput$summary

save(asso.l.test.a, file = "out/asso_l_test_a.RData")

##data point 3: baro
##Data to pass to the model
dat = list(d13C.mu.est = Afr.d13C.mu, d13C.sigma.est = Afr.d13C.vcov,
           prod.mu.est = Afr.prod.mu, prod.sigma.est = Afr.prod.vcov,  
           I = I, N = N, K = K, 
           RA.mix = RA.baro.l, d13C.mix = d13C.baro.l, d13C.mea.sd = d13C.sd.baro.l)

#Start time
t1 = proc.time()

set.seed(t1[3])
#Run it with the model that has the RA likelihood function masked
baro.l.test.a = do.call(jags.parallel,list(model.file = "code/LipidMM-JAGS-Multinorm no RA.R", 
                                          parameters.to.save = parameters, 
                                          data = dat, n.chains = 3, n.iter = n.iter, 
                                          n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1

#use rhat to check convergence
baro.l.test.a$BUGSoutput$summary

save(baro.l.test.a, file = "out/baro_l_test_a.RData")

### b) RA likelihood functions with inflated uncertainty
##data point 1: rhum
##Data to pass to the model
dat = list(d13C.mu.est = Afr.d13C.mu, d13C.sigma.est = Afr.d13C.vcov,
           prod.mu.est = Afr.prod.mu, prod.sigma.est = Afr.prod.vcov, 
           I = I, N = N, K = K, 
           RA.mix = RA.rhum.l, d13C.mix = d13C.rhum.l, d13C.mea.sd = d13C.sd.rhum.l)

#Start time
t1 = proc.time()

set.seed(t1[3])
#Run it with the model that has inflated uncertainty for the RA likelihood function
rhum.l.test.b = do.call(jags.parallel,list(model.file = "code/LipidMM-JAGS-Multinorm inflated RA.R", 
                                           parameters.to.save = parameters, 
                                           data = dat, n.chains = 3, n.iter = n.iter, 
                                           n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1

#use rhat to check convergence
rhum.l.test.b$BUGSoutput$summary

save(rhum.l.test.b, file = "out/rhum_l_test_b.RData")

##data point 2: asso
##Data to pass to the model
dat = list(d13C.mu.est = Afr.d13C.mu, d13C.sigma.est = Afr.d13C.vcov,
           prod.mu.est = Afr.prod.mu, prod.sigma.est = Afr.prod.vcov, 
           I = I, N = N, K = K, 
           RA.mix = RA.asso.l, d13C.mix = d13C.asso.l, d13C.mea.sd = d13C.sd.asso.l)

#Start time
t1 = proc.time()

set.seed(t1[3])
#Run it with the model that has inflated uncertainty for the RA likelihood function
asso.l.test.b = do.call(jags.parallel,list(model.file = "code/LipidMM-JAGS-Multinorm inflated RA.R", 
                                           parameters.to.save = parameters, 
                                           data = dat, n.chains = 3, n.iter = n.iter, 
                                           n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1

#use rhat to check convergence
asso.l.test.b$BUGSoutput$summary

save(asso.l.test.b, file = "out/asso_l_test_b.RData")

##data point 3: baro
##Data to pass to the model
dat = list(d13C.mu.est = Afr.d13C.mu, d13C.sigma.est = Afr.d13C.vcov,
           prod.mu.est = Afr.prod.mu, prod.sigma.est = Afr.prod.vcov, 
           I = I, N = N, K = K, 
           RA.mix = RA.baro.l, d13C.mix = d13C.baro.l, d13C.mea.sd = d13C.sd.baro.l)

#Start time
t1 = proc.time()

set.seed(t1[3])
#Run it with the model that has inflated uncertainty for the RA likelihood function
baro.l.test.b = do.call(jags.parallel,list(model.file = "code/LipidMM-JAGS-Multinorm inflated RA.R", 
                                           parameters.to.save = parameters, 
                                           data = dat, n.chains = 3, n.iter = n.iter, 
                                           n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1

#use rhat to check convergence
baro.l.test.b$BUGSoutput$summary

save(baro.l.test.b, file = "out/baro_l_test_b.RData")

### c) d13C likelihood functions completely ignored

##data point 1: rhum
##Data to pass to the model
dat = list(d13C.mu.est = Afr.d13C.mu, d13C.sigma.est = Afr.d13C.vcov,
           prod.mu.est = Afr.prod.mu, prod.sigma.est = Afr.prod.vcov,
           I = I, N = N, K = K, 
           RA.mix = RA.rhum.l, d13C.mix = d13C.rhum.l, d13C.mea.sd = d13C.sd.rhum.l)

#Start time
t1 = proc.time()

set.seed(t1[3])
#Run it with the model that has the d13C likelihood function masked
rhum.l.test.c = do.call(jags.parallel,list(model.file = "code/LipidMM-JAGS-Multinorm no d13C.R", 
                                           parameters.to.save = parameters, 
                                           data = dat, n.chains = 3, n.iter = n.iter, 
                                           n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1

#use rhat to check convergence
rhum.l.test.c$BUGSoutput$summary

save(rhum.l.test.c, file = "out/rhum_l_test_c.RData")

##data point 2: asso
##Data to pass to the model
dat = list(d13C.mu.est = Afr.d13C.mu, d13C.sigma.est = Afr.d13C.vcov,
           prod.mu.est = Afr.prod.mu, prod.sigma.est = Afr.prod.vcov,
           I = I, N = N, K = K, 
           RA.mix = RA.asso.l, d13C.mix = d13C.asso.l, d13C.mea.sd = d13C.sd.asso.l)

#Start time
t1 = proc.time()

set.seed(t1[3])
#Run it with the model that has the d13C likelihood function masked
asso.l.test.c = do.call(jags.parallel,list(model.file = "code/LipidMM-JAGS-Multinorm no d13C.R", 
                                           parameters.to.save = parameters, 
                                           data = dat, n.chains = 3, n.iter = n.iter, 
                                           n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1

#use rhat to check convergence
asso.l.test.c$BUGSoutput$summary

save(asso.l.test.c, file = "out/asso_l_test_c.RData")

##data point 3: baro
##Data to pass to the model
dat = list(d13C.mu.est = Afr.d13C.mu, d13C.sigma.est = Afr.d13C.vcov,
           prod.mu.est = Afr.prod.mu, prod.sigma.est = Afr.prod.vcov, 
           I = I, N = N, K = K, 
           RA.mix = RA.baro.l, d13C.mix = d13C.baro.l, d13C.mea.sd = d13C.sd.baro.l)

#Start time
t1 = proc.time()

set.seed(t1[3])
#Run it with the model that has the d13C likelihood function masked
baro.l.test.c = do.call(jags.parallel,list(model.file = "code/LipidMM-JAGS-Multinorm no d13C.R", 
                                           parameters.to.save = parameters, 
                                           data = dat, n.chains = 3, n.iter = n.iter, 
                                           n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1

#use rhat to check convergence
baro.l.test.c$BUGSoutput$summary

save(baro.l.test.c, file = "out/baro_l_test_c.RData")

### d) d13C likelihood functions with inflated uncertainty
#all models will use this uncertainty term
d13C.sd.inflated <- c(3, 3, 3)

##data point 1: rhum
##Data to pass to the model
dat = list(d13C.mu.est = Afr.d13C.mu, d13C.sigma.est = Afr.d13C.vcov,
           prod.mu.est = Afr.prod.mu, prod.sigma.est = Afr.prod.vcov,
           I = I, N = N, K = K, 
           RA.mix = RA.rhum.l, d13C.mix = d13C.rhum.l, d13C.mea.sd = d13C.sd.inflated)

#Start time
t1 = proc.time()

set.seed(t1[3])
#Run it with the default model 
rhum.l.test.d = do.call(jags.parallel,list(model.file = "code/LipidMM-JAGS-Multinorm.R", 
                                           parameters.to.save = parameters, 
                                           data = dat, n.chains = 3, n.iter = n.iter, 
                                           n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1

#use rhat to check convergence
rhum.l.test.d$BUGSoutput$summary

save(rhum.l.test.d, file = "out/rhum_l_test_d.RData")

##data point 2: asso
##Data to pass to the model
dat = list(d13C.mu.est = Afr.d13C.mu, d13C.sigma.est = Afr.d13C.vcov,
           prod.mu.est = Afr.prod.mu, prod.sigma.est = Afr.prod.vcov,
           I = I, N = N, K = K, 
           RA.mix = RA.asso.l, d13C.mix = d13C.asso.l, d13C.mea.sd = d13C.sd.inflated)

#Start time
t1 = proc.time()

set.seed(t1[3])
#Run it with the default model 
asso.l.test.d = do.call(jags.parallel,list(model.file = "code/LipidMM-JAGS-Multinorm.R", 
                                           parameters.to.save = parameters, 
                                           data = dat, n.chains = 3, n.iter = n.iter, 
                                           n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1

#use rhat to check convergence
asso.l.test.d$BUGSoutput$summary

save(asso.l.test.d, file = "out/asso_l_test_d.RData")

##data point 3: baro
##Data to pass to the model
dat = list(d13C.mu.est = Afr.d13C.mu, d13C.sigma.est = Afr.d13C.vcov,
           prod.mu.est = Afr.prod.mu, prod.sigma.est = Afr.prod.vcov,
           I = I, N = N, K = K, 
           RA.mix = RA.baro.l, d13C.mix = d13C.baro.l, d13C.mea.sd = d13C.sd.inflated)

#Start time
t1 = proc.time()

set.seed(t1[3])
#Run it with the default model 
baro.l.test.d = do.call(jags.parallel,list(model.file = "code/LipidMM-JAGS-Multinorm.R", 
                                           parameters.to.save = parameters, 
                                           data = dat, n.chains = 3, n.iter = n.iter, 
                                           n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1

#use rhat to check convergence
baro.l.test.d$BUGSoutput$summary

save(baro.l.test.d, file = "out/baro_l_test_d.RData")