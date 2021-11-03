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

Q.Tibet <- read.csv("data/EA-2 data qtp.csv")

#subset tje data
Mac <- Q.Tibet[which(Q.Tibet$Type == "macrophyte"),]
Ter <- Q.Tibet[which(Q.Tibet$Type == "land"),]
Alg <- Q.Tibet[which(Q.Tibet$Type == "algae"),]

###compile prior parameters######

##d13C
###means and vcovs
QTP.Ter.d13C.mean.est <- c(mean(Ter$d.n.C27, na.rm = T),
                         mean(Ter$d.n.C29, na.rm = T), mean(Ter$d.n.C31, na.rm = T))
QTP.Ter.d13C.vcov <- var(data.frame(Ter$d.n.C27, Ter$d.n.C29, Ter$d.n.C31), use = "complete.obs")

QTP.Mac.d13C.mean.est <- c(mean(Mac$d.n.C27, na.rm = T),
                         mean(Mac$d.n.C29, na.rm = T), mean(Mac$d.n.C31, na.rm = T))
QTP.Mac.d13C.vcov <- var(data.frame(Mac$d.n.C27, Mac$d.n.C29, Mac$d.n.C31),use = "complete.obs")

QTP.Alg.d13C.mean.est <- c(mean(Alg$d.n.C27, na.rm = T),
                         mean(Alg$d.n.C29, na.rm = T), mean(Alg$d.n.C31, na.rm = T))
QTP.Alg.d13C.vcov <- var(data.frame(Alg$d.n.C27, Alg$d.n.C29, Alg$d.n.C31),use="complete.obs")

###rows are sources i
###columns are chain n
QTP.d13C.mu <- rbind(QTP.Ter.d13C.mean.est, QTP.Mac.d13C.mean.est, QTP.Alg.d13C.mean.est)

QTP.d13C.vcov <- rbind(QTP.Ter.d13C.vcov, QTP.Mac.d13C.vcov, QTP.Alg.d13C.vcov)

#the compiled prior parameters
QTP.d13C.mu
QTP.d13C.vcov

###concentrations are natrual log transformed
###means and vcovs
QTP.Ter.prod.mean.est <- c(mean(log(Ter$c.n.C27), na.rm = T),
                         mean(log(Ter$c.n.C29), na.rm = T), mean(log(Ter$c.n.C31), na.rm = T))
QTP.Ter.prod.vcov <- var(data.frame(log(Ter$c.n.C27), log(Ter$c.n.C29),
                                  log(Ter$c.n.C31)), use = "complete.obs")

QTP.Mac.prod.mean.est <- c(mean(log(Mac$c.n.C27), na.rm = T),
                         mean(log(Mac$c.n.C29),na.rm = T), mean(log(Mac$c.n.C31), na.rm = T))
QTP.Mac.prod.vcov <- var(data.frame(log(Mac$c.n.C27), log(Mac$c.n.C29),
                                  log(Mac$c.n.C31)), use = "complete.obs")

QTP.Alg.prod.mean.est <- c(mean(log(Alg$c.n.C27), na.rm = T),
                         mean(log(Alg$c.n.C29), na.rm = T), mean(log(Alg$c.n.C31), na.rm = T))
QTP.Alg.prod.vcov <- var(data.frame(log(Alg$c.n.C27), log(Alg$c.n.C29),
                                  log(Alg$c.n.C31)), use = "complete.obs")

###rows are sources i
###columns are chain n
QTP.prod.mu <- rbind(QTP.Ter.prod.mean.est, QTP.Mac.prod.mean.est, QTP.Alg.prod.mean.est)

QTP.prod.vcov <- rbind(QTP.Ter.prod.vcov, QTP.Mac.prod.vcov, QTP.Alg.prod.vcov)

#the compiled prior parameters
QTP.prod.mu
QTP.prod.vcov

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

####first data point: Liu et al. 2015 QHS13-05S######

#relative abundance among the tree chains from raw concentration values
#in the order of n-C27, n-C29, n-C31
RA.QHS13_5S <- c(570, 930, 1125) / (570 + 930 + 1125)

#d13C
d13C.QHS13_5S <- c(-28.5, -31.5, -32.1) 

#analytical precision
d13C.sd.Liu <- c(0.3, 0.3, 0.3) 

##Data to pass to the model
#prior parameters in the first two lines
#model parameters in the third
#data in the fourth
dat = list(d13C.mu.est = QTP.d13C.mu, d13C.sigma.est = QTP.d13C.vcov,
           prod.mu.est = QTP.prod.mu, prod.sigma.est = QTP.prod.vcov, 
           I = I, N = N, K = K, 
           RA.mix = RA.QHS13_5S, d13C.mix = d13C.QHS13_5S, d13C.mea.sd = d13C.sd.Liu)

#Start time
t1 = proc.time()

set.seed(t1[3])
#Run it
QHS13_5S.mix = do.call(jags.parallel,list(model.file = "code/LipidMM-JAGS-Multinorm.R", 
                                          parameters.to.save = parameters, 
                                          data = dat, n.chains = 3, n.iter = n.iter, 
                                          n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1

#use rhat to check convergence
QHS13_5S.mix$BUGSoutput$summary

save(QHS13_5S.mix, file = "out/QHS13_5S_results.RData")

####Second data point: Liu et al. 2015 QHS13-07S######

#relative abundance among the tree chains from raw concentration values
#in the order of n-C27, n-C29, n-C31
RA.QHS13_7S <- c(420, 825, 840) / (420 + 825 + 840)

#d13C
d13C.QHS13_7S <- c(-34.0, -33.1, -32.7)

#analytical precision
d13C.sd.Liu <- c(0.3, 0.3, 0.3) 

##Data to pass to the model
#prior parameters in the first two lines
#model parameters in the third
#data in the fourth
dat = list(d13C.mu.est = QTP.d13C.mu, d13C.sigma.est = QTP.d13C.vcov,
           prod.mu.est = QTP.prod.mu, prod.sigma.est = QTP.prod.vcov, 
           I = I, N = N, K = K, 
           RA.mix = RA.QHS13_7S, d13C.mix = d13C.QHS13_7S, d13C.mea.sd = d13C.sd.Liu)

#Start time
t1 = proc.time()

set.seed(t1[3])
#Run it
QHS13_7S.mix = do.call(jags.parallel,list(model.file = "code/LipidMM-JAGS-Multinorm.R", 
                                          parameters.to.save = parameters, 
                                          data = dat, n.chains = 3, n.iter = n.iter, 
                                          n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1

#use rhat to check convergence
QHS13_7S.mix$BUGSoutput$summary

save(QHS13_7S.mix, file = "out/QHS13_7S_results.RData")

####third data point: Liu et al. 2015 QHS13-09S######

#relative abundance among the tree chains from raw concentration values
#in the order of n-C27, n-C29, n-C31

RA.QHS13_9S <- c(885, 1920, 2220) / (885 + 1920 + 2220)

#d13C
d13C.QHS13_9S <- c(-32.0, -32.6, -32.5) 

#analytical precision
d13C.sd.Liu <- c(0.3, 0.3, 0.3) 

##Data to pass to the model
#prior parameters in the first two lines
#model parameters in the third
#data in the fourth
dat = list(d13C.mu.est = QTP.d13C.mu, d13C.sigma.est = QTP.d13C.vcov,
           prod.mu.est = QTP.prod.mu, prod.sigma.est = QTP.prod.vcov, 
           I = I, N = N, K = K, 
           RA.mix = RA.QHS13_9S, d13C.mix = d13C.QHS13_9S, d13C.mea.sd = d13C.sd.Liu)

#Start time
t1 = proc.time()

set.seed(t1[3])
#Run it
QHS13_9S.mix = do.call(jags.parallel,list(model.file = "code/LipidMM-JAGS-Multinorm.R", 
                                          parameters.to.save = parameters, 
                                          data = dat, n.chains = 3, n.iter = n.iter, 
                                          n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1

#use rhat to check convergence
QHS13_9S.mix$BUGSoutput$summary

save(QHS13_9S.mix, file = "out/QHS13_9S_results.RData")

#### MAPs, medians, and 89% HDIs ####
load("out/QHS13_5S_results.RData")
load("out/QHS13_7S_results.RData")
load("out/QHS13_9S_results.RData")

QHS13_5S.ter.map <- map_estimate(QHS13_5S.mix$BUGSoutput$sims.list$f[,1], method = "KernSmooth")
QHS13_5S.ter.median <- median(QHS13_5S.mix$BUGSoutput$sims.list$f[,1])
QHS13_5S.ter.hdi <- hdi(QHS13_5S.mix$BUGSoutput$sims.list$f[,1], ci = .89)

QHS13_5S.mac.map <- map_estimate(QHS13_5S.mix$BUGSoutput$sims.list$f[,2], method = "KernSmooth")
QHS13_5S.mac.median <- median(QHS13_5S.mix$BUGSoutput$sims.list$f[,2])
QHS13_5S.mac.hdi <- hdi(QHS13_5S.mix$BUGSoutput$sims.list$f[,2], ci = .89)

QHS13_5S.alg.map <- map_estimate(QHS13_5S.mix$BUGSoutput$sims.list$f[,3], method = "KernSmooth")
QHS13_5S.alg.median <- median(QHS13_5S.mix$BUGSoutput$sims.list$f[,3])
QHS13_5S.alg.hdi <- hdi(QHS13_5S.mix$BUGSoutput$sims.list$f[,3], ci = .89)

QHS13_7S.ter.map <- map_estimate(QHS13_7S.mix$BUGSoutput$sims.list$f[,1], method = "KernSmooth")
QHS13_7S.ter.median <- median(QHS13_7S.mix$BUGSoutput$sims.list$f[,1])
QHS13_7S.ter.hdi <- hdi(QHS13_7S.mix$BUGSoutput$sims.list$f[,1], ci = .89)

QHS13_7S.mac.map <- map_estimate(QHS13_7S.mix$BUGSoutput$sims.list$f[,2], method = "KernSmooth")
QHS13_7S.mac.median <- median(QHS13_7S.mix$BUGSoutput$sims.list$f[,2])
QHS13_7S.mac.hdi <- hdi(QHS13_7S.mix$BUGSoutput$sims.list$f[,2], ci = .89)

QHS13_7S.alg.map <- map_estimate(QHS13_7S.mix$BUGSoutput$sims.list$f[,3], method = "KernSmooth")
QHS13_7S.alg.median <- median(QHS13_7S.mix$BUGSoutput$sims.list$f[,3])
QHS13_7S.alg.hdi <- hdi(QHS13_7S.mix$BUGSoutput$sims.list$f[,3], ci = .89)

QHS13_9S.ter.map <- map_estimate(QHS13_9S.mix$BUGSoutput$sims.list$f[,1], method = "KernSmooth")
QHS13_9S.ter.median <- median(QHS13_9S.mix$BUGSoutput$sims.list$f[,1])
QHS13_9S.ter.hdi <- hdi(QHS13_9S.mix$BUGSoutput$sims.list$f[,1], ci = .89)

QHS13_9S.mac.map <- map_estimate(QHS13_9S.mix$BUGSoutput$sims.list$f[,2], method = "KernSmooth")
QHS13_9S.mac.median <- median(QHS13_9S.mix$BUGSoutput$sims.list$f[,2])
QHS13_9S.mac.hdi <- hdi(QHS13_9S.mix$BUGSoutput$sims.list$f[,2], ci = .89)

QHS13_9S.alg.map <- map_estimate(QHS13_9S.mix$BUGSoutput$sims.list$f[,3], method = "KernSmooth")
QHS13_9S.alg.median <- median(QHS13_9S.mix$BUGSoutput$sims.list$f[,3])
QHS13_9S.alg.hdi <- hdi(QHS13_9S.mix$BUGSoutput$sims.list$f[,3], ci = .89)
