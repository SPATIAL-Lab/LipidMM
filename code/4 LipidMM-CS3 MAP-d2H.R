library(coda)
library(rjags)
library(R2jags)
library(mcmcplots)
library(MASS)
library(viridisLite)
library(EnvStats)
library(bayestestR)
library(dplyr)

plot.col<-viridis(7)

#setwd("C:/Users/ydmag/Google Drive/U of U/Proxy project/LipidMM")

source("code/1 LipidMM-helper functions.R")

map_rec<- read.csv("data/EA-5 data map_rec.csv")

map_rec.GR <- subset(map_rec, map_rec$Source == "GR")
map_rec.SV <- subset(map_rec, map_rec$Source == "SV")
map_rec.RF <- subset(map_rec, map_rec$Source == "RF")

#the three chains used in this case study are Cn-27, Cn-29, and Cn-31
###compile prior parameters######

#epsilon alkane-MAP apparent fractionation factor
###means and vcovs

map_rec.GR.e.mean.est <- c(mean(map_rec.GR$eC27.MAP_OIPC, na.rm = T),
                           mean(map_rec.GR$eC29.MAP_OIPC, na.rm = T), 
                           mean(map_rec.GR$eC31.MAP_OIPC, na.rm = T),
                           mean(map_rec.GR$eC33.MAP_OIPC, na.rm = T))
map_rec.GR.e.vcov <- var(data.frame(map_rec.GR$eC27.MAP_OIPC, 
                                    map_rec.GR$eC29.MAP_OIPC, 
                                    map_rec.GR$eC31.MAP_OIPC,
                                    map_rec.GR$eC33.MAP_OIPC), use = "complete.obs")

map_rec.SV.e.mean.est <- c(mean(map_rec.SV$eC27.MAP_OIPC, na.rm = T),
                           mean(map_rec.SV$eC29.MAP_OIPC, na.rm = T), 
                           mean(map_rec.SV$eC31.MAP_OIPC, na.rm = T),
                           mean(map_rec.SV$eC33.MAP_OIPC, na.rm = T))
map_rec.SV.e.vcov <- var(data.frame(map_rec.SV$eC27.MAP_OIPC, 
                                    map_rec.SV$eC29.MAP_OIPC, 
                                    map_rec.SV$eC31.MAP_OIPC,
                                    map_rec.SV$eC33.MAP_OIPC), use = "complete.obs")

map_rec.RF.e.mean.est <- c(mean(map_rec.RF$eC27.MAP_OIPC, na.rm = T),
                           mean(map_rec.RF$eC29.MAP_OIPC, na.rm = T), 
                           mean(map_rec.RF$eC31.MAP_OIPC, na.rm = T),
                           mean(map_rec.RF$eC33.MAP_OIPC, na.rm = T))
map_rec.RF.e.vcov <- var(data.frame(map_rec.RF$eC27.MAP_OIPC, 
                                    map_rec.RF$eC29.MAP_OIPC, 
                                    map_rec.RF$eC31.MAP_OIPC,
                                    map_rec.RF$eC33.MAP_OIPC), use = "complete.obs")

###rows are sources i
###columns are chain n
map_rec.epsilon.app.mu <- rbind(map_rec.GR.e.mean.est, 
                                map_rec.SV.e.mean.est, map_rec.RF.e.mean.est)

map_rec.epsilon.app.vcov <- rbind(map_rec.GR.e.vcov, map_rec.SV.e.vcov, map_rec.RF.e.vcov)

#the compiled prior parameters
map_rec.epsilon.app.mu
map_rec.epsilon.app.vcov

####Adopt the d13C and n-alkane concentrations from CS2
African <- read.csv("data/EA-3 data afr.csv")

#subset the data
GR <- subset(African, African$Source == "GR")
SV <- subset(African, African$Source == "SV")
RF <- subset(African, African$Source == "RF")

###compile prior parameters######

##d13C
###means and vcovs
Afr.GR.d13C.mean.est <- c(mean(GR$d.n.C27, na.rm = T), mean(GR$d.n.C29, na.rm = T),
                          mean(GR$d.n.C31, na.rm = T), mean(GR$d.n.C33, na.rm = T))
Afr.GR.d13C.vcov <- var(data.frame(GR$d.n.C27, GR$d.n.C29, 
                                   GR$d.n.C31, GR$d.n.C33), use = "complete.obs")

Afr.SV.d13C.mean.est <- c(mean(SV$d.n.C27, na.rm = T), mean(SV$d.n.C29, na.rm = T),
                          mean(SV$d.n.C31, na.rm = T), mean(SV$d.n.C33, na.rm = T))
Afr.SV.d13C.vcov <- var(data.frame(SV$d.n.C27, SV$d.n.C29, 
                                   SV$d.n.C31, SV$d.n.C33), use = "complete.obs")

Afr.RF.d13C.mean.est <- c(mean(RF$d.n.C27, na.rm = T), mean(RF$d.n.C29, na.rm = T),
                          mean(RF$d.n.C31, na.rm = T), mean(RF$d.n.C33, na.rm = T))
Afr.RF.d13C.vcov <- var(data.frame(RF$d.n.C27, RF$d.n.C29, 
                                   RF$d.n.C31, RF$d.n.C33), use = "complete.obs")

###rows are sources i
###columns are chain n
Afr.d13C.mu <- rbind(Afr.GR.d13C.mean.est, Afr.SV.d13C.mean.est, Afr.RF.d13C.mean.est)

Afr.d13C.vcov <- rbind(Afr.GR.d13C.vcov, Afr.SV.d13C.vcov, Afr.RF.d13C.vcov)

#the compiled prior parameters
Afr.d13C.mu
Afr.d13C.vcov

###concentrations are natrual log transformed
###means and vcovs
Afr.GR.conc.mean.est <- c(mean(log(GR$c.n.C27), na.rm = T), mean(log(GR$c.n.C29), na.rm = T), 
                          mean(log(GR$c.n.C31), na.rm = T), mean(log(GR$c.n.C33), na.rm = T))
Afr.GR.conc.vcov <- var(data.frame(log(GR$c.n.C27), log(GR$c.n.C29), 
                                   log(GR$c.n.C31), log(GR$c.n.C33)), use = "complete.obs")

Afr.SV.conc.mean.est <- c(mean(log(SV$c.n.C27), na.rm = T), mean(log(SV$c.n.C29), na.rm = T),
                          mean(log(SV$c.n.C31), na.rm = T), mean(log(SV$c.n.C33), na.rm = T))
Afr.SV.conc.vcov <- var(data.frame(log(SV$c.n.C27), log(SV$c.n.C29), 
                                   log(SV$c.n.C31), log(SV$c.n.C33)), use = "complete.obs")

Afr.RF.conc.mean.est <- c(mean(log(RF$c.n.C27), na.rm = T), mean(log(RF$c.n.C29), na.rm = T), 
                          mean(log(RF$c.n.C31), na.rm = T), mean(log(RF$c.n.C33), na.rm = T))
Afr.RF.conc.vcov <- var(data.frame(log(RF$c.n.C27), log(RF$c.n.C29), 
                                   log(RF$c.n.C31), log(RF$c.n.C33)), use = "complete.obs")

###rows are sources i
###columns are chain n
Afr.conc.mu <- rbind(Afr.GR.conc.mean.est, Afr.SV.conc.mean.est, Afr.RF.conc.mean.est)

Afr.conc.vcov <- rbind(Afr.GR.conc.vcov, Afr.SV.conc.vcov, Afr.RF.conc.vcov)

#the compiled prior parameters
Afr.conc.mu
Afr.conc.vcov

#######model parameters#####
#initialize parameters
I <- 3  #number of sources
N <- 4  #number of chains
K <- 50 #number of grams of leaves to integrate per source

#common MCMC parameters, 5 parallel chains are used
n.iter = 1e6
n.burnin = 2e5
n.thin = floor(n.iter-n.burnin)/1500
#average runtime is ~8 hours/sample


######1st data point: 33.5 Ka#####
#relative abundance among the three chains from raw concentration values
#in the order of n-C27, n-C29, n-C31, c-C33

RA.33.5ka <- c(226, 429.6, 539.6, 452.9)/(226 + 429.6 + 539.6 + 452.9)

#d2H
d2H.33.5ka <- c(-145.3, -143.4, -154.6, -161)
#use ice volume correcte values
d2H.33.5ka.ivc <- c(-150.5, -148.7, -159.8, -166.2)

#analytical precision
d2H.sd.33.5ka <- c(0.9, 0.4, 0.3, 0.5)

#d13C
d13C.33.5ka <- c(-26.5, -28, -26.2, -23.6)

#analytical precision
d13C.sd.33.5ka <- c(0.1, 0.1, 0.1, 0.1)

###correction for atm CO2 d13C, should set it as a parameter?
#-8.3 year 2010 (Graven et al 2017)
#no data on 33.5 Ka, so -6.4 is used here as an approximate value
#according to Schmitt et al 2012, the values don't vary too much, at +- 0.2 per mil
#so this should have a very minor effect on the reconstructed vegetation
Afr.d13C.mu.cor <- Afr.d13C.mu + 8.3 - 6.4 

##Data to pass to the model
#prior parameters in the first three lines
#model parameters in the fourth
#data in the fifth and sixth lines
dat = list(d13C.mu.est = Afr.d13C.mu.cor, d13C.omega.est = Afr.d13C.vcov,
           conc.mu.est = Afr.conc.mu, conc.omega.est = Afr.conc.vcov, 
           epsilon.app.mu.est = map_rec.epsilon.app.mu, epsilon.app.omega.est = map_rec.epsilon.app.vcov,
           I = I, N = N, K = K, 
           RA.mix = RA.33.5ka, d13C.mix = d13C.33.5ka, d13C.mea.sd = d13C.sd.33.5ka,
           d2H.mix = d2H.33.5ka.ivc, d2H.mea.sd = d2H.sd.33.5ka)

#model parameters to save
parameters <- c("d13C.mix.m","RA.mix.m", "FLMC","f.sum.conc_n_i","exp.conc_k","d13C.k",
                "d2H.mix.m", "d2H.k", "d2H.MAP", "epsilon.app.k")

#Start time
t1 = proc.time()

set.seed(t1[3])
#Run it
Wang.33.5ka.4ch = do.call(jags.parallel,list(model.file = "code/LipidMM-JAGS-Multinorm-plus-d2H.R", 
                                             parameters.to.save = parameters, 
                                             data = dat, n.chains = 5, n.iter = n.iter, 
                                             n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1

#use rhat to check convergence
Wang.33.5ka.4ch$BUGSoutput$summary[1:3,]
#traceplots
traplot(Wang.33.5ka.4ch, parms = "FLMC")

save(Wang.33.5ka.4ch, file = "out/Wang_33.5ka_results.RData")

######2nd data point: 30.0 Ka#####
#relative abundance among the three chains from raw concentration values
#in the order of n-C27, n-C29, n-C31, c-C33

RA.30.0ka <- c(194.9, 409.8, 481.9, 386.3)/(194.9 + 409.8 + 481.9 + 386.3)

#d2H
d2H.30.0ka <- c(-141.9, -141, -152.2, -159.7)
#use ice volume correcte values
d2H.30.0ka.ivc <- c(-147.5, -146.7, -157.8, -165.3)

#analytical precision
d2H.sd.30.0ka <- c(0.7, 0.5, 0.5, 0.5)

#d13C
d13C.30.0ka <- c(-26.5, -28.4, -26.2, -24)

#analytical precision
d13C.sd.30.0ka <- c(0.2, 0.04, 0.2, 0.1)

###correction for atm CO2 d13C, should set it as a parameter?
#-8.3 year 2010 (Graven et al 2017)
#no data on 30.0 Ka, so -6.4 is used here as an approximate value
#according to Schmitt et al 2012, the values don't vary too much, at +- 0.2 per mil
#so this should have a very minor effect on the reconstructed vegetation
Afr.d13C.mu.cor <- Afr.d13C.mu + 8.3 - 6.4 

##Data to pass to the model
#prior parameters in the first three lines
#model parameters in the fourth
#data in the fifth and sixth lines
dat = list(d13C.mu.est = Afr.d13C.mu.cor, d13C.omega.est = Afr.d13C.vcov,
           conc.mu.est = Afr.conc.mu, conc.omega.est = Afr.conc.vcov, 
           epsilon.app.mu.est = map_rec.epsilon.app.mu, epsilon.app.omega.est = map_rec.epsilon.app.vcov,
           I = I, N = N, K = K, 
           RA.mix = RA.30.0ka, d13C.mix = d13C.30.0ka, d13C.mea.sd = d13C.sd.30.0ka,
           d2H.mix = d2H.30.0ka.ivc, d2H.mea.sd = d2H.sd.30.0ka)

#model parameters to save
parameters <- c("d13C.mix.m","RA.mix.m", "FLMC","f.sum.conc_n_i","exp.conc_k","d13C.k",
                "d2H.mix.m", "d2H.k", "d2H.MAP", "epsilon.app.k")

#Start time
t1 = proc.time()

set.seed(t1[3])
#Run it
Wang.30.0ka.4ch = do.call(jags.parallel,list(model.file = "code/LipidMM-JAGS-Multinorm-plus-d2H.R", 
                                             parameters.to.save = parameters, 
                                             data = dat, n.chains = 5, n.iter = n.iter, 
                                             n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1

#use rhat to check convergence
Wang.30.0ka.4ch$BUGSoutput$summary[1:3,]
#traceplots
traplot(Wang.30.0ka.4ch, parms = "FLMC")

save(Wang.30.0ka.4ch, file = "out/Wang_30.0ka_results.RData")

######3rd data point: 22.7 Ka, LGM#####
#relative abundance among the three chains from raw concentration values
#in the order of n-C27, n-C29, n-C31, c-C33

RA.22.7ka <- c(410.4, 927.1, 1033, 675.3)/(410.4 + 927.1 + 1033 + 675.3)

#d2H
d2H.22.7ka <- c(-146.9, -144.2, -157.1, -162.1)
#use ice volume correcte values
d2H.22.7ka.ivc <- c(-153.7, -150.9, -163.8, -168.7)

#analytical precision
d2H.sd.22.7ka <- c(1.1, 0.3, 0.3, 0.5)

#d13C
d13C.22.7ka <- c(-26.2, -29.5, -27, -23.9)

#analytical precision
d13C.sd.22.7ka <- c(0.2, 0.1, 0.1, 0.2)

###correction for atm CO2 d13C, should set it as a parameter?
#-8.3 year 2010 (Graven et al 2017)
#-6.45 22.7 Ka (Schmitt et al 2012)
Afr.d13C.mu.cor <- Afr.d13C.mu + 8.3 - 6.45 

##Data to pass to the model
#prior parameters in the first three lines
#model parameters in the fourth
#data in the fifth and sixth lines
dat = list(d13C.mu.est = Afr.d13C.mu.cor, d13C.omega.est = Afr.d13C.vcov,
           conc.mu.est = Afr.conc.mu, conc.omega.est = Afr.conc.vcov, 
           epsilon.app.mu.est = map_rec.epsilon.app.mu, epsilon.app.omega.est = map_rec.epsilon.app.vcov,
           I = I, N = N, K = K, 
           RA.mix = RA.22.7ka, d13C.mix = d13C.22.7ka, d13C.mea.sd = d13C.sd.22.7ka,
           d2H.mix = d2H.22.7ka.ivc, d2H.mea.sd = d2H.sd.22.7ka)

#model parameters to save
parameters <- c("d13C.mix.m","RA.mix.m", "FLMC","f.sum.conc_n_i","exp.conc_k","d13C.k",
                "d2H.mix.m", "d2H.k", "d2H.MAP", "epsilon.app.k")

#Start time
t1 = proc.time()

set.seed(t1[3])
#Run it
Wang.22.7ka.4ch = do.call(jags.parallel,list(model.file = "code/LipidMM-JAGS-Multinorm-plus-d2H.R", 
                                             parameters.to.save = parameters, 
                                             data = dat, n.chains = 5, n.iter = n.iter, 
                                             n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1

#use rhat to check convergence
Wang.22.7ka.4ch$BUGSoutput$summary[1:3,]
#traceplots
traplot(Wang.22.7ka.4ch, parms = "FLMC")

save(Wang.22.7ka.4ch, file = "out/Wang_22.7ka_results.RData")

######4th data point: 18.6 Ka#####
#relative abundance among the three chains from raw concentration values
#in the order of n-C27, n-C29, n-C31, c-C33

RA.18.6ka <- c(287.6, 693, 635.4, 364.5)/(287.6 + 693 + 635.4 + 364.5)

#d2H
d2H.18.6ka <- c(-147, -139.6, -150.4, -158.6)
#use ice volume correcte values
d2H.18.6ka.ivc <- c(-154, -146.6, -157.4, -165.5)

#analytical precision
d2H.sd.18.6ka <- c(1, 0.5, 0.5, 0.5)

#d13C
d13C.18.6ka <- c(-27.6, -30.9, -28.6, -24.9)

#analytical precision
d13C.sd.18.6ka <- c(0.1, 0.04, 0.1, 0.1)

###correction for atm CO2 d13C, should set it as a parameter?
#-8.3 year 2010 (Graven et al 2017)
#-6.43 18.6 Ka (Schmitt et al 2012)
Afr.d13C.mu.cor <- Afr.d13C.mu + 8.3 - 6.43 

##Data to pass to the model
#prior parameters in the first three lines
#model parameters in the fourth
#data in the fifth and sixth lines
dat = list(d13C.mu.est = Afr.d13C.mu.cor, d13C.omega.est = Afr.d13C.vcov,
           conc.mu.est = Afr.conc.mu, conc.omega.est = Afr.conc.vcov, 
           epsilon.app.mu.est = map_rec.epsilon.app.mu, epsilon.app.omega.est = map_rec.epsilon.app.vcov,
           I = I, N = N, K = K, 
           RA.mix = RA.18.6ka, d13C.mix = d13C.18.6ka, d13C.mea.sd = d13C.sd.18.6ka,
           d2H.mix = d2H.18.6ka.ivc, d2H.mea.sd = d2H.sd.18.6ka)

#model parameters to save
parameters <- c("d13C.mix.m","RA.mix.m", "FLMC","f.sum.conc_n_i","exp.conc_k","d13C.k",
                "d2H.mix.m", "d2H.k", "d2H.MAP", "epsilon.app.k")

#Start time
t1 = proc.time()

set.seed(t1[3])
#Run it
Wang.18.6ka.4ch = do.call(jags.parallel,list(model.file = "code/LipidMM-JAGS-Multinorm-plus-d2H.R", 
                                             parameters.to.save = parameters, 
                                             data = dat, n.chains = 5, n.iter = n.iter, 
                                             n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1

#use rhat to check convergence
Wang.18.6ka.4ch$BUGSoutput$summary[1:3,]
#traceplots
traplot(Wang.18.6ka.4ch, parms = "FLMC")

save(Wang.18.6ka.4ch, file = "out/Wang_18.6ka_results.RData")

######4th data point: 15.8 Ka, HS1#####
#relative abundance among the three chains from raw concentration values
#in the order of n-C27, n-C29, n-C31, c-C33

RA.15.8ka <- c(181.6, 588.7, 541.8, 297.7)/(181.6 + 588.7 + 541.8 + 297.7)

#d2H
d2H.15.8ka <- c(-137.8, -134.9, -145, -148.3)
#use ice volume correcte values
d2H.15.8ka.ivc <- c(-143.8, -141, -151.1, -154.3)

#analytical precision
d2H.sd.15.8ka <- c(1.3, 0.6, 0.3, 1.1)

#d13C
d13C.15.8ka <- c(-28.5, -31.8, -30.4, -26.6)

#analytical precision
d13C.sd.15.8ka <- c(0.2, 0.04, 0.1, 0.1)

###correction for atm CO2 d13C, should set it as a parameter?
#-8.3 year 2010 (Graven et al 2017)
#-6.7 15.8 Ka (Schmitt et al 2012)
Afr.d13C.mu.cor <- Afr.d13C.mu + 8.3 - 6.7 

##Data to pass to the model
#prior parameters in the first three lines
#model parameters in the fourth
#data in the fifth and sixth lines
dat = list(d13C.mu.est = Afr.d13C.mu.cor, d13C.omega.est = Afr.d13C.vcov,
           conc.mu.est = Afr.conc.mu, conc.omega.est = Afr.conc.vcov, 
           epsilon.app.mu.est = map_rec.epsilon.app.mu, epsilon.app.omega.est = map_rec.epsilon.app.vcov,
           I = I, N = N, K = K, 
           RA.mix = RA.15.8ka, d13C.mix = d13C.15.8ka, d13C.mea.sd = d13C.sd.15.8ka,
           d2H.mix = d2H.15.8ka.ivc, d2H.mea.sd = d2H.sd.15.8ka)

#model parameters to save
parameters <- c("d13C.mix.m","RA.mix.m", "FLMC","f.sum.conc_n_i","exp.conc_k","d13C.k",
                "d2H.mix.m", "d2H.k", "d2H.MAP", "epsilon.app.k")

#Start time
t1 = proc.time()

set.seed(t1[3])
#Run it
Wang.15.8ka.4ch = do.call(jags.parallel,list(model.file = "code/LipidMM-JAGS-Multinorm-plus-d2H.R", 
                                             parameters.to.save = parameters, 
                                             data = dat, n.chains = 5, n.iter = n.iter, 
                                             n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1

#use rhat to check convergence
Wang.15.8ka.4ch$BUGSoutput$summary[1:3,]
#traceplots
traplot(Wang.15.8ka.4ch, parms = "FLMC")

save(Wang.15.8ka.4ch, file = "out/Wang_15.8ka_results.RData")

######5th data point: ~13.8 Ka, #####
#relative abundance among the three chains from raw concentration values
#in the order of n-C27, n-C29, n-C31, c-C33

RA.13.8ka <- c(182, 358.3, 373.3, 284.9)/(182 + 358.3 + 373.3 + 284.9)

#d2H
d2H.13.8ka <- c(-142.1, -139.2, -147.5, -154.9)
#use ice volume correcte values
d2H.13.8ka.ivc <- c(-146.6, -143.7, -152, -159.3)

#analytical precision
d2H.sd.13.8ka <- c(0.9, 0.9, 0.4, 0.5)

#d13C
d13C.13.8ka <- c(-29.3, -30, -29.2, -26.7)

#analytical precision
d13C.sd.13.8ka <- c(0.1, 0.04, 0.3, 0.1)

###correction for atm CO2 d13C, should set it as a parameter?
#-8.3 year 2010 (Graven et al 2017)
#-6.63 13.8 Ka (Schmitt et al 2012)
Afr.d13C.mu.cor <- Afr.d13C.mu + 8.3 - 6.63 

##Data to pass to the model
#prior parameters in the first three lines
#model parameters in the fourth
#data in the fifth and sixth lines
dat = list(d13C.mu.est = Afr.d13C.mu.cor, d13C.omega.est = Afr.d13C.vcov,
           conc.mu.est = Afr.conc.mu, conc.omega.est = Afr.conc.vcov, 
           epsilon.app.mu.est = map_rec.epsilon.app.mu, epsilon.app.omega.est = map_rec.epsilon.app.vcov,
           I = I, N = N, K = K, 
           RA.mix = RA.13.8ka, d13C.mix = d13C.13.8ka, d13C.mea.sd = d13C.sd.13.8ka,
           d2H.mix = d2H.13.8ka.ivc, d2H.mea.sd = d2H.sd.13.8ka)

#model parameters to save
parameters <- c("d13C.mix.m","RA.mix.m", "FLMC","f.sum.conc_n_i","exp.conc_k","d13C.k",
                "d2H.mix.m", "d2H.k", "d2H.MAP", "epsilon.app.k")

#Start time
t1 = proc.time()

set.seed(t1[3])
#Run it
Wang.13.8ka.4ch = do.call(jags.parallel,list(model.file = "code/LipidMM-JAGS-Multinorm-plus-d2H.R", 
                                             parameters.to.save = parameters, 
                                             data = dat, n.chains = 5, n.iter = n.iter, 
                                             n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1

#use rhat to check convergence
Wang.13.8ka.4ch$BUGSoutput$summary[1:3,]
#traceplots
traplot(Wang.13.8ka.4ch, parms = "FLMC")

save(Wang.13.8ka.4ch, file = "out/Wang_13.8ka_results.RData")

######6th data point: ~12.1 Ka, YD #####
#relative abundance among the three chains from raw concentration values
#in the order of n-C27, n-C29, n-C31, c-C33

RA.12.1ka <- c(146, 308, 343.2, 246.5)/(146 + 308 + 343.2 + 246.5)

#d2H
d2H.12.1ka <- c(-139.1, -138.8, -145.6, -147)
#use ice volume correcte values
d2H.12.1ka.ivc <- c(-142, -141.8, -148.5, -149.9)

#analytical precision
d2H.sd.12.1ka <- c(1.3, 0.7, 0.7, 1)

#d13C
d13C.12.1ka <- c(-29.4, -30.1, -29.8, -27.9)

#analytical precision
d13C.sd.12.1ka <- c(0.2, 0.1, 0.1, 0.2)

###correction for atm CO2 d13C, should set it as a parameter?
#-8.3 year 2010 (Graven et al 2017)
#-6.69 12.1 Ka (Schmitt et al 2012)
Afr.d13C.mu.cor <- Afr.d13C.mu + 8.3 - 6.69 

##Data to pass to the model
#prior parameters in the first three lines
#model parameters in the fourth
#data in the fifth and sixth lines
dat = list(d13C.mu.est = Afr.d13C.mu.cor, d13C.omega.est = Afr.d13C.vcov,
           conc.mu.est = Afr.conc.mu, conc.omega.est = Afr.conc.vcov, 
           epsilon.app.mu.est = map_rec.epsilon.app.mu, epsilon.app.omega.est = map_rec.epsilon.app.vcov,
           I = I, N = N, K = K, 
           RA.mix = RA.12.1ka, d13C.mix = d13C.12.1ka, d13C.mea.sd = d13C.sd.12.1ka,
           d2H.mix = d2H.12.1ka.ivc, d2H.mea.sd = d2H.sd.12.1ka)

#model parameters to save
parameters <- c("d13C.mix.m","RA.mix.m", "FLMC","f.sum.conc_n_i","exp.conc_k","d13C.k",
                "d2H.mix.m", "d2H.k", "d2H.MAP", "epsilon.app.k")

#Start time
t1 = proc.time()

set.seed(t1[3])
#Run it
Wang.12.1ka.4ch = do.call(jags.parallel,list(model.file = "code/LipidMM-JAGS-Multinorm-plus-d2H.R", 
                                             parameters.to.save = parameters, 
                                             data = dat, n.chains = 5, n.iter = n.iter, 
                                             n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1

#use rhat to check convergence
Wang.12.1ka.4ch$BUGSoutput$summary[1:3,]
#traceplots
traplot(Wang.12.1ka.4ch, parms = "FLMC")

save(Wang.12.1ka.4ch, file = "out/Wang_12.1ka_results.RData")

######6th data point: ~10.6 Ka, #####
#relative abundance among the three chains from raw concentration values
#in the order of n-C27, n-C29, n-C31, c-C33

RA.10.6ka <- c(143, 280.3, 366.1, 281.3)/(143 + 280.3 + 366.1 + 281.3)

#d2H
d2H.10.6ka <- c(-148.6, -143.2, -145.2, -146.7)
#use ice volume correcte values
d2H.10.6ka.ivc <- c(-150.3, -145, -147, -148.4)

#analytical precision
d2H.sd.10.6ka <- c(1.1, 0.4, 0.3, 0.5)

#d13C
d13C.10.6ka <- c(-28.7, -29.3, -29.5, -28.2)

#analytical precision
d13C.sd.10.6ka <- c(0.3, 0.04, 0.1, 0.3)

###correction for atm CO2 d13C, should set it as a parameter?
#-8.3 year 2010 (Graven et al 2017)
#-6.58 10.6 Ka (Schmitt et al 2012)
Afr.d13C.mu.cor <- Afr.d13C.mu + 8.3 - 6.58 

##Data to pass to the model
#prior parameters in the first three lines
#model parameters in the fourth
#data in the fifth and sixth lines
dat = list(d13C.mu.est = Afr.d13C.mu.cor, d13C.omega.est = Afr.d13C.vcov,
           conc.mu.est = Afr.conc.mu, conc.omega.est = Afr.conc.vcov, 
           epsilon.app.mu.est = map_rec.epsilon.app.mu, epsilon.app.omega.est = map_rec.epsilon.app.vcov,
           I = I, N = N, K = K, 
           RA.mix = RA.10.6ka, d13C.mix = d13C.10.6ka, d13C.mea.sd = d13C.sd.10.6ka,
           d2H.mix = d2H.10.6ka.ivc, d2H.mea.sd = d2H.sd.10.6ka)

#model parameters to save
parameters <- c("d13C.mix.m","RA.mix.m", "FLMC","f.sum.conc_n_i","exp.conc_k","d13C.k",
                "d2H.mix.m", "d2H.k", "d2H.MAP", "epsilon.app.k")

#Start time
t1 = proc.time()

set.seed(t1[3])
#Run it
Wang.10.6ka.4ch = do.call(jags.parallel,list(model.file = "code/LipidMM-JAGS-Multinorm-plus-d2H.R", 
                                             parameters.to.save = parameters, 
                                             data = dat, n.chains = 5, n.iter = n.iter, 
                                             n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1

#use rhat to check convergence
Wang.10.6ka.4ch$BUGSoutput$summary[1:3,]
#traceplots
traplot(Wang.10.6ka.4ch, parms = "FLMC")

save(Wang.10.6ka.4ch, file = "out/Wang_10.6ka_results.RData")

######7th data point: 7.8 Ka, #####
#relative abundance among the three chains from raw concentration values
#in the order of n-C27, n-C29, n-C31, c-C33

RA.7.8ka <- c(84.3, 142.1, 174.6, 142.7)/(84.3 + 142.1 + 174.6 + 142.7)

#d2H
d2H.7.8ka <- c(-142.6, -140.4, -146.2, -145.4)
#use ice volume correcte values
d2H.7.8ka.ivc <- c(-143.3, -141.1, -146.9, -146.1)

#analytical precision
d2H.sd.7.8ka <- c(1.2, 1, 0.7, 1)

#d13C
d13C.7.8ka <- c(-28.3, -29.7, -28.9, -27.6)

#analytical precision
d13C.sd.7.8ka <- c(0.7, 0.04, 0.1, 0.5)

###correction for atm CO2 d13C, should set it as a parameter?
#-8.3 year 2010 (Graven et al 2017)
#-6.39 7.8 Ka (Schmitt et al 2012)
Afr.d13C.mu.cor <- Afr.d13C.mu + 8.3 - 6.39 

##Data to pass to the model
#prior parameters in the first three lines
#model parameters in the fourth
#data in the fifth and sixth lines
dat = list(d13C.mu.est = Afr.d13C.mu.cor, d13C.omega.est = Afr.d13C.vcov,
           conc.mu.est = Afr.conc.mu, conc.omega.est = Afr.conc.vcov, 
           epsilon.app.mu.est = map_rec.epsilon.app.mu, epsilon.app.omega.est = map_rec.epsilon.app.vcov,
           I = I, N = N, K = K, 
           RA.mix = RA.7.8ka, d13C.mix = d13C.7.8ka, d13C.mea.sd = d13C.sd.7.8ka,
           d2H.mix = d2H.7.8ka.ivc, d2H.mea.sd = d2H.sd.7.8ka)

#model parameters to save
parameters <- c("d13C.mix.m","RA.mix.m", "FLMC","f.sum.conc_n_i","exp.conc_k","d13C.k",
                "d2H.mix.m", "d2H.k", "d2H.MAP", "epsilon.app.k")

#Start time
t1 = proc.time()

set.seed(t1[3])
#Run it
Wang.7.8ka.4ch = do.call(jags.parallel,list(model.file = "code/LipidMM-JAGS-Multinorm-plus-d2H.R", 
                                            parameters.to.save = parameters, 
                                            data = dat, n.chains = 5, n.iter = n.iter, 
                                            n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1

#use rhat to check convergence
Wang.7.8ka.4ch$BUGSoutput$summary[1:3,]
#traceplots
traplot(Wang.7.8ka.4ch, parms = "FLMC")

save(Wang.7.8ka.4ch, file = "out/Wang_7.8ka_results.RData")

######8th data point: 4.0 Ka, #####
#relative abundance among the three chains from raw concentration values
#in the order of n-C27, n-C29, n-C31, c-C33

RA.4.0ka <- c(144.5, 268.3, 302, 220.1)/(144.5 + 268.3 + 302 + 220.1)

#d2H
d2H.4.0ka <- c(-144.8, -143.1, -151.1, -152.8)
#use ice volume correcte values
d2H.4.0ka.ivc <- c(-145.1, -143.3, -151.3, -153)

#analytical precision
d2H.sd.4.0ka <- c(0.6, 0.4, 0.4, 0.7)

#d13C
d13C.4.0ka <- c(-28.3, -28.9, -29.1, -28)

#analytical precision
d13C.sd.4.0ka <- c(0.1, 0.2, 0.1, 0.2)

###correction for atm CO2 d13C, should set it as a parameter?
#-8.3 year 2010 (Graven et al 2017)
#-6.34 4.0 Ka (Schmitt et al 2012)
Afr.d13C.mu.cor <- Afr.d13C.mu + 8.3 - 6.34 

##Data to pass to the model
#prior parameters in the first three lines
#model parameters in the fourth
#data in the fifth and sixth lines
dat = list(d13C.mu.est = Afr.d13C.mu.cor, d13C.omega.est = Afr.d13C.vcov,
           conc.mu.est = Afr.conc.mu, conc.omega.est = Afr.conc.vcov, 
           epsilon.app.mu.est = map_rec.epsilon.app.mu, epsilon.app.omega.est = map_rec.epsilon.app.vcov,
           I = I, N = N, K = K, 
           RA.mix = RA.4.0ka, d13C.mix = d13C.4.0ka, d13C.mea.sd = d13C.sd.4.0ka,
           d2H.mix = d2H.4.0ka.ivc, d2H.mea.sd = d2H.sd.4.0ka)

#model parameters to save
parameters <- c("d13C.mix.m","RA.mix.m", "FLMC","f.sum.conc_n_i","exp.conc_k","d13C.k",
                "d2H.mix.m", "d2H.k", "d2H.MAP", "epsilon.app.k")

#Start time
t1 = proc.time()

set.seed(t1[3])
#Run it
Wang.4.0ka.4ch = do.call(jags.parallel,list(model.file = "code/LipidMM-JAGS-Multinorm-plus-d2H.R", 
                                            parameters.to.save = parameters, 
                                            data = dat, n.chains = 5, n.iter = n.iter, 
                                            n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1

#use rhat to check convergence
Wang.4.0ka.4ch$BUGSoutput$summary[1:3,]
#traceplots
traplot(Wang.4.0ka.4ch, parms = "FLMC")

save(Wang.4.0ka.4ch, file = "out/Wang_4.0ka_results.RData")

###summary of results###
load("out/Wang_33.5ka_results.RData")
load("out/Wang_30.0ka_results.RData")
load("out/Wang_22.7ka_results.RData")
load("out/Wang_18.6ka_results.RData")
load("out/Wang_15.8ka_results.RData")
load("out/Wang_13.8ka_results.RData")
load("out/Wang_12.1ka_results.RData")
load("out/Wang_10.6ka_results.RData")
load("out/Wang_7.8ka_results.RData")
load("out/Wang_4.0ka_results.RData")

#### MAPs, medians, and 89% HDIs ####
Wang.33.5ka.GR <- MAP_HDI(Wang.33.5ka.4ch$BUGSoutput$sims.list$FLMC[,1],.89)
Wang.33.5ka.SV <- MAP_HDI(Wang.33.5ka.4ch$BUGSoutput$sims.list$FLMC[,2],.89)
Wang.33.5ka.RF <- MAP_HDI(Wang.33.5ka.4ch$BUGSoutput$sims.list$FLMC[,3],.89)
Wang.33.5ka.d2H.MAP <- MAP_HDI(Wang.33.5ka.4ch$BUGSoutput$sims.list$d2H.MAP,.89)

Wang.30.0ka.GR <- MAP_HDI(Wang.30.0ka.4ch$BUGSoutput$sims.list$FLMC[,1],.89)
Wang.30.0ka.SV <- MAP_HDI(Wang.30.0ka.4ch$BUGSoutput$sims.list$FLMC[,2],.89)
Wang.30.0ka.RF <- MAP_HDI(Wang.30.0ka.4ch$BUGSoutput$sims.list$FLMC[,3],.89)
Wang.30.0ka.d2H.MAP <- MAP_HDI(Wang.30.0ka.4ch$BUGSoutput$sims.list$d2H.MAP,.89)

Wang.22.7ka.GR <- MAP_HDI(Wang.22.7ka.4ch$BUGSoutput$sims.list$FLMC[,1],.89)
Wang.22.7ka.SV <- MAP_HDI(Wang.22.7ka.4ch$BUGSoutput$sims.list$FLMC[,2],.89)
Wang.22.7ka.RF <- MAP_HDI(Wang.22.7ka.4ch$BUGSoutput$sims.list$FLMC[,3],.89)
Wang.22.7ka.d2H.MAP <- MAP_HDI(Wang.22.7ka.4ch$BUGSoutput$sims.list$d2H.MAP,.89)

Wang.18.6ka.GR <- MAP_HDI(Wang.18.6ka.4ch$BUGSoutput$sims.list$FLMC[,1],.89)
Wang.18.6ka.SV <- MAP_HDI(Wang.18.6ka.4ch$BUGSoutput$sims.list$FLMC[,2],.89)
Wang.18.6ka.RF <- MAP_HDI(Wang.18.6ka.4ch$BUGSoutput$sims.list$FLMC[,3],.89)
Wang.18.6ka.d2H.MAP <- MAP_HDI(Wang.18.6ka.4ch$BUGSoutput$sims.list$d2H.MAP,.89)

Wang.15.8ka.GR <- MAP_HDI(Wang.15.8ka.4ch$BUGSoutput$sims.list$FLMC[,1],.89)
Wang.15.8ka.SV <- MAP_HDI(Wang.15.8ka.4ch$BUGSoutput$sims.list$FLMC[,2],.89)
Wang.15.8ka.RF <- MAP_HDI(Wang.15.8ka.4ch$BUGSoutput$sims.list$FLMC[,3],.89)
Wang.15.8ka.d2H.MAP <- MAP_HDI(Wang.15.8ka.4ch$BUGSoutput$sims.list$d2H.MAP,.89)

Wang.13.8ka.GR <- MAP_HDI(Wang.13.8ka.4ch$BUGSoutput$sims.list$FLMC[,1],.89)
Wang.13.8ka.SV <- MAP_HDI(Wang.13.8ka.4ch$BUGSoutput$sims.list$FLMC[,2],.89)
Wang.13.8ka.RF <- MAP_HDI(Wang.13.8ka.4ch$BUGSoutput$sims.list$FLMC[,3],.89)
Wang.13.8ka.d2H.MAP <- MAP_HDI(Wang.13.8ka.4ch$BUGSoutput$sims.list$d2H.MAP,.89)

Wang.12.1ka.GR <- MAP_HDI(Wang.12.1ka.4ch$BUGSoutput$sims.list$FLMC[,1],.89)
Wang.12.1ka.SV <- MAP_HDI(Wang.12.1ka.4ch$BUGSoutput$sims.list$FLMC[,2],.89)
Wang.12.1ka.RF <- MAP_HDI(Wang.12.1ka.4ch$BUGSoutput$sims.list$FLMC[,3],.89)
Wang.12.1ka.d2H.MAP <- MAP_HDI(Wang.12.1ka.4ch$BUGSoutput$sims.list$d2H.MAP,.89)

Wang.10.6ka.GR <- MAP_HDI(Wang.10.6ka.4ch$BUGSoutput$sims.list$FLMC[,1],.89)
Wang.10.6ka.SV <- MAP_HDI(Wang.10.6ka.4ch$BUGSoutput$sims.list$FLMC[,2],.89)
Wang.10.6ka.RF <- MAP_HDI(Wang.10.6ka.4ch$BUGSoutput$sims.list$FLMC[,3],.89)
Wang.10.6ka.d2H.MAP <- MAP_HDI(Wang.10.6ka.4ch$BUGSoutput$sims.list$d2H.MAP,.89)

Wang.7.8ka.GR <- MAP_HDI(Wang.7.8ka.4ch$BUGSoutput$sims.list$FLMC[,1],.89)
Wang.7.8ka.SV <- MAP_HDI(Wang.7.8ka.4ch$BUGSoutput$sims.list$FLMC[,2],.89)
Wang.7.8ka.RF <- MAP_HDI(Wang.7.8ka.4ch$BUGSoutput$sims.list$FLMC[,3],.89)
Wang.7.8ka.d2H.MAP <- MAP_HDI(Wang.7.8ka.4ch$BUGSoutput$sims.list$d2H.MAP,.89)

Wang.4.0ka.GR <- MAP_HDI(Wang.4.0ka.4ch$BUGSoutput$sims.list$FLMC[,1],.89)
Wang.4.0ka.SV <- MAP_HDI(Wang.4.0ka.4ch$BUGSoutput$sims.list$FLMC[,2],.89)
Wang.4.0ka.RF <- MAP_HDI(Wang.4.0ka.4ch$BUGSoutput$sims.list$FLMC[,3],.89)
Wang.4.0ka.d2H.MAP <- MAP_HDI(Wang.4.0ka.4ch$BUGSoutput$sims.list$d2H.MAP,.89)

### organize data for analysis and plots ###
Wang.GR <- rbind(Wang.33.5ka.GR, Wang.30.0ka.GR,Wang.22.7ka.GR, 
                 Wang.18.6ka.GR, Wang.15.8ka.GR, Wang.13.8ka.GR, 
                 Wang.12.1ka.GR, Wang.10.6ka.GR,
                 Wang.7.8ka.GR, Wang.4.0ka.GR)
colnames(Wang.GR) <- c("Wang.GR.MAP","Wang.GR.hdiL","Wang.GR.hdiH")

Wang.SV <- rbind(Wang.33.5ka.SV, Wang.30.0ka.SV, Wang.22.7ka.SV, 
                 Wang.18.6ka.SV, Wang.15.8ka.SV, Wang.13.8ka.SV, 
                 Wang.12.1ka.SV, Wang.10.6ka.SV,
                 Wang.7.8ka.SV, Wang.4.0ka.SV)
colnames(Wang.SV) <- c("Wang.SV.MAP","Wang.SV.hdiL","Wang.SV.hdiH")

Wang.RF <- rbind(Wang.33.5ka.RF, Wang.30.0ka.RF, Wang.22.7ka.RF, 
                 Wang.18.6ka.RF, Wang.15.8ka.RF, Wang.13.8ka.RF, 
                 Wang.12.1ka.RF, Wang.10.6ka.RF,
                 Wang.7.8ka.RF, Wang.4.0ka.RF)
colnames(Wang.RF) <- c("Wang.RF.MAP","Wang.RF.hdiL","Wang.RF.hdiH")

Wang.d2H.MAP <- rbind(Wang.33.5ka.d2H.MAP, Wang.30.0ka.d2H.MAP, Wang.22.7ka.d2H.MAP, 
                      Wang.18.6ka.d2H.MAP, Wang.15.8ka.d2H.MAP, Wang.13.8ka.d2H.MAP,
                      Wang.12.1ka.d2H.MAP, Wang.10.6ka.d2H.MAP,
                      Wang.7.8ka.d2H.MAP, Wang.4.0ka.d2H.MAP)
colnames(Wang.d2H.MAP) <- c("Wang.d2H_MAP.MAP","Wang.d2H_MAP.hdiL","Wang.d2H_MAP.hdiH")

Wang.RA <- rbind(RA.33.5ka, RA.30.0ka, RA.22.7ka, RA.18.6ka, RA.15.8ka,
                 RA.13.8ka, RA.12.1ka, RA.10.6ka, RA.7.8ka, RA.4.0ka)
colnames(Wang.RA) <- c("RA.n-C27", "RA.n-C29", "RA.n-C31", "RA.n-C33")

Wang.d13C <- rbind(d13C.33.5ka, d13C.30.0ka, d13C.22.7ka, d13C.18.6ka, d13C.15.8ka,
                   d13C.13.8ka, d13C.12.1ka, d13C.10.6ka, d13C.7.8ka, d13C.4.0ka)
colnames(Wang.d13C) <- c("d13C.n-C27", "d13C.n-C29", "d13C.n-C31", "d13C.n-C33")

Wang.age <- c(33.375, 29.889, 22.691, 18.598, 15.785, 13.744, 12.063, 10.556, 7.985, 3.949)

Wang.d2H.C29.ivc <- c(-148.7, -146.7, -150.9, -146.6, -141, -143.7, -141.8, -145, -141.1, -143.3)
Wang.d2H.C31.ivc <- c(-159.8, -157.8, -163.8, -157.4, -151.1, -152, -148.5, -147, -146.9, -151.3)
Wang.veg <- data.frame(Wang.age, Wang.RA, Wang.d13C, Wang.GR, Wang.SV, Wang.RF, Wang.d2H.MAP, 
                       Wang.d2H.C29.ivc, Wang.d2H.C31.ivc)

#correlatin tests between estimated MAP d2H and d2H of two n-alkane chains
#n-C29
cor.test(Wang.veg$Wang.d2H.C29.ivc,Wang.veg$Wang.d2H_MAP.MAP)

#Pearson's product-moment correlation
#data:  Wang.veg$Wang.d2H.C29.ivc and Wang.veg$Wang.d2H_MAP.MAP
#t = 2.8786, df = 8, p-value = 0.02056
#r = 0.7132981

lm.C29.ivc <- lm(Wang.veg$Wang.d2H_MAP.MAP ~ Wang.veg$Wang.d2H.C29.ivc)
summary(lm.C29.ivc)
#                           Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                86.2794    42.7030   2.020   0.0780 .
#Wang.veg$Wang.d2H.C29.ivc   0.8483     0.2947   2.879   0.0206 * 

#n-C31
cor.test(Wang.veg$Wang.d2H.C31.ivc,Wang.veg$Wang.d2H_MAP.MAP)

#Pearson's product-moment correlation
#data:  Wang.veg$Wang.d2H.C29.ivc and Wang.veg$Wang.d2H_MAP.MAP
#t = 2.4703, df = 8, p-value = 0.03869
#r = 0.6578097 

#plot(Wang.veg$Wang.d2H_MAP.MAP ~ Wang.veg$Wang.d2H.C29.ivc)

lm.C31.ivc <- lm(Wang.veg$Wang.d2H_MAP.MAP ~ Wang.veg$Wang.d2H.C31.ivc)
summary(lm.C31.ivc)
#                           Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                32.3252    27.9268   1.157   0.2805  
#Wang.veg$Wang.d2H.C31.ivc   0.4490     0.1817   2.470   0.0387 * 

#plot(Wang.veg$Wang.d2H_MAP.MAP ~ Wang.veg$Wang.d2H.C31.ivc)

###n-C29 has a higher correlation coefficient


####Qualitative comparison with the GeoB9311 pollen record####
####loading GeoB9311 pollen record
pollen_rec<- read.csv("data/EA-6 GeoB9311 pollen.csv")
#note that the data are at a much lower temporal resolutoin
#12.1 ka data point not available

#use taxon groups by Dupont and Kuhlmann, 2017
GeoB9311.RF <- subset(pollen_rec, Habitat == "RF")
GeoB9311.SV <- subset(pollen_rec, Habitat == "SV")
GeoB9311.C4 <- subset(pollen_rec, Habitat == "C4")
GeoB9311.MT <- subset(pollen_rec, Habitat == "MT")
GeoB9311.SW <- subset(pollen_rec, Habitat == "SW")
GeoB9311.Fern <- subset(pollen_rec, Habitat == "Fern")

GeoB9311.RF.sum <- dplyr::summarise_each(GeoB9311.RF[4:11], list(total = sum))
GeoB9311.SV.sum <- dplyr::summarise_each(GeoB9311.SV[4:11], list(total = sum))
GeoB9311.C4.sum <- dplyr::summarise_each(GeoB9311.C4[4:11], list(total = sum))
GeoB9311.MT.sum <- dplyr::summarise_each(GeoB9311.MT[4:11], list(total = sum))
GeoB9311.SW.sum <- dplyr::summarise_each(GeoB9311.SW[4:11], list(total = sum))
GeoB9311.Fern.sum <- dplyr::summarise_each(GeoB9311.Fern[4:11], list(total = sum))

#fern spores excluded in the analysis#
GeoB9311.pollen.sum <- GeoB9311.RF.sum + GeoB9311.SV.sum + GeoB9311.C4.sum + GeoB9311.MT.sum + GeoB9311.SW.sum

GeoB9311.RF.pf <- as.numeric(GeoB9311.RF.sum/GeoB9311.pollen.sum)
GeoB9311.SV.pf <- as.numeric(GeoB9311.SV.sum/GeoB9311.pollen.sum)
GeoB9311.C4.pf <- as.numeric(GeoB9311.C4.sum/GeoB9311.pollen.sum)
GeoB9311.MT.pf <- as.numeric(GeoB9311.MT.sum/GeoB9311.pollen.sum)
GeoB9311.SW.pf <- as.numeric(GeoB9311.SW.sum/GeoB9311.pollen.sum)

##grouping Afromontane and swamp taxa into RF for comparison
GeoB9311.forest.pf<- GeoB9311.RF.pf + GeoB9311.MT.pf + GeoB9311.SW.pf

###pull out the sequence that is comparable between the two records
Wang.indeces <- c(1:3,5, 6, 8:10)
Wang.veg.compare <- Wang.veg[Wang.indeces,]

####qualitative comparison with non-parametric tests#####
#spearman's rank correlation#

cor.test(Wang.veg.compare$Wang.GR.MAP, rev(GeoB9311.C4.pf), method = "spearman")
#S = 24, p-value = 0.05759
#      rho 
#0.7142857
cor.test(Wang.veg.compare$Wang.SV.MAP, rev(GeoB9311.SV.pf), method = "spearman")
#S = 32, p-value = 0.115
#rho 
#0.6190476 
cor.test(Wang.veg.compare$Wang.RF.MAP, rev(GeoB9311.forest.pf), method = "spearman")
#S = 116, p-value = 0.3599
#       rho 
#-0.3809524

cor.test(Wang.veg.compare$Wang.GR.MAP, rev(GeoB9311.C4.pf), method = "kendall")
#T = 22, p-value = 0.06101
#      tau 
#0.5714286  
cor.test(Wang.veg.compare$Wang.SV.MAP, rev(GeoB9311.SV.pf), method = "kendall")
#T = 20, p-value = 0.1789
#      tau 
#0.4285714
cor.test(Wang.veg.compare$Wang.RF.MAP, rev(GeoB9311.forest.pf), method = "kendall")
#T = 10, p-value = 0.3988
#       tau 
#-0.2857143