library(coda)
library(rjags)
library(R2jags)
library(mcmcplots)
library(MASS)
library(viridisLite)
library(EnvStats)
library(bayestestR)

setwd("C:/Users/ydmag/Google Drive/U of U/Proxy project/LipidMM")

####make sure to run all the helper functions in "code/LipidMM-helper functions.R"
map_rec<- read.csv("data/EA-5 data map_rec.csv")

map_rec.GR <- map_rec[which(map_rec$Source == "GR"),]
map_rec.SV <- map_rec[which(map_rec$Source == "SV"),]
map_rec.RF <- map_rec[which(map_rec$Source == "RF"),]

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
GR <- African[which(African$Source == "GR"),]
SV <- African[which(African$Source == "SV"),]
RF <- African[which(African$Source == "RF"),]

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

#common MCMC parameters
n.iter = 8e5
n.burnin = 2e5
n.thin = floor(n.iter-n.burnin)/2500
#average runtime is ~4 hours/sample

#MAP reconstruction based on Core GIK16160-3, Zambezi River mouth, Wang et al. 2013 
####first data point: LGM (20.1 ka), HS1 (16.1 ka), YD (12.1 Ka), Holocene (2.15 Ka)


#relative abundance among the three chains from raw concentration values
#in the order of n-C27, n-C29, n-C31, n-C33
RA.LGM <- c(195.9, 430.3, 465, 313.3)/(195.9 + 430.3 + 465 + 313.3)

#d2H, 
d2H.LGM <- c(-146.2, -141.3, -154.3, -158.9)
#use ice volume correcte values
d2H.LGM.ivc <- c(-153.1, -148.3, -161.2, -165.7)

#analytical precision
d2H.sd.LGM <- c(1.2, 0.6, 1, 0.9)

#d13C
d13C.LGM <- c(-26.9, -29.8, -27.8, -24.9)

#analytical precision
d13C.sd.LGM <- c(0.1, 0.1, 0.1, 0.2)

###correction for atm CO2 d13C, should set it as a parameter?
#-8.3 year 2010 (Graven et al 2017)
#-6.5 LGM (Schmitt et al 2012)
Afr.d13C.mu.cor <- Afr.d13C.mu + 8.3 - 6.5 

##Data to pass to the model
#prior parameters in the first three lines
#model parameters in the fourth
#data in the fifth and sixth lines
dat = list(d13C.mu.est = Afr.d13C.mu.cor, d13C.omega.est = Afr.d13C.vcov,
           conc.mu.est = Afr.conc.mu, conc.omega.est = Afr.conc.vcov, 
           epsilon.app.mu.est = map_rec.epsilon.app.mu, epsilon.app.omega.est = map_rec.epsilon.app.vcov,
           I = I, N = N, K = K, 
           RA.mix = RA.LGM, d13C.mix = d13C.LGM, d13C.mea.sd = d13C.sd.LGM,
           d2H.mix = d2H.LGM.ivc, d2H.mea.sd = d2H.sd.LGM)

#model parameters to save
parameters <- c("d13C.mix.m","RA.mix.m", "FLMC","f.sum.conc_n_i","exp.conc_k","d13C.k",
                "d2H.mix.m", "d2H.k", "d2H.MAP", "epsilon.app.k")

#Start time
t1 = proc.time()

set.seed(t1[3])
#Run it
Wang.LGM.4ch = do.call(jags.parallel,list(model.file = "code/LipidMM-JAGS-Multinorm-plus-d2H.R", 
                                          parameters.to.save = parameters, 
                                          data = dat, n.chains = 3, n.iter = n.iter, 
                                          n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~5 hours

#use rhat to check convergence
Wang.LGM.4ch$BUGSoutput$summary[1:3,]
#traceplots
traplot(Wang.LGM.4ch, parms = "FLMC")

save(Wang.LGM.4ch, file = "out/Wang_LGM_results.RData")

#relative abundance among the three chains from raw concentration values
#in the order of n-C27, n-C29, n-C31, n-C33

RA.HS1 <- c(305.8, 962.3, 885.9, 473.2)/(305.8 + 962.3 + 885.9 + 473.2)

#d2H
d2H.HS1 <- c(-142.3, -137.3, -148.2, -152.5)
#use ice volume correcte values
d2H.HS1.ivc <- c(-148.6, -143.6, -154.4, -158.7)

#analytical precision
d2H.sd.HS1 <- c(0.7, 0.2, 0.4, 0.6)

#d13C
d13C.HS1 <- c(-28.6, -31.8, -30.3, -27)

#analytical precision
d13C.sd.HS1 <- c(0.2, 0.1, 0.1, 0.1)

###correction for atm CO2 d13C, should set it as a parameter?
#-8.3 year 2010 (Graven et al 2017)
#-6.7 HS1 (Schmitt et al 2012)
Afr.d13C.mu.cor <- Afr.d13C.mu + 8.3 - 6.7 

##Data to pass to the model
#prior parameters in the first three lines
#model parameters in the fourth
#data in the fifth and sixth lines
dat = list(d13C.mu.est = Afr.d13C.mu.cor, d13C.omega.est = Afr.d13C.vcov,
           conc.mu.est = Afr.conc.mu, conc.omega.est = Afr.conc.vcov, 
           epsilon.app.mu.est = map_rec.epsilon.app.mu, epsilon.app.omega.est = map_rec.epsilon.app.vcov,
           I = I, N = N, K = K, 
           RA.mix = RA.HS1, d13C.mix = d13C.HS1, d13C.mea.sd = d13C.sd.HS1,
           d2H.mix = d2H.HS1.ivc, d2H.mea.sd = d2H.sd.HS1)

#model parameters to save
parameters <- c("d13C.mix.m","RA.mix.m", "FLMC","f.sum.conc_n_i","exp.conc_k","d13C.k",
                "d2H.mix.m", "d2H.k", "d2H.MAP", "epsilon.app.k")

#Start time
t1 = proc.time()

set.seed(t1[3])
#Run it
Wang.HS1.4ch = do.call(jags.parallel,list(model.file = "code/LipidMM-JAGS-Multinorm-plus-d2H.R", 
                                          parameters.to.save = parameters, 
                                          data = dat, n.chains = 3, n.iter = n.iter, 
                                          n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1

#use rhat to check convergence
Wang.HS1.4ch$BUGSoutput$summary[1:3,]
#traceplots
traplot(Wang.HS1.4ch, parms = "FLMC")

save(Wang.HS1.4ch, file = "out/Wang_HS1_results.RData")

#relative abundance among the three chains from raw concentration values
#in the order of n-C27, n-C29, n-C31, n-C33

RA.YD <- c(146, 308, 343.2, 246.5)/(146 + 308 + 343.2 + 246.5)

#d2H 
d2H.YD <- c(-139.1, -138.8, -145.6, -147)
#use ice volume correcte values
d2H.YD.ivc <- c(-142, -141.8, -148.5, -149.9)

#analytical precision
d2H.sd.YD <- c(1.3, 0.7, 0.7, 1)

#d13C
d13C.YD <- c(-29.4, -30.1, -29.8, -27.9)

#analytical precision
d13C.sd.YD <- c(0.2, 0.1, 0.1, 0.2)

###correction for atm CO2 d13C, should set it as a parameter?
#-8.3 year 2010 (Graven et al 2017)
#-6.7 YD (Schmitt et al 2012)
Afr.d13C.mu.cor <- Afr.d13C.mu + 8.3 - 6.7 

##Data to pass to the model
#prior parameters in the first three lines
#model parameters in the fourth
#data in the fifth and sixth lines
dat = list(d13C.mu.est = Afr.d13C.mu.cor, d13C.omega.est = Afr.d13C.vcov,
           conc.mu.est = Afr.conc.mu, conc.omega.est = Afr.conc.vcov, 
           epsilon.app.mu.est = map_rec.epsilon.app.mu, epsilon.app.omega.est = map_rec.epsilon.app.vcov,
           I = I, N = N, K = K, 
           RA.mix = RA.YD, d13C.mix = d13C.YD, d13C.mea.sd = d13C.sd.YD,
           d2H.mix = d2H.YD.ivc, d2H.mea.sd = d2H.sd.YD)

#model parameters to save
parameters <- c("d13C.mix.m","RA.mix.m", "FLMC","f.sum.conc_n_i","exp.conc_k","d13C.k",
                "d2H.mix.m", "d2H.k", "d2H.MAP", "epsilon.app.k")

#Start time
t1 = proc.time()

set.seed(t1[3])
#Run it
Wang.YD.4ch = do.call(jags.parallel,list(model.file = "code/LipidMM-JAGS-Multinorm-plus-d2H.R", 
                                          parameters.to.save = parameters, 
                                          data = dat, n.chains = 3, n.iter = n.iter, 
                                          n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1

#use rhat to check convergence
Wang.YD.4ch$BUGSoutput$summary[1:3,]
#traceplots
traplot(Wang.YD.4ch, parms = "FLMC")

save(Wang.YD.4ch, file = "out/Wang_YD_results.RData")

######2.1 Ka#####
#relative abundance among the three chains from raw concentration values
#in the order of n-C27, n-C29, n-C31, c-C33

RA.Ho2.1ka <- c(102.6, 199.8, 223.1, 162.3)/(102.6 + 199.8 + 223.1 + 162.3)

#d2H
d2H.Ho2.1ka <- c(-148.1, -146.7, -154.7, -155.7)
#use ice volume correcte values
d2H.Ho2.1ka.ivc <- c(-148.3, -146.9, -154.9, -155.8)

#analytical precision
d2H.sd.Ho2.1ka <- c(0.6, 0.5, 0.3, 0.5)

#d13C
d13C.Ho2.1ka <- c(-29.1, -30.2, -30.1, -27.8)

#analytical precision
d13C.sd.Ho2.1ka <- c(0.2, 0.2, 0.2, 0.2)

###correction for atm CO2 d13C, should set it as a parameter?
#-8.3 year 2010 (Graven et al 2017)
#-6.4 Holocene 2.1 Ka (Schmitt et al 2012)
Afr.d13C.mu.cor <- Afr.d13C.mu + 8.3 - 6.4 

##Data to pass to the model
#prior parameters in the first three lines
#model parameters in the fourth
#data in the fifth and sixth lines
dat = list(d13C.mu.est = Afr.d13C.mu.cor, d13C.omega.est = Afr.d13C.vcov,
           conc.mu.est = Afr.conc.mu, conc.omega.est = Afr.conc.vcov, 
           epsilon.app.mu.est = map_rec.epsilon.app.mu, epsilon.app.omega.est = map_rec.epsilon.app.vcov,
           I = I, N = N, K = K, 
           RA.mix = RA.Ho2.1ka, d13C.mix = d13C.Ho2.1ka, d13C.mea.sd = d13C.sd.Ho2.1ka,
           d2H.mix = d2H.Ho2.1ka.ivc, d2H.mea.sd = d2H.sd.Ho2.1ka)

#model parameters to save
parameters <- c("d13C.mix.m","RA.mix.m", "FLMC","f.sum.conc_n_i","exp.conc_k","d13C.k",
                "d2H.mix.m", "d2H.k", "d2H.MAP", "epsilon.app.k")

#Start time
t1 = proc.time()

set.seed(t1[3])
#Run it
Wang.Ho2.1ka.4ch = do.call(jags.parallel,list(model.file = "code/LipidMM-JAGS-Multinorm-plus-d2H.R", 
                                              parameters.to.save = parameters, 
                                              data = dat, n.chains = 3, n.iter = n.iter, 
                                              n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1

#use rhat to check convergence
Wang.Ho2.1ka.4ch$BUGSoutput$summary[1:3,]
#traceplots
traplot(Wang.Ho2.1ka.4ch, parms = "FLMC")

save(Wang.Ho2.1ka.4ch, file = "out/Wang_Ho2.1ka_results.RData")

######22.7 Ka, LGM#####
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
                                         data = dat, n.chains = 3, n.iter = n.iter, 
                                         n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1

#use rhat to check convergence
Wang.22.7ka.4ch$BUGSoutput$summary[1:3,]
#traceplots
traplot(Wang.22.7ka.4ch, parms = "FLMC")

save(Wang.22.7ka.4ch, file = "out/Wang_22.7ka_results.RData")

######15.8 Ka, LGM#####
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
                                             data = dat, n.chains = 3, n.iter = n.iter, 
                                             n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1

#use rhat to check convergence
Wang.15.8ka.4ch$BUGSoutput$summary[1:3,]
#traceplots
traplot(Wang.15.8ka.4ch, parms = "FLMC")

save(Wang.15.8ka.4ch, file = "out/Wang_15.8ka_results.RData")

######7.8 Ka, #####
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
                                             data = dat, n.chains = 3, n.iter = n.iter, 
                                             n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1

#use rhat to check convergence
Wang.7.8ka.4ch$BUGSoutput$summary[1:3,]
#traceplots
traplot(Wang.7.8ka.4ch, parms = "FLMC")

save(Wang.7.8ka.4ch, file = "out/Wang_7.8ka_results.RData")

####plots###
#calculate prior densities by using the helper function "pri.multi.norm.den"
map_rec.epsilon.app.prior.GR<-pri.multi.norm.den(-200,0,map_rec.epsilon.app.mu[1,],
                                      map_rec.epsilon.app.vcov[1:4,1:4])
map_rec.epsilon.app.prior.SV<-pri.multi.norm.den(-200,0,map_rec.epsilon.app.mu[2,],
                                      map_rec.epsilon.app.vcov[5:8,1:4])
map_rec.epsilon.app.prior.RF<-pri.multi.norm.den(-200,0,map_rec.epsilon.app.mu[3,],
                                      map_rec.epsilon.app.vcov[9:12,1:4])

par(mfrow=c(3,4)) #900*800
#C4 grasses 
hist(map_rec.GR$eC27.MAP_OIPC, xlim = c(-200,0),freq = F,breaks=11, ylim = c(0,0.08), 
     main = "C4 n-C27 epsilon")
lines(map_rec.epsilon.app.prior.GR$x,map_rec.epsilon.app.prior.GR$y[1,], col = "blue", lwd = 2)

hist(map_rec.GR$eC29.MAP_OIPC, xlim = c(-200,0),freq = F,breaks=11, ylim = c(0,0.08), 
     main = "C4 n-C29 epsilon")
lines(map_rec.epsilon.app.prior.GR$x,map_rec.epsilon.app.prior.GR$y[2,], col = "blue", lwd = 2)

hist(map_rec.GR$eC31.MAP_OIPC, xlim = c(-200,0),freq = F,breaks=11, ylim = c(0,0.08), 
     main = "C4 n-C31 epsilon")
lines(map_rec.epsilon.app.prior.GR$x,map_rec.epsilon.app.prior.GR$y[3,], col = "blue", lwd = 2)

hist(map_rec.GR$eC33.MAP_OIPC, xlim = c(-200,0),freq = F,breaks=11, ylim = c(0,0.08), 
     main = "C4 n-C33 epsilon")
lines(map_rec.epsilon.app.prior.GR$x,map_rec.epsilon.app.prior.GR$y[4,], col = "blue", lwd = 2)

#savanna C3 plants
hist(map_rec.SV$eC27.MAP_OIPC, xlim = c(-200,0),freq = F,breaks=11, ylim = c(0,0.04), 
     main = "SV n-C27 epsilon")
lines(map_rec.epsilon.app.prior.SV$x,map_rec.epsilon.app.prior.SV$y[1,], col = "blue", lwd = 2)

hist(map_rec.SV$eC29.MAP_OIPC, xlim = c(-200,0),freq = F,breaks=11, ylim = c(0,0.04), 
     main = "SV n-C29 epsilon")
lines(map_rec.epsilon.app.prior.SV$x,map_rec.epsilon.app.prior.SV$y[2,], col = "blue", lwd = 2)

hist(map_rec.SV$eC31.MAP_OIPC, xlim = c(-200,0),freq = F,breaks=11, ylim = c(0,0.04), 
     main = "SV n-C31 epsilon")
lines(map_rec.epsilon.app.prior.SV$x,map_rec.epsilon.app.prior.SV$y[3,], col = "blue", lwd = 2)

hist(map_rec.SV$eC33.MAP_OIPC, xlim = c(-200,0),freq = F,breaks=11, ylim = c(0,0.04), 
     main = "SV n-C33 epsilon")
lines(map_rec.epsilon.app.prior.SV$x,map_rec.epsilon.app.prior.SV$y[4,], col = "blue", lwd = 2)

#rainforest C3 plants
hist(map_rec.RF$eC27.MAP_OIPC, xlim = c(-200,0),freq = F,breaks=11, ylim = c(0,0.03), 
     main = "RF n-C27 epsilon")
lines(map_rec.epsilon.app.prior.RF$x,map_rec.epsilon.app.prior.RF$y[1,], col = "blue", lwd = 2)

hist(map_rec.RF$eC29.MAP_OIPC, xlim = c(-200,0),freq = F,breaks=11, ylim = c(0,0.03), 
     main = "RF n-C29 epsilon")
lines(map_rec.epsilon.app.prior.RF$x,map_rec.epsilon.app.prior.RF$y[2,], col = "blue", lwd = 2)

hist(map_rec.RF$eC31.MAP_OIPC, xlim = c(-200,0),freq = F,breaks=11, ylim = c(0,0.03), 
     main = "RF n-C31 epsilon")
lines(map_rec.epsilon.app.prior.RF$x,map_rec.epsilon.app.prior.RF$y[3,], col = "blue", lwd = 2)

hist(map_rec.RF$eC33.MAP_OIPC, xlim = c(-200,0),freq = F,breaks=11, ylim = c(0,0.03), 
     main = "RF n-C33 epsilon")
lines(map_rec.epsilon.app.prior.RF$x,map_rec.epsilon.app.prior.RF$y[4,], col = "blue", lwd = 2)

#posterior densities of mixing fractions

plot.col<-viridis(7)

par(mfrow=c(2,2)) #900*800
# plot(density(Wang.LGM.mix$BUGSoutput$sims.list$FLMC[,1],from=0,to=1), 
#      main = "LGM", xlab="FLMC", xlim=c(0,1), ylim =c(0,10),lwd=2,col = plot.col[6])
# lines(density(Wang.LGM.mix$BUGSoutput$sims.list$FLMC[,2],from=0,to=1),lwd=2,col = plot.col[4])
# lines(density(Wang.LGM.mix$BUGSoutput$sims.list$FLMC[,3],from=0,to=1),lwd=2,col = plot.col[2])
# legend(0.4,10,c("C4","Savanna C3","RF C3"),lwd=c(2,2,2),
#        col=plot.col[c(6,4,2)])

plot(density(Wang.LGM.4ch$BUGSoutput$sims.list$FLMC[,1],from=0,to=1), 
     main = "LGM", xlab="FLMC", xlim=c(0,1), ylim =c(0,10),lwd=2,col = plot.col[6])
lines(density(Wang.LGM.4ch$BUGSoutput$sims.list$FLMC[,2],from=0,to=1),lwd=2,col = plot.col[4])
lines(density(Wang.LGM.4ch$BUGSoutput$sims.list$FLMC[,3],from=0,to=1),lwd=2,col = plot.col[2])
legend(0.4,10,c("C4","Savanna C3","RF C3"),lwd=c(2,2,2),
       col=plot.col[c(6,4,2)])

plot(density(Wang.HS1.4ch$BUGSoutput$sims.list$FLMC[,1],from=0,to=1), 
     main = "HS1", xlab="FLMC", xlim=c(0,1), ylim =c(0,10),lwd=2,col = plot.col[6])
lines(density(Wang.HS1.4ch$BUGSoutput$sims.list$FLMC[,2],from=0,to=1),lwd=2,col = plot.col[4])
lines(density(Wang.HS1.4ch$BUGSoutput$sims.list$FLMC[,3],from=0,to=1),lwd=2,col = plot.col[2])
legend(0.4,10,c("C4","Savanna C3","RF C3"),lwd=c(2,2,2),
       col=plot.col[c(6,4,2)])

plot(density(Wang.YD.4ch$BUGSoutput$sims.list$FLMC[,1],from=0,to=1), 
     main = "YD", xlab="FLMC", xlim=c(0,1), ylim =c(0,10),lwd=2,col = plot.col[6])
lines(density(Wang.YD.4ch$BUGSoutput$sims.list$FLMC[,2],from=0,to=1),lwd=2,col = plot.col[4])
lines(density(Wang.YD.4ch$BUGSoutput$sims.list$FLMC[,3],from=0,to=1),lwd=2,col = plot.col[2])
legend(0.4,10,c("C4","Savanna C3","RF C3"),lwd=c(2,2,2),
       col=plot.col[c(6,4,2)])

plot(density(Wang.22.7ka.4ch$BUGSoutput$sims.list$FLMC[,1],from=0,to=1), 
     main = "22.7 Ka", xlab="FLMC", xlim=c(0,1), ylim =c(0,10),lwd=2,col = plot.col[6])
lines(density(Wang.22.7ka.4ch$BUGSoutput$sims.list$FLMC[,2],from=0,to=1),lwd=2,col = plot.col[4])
lines(density(Wang.22.7ka.4ch$BUGSoutput$sims.list$FLMC[,3],from=0,to=1),lwd=2,col = plot.col[2])
legend(0.4,10,c("C4","Savanna C3","RF C3"),lwd=c(2,2,2),
       col=plot.col[c(6,4,2)])


# plot(density(Wang.LGM.mix$BUGSoutput$sims.list$d2H.MAP), xlim = c(-60,-10), ylim = c(0,0.1),
#      main = "", xlab="d2H",lwd=2,col = plot.col[7])
par(mfrow=c(1,1))
plot(density(Wang.LGM.4ch$BUGSoutput$sims.list$d2H.MAP), xlim = c(-60,-10), ylim = c(0,0.12),
     main = "", xlab="d2H",lwd=2,col = plot.col[7])

lines(density(Wang.HS1.4ch$BUGSoutput$sims.list$d2H.MAP), 
     lwd=2,col = plot.col[5])

lines(density(Wang.YD.4ch$BUGSoutput$sims.list$d2H.MAP), 
      lwd=2,col = plot.col[3])

lines(density(Wang.Ho2.1ka.4ch$BUGSoutput$sims.list$d2H.MAP), 
      lwd=2,col = plot.col[1])

legend(-60,0.1,c("LGM","HS1","DY","Ho 2.1 Ka"),lwd=c(2,2,2,2),
       col=plot.col[c(7,5,3,1)])
# ####Try C29, C31, C33 only####
# ####make sure to run all the helper functions in "code/LipidMM-helper functions.R"
# map_rec<- read.csv("data/EA-5 data map_rec.csv")
# 
# map_rec.GR <- map_rec[which(map_rec$Source == "GR"),]
# map_rec.SV <- map_rec[which(map_rec$Source == "SV"),]
# map_rec.RF <- map_rec[which(map_rec$Source == "RF"),]
# 
# #the three chains used in this case study are Cn-27, Cn-29, and Cn-31
# ###compile prior parameters######
# 
# #epsilon alkane-MAP apparent fractionation factor
# ###means and vcovs
# 
# map_rec.GR.e.mean.est <- c(mean(map_rec.GR$eC29.MAP_OIPC, na.rm = T), 
#                            mean(map_rec.GR$eC31.MAP_OIPC, na.rm = T),
#                            mean(map_rec.GR$eC33.MAP_OIPC, na.rm = T))
# map_rec.GR.e.vcov <- var(data.frame(map_rec.GR$eC29.MAP_OIPC, 
#                                     map_rec.GR$eC31.MAP_OIPC,
#                                     map_rec.GR$eC33.MAP_OIPC), use = "complete.obs")
# 
# map_rec.SV.e.mean.est <- c(mean(map_rec.SV$eC29.MAP_OIPC, na.rm = T), 
#                            mean(map_rec.SV$eC31.MAP_OIPC, na.rm = T),
#                            mean(map_rec.SV$eC33.MAP_OIPC, na.rm = T))
# map_rec.SV.e.vcov <- var(data.frame(map_rec.SV$eC29.MAP_OIPC, 
#                                     map_rec.SV$eC31.MAP_OIPC,
#                                     map_rec.SV$eC33.MAP_OIPC), use = "complete.obs")
# 
# map_rec.RF.e.mean.est <- c(mean(map_rec.RF$eC29.MAP_OIPC, na.rm = T), 
#                            mean(map_rec.RF$eC31.MAP_OIPC, na.rm = T),
#                            mean(map_rec.RF$eC33.MAP_OIPC, na.rm = T))
# map_rec.RF.e.vcov <- var(data.frame(map_rec.RF$eC29.MAP_OIPC, 
#                                     map_rec.RF$eC31.MAP_OIPC,
#                                     map_rec.RF$eC33.MAP_OIPC), use = "complete.obs")
# 
# ###rows are sources i
# ###columns are chain n
# map_rec.epsilon.app.mu <- rbind(map_rec.GR.e.mean.est, 
#                                 map_rec.SV.e.mean.est, map_rec.RF.e.mean.est)
# 
# map_rec.epsilon.app.vcov <- rbind(map_rec.GR.e.vcov, map_rec.SV.e.vcov, map_rec.RF.e.vcov)
# 
# #the compiled prior parameters
# map_rec.epsilon.app.mu
# map_rec.epsilon.app.vcov
# 
# ####Adopt the d13C and n-alkane concentrations from CS2
# African <- read.csv("data/EA-3 data afr.csv")
# 
# #subset the data
# GR <- African[which(African$Source == "GR"),]
# SV <- African[which(African$Source == "SV"),]
# RF <- African[which(African$Source == "RF"),]
# 
# ###compile prior parameters######
# 
# ##d13C
# ###means and vcovs
# Afr.GR.d13C.mean.est <- c(mean(GR$d.n.C29, na.rm = T),
#                           mean(GR$d.n.C31, na.rm = T), mean(GR$d.n.C33, na.rm = T))
# Afr.GR.d13C.vcov <- var(data.frame(GR$d.n.C29, 
#                                    GR$d.n.C31, GR$d.n.C33), use = "complete.obs")
# 
# Afr.SV.d13C.mean.est <- c(mean(SV$d.n.C29, na.rm = T),
#                           mean(SV$d.n.C31, na.rm = T), mean(SV$d.n.C33, na.rm = T))
# Afr.SV.d13C.vcov <- var(data.frame(SV$d.n.C29, 
#                                    SV$d.n.C31, SV$d.n.C33), use = "complete.obs")
# 
# Afr.RF.d13C.mean.est <- c(mean(RF$d.n.C29, na.rm = T),
#                           mean(RF$d.n.C31, na.rm = T), mean(RF$d.n.C33, na.rm = T))
# Afr.RF.d13C.vcov <- var(data.frame(RF$d.n.C29, 
#                                    RF$d.n.C31, RF$d.n.C33), use = "complete.obs")
# 
# ###rows are sources i
# ###columns are chain n
# Afr.d13C.mu <- rbind(Afr.GR.d13C.mean.est, Afr.SV.d13C.mean.est, Afr.RF.d13C.mean.est)
# 
# Afr.d13C.vcov <- rbind(Afr.GR.d13C.vcov, Afr.SV.d13C.vcov, Afr.RF.d13C.vcov)
# 
# #the compiled prior parameters
# Afr.d13C.mu
# Afr.d13C.vcov
# 
# ###concentrations are natrual log transformed
# ###means and vcovs
# Afr.GR.conc.mean.est <- c(mean(log(GR$c.n.C29), na.rm = T), 
#                           mean(log(GR$c.n.C31), na.rm = T), mean(log(GR$c.n.C33), na.rm = T))
# Afr.GR.conc.vcov <- var(data.frame(log(GR$c.n.C29), 
#                                    log(GR$c.n.C31), log(GR$c.n.C33)), use = "complete.obs")
# 
# Afr.SV.conc.mean.est <- c(mean(log(SV$c.n.C29), na.rm = T),
#                           mean(log(SV$c.n.C31), na.rm = T), mean(log(SV$c.n.C33), na.rm = T))
# Afr.SV.conc.vcov <- var(data.frame(log(SV$c.n.C29), 
#                                    log(SV$c.n.C31), log(SV$c.n.C33)), use = "complete.obs")
# 
# Afr.RF.conc.mean.est <- c(mean(log(RF$c.n.C29), na.rm = T), 
#                           mean(log(RF$c.n.C31), na.rm = T), mean(log(RF$c.n.C33), na.rm = T))
# Afr.RF.conc.vcov <- var(data.frame(log(RF$c.n.C29), 
#                                    log(RF$c.n.C31), log(RF$c.n.C33)), use = "complete.obs")
# 
# ###rows are sources i
# ###columns are chain n
# Afr.conc.mu <- rbind(Afr.GR.conc.mean.est, Afr.SV.conc.mean.est, Afr.RF.conc.mean.est)
# 
# Afr.conc.vcov <- rbind(Afr.GR.conc.vcov, Afr.SV.conc.vcov, Afr.RF.conc.vcov)
# 
# #the compiled prior parameters
# Afr.conc.mu
# Afr.conc.vcov
# 
# #######model parameters#####
# #initialize parameters
# I <- 3  #number of sources
# N <- 3  #number of chains
# K <- 50 #number of grams of leaves to integrate per source
# 
# #common MCMC parameters
# n.iter = 8e5
# n.burnin = 2e5
# n.thin = floor(n.iter-n.burnin)/2500
# #average runtime is ~4 hours/sample
# 
# #MAP reconstruction based on Core GIK16160-3, Zambezi River mouth, Wang et al. 2013 
# ####first data point: LGM (20.1 ka), HS1 (16.1 ka), YD (12.1 Ka), Holocene (2.15 Ka)
# 
# 
# #relative abundance among the three chains from raw concentration values
# #in the order of n-C27, n-C29, n-C31
# RA.LGM <- c(430.3, 465, 313.3)/(430.3 + 465 + 313.3)
# 
# #d2H, 
# d2H.LGM <- c(-141.3, -154.3, -158.9)
# #use ice volume correcte values
# d2H.LGM.ivc <- c(-148.3, -161.2, -165.7)
# 
# #analytical precision
# d2H.sd.LGM <- c(0.6, 1, 0.9)
# 
# #d13C
# d13C.LGM <- c(-29.8, -27.8, -24.9)
# 
# #analytical precision
# d13C.sd.LGM <- c(0.1, 0.1, 0.2)
# 
# ###correction for atm CO2 d13C, should set it as a parameter?
# #-8.3 year 2010 (Graven et al 2017)
# #-6.5 LGM (Schmitt et al 2012)
# Afr.d13C.mu.cor <- Afr.d13C.mu + 8.3 - 6.5 
# 
# ##Data to pass to the model
# #prior parameters in the first three lines
# #model parameters in the fourth
# #data in the fifth and sixth lines
# dat = list(d13C.mu.est = Afr.d13C.mu.cor, d13C.omega.est = Afr.d13C.vcov,
#            conc.mu.est = Afr.conc.mu, conc.omega.est = Afr.conc.vcov, 
#            epsilon.app.mu.est = map_rec.epsilon.app.mu, epsilon.app.omega.est = map_rec.epsilon.app.vcov,
#            I = I, N = N, K = K, 
#            RA.mix = RA.LGM, d13C.mix = d13C.LGM, d13C.mea.sd = d13C.sd.LGM,
#            d2H.mix = d2H.LGM.ivc, d2H.mea.sd = d2H.sd.LGM)
# 
# #model parameters to save
# parameters <- c("d13C.mix.m","RA.mix.m", "FLMC","f.sum.conc_n_i","exp.conc_k","d13C.k",
#                 "d2H.mix.m", "d2H.k", "d2H.MAP", "epsilon.app.k")
# 
# #Start time
# t1 = proc.time()
# 
# set.seed(t1[3])
# #Run it
# Wang.LGM.mix = do.call(jags.parallel,list(model.file = "code/LipidMM-JAGS-Multinorm-plus-d2H.R", 
#                                           parameters.to.save = parameters, 
#                                           data = dat, n.chains = 3, n.iter = n.iter, 
#                                           n.burnin = n.burnin, n.thin = n.thin))
# 
# #Time taken
# proc.time() - t1 #~3 hours
# 
# #use rhat to check convergence
# Wang.LGM.mix$BUGSoutput$summary[1:3,]
# #traceplots
# traplot(Wang.LGM.mix, parms = "FLMC")