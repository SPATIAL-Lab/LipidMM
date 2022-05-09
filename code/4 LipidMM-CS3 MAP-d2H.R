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
                          mean(map_rec.GR$eC31.MAP_OIPC, na.rm = T))
map_rec.GR.e.vcov <- var(data.frame(map_rec.GR$eC27.MAP_OIPC, 
                                    map_rec.GR$eC29.MAP_OIPC, 
                                    map_rec.GR$eC31.MAP_OIPC), use = "complete.obs")

map_rec.SV.e.mean.est <- c(mean(map_rec.SV$eC27.MAP_OIPC, na.rm = T),
                           mean(map_rec.SV$eC29.MAP_OIPC, na.rm = T), 
                           mean(map_rec.SV$eC31.MAP_OIPC, na.rm = T))
map_rec.SV.e.vcov <- var(data.frame(map_rec.SV$eC27.MAP_OIPC, 
                                    map_rec.SV$eC29.MAP_OIPC, 
                                    map_rec.SV$eC31.MAP_OIPC), use = "complete.obs")

map_rec.RF.e.mean.est <- c(mean(map_rec.RF$eC27.MAP_OIPC, na.rm = T),
                           mean(map_rec.RF$eC29.MAP_OIPC, na.rm = T), 
                           mean(map_rec.RF$eC31.MAP_OIPC, na.rm = T))
map_rec.RF.e.vcov <- var(data.frame(map_rec.RF$eC27.MAP_OIPC, 
                                    map_rec.RF$eC29.MAP_OIPC, 
                                    map_rec.RF$eC31.MAP_OIPC), use = "complete.obs")

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
Afr.GR.d13C.mean.est <- c(mean(GR$d.n.C27, na.rm = T),
                          mean(GR$d.n.C29, na.rm = T), mean(GR$d.n.C31, na.rm = T))
Afr.GR.d13C.vcov <- var(data.frame(GR$d.n.C27, GR$d.n.C29, GR$d.n.C31), use = "complete.obs")

Afr.SV.d13C.mean.est <- c(mean(SV$d.n.C27, na.rm = T),
                          mean(SV$d.n.C29, na.rm = T), mean(SV$d.n.C31, na.rm = T))
Afr.SV.d13C.vcov <- var(data.frame(SV$d.n.C27, SV$d.n.C29, SV$d.n.C31), use = "complete.obs")

Afr.RF.d13C.mean.est <- c(mean(RF$d.n.C27, na.rm = T),
                          mean(RF$d.n.C29, na.rm = T), mean(RF$d.n.C31, na.rm = T))
Afr.RF.d13C.vcov <- var(data.frame(RF$d.n.C27, RF$d.n.C29, RF$d.n.C31), use = "complete.obs")

###rows are sources i
###columns are chain n
Afr.d13C.mu <- rbind(Afr.GR.d13C.mean.est, Afr.SV.d13C.mean.est, Afr.RF.d13C.mean.est)

Afr.d13C.vcov <- rbind(Afr.GR.d13C.vcov, Afr.SV.d13C.vcov, Afr.RF.d13C.vcov)

#the compiled prior parameters
Afr.d13C.mu
Afr.d13C.vcov

###concentrations are natrual log transformed
###means and vcovs
Afr.GR.conc.mean.est <- c(mean(log(GR$c.n.C27), na.rm = T),
                          mean(log(GR$c.n.C29), na.rm = T), mean(log(GR$c.n.C31), na.rm = T))
Afr.GR.conc.vcov <- var(data.frame(log(GR$c.n.C27), log(GR$c.n.C29), log(GR$c.n.C31)), use = "complete.obs")

Afr.SV.conc.mean.est <- c(mean(log(SV$c.n.C27), na.rm = T),
                          mean(log(SV$c.n.C29), na.rm = T), mean(log(SV$c.n.C31), na.rm = T))
Afr.SV.conc.vcov <- var(data.frame(log(SV$c.n.C27), log(SV$c.n.C29), log(SV$c.n.C31)), use = "complete.obs")

Afr.RF.conc.mean.est <- c(mean(log(RF$c.n.C27), na.rm = T),
                          mean(log(RF$c.n.C29), na.rm = T), mean(log(RF$c.n.C31), na.rm = T))
Afr.RF.conc.vcov <- var(data.frame(log(RF$c.n.C27), log(RF$c.n.C29), log(RF$c.n.C31)), use = "complete.obs")

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
N <- 3  #number of chains
K <- 50 #number of grams of leaves to integrate per source

#common MCMC parameters
n.iter = 8e5
n.burnin = 2e5
n.thin = floor(n.iter-n.burnin)/2500
#average runtime is ~4 hours/sample

#model parameters to save
parameters <- c("d13C.mix.m","RA.mix.m", "FLMC","f.sum.conc_n_i","exp.conc_k","d13C.k",
                "d2H.mix.m", "d2H.k", "d2H.MAP")

#MAP reconstruction based on Core GIK16160-3, Zambezi River mouth, Wang et al. 2013 
####first data point: LGM (20.1 ka), HS1 (16.1 ka), YD (12.1 Ka), Holocene (2.15 Ka)


#relative abundance among the tree chains from raw concentration values
#in the order of n-C27, n-C29, n-C31
RA.LGM <- c(0.317, 0.384, 0.299)

#d2H
d2H.LGM <- c(-146.2, -141.3, -154.3)

#analytical precision
d2H.sd.LGM <- c(0.3, 0.3, 0.3)

#d13C
d13C.LGM <- c(-28.6, -29.7, -24.8)

#analytical precision
d13C.sd.LGM <- c(0.3, 0.3, 0.3)

###correction for atm CO2 d13C
Afr.d13C.mu.cor <- Afr.d13C.mu + x

##Data to pass to the model
#prior parameters in the first three lines
#model parameters in the third
#data in the fourth
dat = list(d13C.mu.est = Afr.d13C.mu.cor, d13C.omega.est = Afr.d13C.vcov,
           conc.mu.est = Afr.conc.mu, conc.omega.est = Afr.conc.vcov, 
           epsilon.app.mu.est = map_rec.epsilon.app.mu, epsilon.app.omega = map_rec.epsilon.app.vcov,
           I = I, N = N, K = K, 
           RA.mix = RA.rhum.l, d13C.mix = d13C.rhum.l, d13C.mea.sd = d13C.sd.rhum.l)


#what is ice volume correction?


