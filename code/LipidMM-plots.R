library(coda)
library(lattice)
library(rjags)
library(R2jags)
library(MASS)
library(scales)
library(viridisLite)
library(EnvStats)
library(bayestestR)

plot.col<-viridis(7)
###CS1 QTP mixing####

###CS1 QTP d13C prior distribution#####
#calculate prior densities by using the helper function "pri.multi.norm.den"
QTP.d13C.prior.ter<-pri.multi.norm.den(-45,-5,QTP.d13C.mu[1,],QTP.d13C.vcov[1:3,1:3])
QTP.d13C.prior.mac<-pri.multi.norm.den(-45,-5,QTP.d13C.mu[2,],QTP.d13C.vcov[4:6,1:3])
QTP.d13C.prior.alg<-pri.multi.norm.den(-45,-5,QTP.d13C.mu[3,],QTP.d13C.vcov[7:9,1:3])


#histograms are empirical distributions, blue lines are parameterized distributions
par(mfrow=c(3,3)) #900*800
#terrestrial plants
hist(Ter$d.n.C27, xlim = c(-40,-10),freq = F,breaks=11, ylim = c(0,0.8), main = "n-C27 d13C")
lines(QTP.d13C.prior.ter$x,QTP.d13C.prior.ter$y1, col = "blue", lwd = 2)

hist(Ter$d.n.C29, xlim = c(-40,-10),freq = F,breaks=11, ylim = c(0,0.8), main = "n-C29 d13C")
lines(QTP.d13C.prior.ter$x,QTP.d13C.prior.ter$y2, col = "blue", lwd = 2)

hist(Ter$d.n.C31, xlim = c(-40,-10),freq = F,breaks=11, ylim = c(0,0.8), main = "n-C31 d13C")
lines(QTP.d13C.prior.ter$x,QTP.d13C.prior.ter$y3, col = "blue", lwd = 2)

#aquatic macrophyte
hist(Mac$d.n.C27, xlim = c(-40,-10),freq = F,breaks=11, ylim = c(0,0.15), main = "n-C27 d13C")
lines(QTP.d13C.prior.mac$x,QTP.d13C.prior.mac$y1, col = "blue", lwd = 2)

hist(Mac$d.n.C29, xlim = c(-40,-10),freq = F,breaks=11, ylim = c(0,0.15), main = "n-C29 d13C")
lines(QTP.d13C.prior.mac$x,QTP.d13C.prior.mac$y2, col = "blue", lwd = 2)

hist(Mac$d.n.C31, xlim = c(-40,-10),freq = F,breaks=11, ylim = c(0,0.15), main = "n-C31 d13C")
lines(QTP.d13C.prior.mac$x,QTP.d13C.prior.mac$y3, col = "blue", lwd = 2)

#algae
hist(Alg$d.n.C27, xlim = c(-40,-10),freq = F,breaks=11, ylim = c(0,0.6), main = "n-C27 d13C")
lines(QTP.d13C.prior.alg$x,QTP.d13C.prior.alg$y1, col = "blue", lwd = 2)

hist(Alg$d.n.C29, xlim = c(-40,-10),freq = F,breaks=11, ylim = c(0,0.6), main = "n-C29 d13C")
lines(QTP.d13C.prior.alg$x,QTP.d13C.prior.alg$y2, col = "blue", lwd = 2)

hist(Alg$d.n.C31, xlim = c(-40,-10),freq = F,breaks=11, ylim = c(0,0.6), main = "n-C31 d13C")
lines(QTP.d13C.prior.alg$x,QTP.d13C.prior.alg$y3, col = "blue", lwd = 2)

####CS1 QTP concentration prior distribution###
#empirical distribution in histograms, parameterized prior distribution in density curves
#calculate prior densities by using the helper function "pri.multi.norm.den"
QTP.prod.prior.ter<-pri.multi.norm.den(-5,10,QTP.prod.mu[1,],QTP.prod.vcov[1:3,1:3])
QTP.prod.prior.mac<-pri.multi.norm.den(-5,10,QTP.prod.mu[2,],QTP.prod.vcov[4:6,1:3])
QTP.prod.prior.alg<-pri.multi.norm.den(-5,10,QTP.prod.mu[3,],QTP.prod.vcov[7:9,1:3])

#histograms are empirical distributions, blue lines are parameterized distributions
par(mfrow=c(3,3)) #900*800
#terrestrial plants

hist(log(Ter$c.n.C27), xlim = c(-5,10),freq = F,breaks=11, ylim = c(0,0.5), main = "n-C27 concentration")
lines(QTP.prod.prior.ter$x,QTP.prod.prior.ter$y1, col = "blue", lwd = 2)

hist(log(Ter$c.n.C29), xlim = c(-5,10),freq = F,breaks=11, ylim = c(0,0.5), main = "n-C29 concentration")
lines(QTP.prod.prior.ter$x,QTP.prod.prior.ter$y2, col = "blue", lwd = 2)

hist(log(Ter$c.n.C31), xlim = c(-5,10),freq = F,breaks=11, ylim = c(0,0.5), main = "n-C31 concentration")
lines(QTP.prod.prior.ter$x,QTP.prod.prior.ter$y3, col = "blue", lwd = 2)

hist(log(Mac$c.n.C27), xlim = c(-5,10),freq = F,breaks=11, ylim = c(0,0.8), main = "n-C27 concentration")
lines(QTP.prod.prior.mac$x,QTP.prod.prior.mac$y1, col = "blue", lwd = 2)

hist(log(Mac$c.n.C29), xlim = c(-5,10),freq = F,breaks=11, ylim = c(0,0.8), main = "n-C29 concentration")
lines(QTP.prod.prior.mac$x,QTP.prod.prior.mac$y2, col = "blue", lwd = 2)

hist(log(Mac$c.n.C31), xlim = c(-5,10),freq = F,breaks=11, ylim = c(0,0.8), main = "n-C31 concentration")
lines(QTP.prod.prior.mac$x,QTP.prod.prior.mac$y3, col = "blue", lwd = 2)

hist(log(Alg$c.n.C27), xlim = c(-5,10),freq = F,breaks=11, ylim = c(0,2.5), main = "n-C27 concentration")
lines(QTP.prod.prior.alg$x,QTP.prod.prior.alg$y1, col = "blue", lwd = 2)

hist(log(Alg$c.n.C29), xlim = c(-5,10),freq = F,breaks=11, ylim = c(0,2.5), main = "n-C29 concentration")
lines(QTP.prod.prior.alg$x,QTP.prod.prior.alg$y2, col = "blue", lwd = 2)

hist(log(Alg$c.n.C31), xlim = c(-5,10),freq = F,breaks=11, ylim = c(0,2.5), 
     main = "n-C31 concentration", axes = F)
lines(QTP.prod.prior.alg$x,QTP.prod.prior.alg$y3, col = "blue", lwd = 2)
#use log scale axes
axis(2,c(0,2.5))
axis(1,log(c(0.01,0.1,1,10,100,1000,10000)))

###CS1 QTP results: mixing ratios####
#run the corresponding code in the CS1 file
load("out/QHS13_5S_results.RData")
load("out/QHS13_7S_results.RData")
load("out/QHS13_9S_results.RData")

par(mfrow=c(3,3)) #800*800
###barplots for RA of data
barplot(RA.QHS13_5S,ylim = c(0,1),names.arg=c("C27","C29","C31"),space=0.5,xlim = c(0,5),axes=F,
        main="QHS13_5S, Lake sediment, Paq = 0.34")
axis(2,c(0,0.2,0.4,0.6))
text(c(1,2.5,4),RA.QHS13_5S,labels=c("0.217","0.354","0.429"),pos=1)

barplot(RA.QHS13_7S,ylim = c(0,1),names.arg=c("C27","C29","C31"),space=0.5,xlim = c(0,5),axes=F,
        main="QHS13_5S, Lake sediment, Paq = 0.22")
axis(2,c(0,0.2,0.4,0.6))
text(c(1,2.5,4),RA.QHS13_7S,labels=c("0.201","0.396","0.403"),pos=1)

barplot(RA.QHS13_9S,ylim = c(0,1),names.arg=c("C27","C29","C31"),space=0.5,xlim = c(0,5),axes=F,
        main="QHS13_9S, Lake sediment, Paq = 0.19")
axis(2,c(0,0.2,0.4,0.6))
text(c(1,2.5,4),RA.QHS13_9S,labels=c("0.176","0.382","0.442"),pos=1)

###line plots for d13C of data
barplot(RA.QHS13_5S,ylim = c(-44,-26),space=0.5,xlim = c(0,5),axes=F)
axis(4,c(-35,-32,-29,-26))
lines(c(1,2.5,4),d13C.QHS13_5S,lwd = 1.5, lty = 2)
arrows(c(1,2.5,4), d13C.QHS13_5S+d13C.sd.Liu, c(1,2.5,4), d13C.QHS13_5S-d13C.sd.Liu, 
       code=3, angle=90, length=0.05)
points(c(1,2.5,4),d13C.QHS13_5S,cex = 1.5)

barplot(RA.QHS13_7S,ylim = c(-44,-26),space=0.5,xlim = c(0,5),axes=F)
axis(4,c(-35,-32,-29,-26))
lines(c(1,2.5,4),d13C.QHS13_7S,lwd = 1.5, lty = 2)
arrows(c(1,2.5,4), d13C.QHS13_7S+d13C.sd.Liu, c(1,2.5,4), d13C.QHS13_7S-d13C.sd.Liu, 
       code=3, angle=90, length=0.05)
points(c(1,2.5,4),d13C.QHS13_7S,cex = 1.5)

barplot(RA.QHS13_9S,ylim = c(-44,-26),space=0.5,xlim = c(0,5),axes=F)
axis(4,c(-35,-32,-29,-26))
lines(c(1,2.5,4),d13C.QHS13_9S,lwd = 1.5, lty = 2)
arrows(c(1,2.5,4), d13C.QHS13_9S+d13C.sd.Liu, c(1,2.5,4), d13C.QHS13_9S-d13C.sd.Liu, 
       code=3, angle=90, length=0.05)
points(c(1,2.5,4),d13C.QHS13_9S,cex = 1.5)

#density plots for posterior distribution of mixing ratios
plot(density(QHS13_5S.mix$BUGSoutput$sims.list$f[,1], kernel = "gaussian", from = 0, to = 1),
     xlab="fraction", xlim = c(0, 1), col = plot.col[6], ylim = c(0, 8),lwd = 2, 
     main="QHS13_5S Paq = 0.34")
lines(density(QHS13_5S.mix$BUGSoutput$sims.list$f[,2], kernel = "gaussian", from = 0, to = 1),
      col = plot.col[4], lwd = 2)
lines(density(QHS13_5S.mix$BUGSoutput$sims.list$f[,3], kernel = "gaussian", from = 0, to = 1),
      col = plot.col[2], lwd = 2)
legend(0.6, 8, c("Terrestrial","Macrophyte","Algae"),lwd = c(2, 2, 2),
       col = plot.col[c(2, 4, 6)])

plot(density(QHS13_7S.mix$BUGSoutput$sims.list$f[,1], kernel = "gaussian", from = 0, to = 1),
     xlab="fraction", xlim = c(0, 1), col = plot.col[6], ylim = c(0, 8),lwd = 2, 
     main="QHS13_7S Paq = 0.22")
lines(density(QHS13_7S.mix$BUGSoutput$sims.list$f[,2], kernel = "gaussian", from = 0, to = 1),
      col = plot.col[4], lwd = 2)
lines(density(QHS13_7S.mix$BUGSoutput$sims.list$f[,3], kernel = "gaussian", from = 0, to = 1),
      col = plot.col[2], lwd = 2)
legend(0.6, 8, c("Terrestrial","Macrophyte","Algae"),lwd = c(2, 2, 2),
       col = plot.col[c(2, 4, 6)])

plot(density(QHS13_9S.mix$BUGSoutput$sims.list$f[,1], kernel = "gaussian", from = 0, to = 1),
     xlab="fraction", xlim = c(0, 1), col = plot.col[6], ylim = c(0, 8),lwd = 2, 
     main="QHS13_9S Paq = 0.19")
lines(density(QHS13_9S.mix$BUGSoutput$sims.list$f[,2], kernel = "gaussian", from = 0, to = 1),
      col = plot.col[4], lwd = 2)
lines(density(QHS13_9S.mix$BUGSoutput$sims.list$f[,3], kernel = "gaussian", from = 0, to = 1),
      col = plot.col[2], lwd = 2)
legend(0.6, 8, c("Terrestrial","Macrophyte","Algae"),lwd = c(2, 2, 2),
       col = plot.col[c(2, 4, 6)])

######CS1 QTP Bivariate density plots####
QHS13_5S.countours<-contour.fn(QHS13_5S.mix$BUGSoutput$sims.list$f)
QHS13_7S.countours<-contour.fn(QHS13_7S.mix$BUGSoutput$sims.list$f)
QHS13_9S.countours<-contour.fn(QHS13_9S.mix$BUGSoutput$sims.list$f)

#plot the data by column, each column is a sample
#from left to right: QHS13_5S, QHS13_7S, QHS13_9S
par(mfrow=c(3,3)) #700*770
image(QHS13_5S.countours[[4]],col=viridis(64),xlab="f terrestrial",ylab="f macrophyte",
      main="QHS13_5S")
contour(QHS13_5S.countours[[4]], lwd = 0.6, add = TRUE, labcex = 0.5)

image(QHS13_7S.countours[[4]],col=viridis(64),xlab="f terrestrial",ylab="f macrophyte",
      main="QHS13_7S")
contour(QHS13_7S.countours[[4]], lwd = 0.6, add = TRUE, labcex = 0.5)

image(QHS13_9S.countours[[4]],col=viridis(64),xlab="f terrestrial",ylab="f macrophyte",
      main="QHS13_9S")
contour(QHS13_9S.countours[[4]], lwd = 0.6, add = TRUE, labcex = 0.5)

image(QHS13_5S.countours[[5]],col=viridis(64),xlab="f terrestrial",ylab="f algae")
contour(QHS13_5S.countours[[5]], lwd = 0.6, add = TRUE, labcex = 0.5)

image(QHS13_7S.countours[[5]],col=viridis(64),xlab="f terrestrial",ylab="f algae")
contour(QHS13_7S.countours[[5]], lwd = 0.6, add = TRUE, labcex = 0.5)

image(QHS13_9S.countours[[5]],col=viridis(64),xlab="f terrestrial",ylab="f algae")
contour(QHS13_9S.countours[[5]], lwd = 0.6, add = TRUE, labcex = 0.5)

image(QHS13_5S.countours[[6]],col=viridis(64),xlab="f macrophyte",ylab="f algae")
contour(QHS13_5S.countours[[6]], lwd = 0.6, add = TRUE, labcex = 0.5)

image(QHS13_7S.countours[[6]],col=viridis(64),xlab="f macrophyte",ylab="f algae")
contour(QHS13_7S.countours[[6]], lwd = 0.6, add = TRUE, labcex = 0.5)

image(QHS13_9S.countours[[6]],col=viridis(64),xlab="f macrophyte",ylab="f algae")
contour(QHS13_9S.countours[[6]], lwd = 0.6, add = TRUE, labcex = 0.5)

###CS1 QTP chain specific mixing ratios####
#calculate prior densities by using the helper function "pri.multi.norm.den"
QHS13_5S.c.spec.ratio<-chain.spec.ratio(QHS13_5S.mix$BUGSoutput$sims.list$f.sum.prod_n_i)
QHS13_7S.c.spec.ratio<-chain.spec.ratio(QHS13_7S.mix$BUGSoutput$sims.list$f.sum.prod_n_i)
QHS13_9S.c.spec.ratio<-chain.spec.ratio(QHS13_9S.mix$BUGSoutput$sims.list$f.sum.prod_n_i)

#plot the data by column, each column is a sample
#from left to right: QHS13_5S, QHS13_7S, QHS13_9S
#the rows are the chains: 27, 29, 31
par(mfrow=c(3,3))
plot(density(QHS13_5S.c.spec.ratio[,1,1],from=0,to=1),main="n-alkane contribution C27",
     xlim=c(0,1),ylim=c(0,30),lwd=2,col=plot.col[6], xlab="fraction")
lines(density(QHS13_5S.c.spec.ratio[,2,1],from=0,to=1),lwd=2,col=plot.col[4])
lines(density(QHS13_5S.c.spec.ratio[,3,1],from=0,to=1),lwd=2,col=plot.col[2])

plot(density(QHS13_7S.c.spec.ratio[,1,1],from=0,to=1),main="n-alkane contribution C27",
     xlim=c(0,1),ylim=c(0,30),lwd=2,col=plot.col[6], xlab="fraction")
lines(density(QHS13_7S.c.spec.ratio[,2,1],from=0,to=1),lwd=2,col=plot.col[4])
lines(density(QHS13_7S.c.spec.ratio[,3,1],from=0,to=1),lwd=2,col=plot.col[2])

plot(density(QHS13_9S.c.spec.ratio[,1,1],from=0,to=1),main="n-alkane contribution C27",
     xlim=c(0,1),ylim=c(0,30),lwd=2,col=plot.col[6], xlab="fraction")
lines(density(QHS13_9S.c.spec.ratio[,2,1],from=0,to=1),lwd=2,col=plot.col[4])
lines(density(QHS13_9S.c.spec.ratio[,3,1],from=0,to=1),lwd=2,col=plot.col[2])

plot(density(QHS13_5S.c.spec.ratio[,1,2],from=0,to=1),main="n-alkane contribution C29",
     xlim=c(0,1),ylim=c(0,30),lwd=2,col=plot.col[6], xlab="fraction")
lines(density(QHS13_5S.c.spec.ratio[,2,2],from=0,to=1),lwd=2,col=plot.col[4])
lines(density(QHS13_5S.c.spec.ratio[,3,2],from=0,to=1),lwd=2,col=plot.col[2])

plot(density(QHS13_7S.c.spec.ratio[,1,2],from=0,to=1),main="n-alkane contribution C29",
     xlim=c(0,1),ylim=c(0,30),lwd=2,col=plot.col[6], xlab="fraction")
lines(density(QHS13_7S.c.spec.ratio[,2,2],from=0,to=1),lwd=2,col=plot.col[4])
lines(density(QHS13_7S.c.spec.ratio[,3,2],from=0,to=1),lwd=2,col=plot.col[2])

plot(density(QHS13_9S.c.spec.ratio[,1,2],from=0,to=1),main="n-alkane contribution C29",
     xlim=c(0,1),ylim=c(0,30),lwd=2,col=plot.col[6], xlab="fraction")
lines(density(QHS13_9S.c.spec.ratio[,2,2],from=0,to=1),lwd=2,col=plot.col[4])
lines(density(QHS13_9S.c.spec.ratio[,3,2],from=0,to=1),lwd=2,col=plot.col[2])

plot(density(QHS13_5S.c.spec.ratio[,1,3],from=0,to=1),main="n-alkane contribution C31",
     xlim=c(0,1),ylim=c(0,30),lwd=2,col=plot.col[6], xlab="fraction")
lines(density(QHS13_5S.c.spec.ratio[,2,3],from=0,to=1),lwd=2,col=plot.col[4])
lines(density(QHS13_5S.c.spec.ratio[,3,3],from=0,to=1),lwd=2,col=plot.col[2])

plot(density(QHS13_7S.c.spec.ratio[,1,3],from=0,to=1),main="n-alkane contribution C31",
     xlim=c(0,1),ylim=c(0,30),lwd=2,col=plot.col[6], xlab="fraction")
lines(density(QHS13_7S.c.spec.ratio[,2,3],from=0,to=1),lwd=2,col=plot.col[4])
lines(density(QHS13_7S.c.spec.ratio[,3,3],from=0,to=1),lwd=2,col=plot.col[2])

plot(density(QHS13_9S.c.spec.ratio[,1,3],from=0,to=1),main="n-alkane contribution C31",
     xlim=c(0,1),ylim=c(0,30),lwd=2,col=plot.col[6], xlab="fraction")
lines(density(QHS13_9S.c.spec.ratio[,2,3],from=0,to=1),lwd=2,col=plot.col[4])
lines(density(QHS13_9S.c.spec.ratio[,3,3],from=0,to=1),lwd=2,col=plot.col[2])
legend(0.3,30,c("Terrestrial","Macrophyte","Algae"),lwd=c(2,2,2),
       col=c(plot.col[c(6,4,2)]))

######CS2 Western Africa transect####

###CS2 n-alkane d13C prior distributions####
#calculate prior densities by using the helper function "pri.multi.norm.den"
Afr.d13C.prior.GR<-pri.multi.norm.den(-50,-15,Afr.d13C.mu[1,],Afr.d13C.vcov[1:3,1:3])
Afr.d13C.prior.SV<-pri.multi.norm.den(-50,-15,Afr.d13C.mu[2,],Afr.d13C.vcov[4:6,1:3])
Afr.d13C.prior.RF<-pri.multi.norm.den(-50,-15,Afr.d13C.mu[3,],Afr.d13C.vcov[7:9,1:3])

#histograms are empirical distributions, blue lines are parameterized distributions
par(mfrow=c(3,3)) #900*800
#C4 grasses 
hist(GR$d.n.C29, xlim = c(-50,-15),freq = F,breaks=11, ylim = c(0,0.3), main = "C4 n-C29 d13C")
lines(Afr.d13C.prior.GR$x,Afr.d13C.prior.GR$y1, col = "blue", lwd = 2)

hist(GR$d.n.C31, xlim = c(-50,-15),freq = F,breaks=11, ylim = c(0,0.3), main = "C4 n-C31 d13C")
lines(Afr.d13C.prior.GR$x,Afr.d13C.prior.GR$y2, col = "blue", lwd = 2)

hist(GR$d.n.C33, xlim = c(-50,-15),freq = F,breaks=11, ylim = c(0,0.3), main = "C4 n-C33 d13C")
lines(Afr.d13C.prior.GR$x,Afr.d13C.prior.GR$y3, col = "blue", lwd = 2)

#savanna C3 plants
hist(SV$d.n.C29, xlim = c(-50,-15),freq = F,breaks=11, ylim = c(0,0.3), main = "SV n-C29 d13C")
lines(Afr.d13C.prior.SV$x,Afr.d13C.prior.SV$y1, col = "blue", lwd = 2)

hist(SV$d.n.C31, xlim = c(-50,-15),freq = F,breaks=11, ylim = c(0,0.3), main = "SV n-C31 d13C")
lines(Afr.d13C.prior.SV$x,Afr.d13C.prior.SV$y2, col = "blue", lwd = 2)

hist(SV$d.n.C33, xlim = c(-50,-15),freq = F,breaks=11, ylim = c(0,0.3), main = "SV n-C33 d13C")
lines(Afr.d13C.prior.SV$x,Afr.d13C.prior.SV$y3, col = "blue", lwd = 2)

#Rainforest plants
hist(RF$d.n.C29, xlim = c(-50,-15),freq = F,breaks=11, ylim = c(0,0.3), main = "RF n-C29 d13C")
lines(Afr.d13C.prior.RF$x,Afr.d13C.prior.RF$y1, col = "blue", lwd = 2)

hist(RF$d.n.C31, xlim = c(-50,-15),freq = F,breaks=11, ylim = c(0,0.3), main = "RF n-C31 d13C")
lines(Afr.d13C.prior.RF$x,Afr.d13C.prior.RF$y2, col = "blue", lwd = 2)

hist(RF$d.n.C33, xlim = c(-50,-15),freq = F,breaks=11, ylim = c(0,0.3), main = "RF n-C33 d13C")
lines(Afr.d13C.prior.RF$x,Afr.d13C.prior.RF$y3, col = "blue", lwd = 2)


###CS2 n-alkane concentration prior distributions####
Afr.prod.prior.GR<-pri.multi.norm.den(-4,13,Afr.prod.mu[1,],Afr.prod.vcov[1:3,1:3])
Afr.prod.prior.SV<-pri.multi.norm.den(-4,13,Afr.prod.mu[2,],Afr.prod.vcov[4:6,1:3])
Afr.prod.prior.RF<-pri.multi.norm.den(-4,13,Afr.prod.mu[3,],Afr.prod.vcov[7:9,1:3])

par(mfrow=c(3,3))
hist(log(GR$c.n.C29), xlim = c(-4,13),freq = F,breaks=11, ylim = c(0,0.5), main = "n-C29 concentration")
lines(Afr.prod.prior.GR$x,Afr.prod.prior.GR$y1, col = "blue", lwd = 2)

hist(log(GR$c.n.C31), xlim = c(-4,13),freq = F,breaks=11, ylim = c(0,0.5), main = "n-C31 concentration")
lines(Afr.prod.prior.GR$x,Afr.prod.prior.GR$y2, col = "blue", lwd = 2)

hist(log(GR$c.n.C33), xlim = c(-4,13),freq = F,breaks=11, ylim = c(0,0.5), main = "n-C33 concentration")
lines(Afr.prod.prior.GR$x,Afr.prod.prior.GR$y3, col = "blue", lwd = 2)

hist(log(SV$c.n.C29), xlim = c(-4,13),freq = F,breaks=11, ylim = c(0,0.5), main = "n-C29 concentration")
lines(Afr.prod.prior.SV$x,Afr.prod.prior.SV$y1, col = "blue", lwd = 2)

hist(log(SV$c.n.C31), xlim = c(-4,13),freq = F,breaks=11, ylim = c(0,0.5), main = "n-C31 concentration")
lines(Afr.prod.prior.SV$x,Afr.prod.prior.SV$y2, col = "blue", lwd = 2)

hist(log(SV$c.n.C33), xlim = c(-4,13),freq = F,breaks=11, ylim = c(0,0.5), main = "n-C33 concentration")
lines(Afr.prod.prior.SV$x,Afr.prod.prior.SV$y3, col = "blue", lwd = 2)

hist(log(RF$c.n.C29), xlim = c(-4,13),freq = F,breaks=11, ylim = c(0,0.5), main = "n-C29 concentration")
lines(Afr.prod.prior.RF$x,Afr.prod.prior.RF$y1, col = "blue", lwd = 2)

hist(log(RF$c.n.C31), xlim = c(-4,13),freq = F,breaks=11, ylim = c(0,0.5), main = "n-C31 concentration")
lines(Afr.prod.prior.RF$x,Afr.prod.prior.RF$y2, col = "blue", lwd = 2)

hist(log(RF$c.n.C33), xlim = c(-4,13),freq = F,breaks=11, ylim = c(0,0.5), 
     main = "n-C33 concentration", axes=F)
lines(Afr.prod.prior.RF$x,Afr.prod.prior.RF$y3, col = "blue", lwd = 2)
#use log scale axes
axis(2,c(0,0.5))
axis(1,log(c(0.01,0.1,1,10,100,1e3,1e4,1e5)))

###CS2 mixing results####
#run the corresponding code in the CS2 file
load("out/rhum_l_results.RData")
load("out/asso_l_results.RData")
load("out/baro_l_results.RData")

#bar plots for relative abundance
#lines for d13C values
par(mfrow=c(3,3)) #900*800
barplot(RA.rhum.l,ylim = c(0,1),names.arg=c("C29","C31","C33"),space=0.5,xlim = c(0,5),axes=F,
        main="Rhum, fC4 = 0.72")
axis(2,c(0,0.2,0.4,0.6,0.8))
text(c(1,2.5,4),RA.rhum.l,labels=c("0.317","0.384","0.299"),pos=1)

barplot(RA.asso.l,ylim = c(0,1),names.arg=c("C29","C31","C33"),space=0.5,xlim = c(0,5),axes=F,
        main="Asso, fC4 = 0.31")
axis(2,c(0,0.2,0.4,0.6,0.8))
text(c(1,2.5,4),RA.asso.l,labels=c("0.502","0.289","0.209"),pos=1)

barplot(RA.baro.l,ylim = c(0,1),names.arg=c("C29","C31","C33"),space=0.5,xlim = c(0,5),axes=F,
        main="Baro, fC4 = 0.05")
axis(2,c(0,0.2,0.4,0.6,0.8))
text(c(1,2.5,4),RA.baro.l,labels=c("0.317","0.462","0.167"),pos=1)

barplot(RA.rhum.l,ylim = c(-40,-20),space=0.5,xlim = c(0,5),axes=F)
axis(4,c(-38,-35,-32,-29,-26,-23))
lines(c(1,2.5,4),d13C.rhum.l,lwd = 1.5, lty = 2)
arrows(c(1,2.5,4), d13C.rhum.l+d13C.sd.rhum.l, c(1,2.5,4), d13C.rhum.l-d13C.sd.rhum.l, 
       code=3, angle=90, length=0.05)
points(c(1,2.5,4),d13C.rhum.l,cex = 1.5)

barplot(RA.asso.l,ylim = c(-40,-20),space=0.5,xlim = c(0,5),axes=F)
axis(4,c(-38,-35,-32,-29,-26,-23))
lines(c(1,2.5,4),d13C.asso.l,lwd = 1.5, lty = 2)
arrows(c(1,2.5,4), d13C.asso.l+d13C.sd.asso.l, c(1,2.5,4), d13C.asso.l-d13C.sd.asso.l, 
       code=3, angle=90, length=0.05)
points(c(1,2.5,4),d13C.asso.l,cex = 1.5)

barplot(RA.baro.l,ylim = c(-40,-20),space=0.5,xlim = c(0,5),axes=F)
axis(4,c(-38,-35,-32,-29,-26,-23))
lines(c(1,2.5,4),d13C.baro.l,lwd = 1.5, lty = 2)
arrows(c(1,2.5,4), d13C.baro.l+d13C.sd.asso.l, c(1,2.5,4), d13C.baro.l-d13C.sd.asso.l, 
       code=3, angle=90, length=0.05)
points(c(1,2.5,4),d13C.baro.l,cex = 1.5)

#posterior densites
plot(density(rhum.l.mix$BUGSoutput$sims.list$f[,1],from=0,to=1), 
     main = "", xlab="fraction", xlim=c(0,1), ylim =c(0,8),lwd=2,col = plot.col[6])
lines(density(rhum.l.mix$BUGSoutput$sims.list$f[,2],from=0,to=1),lwd=2,col = plot.col[4])
lines(density(rhum.l.mix$BUGSoutput$sims.list$f[,3],from=0,to=1),lwd=2,col = plot.col[2])

plot(density(asso.l.mix$BUGSoutput$sims.list$f[,1],from=0,to=1), 
     main = "", xlab="fraction", xlim=c(0,1), ylim =c(0,8),lwd=2,col = plot.col[6])
lines(density(asso.l.mix$BUGSoutput$sims.list$f[,2],from=0,to=1),lwd=2,col = plot.col[4])
lines(density(asso.l.mix$BUGSoutput$sims.list$f[,3],from=0,to=1),lwd=2,col = plot.col[2])

plot(density(baro.l.mix$BUGSoutput$sims.list$f[,1],from=0,to=1), 
     main = "", xlab="fraction", xlim=c(0,1), ylim =c(0,8),lwd=2,col = plot.col[6])
lines(density(baro.l.mix$BUGSoutput$sims.list$f[,2],from=0,to=1),lwd=2,col = plot.col[4])
lines(density(baro.l.mix$BUGSoutput$sims.list$f[,3],from=0,to=1),lwd=2,col = plot.col[2])
legend(0.4,8,c("C4","Savanna C3","RF C3"),lwd=c(2,2,2),
       col=plot.col[c(6,4,2)])

######CS2 Bivariate density plots####
rhum.l.countours<-contour.fn(rhum.l.mix$BUGSoutput$sims.list$f)
asso.l.countours<-contour.fn(asso.l.mix$BUGSoutput$sims.list$f)
baro.l.countours<-contour.fn(baro.l.mix$BUGSoutput$sims.list$f)

#plot the data by column, each column is a sample
#from left to right: Rhum, Asso, Baro
par(mfrow=c(3,3)) #700*760
image(rhum.l.countours[[4]],col=viridis(64),xlab="f C4 GR",ylab="f SV C3",
      main="Rhum")
contour(rhum.l.countours[[4]], lwd = 0.6, add = TRUE, labcex = 0.5)

image(asso.l.countours[[4]],col=viridis(64),xlab="f C4 GR",ylab="f SV C3",
      main="Asso")
contour(asso.l.countours[[4]], lwd = 0.6, add = TRUE, labcex = 0.5)

image(baro.l.countours[[4]],col=viridis(64),xlab="f C4 GR",ylab="f SV C3",
      main="Baro")
contour(baro.l.countours[[4]], lwd = 0.6, add = TRUE, labcex = 0.5)

image(rhum.l.countours[[5]],col=viridis(64),xlab="f C4 GR",ylab="f RF C3")
contour(rhum.l.countours[[5]], lwd = 0.6, add = TRUE, labcex = 0.5)

image(asso.l.countours[[5]],col=viridis(64),xlab="f C4 GR",ylab="f RF C3")
contour(asso.l.countours[[5]], lwd = 0.6, add = TRUE, labcex = 0.5)

image(baro.l.countours[[5]],col=viridis(64),xlab="f C4 GR",ylab="f RF C3")
contour(baro.l.countours[[5]], lwd = 0.6, add = TRUE, labcex = 0.5)

image(rhum.l.countours[[6]],col=viridis(64),xlab="f SV C3",ylab="f RF C3")
contour(rhum.l.countours[[6]], lwd = 0.6, add = TRUE, labcex = 0.5)

image(asso.l.countours[[6]],col=viridis(64),xlab="f SV C3",ylab="f RF C3")
contour(asso.l.countours[[6]], lwd = 0.6, add = TRUE, labcex = 0.5)

image(baro.l.countours[[6]],col=viridis(64),xlab="f SV C3",ylab="f RF C3")
contour(baro.l.countours[[6]], lwd = 0.6, add = TRUE, labcex = 0.5)

###CS2: chain specific mixing ratios can be used to help to interpret d2H####
rhum.l.c.spec.ratio<-chain.spec.ratio(rhum.l.mix$BUGSoutput$sims.list$f.sum.prod_n_i)
asso.l.c.spec.ratio<-chain.spec.ratio(asso.l.mix$BUGSoutput$sims.list$f.sum.prod_n_i)
baro.l.c.spec.ratio<-chain.spec.ratio(baro.l.mix$BUGSoutput$sims.list$f.sum.prod_n_i)

#plot the data by column, each column is a sample
#from left to right: Rhum, Asso, Baro
#rows are by chain: 29, 31, 33
par(mfrow=c(3,3))
plot(density(rhum.l.c.spec.ratio[,1,1],from=0,to=1),main="n-alkane contribution C29",
     xlim=c(0,1),ylim=c(0,15),lwd=2,col=plot.col[6], xlab="fraction")
lines(density(rhum.l.c.spec.ratio[,2,1],from=0,to=1),lwd=2,col=plot.col[4])
lines(density(rhum.l.c.spec.ratio[,3,1],from=0,to=1),lwd=2,col=plot.col[2])

plot(density(asso.l.c.spec.ratio[,1,1],from=0,to=1),main="n-alkane contribution C29",
     xlim=c(0,1),ylim=c(0,15),lwd=2,col=plot.col[6], xlab="fraction")
lines(density(asso.l.c.spec.ratio[,2,1],from=0,to=1),lwd=2,col=plot.col[4])
lines(density(asso.l.c.spec.ratio[,3,1],from=0,to=1),lwd=2,col=plot.col[2])

plot(density(baro.l.c.spec.ratio[,1,1],from=0,to=1),main="n-alkane contribution C29",
     xlim=c(0,1),ylim=c(0,15),lwd=2,col=plot.col[6], xlab="fraction")
lines(density(baro.l.c.spec.ratio[,2,1],from=0,to=1),lwd=2,col=plot.col[4])
lines(density(baro.l.c.spec.ratio[,3,1],from=0,to=1),lwd=2,col=plot.col[2])

plot(density(rhum.l.c.spec.ratio[,1,2],from=0,to=1),main="n-alkane contribution C31",
     xlim=c(0,1),ylim=c(0,15),lwd=2,col=plot.col[6], xlab="fraction")
lines(density(rhum.l.c.spec.ratio[,2,2],from=0,to=1),lwd=2,col=plot.col[4])
lines(density(rhum.l.c.spec.ratio[,3,2],from=0,to=1),lwd=2,col=plot.col[2])

plot(density(asso.l.c.spec.ratio[,1,2],from=0,to=1),main="n-alkane contribution C31",
     xlim=c(0,1),ylim=c(0,15),lwd=2,col=plot.col[6], xlab="fraction")
lines(density(asso.l.c.spec.ratio[,2,2],from=0,to=1),lwd=2,col=plot.col[4])
lines(density(asso.l.c.spec.ratio[,3,2],from=0,to=1),lwd=2,col=plot.col[2])

plot(density(baro.l.c.spec.ratio[,1,2],from=0,to=1),main="n-alkane contribution C31",
     xlim=c(0,1),ylim=c(0,15),lwd=2,col=plot.col[6], xlab="fraction")
lines(density(baro.l.c.spec.ratio[,2,2],from=0,to=1),lwd=2,col=plot.col[4])
lines(density(baro.l.c.spec.ratio[,3,2],from=0,to=1),lwd=2,col=plot.col[2])

plot(density(rhum.l.c.spec.ratio[,1,3],from=0,to=1),main="n-alkane contribution C33",
     xlim=c(0,1),ylim=c(0,15),lwd=2,col=plot.col[6], xlab="fraction")
lines(density(rhum.l.c.spec.ratio[,2,3],from=0,to=1),lwd=2,col=plot.col[4])
lines(density(rhum.l.c.spec.ratio[,3,3],from=0,to=1),lwd=2,col=plot.col[2])

plot(density(asso.l.c.spec.ratio[,1,3],from=0,to=1),main="n-alkane contribution C33",
     xlim=c(0,1),ylim=c(0,15),lwd=2,col=plot.col[6], xlab="fraction")
lines(density(asso.l.c.spec.ratio[,2,3],from=0,to=1),lwd=2,col=plot.col[4])
lines(density(asso.l.c.spec.ratio[,3,3],from=0,to=1),lwd=2,col=plot.col[2])

plot(density(baro.l.c.spec.ratio[,1,3],from=0,to=1),main="n-alkane contribution C33",
     xlim=c(0,1),ylim=c(0,15),lwd=2,col=plot.col[6], xlab="fraction")
lines(density(baro.l.c.spec.ratio[,2,3],from=0,to=1),lwd=2,col=plot.col[4])
lines(density(baro.l.c.spec.ratio[,3,3],from=0,to=1),lwd=2,col=plot.col[2])
legend(0.3,15,c("C4","Savanna C3","RF C3"),lwd=c(2,2,2),
       col=c(plot.col[c(6,4,2)]))



###sensitivity test1: different priors####

###prior2: w.aftrica distributions####
#use the same functions to calculate prior densities
#concentration
W.prod.prior.GR<-pri.multi.norm.den(-4,13,W.prod.mu[1,],W.prod.vcov[1:3,1:3])
W.prod.prior.SV<-pri.multi.norm.den(-4,13,W.prod.mu[2,],W.prod.vcov[4:6,1:3])
W.prod.prior.RF<-pri.multi.norm.den(-4,13,W.prod.mu[3,],W.prod.vcov[7:9,1:3])

#comparing two priors
#default prior is blue
#western Africa only prior is red
par(mfrow=c(3,3))
plot(W.prod.prior.GR$x,W.prod.prior.GR$y1, type = "l", lwd=2, col = "red",
     xlim = c(-4,13), ylim = c(0,0.4), main = "n-C29 concentration")
lines(Afr.prod.prior.GR$x,Afr.prod.prior.GR$y1, col = "blue", lwd = 1)

plot(W.prod.prior.GR$x,W.prod.prior.GR$y2, type = "l", lwd=2, col = "red",
     xlim = c(-4,13), ylim = c(0,0.4), main = "n-C31 concentration")
lines(Afr.prod.prior.GR$x,Afr.prod.prior.GR$y2, col = "blue", lwd = 1)

plot(W.prod.prior.GR$x,W.prod.prior.GR$y3, type = "l", lwd=2, col = "red",
     xlim = c(-4,13),ylim = c(0,0.4), main = "n-C33 concentration")
lines(Afr.prod.prior.GR$x,Afr.prod.prior.GR$y3, col = "blue", lwd = 1)

plot(W.prod.prior.SV$x,W.prod.prior.SV$y1, type = "l", lwd=2, col = "red",
     xlim = c(-4,13),ylim = c(0,0.3), main = "n-C29 concentration")
lines(Afr.prod.prior.SV$x,Afr.prod.prior.SV$y1, col = "blue", lwd = 1)

plot(W.prod.prior.SV$x,W.prod.prior.SV$y2, type = "l", lwd=2, col = "red",
     xlim = c(-4,13),ylim = c(0,0.3), main = "n-C31 concentration")
lines(Afr.prod.prior.SV$x,Afr.prod.prior.SV$y2, col = "blue", lwd = 1)

plot(W.prod.prior.SV$x,W.prod.prior.SV$y3, type = "l", lwd=2, col = "red",
     xlim = c(-4,13),ylim = c(0,0.3), main = "n-C33 concentration")
lines(Afr.prod.prior.SV$x,Afr.prod.prior.SV$y3, col = "blue", lwd = 1)

plot(W.prod.prior.RF$x,W.prod.prior.RF$y1, type = "l", lwd=2, col = "red",
     xlim = c(-4,13),ylim = c(0,0.3), main = "n-C29 concentration")
lines(Afr.prod.prior.RF$x,Afr.prod.prior.RF$y1, col = "blue", lwd = 1)

plot(W.prod.prior.RF$x,W.prod.prior.RF$y2, type = "l", lwd=2, col = "red",
     xlim = c(-4,13),ylim = c(0,0.3), main = "n-C31 concentration")
lines(Afr.prod.prior.RF$x,Afr.prod.prior.RF$y2, col = "blue", lwd = 1)

plot(W.prod.prior.RF$x,W.prod.prior.RF$y3, type = "l", lwd=2, col = "red",
     xlim = c(-4,13),ylim = c(0,0.3), main = "n-C33 concentration", axes=F)
lines(Afr.prod.prior.RF$x,Afr.prod.prior.RF$y3, col = "blue", lwd = 1)
#log-scale axis
axis(2,c(0,0.5))
axis(1,log(c(0.01,0.1,1,10,100,1e3,1e4,1e5)))

#d13C
W.d13C.prior.GR<-pri.multi.norm.den(-50,-15,W.d13C.mu[1,],W.d13C.vcov[1:3,1:3])
W.d13C.prior.SV<-pri.multi.norm.den(-50,-15,W.d13C.mu[2,],W.d13C.vcov[4:6,1:3])
W.d13C.prior.RF<-pri.multi.norm.den(-50,-15,W.d13C.mu[3,],W.d13C.vcov[7:9,1:3])

#comparing two priors
#default prior is blue
#western Africa only prior is red
par(mfrow=c(3,3)) #900*800
#C4 grasses 
plot(W.d13C.prior.GR$x,W.d13C.prior.GR$y1, type = "l", lwd=2, col = "red",
     xlim = c(-50,-15), ylim = c(0,0.3), main = "C4 n-C29 d13C")
lines(Afr.d13C.prior.GR$x,Afr.d13C.prior.GR$y1, col = "blue", lwd = 1)

plot(W.d13C.prior.GR$x,W.d13C.prior.GR$y2, type = "l", lwd=2, col = "red",
     xlim = c(-50,-15),ylim = c(0,0.3), main = "C4 n-C31 d13C")
lines(Afr.d13C.prior.GR$x,Afr.d13C.prior.GR$y2, col = "blue", lwd = 1)

plot(W.d13C.prior.GR$x,W.d13C.prior.GR$y3, type = "l", lwd=2, col = "red",
     xlim = c(-50,-15), ylim = c(0,0.3), main = "C4 n-C33 d13C")
lines(Afr.d13C.prior.GR$x,Afr.d13C.prior.GR$y3, col = "blue", lwd = 1)

#savanna C3 plants
plot(W.d13C.prior.SV$x,W.d13C.prior.SV$y1, type = "l", lwd=2, col = "red",
     xlim = c(-50,-15),ylim = c(0,0.3), main = "SV n-C29 d13C")
lines(Afr.d13C.prior.SV$x,Afr.d13C.prior.SV$y1, col = "blue", lwd = 1)

plot(W.d13C.prior.SV$x,W.d13C.prior.SV$y2, type = "l", lwd=2, col = "red",
     xlim = c(-50,-15), ylim = c(0,0.3), main = "SV n-C31 d13C")
lines(Afr.d13C.prior.SV$x,Afr.d13C.prior.SV$y2, col = "blue", lwd = 1)

plot(W.d13C.prior.SV$x,W.d13C.prior.SV$y3, type = "l", lwd=2, col = "red",
     xlim = c(-50,-15), ylim = c(0,0.3), main = "SV n-C33 d13C")
lines(Afr.d13C.prior.SV$x,Afr.d13C.prior.SV$y3, col = "blue", lwd = 1)

#Rainforest plants
plot(W.d13C.prior.RF$x,W.d13C.prior.RF$y1, type = "l", lwd=2, col = "red",
     xlim = c(-50,-15),ylim = c(0,0.3), main = "RF n-C29 d13C")
lines(Afr.d13C.prior.RF$x,Afr.d13C.prior.RF$y1, col = "blue", lwd = 1)

plot(W.d13C.prior.RF$x,W.d13C.prior.RF$y2, type = "l", lwd=2, col = "red",
     xlim = c(-50,-15),ylim = c(0,0.3), main = "RF n-C31 d13C")
lines(Afr.d13C.prior.RF$x,Afr.d13C.prior.RF$y2, col = "blue", lwd = 1)

plot(W.d13C.prior.RF$x,W.d13C.prior.RF$y3, type = "l", lwd=2, col = "red",
     xlim = c(-50,-15), ylim = c(0,0.3), main = "RF n-C33 d13C")
lines(Afr.d13C.prior.RF$x,Afr.d13C.prior.RF$y3, col = "blue", lwd = 1)

###Sensitivity test1 vs CS2 mixing ratios####
#run the corresponding code in the CS2 file
load("out/W_rhum_l_results.RData")
load("out/W_asso_l_results.RData")
load("out/W_baro_l_results.RData")

#comparing results using two priors
#first row: western Africa prior
#second row: CS2 "subsaharan Africa" prior
par(mfrow=c(2,3))
plot(density(W.rhum.l.mix$BUGSoutput$sims.list$f[,1],from=0,to=1), 
     main="Rhum fC4 = 0.72",xlab="fraction",xlim=c(0,1),col=plot.col[6], ylim =c(0,8),lwd=2)
lines(density(W.rhum.l.mix$BUGSoutput$sims.list$f[,2],from=0,to=1),col=plot.col[4],lwd=2)
lines(density(W.rhum.l.mix$BUGSoutput$sims.list$f[,3],from=0,to=1),col=plot.col[2],lwd=2)

plot(density(W.asso.l.mix$BUGSoutput$sims.list$f[,1],from=0,to=1), 
     main="Asso fC4 = 0.31",xlab="fraction",xlim=c(0,1),col=plot.col[6], ylim =c(0,8),lwd=2)
lines(density(W.asso.l.mix$BUGSoutput$sims.list$f[,2],from=0,to=1),col=plot.col[4],lwd=2)
lines(density(W.asso.l.mix$BUGSoutput$sims.list$f[,3],from=0,to=1),col=plot.col[2],lwd=2)

plot(density(W.baro.l.mix$BUGSoutput$sims.list$f[,1],from=0,to=1), 
     main="Baro fC4 = 0.05",xlab="fraction",xlim=c(0,1),col=plot.col[6], ylim =c(0,8),lwd=2)
lines(density(W.baro.l.mix$BUGSoutput$sims.list$f[,2],from=0,to=1),col=plot.col[4],lwd=2)
lines(density(W.baro.l.mix$BUGSoutput$sims.list$f[,3],from=0,to=1),col=plot.col[2],lwd=2)

plot(density(rhum.l.mix$BUGSoutput$sims.list$f[,1],from=0,to=1), 
     main = "", xlab="fraction", xlim=c(0,1), ylim =c(0,8),lwd=2,col = plot.col[6])
lines(density(rhum.l.mix$BUGSoutput$sims.list$f[,2],from=0,to=1),lwd=2,col = plot.col[4])
lines(density(rhum.l.mix$BUGSoutput$sims.list$f[,3],from=0,to=1),lwd=2,col = plot.col[2])

plot(density(asso.l.mix$BUGSoutput$sims.list$f[,1],from=0,to=1), 
     main = "", xlab="fraction", xlim=c(0,1), ylim =c(0,8),lwd=2,col = plot.col[6])
lines(density(asso.l.mix$BUGSoutput$sims.list$f[,2],from=0,to=1),lwd=2,col = plot.col[4])
lines(density(asso.l.mix$BUGSoutput$sims.list$f[,3],from=0,to=1),lwd=2,col = plot.col[2])

plot(density(baro.l.mix$BUGSoutput$sims.list$f[,1],from=0,to=1), 
     main = "", xlab="fraction", xlim=c(0,1), ylim =c(0,8),lwd=2,col = plot.col[6])
lines(density(baro.l.mix$BUGSoutput$sims.list$f[,2],from=0,to=1),lwd=2,col = plot.col[4])
lines(density(baro.l.mix$BUGSoutput$sims.list$f[,3],from=0,to=1),lwd=2,col = plot.col[2])
legend(0.4,8,c("C4 Grass","C3 Savanna","C3 Forest"),lwd=c(2,2,2),
       col=plot.col[c(6,4,2)])

###Sensitivity test2: sensitivity to likelihood functions####
#run the corresponding code in the CS2 file
load("out/rhum_l_test_a.RData")
load("out/asso_l_test_a.RData")
load("out/baro_l_test_a.RData")
load("out/rhum_l_test_b.RData")
load("out/asso_l_test_b.RData")
load("out/baro_l_test_b.RData")
load("out/rhum_l_test_c.RData")
load("out/asso_l_test_c.RData")
load("out/baro_l_test_c.RData")
load("out/rhum_l_test_d.RData")
load("out/asso_l_test_d.RData")
load("out/baro_l_test_d.RData")

#plotting order
#rows: test a|  test b|  default| test d| test c|
#columns: sample
par(mfrow=c(3,5))#1200*600
#rhum
plot(density(rhum.l.test.a$BUGSoutput$sims.list$f[,1],from=0,to=1), 
     main="No RA",xlab="fraction",xlim=c(0,1),col=plot.col[6], ylim =c(0,6),lwd=2)
lines(density(rhum.l.test.a$BUGSoutput$sims.list$f[,2],from=0,to=1),col=plot.col[4],lwd=2)
lines(density(rhum.l.test.a$BUGSoutput$sims.list$f[,3],from=0,to=1),col=plot.col[2],lwd=2)

plot(density(rhum.l.test.b$BUGSoutput$sims.list$f[,1],from=0,to=1),
     main="inflated RA",xlab="fraction",xlim=c(0,1),col=plot.col[6], ylim =c(0,6),lwd=2)
lines(density(rhum.l.test.b$BUGSoutput$sims.list$f[,2],from=0,to=1),col=plot.col[4],lwd=2)
lines(density(rhum.l.test.b$BUGSoutput$sims.list$f[,3],from=0,to=1),col=plot.col[2],lwd=2)

plot(density(rhum.l.mix$BUGSoutput$sims.list$f[,1],from=0,to=1), 
     main = "Rhum default", xlab="fraction", xlim=c(0,1), ylim =c(0,6),lwd=2,col = plot.col[6])
lines(density(rhum.l.mix$BUGSoutput$sims.list$f[,2],from=0,to=1),lwd=2,col = plot.col[4])
lines(density(rhum.l.mix$BUGSoutput$sims.list$f[,3],from=0,to=1),lwd=2,col = plot.col[2])

plot(density(rhum.l.test.d$BUGSoutput$sims.list$f[,1],from=0,to=1),
     main="inflated d13C",xlab="fraction",xlim=c(0,1),col=plot.col[6], ylim =c(0,6),lwd=2)
lines(density(rhum.l.test.d$BUGSoutput$sims.list$f[,2],from=0,to=1),col=plot.col[4],lwd=2)
lines(density(rhum.l.test.d$BUGSoutput$sims.list$f[,3],from=0,to=1),col=plot.col[2],lwd=2)

plot(density(rhum.l.test.c$BUGSoutput$sims.list$f[,1],from=0,to=1), 
     main="No d13C",xlab="fraction",xlim=c(0,1),col=plot.col[6], ylim =c(0,6),lwd=2)
lines(density(rhum.l.test.c$BUGSoutput$sims.list$f[,2],from=0,to=1),col=plot.col[4],lwd=2)
lines(density(rhum.l.test.c$BUGSoutput$sims.list$f[,3],from=0,to=1),col=plot.col[2],lwd=2)

#asso
plot(density(asso.l.test.a$BUGSoutput$sims.list$f[,1],from=0,to=1), 
     main="No RA",xlab="fraction",xlim=c(0,1),col=plot.col[6], ylim =c(0,8),lwd=2)
lines(density(asso.l.test.a$BUGSoutput$sims.list$f[,2],from=0,to=1),col=plot.col[4],lwd=2)
lines(density(asso.l.test.a$BUGSoutput$sims.list$f[,3],from=0,to=1),col=plot.col[2],lwd=2)

plot(density(asso.l.test.b$BUGSoutput$sims.list$f[,1],from=0,to=1),
     main="inflated RA",xlab="fraction",xlim=c(0,1),col=plot.col[6], ylim =c(0,8),lwd=2)
lines(density(asso.l.test.b$BUGSoutput$sims.list$f[,2],from=0,to=1),col=plot.col[4],lwd=2)
lines(density(asso.l.test.b$BUGSoutput$sims.list$f[,3],from=0,to=1),col=plot.col[2],lwd=2)

plot(density(asso.l.mix$BUGSoutput$sims.list$f[,1],from=0,to=1), 
     main = "Asso default", xlab="fraction", xlim=c(0,1), ylim =c(0,8),lwd=2,col = plot.col[6])
lines(density(asso.l.mix$BUGSoutput$sims.list$f[,2],from=0,to=1),lwd=2,col = plot.col[4])
lines(density(asso.l.mix$BUGSoutput$sims.list$f[,3],from=0,to=1),lwd=2,col = plot.col[2])

plot(density(asso.l.test.d$BUGSoutput$sims.list$f[,1],from=0,to=1),
     main="inflated d13C",xlab="fraction",xlim=c(0,1),col=plot.col[6], ylim =c(0,8),lwd=2)
lines(density(asso.l.test.d$BUGSoutput$sims.list$f[,2],from=0,to=1),col=plot.col[4],lwd=2)
lines(density(asso.l.test.d$BUGSoutput$sims.list$f[,3],from=0,to=1),col=plot.col[2],lwd=2)

plot(density(asso.l.test.c$BUGSoutput$sims.list$f[,1],from=0,to=1), 
     main="No d13C",xlab="fraction",xlim=c(0,1),col=plot.col[6], ylim =c(0,8),lwd=2)
lines(density(asso.l.test.c$BUGSoutput$sims.list$f[,2],from=0,to=1),col=plot.col[4],lwd=2)
lines(density(asso.l.test.c$BUGSoutput$sims.list$f[,3],from=0,to=1),col=plot.col[2],lwd=2)

#Baro
plot(density(baro.l.test.a$BUGSoutput$sims.list$f[,1],from=0,to=1), 
     main="No RA",xlab="fraction",xlim=c(0,1),col=plot.col[6], ylim =c(0,8),lwd=2)
lines(density(baro.l.test.a$BUGSoutput$sims.list$f[,2],from=0,to=1),col=plot.col[4],lwd=2)
lines(density(baro.l.test.a$BUGSoutput$sims.list$f[,3],from=0,to=1),col=plot.col[2],lwd=2)

plot(density(baro.l.test.b$BUGSoutput$sims.list$f[,1],from=0,to=1),
     main="inflated RA",xlab="fraction",xlim=c(0,1),col=plot.col[6], ylim =c(0,8),lwd=2)
lines(density(baro.l.test.b$BUGSoutput$sims.list$f[,2],from=0,to=1),col=plot.col[4],lwd=2)
lines(density(baro.l.test.b$BUGSoutput$sims.list$f[,3],from=0,to=1),col=plot.col[2],lwd=2)

plot(density(baro.l.mix$BUGSoutput$sims.list$f[,1],from=0,to=1), 
     main = "Baro default", xlab="fraction", xlim=c(0,1), ylim =c(0,8),lwd=2,col = plot.col[6])
lines(density(baro.l.mix$BUGSoutput$sims.list$f[,2],from=0,to=1),lwd=2,col = plot.col[4])
lines(density(baro.l.mix$BUGSoutput$sims.list$f[,3],from=0,to=1),lwd=2,col = plot.col[2])

plot(density(baro.l.test.d$BUGSoutput$sims.list$f[,1],from=0,to=1),
     main="inflated d13C",xlab="fraction",xlim=c(0,1),col=plot.col[6], ylim =c(0,8),lwd=2)
lines(density(baro.l.test.d$BUGSoutput$sims.list$f[,2],from=0,to=1),col=plot.col[4],lwd=2)
lines(density(baro.l.test.d$BUGSoutput$sims.list$f[,3],from=0,to=1),col=plot.col[2],lwd=2)

plot(density(baro.l.test.c$BUGSoutput$sims.list$f[,1],from=0,to=1), 
     main="No d13C",xlab="fraction",xlim=c(0,1),col=plot.col[6], ylim =c(0,8),lwd=2)
lines(density(baro.l.test.c$BUGSoutput$sims.list$f[,2],from=0,to=1),col=plot.col[4],lwd=2)
lines(density(baro.l.test.c$BUGSoutput$sims.list$f[,3],from=0,to=1),col=plot.col[2],lwd=2)
legend(0.2,8,c("C4 Grass","C3 Savanna","C3 Forest"),lwd=c(2,2,2),
       col=plot.col[c(6,4,2)])

#####Supplementary plots#####
#Q-Q plot of QTP sources#
####parameter estimates leaf wax production, log-normal distribution, 
##using the elnorm function in the EnvStat package
#then use the estimates to make Q-Q plots

#terrestrial source
elnorm(Ter$c.n.C27,method = "mle")
Ter.C27.mean<-3.678073
Ter.C27.sd<-1.041735
elnorm(Ter$c.n.C29,method = "mle")
Ter.C29.mean<-4.307679
Ter.C29.sd<-1.344295
elnorm(Ter$c.n.C31,method = "mle")
Ter.C31.mean<-3.522836
Ter.C31.sd<-1.616805

#macrophyte source
elnorm(Mac$c.n.C27,method = "mle")
Mac.C27.mean<-2.747748
Mac.C27.sd<-1.083799
elnorm(Mac$c.n.C29,method = "mle")
Mac.C29.mean<-2.2110273
Mac.C29.sd<-0.9335903
elnorm(Mac$c.n.C31,method = "mle")
Mac.C31.mean<--0.0476248
Mac.C31.sd<-1.0787635

#algal source
elnorm(Alg$c.n.C27,method = "mle")
Alg.C27.mean<--0.8769988
Alg.C27.sd<-1.1243304
elnorm(Alg$c.n.C29,method = "mle")
Alg.C29.mean<--1.045367
Alg.C29.sd<-1.004963
elnorm(Alg$c.n.C31,method = "mle")
Alg.C31.mean<--1.6010903
Alg.C31.sd<-0.7567727

par(mfrow=c(3,3))
qqPlot(Ter$c.n.C27,distribution = "lnorm",
       param.list=list(mean=Ter.C27.mean,sd=Ter.C27.sd),add.line=T)
qqPlot(Ter$c.n.C29,distribution = "lnorm",
       param.list=list(mean=Ter.C29.mean,sd=Ter.C29.sd),add.line=T)
qqPlot(Ter$c.n.C31,distribution = "lnorm",
       param.list=list(mean=Ter.C31.mean,sd=Ter.C31.sd),add.line=T)

qqPlot(Mac$c.n.C27,distribution = "lnorm",
       param.list=list(mean=Mac.C27.mean,sd=Mac.C27.sd),add.line=T)
qqPlot(Mac$c.n.C29,distribution = "lnorm",
       param.list=list(mean=Mac.C29.mean,sd=Mac.C29.sd),add.line=T)
qqPlot(Mac$c.n.C31,distribution = "lnorm",
       param.list=list(mean=Mac.C31.mean,sd=Mac.C31.sd),add.line=T)

qqPlot(Alg$c.n.C27,distribution = "lnorm",
       param.list=list(mean=Alg.C27.mean,sd=Alg.C27.sd),add.line=T)
qqPlot(Alg$c.n.C29,distribution = "lnorm",
       param.list=list(mean=Alg.C29.mean,sd=Alg.C29.sd),add.line=T)
qqPlot(Alg$c.n.C31,distribution = "lnorm",
       param.list=list(mean=Alg.C31.mean,sd=Alg.C31.sd),add.line=T)

#Q-Q plot of Afr sources#
####parameter estimates leaf wax production, log-normal distribution, 
##using the elnorm function in the EnvStat package
#then use the estimates to make Q-Q plots

#C4 source
elnorm(GR$c.n.C29,method = "mle")
GR.c.n.C29.mean<-3.443274
GR.c.n.C29.sd<-1.161768
elnorm(GR$c.n.C31,method = "mle")
GR.c.n.C31.mean<-4.428663
GR.c.n.C31.sd<-1.244289
elnorm(GR$c.n.C33,method = "mle")
GR.c.n.C33.mean<-3.749900
GR.c.n.C33.sd<-1.259917

#savanna C3
elnorm(SV$c.n.C29,method = "mle")
SV.c.n.C29.mean<-3.486172
SV.c.n.C29.sd<-1.428710
elnorm(SV$c.n.C31,method = "mle")
SV.c.n.C31.mean<-3.697397
SV.c.n.C31.sd<-1.555444 
elnorm(SV$c.n.C33,method = "mle")
SV.c.n.C33.mean<-2.494557
SV.c.n.C33.sd<-2.009675

#rainforest C3
elnorm(RF$c.n.C29,method = "mle")
RF.c.n.C29.mean<-3.977994
RF.c.n.C29.sd<-1.797465
elnorm(RF$c.n.C31,method = "mle")
RF.c.n.C31.mean<-3.827783
RF.c.n.C31.sd<-1.340019
elnorm(RF$c.n.C33,method = "mle")
RF.c.n.C33.mean<-2.058397
RF.c.n.C33.sd<-1.469450

par(mfrow=c(3,3))
qqPlot(GR$c.n.C29,distribution = "lnorm",
       param.list=list(mean=GR.c.n.C29.mean,sd=GR.c.n.C29.sd),add.line=T)
qqPlot(GR$c.n.C31,distribution = "lnorm",
       param.list=list(mean=GR.c.n.C31.mean,sd=GR.c.n.C31.sd),add.line=T)
qqPlot(GR$c.n.C33,distribution = "lnorm",
       param.list=list(mean=GR.c.n.C33.mean,sd=GR.c.n.C33.sd),add.line=T)

qqPlot(SV$c.n.C29,distribution = "lnorm",
       param.list=list(mean=SV.c.n.C29.mean,sd=SV.c.n.C29.sd),add.line=T)
qqPlot(SV$c.n.C31,distribution = "lnorm",
       param.list=list(mean=SV.c.n.C31.mean,sd=SV.c.n.C31.sd),add.line=T)
qqPlot(SV$c.n.C33,distribution = "lnorm",
       param.list=list(mean=SV.c.n.C33.mean,sd=SV.c.n.C33.sd),add.line=T)

qqPlot(RF$c.n.C29,distribution = "lnorm",
       param.list=list(mean=RF.c.n.C29.mean,sd=RF.c.n.C29.sd),add.line=T)
qqPlot(RF$c.n.C31,distribution = "lnorm",
       param.list=list(mean=RF.c.n.C31.mean,sd=RF.c.n.C31.sd),add.line=T)
qqPlot(RF$c.n.C33,distribution = "lnorm",
       param.list=list(mean=RF.c.n.C33.mean,sd=RF.c.n.C33.sd),add.line=T)

####CS1 QTP prior vs posterior density of n-alkane concentration###
#blue is prior, red is posterior
par(mfrow=c(3,3)) #900*800
plot(density(log(QHS13_5S.mix$BUGSoutput$sims.list$exp.prod_k[,,1,1])),
     col = "red", lwd = 2,type="l",main="n-C27 Terrestrial",xlim=c(-5,10),ylim=c(0,0.4))
lines(QTP.prod.prior.ter$x,QTP.prod.prior.ter$y1, col = "blue", lwd = 1)

plot(density(log(QHS13_5S.mix$BUGSoutput$sims.list$exp.prod_k[,,1,2])),
     col = "red", lwd = 2,type="l",main="n-C29 Terrestrial",xlim=c(-5,10),ylim=c(0,0.4))
lines(QTP.prod.prior.ter$x,QTP.prod.prior.ter$y2, col = "blue", lwd = 1)

plot(density(log(QHS13_5S.mix$BUGSoutput$sims.list$exp.prod_k[,,1,3])),
     col = "red", lwd = 2,type="l",main="n-C33 Terrestrial",xlim=c(-5,10),ylim=c(0,0.4))
lines(QTP.prod.prior.ter$x,QTP.prod.prior.ter$y3, col = "blue", lwd = 1)

#macrophyte
plot(density(log(QHS13_5S.mix$BUGSoutput$sims.list$exp.prod_k[,,2,1])),
     col = "red", lwd = 2,type="l",main="n-C27 Macrophyte",xlim=c(-5,10),ylim = c(0,0.5))
lines(QTP.prod.prior.mac$x,QTP.prod.prior.mac$y1, col = "blue", lwd = 1)

plot(density(log(QHS13_5S.mix$BUGSoutput$sims.list$exp.prod_k[,,2,2])),
     col = "red", lwd = 2,type="l",main="n-C29 Macrophyte",xlim=c(-5,10),ylim = c(0,0.5))
lines(QTP.prod.prior.mac$x,QTP.prod.prior.mac$y2, col = "blue", lwd = 1)

plot(density(log(QHS13_5S.mix$BUGSoutput$sims.list$exp.prod_k[,,2,3])),
     col = "red", lwd = 2,type="l",main="n-C31 Macrophyte",xlim=c(-5,10),ylim = c(0,0.5))
lines(QTP.prod.prior.mac$x,QTP.prod.prior.mac$y3, col = "blue", lwd = 1)

#algae
plot(density(log(QHS13_5S.mix$BUGSoutput$sims.list$exp.prod_k[,,3,1])),
     col = "red", lwd = 2,type="l",main="n-C27 Algae",xlim=c(-5,10),ylim = c(0,0.6))
lines(QTP.prod.prior.alg$x,QTP.prod.prior.alg$y1, col = "blue", lwd = 1)

plot(density(log(QHS13_5S.mix$BUGSoutput$sims.list$exp.prod_k[,,3,2])),
     col = "red", lwd = 2,type="l",main="n-C29 Algae",xlim=c(-5,10),ylim = c(0,0.6))
lines(QTP.prod.prior.alg$x,QTP.prod.prior.alg$y2, col = "blue", lwd = 1)

plot(density(log(QHS13_5S.mix$BUGSoutput$sims.list$exp.prod_k[,,3,3])),
     col = "red", lwd = 2,type="l",main="n-C31 Algae",xlim=c(-5,10),ylim = c(0,0.6),axes=F)
lines(QTP.prod.prior.alg$x,QTP.prod.prior.alg$y3, col = "blue", lwd = 1)

#make a log-scale axis
axis(2,c(0,0.6))
axis(1,log(c(0.01,0.1,1,10,100,1000,10000)))

#####CS1 QTP prior vs posterior density of n-alkane d13C###
#blue is prior, red is posterior
plot(density(QHS13_5S.mix$BUGSoutput$sims.list$d13C.k[,,1,1]),
     col = "red", lwd = 2,type="l",main="n-C27 Terrestrial",xlim=c(-40,-10),ylim = c(0,0.4))
lines(QTP.d13C.prior.ter$x,QTP.d13C.prior.ter$y1, col = "blue", lwd = 1)

plot(density(QHS13_5S.mix$BUGSoutput$sims.list$d13C.k[,,1,2]),
     col = "red", lwd = 2,type="l",main="n-C29 Terrestrial",xlim=c(-40,-10),ylim = c(0,0.4))
lines(QTP.d13C.prior.ter$x,QTP.d13C.prior.ter$y2, col = "blue", lwd = 1)

plot(density(QHS13_5S.mix$BUGSoutput$sims.list$d13C.k[,,1,3]),
     col = "red", lwd = 2,type="l",main="n-C31 Terrestrial",xlim=c(-40,-10),ylim = c(0,0.4))
lines(QTP.d13C.prior.ter$x,QTP.d13C.prior.ter$y3, col = "blue", lwd = 1)

plot(density(QHS13_5S.mix$BUGSoutput$sims.list$d13C.k[,,2,1]),
     col = "red", lwd = 2,type="l",main="n-C27 Macrophyte",xlim=c(-40,-10),ylim = c(0,0.12))
lines(QTP.d13C.prior.mac$x,QTP.d13C.prior.mac$y1, col = "blue", lwd = 1)

plot(density(QHS13_5S.mix$BUGSoutput$sims.list$d13C.k[,,2,2]),
     col = "red", lwd = 2,type="l",main="n-C29 Macrophyte",xlim=c(-40,-10),ylim = c(0,0.12))
lines(QTP.d13C.prior.mac$x,QTP.d13C.prior.mac$y2, col = "blue", lwd = 1)

plot(density(QHS13_5S.mix$BUGSoutput$sims.list$d13C.k[,,2,3]),
     col = "red", lwd = 2,type="l",main="n-C31 Macrophyte",xlim=c(-40,-10),ylim = c(0,0.12))
lines(QTP.d13C.prior.mac$x,QTP.d13C.prior.mac$y3, col = "blue", lwd = 1)

plot(density(QHS13_5S.mix$BUGSoutput$sims.list$d13C.k[,,3,1]),
     col = "red", lwd = 2,type="l",main="n-C27 Algae",xlim=c(-40,-10),ylim = c(0,0.25))
lines(QTP.d13C.prior.alg$x,QTP.d13C.prior.alg$y1, col = "blue", lwd = 1)

plot(density(QHS13_5S.mix$BUGSoutput$sims.list$d13C.k[,,3,2]),
     col = "red", lwd = 2,type="l",main="n-C29 Algae",xlim=c(-40,-10),ylim = c(0,0.25))
lines(QTP.d13C.prior.alg$x,QTP.d13C.prior.alg$y2, col = "blue", lwd = 1)

plot(density(QHS13_5S.mix$BUGSoutput$sims.list$d13C.k[,,3,3]),
     col = "red", lwd = 2,type="l",main="n-C31 Algae",xlim=c(-40,-10),ylim = c(0,0.25))
lines(QTP.d13C.prior.alg$x,QTP.d13C.prior.alg$y3, col = "blue", lwd = 1)

###CS2: prior vs posterior density concentration####
#blue is prior, red is posterior
par(mfrow=c(3,3)) #700*800
plot(density(log(asso.l.mix$BUGSoutput$sims.list$exp.prod_k[,,1,1])),
     col = "red", lwd = 2,type="l",main="n-C29 C4",xlim=c(-4,13),ylim=c(0,0.4))
lines(Afr.prod.prior.GR$x,Afr.prod.prior.GR$y1, col = "blue", lwd = 1)

plot(density(log(asso.l.mix$BUGSoutput$sims.list$exp.prod_k[,,1,2])),
     col = "red", lwd = 2,type="l",main="n-C31 C4",xlim=c(-4,13),ylim=c(0,0.4))
lines(Afr.prod.prior.GR$x,Afr.prod.prior.GR$y2, col = "blue", lwd = 1)

plot(density(log(asso.l.mix$BUGSoutput$sims.list$exp.prod_k[,,1,3])),
     col = "red", lwd = 2,type="l",main="n-C33 C4",xlim=c(-4,13),ylim=c(0,0.4))
lines(Afr.prod.prior.GR$x,Afr.prod.prior.GR$y3, col = "blue", lwd = 1)

#SV
plot(density(log(asso.l.mix$BUGSoutput$sims.list$exp.prod_k[,,2,1])),
     col = "red", lwd = 2,type="l",main="n-C29 SV",xlim=c(-4,13),ylim = c(0,0.3))
lines(Afr.prod.prior.SV$x,Afr.prod.prior.SV$y1, col = "blue", lwd = 1)

plot(density(log(asso.l.mix$BUGSoutput$sims.list$exp.prod_k[,,2,2])),
     col = "red", lwd = 2,type="l",main="n-C31 SV",xlim=c(-4,13),ylim = c(0,0.3))
lines(Afr.prod.prior.SV$x,Afr.prod.prior.SV$y2, col = "blue", lwd = 1)

plot(density(log(asso.l.mix$BUGSoutput$sims.list$exp.prod_k[,,2,3])),
     col = "red", lwd = 2,type="l",main="n-C33 SV",xlim=c(-4,13),ylim = c(0,0.3))
lines(Afr.prod.prior.SV$x,Afr.prod.prior.SV$y3, col = "blue", lwd = 1)

#RF
plot(density(log(asso.l.mix$BUGSoutput$sims.list$exp.prod_k[,,3,1])),
     col = "red", lwd = 2,type="l",main="n-C29 RF",xlim=c(-4,13),ylim = c(0,0.3))
lines(Afr.prod.prior.RF$x,Afr.prod.prior.RF$y1, col = "blue", lwd = 1)

plot(density(log(asso.l.mix$BUGSoutput$sims.list$exp.prod_k[,,3,2])),
     col = "red", lwd = 2,type="l",main="n-C31 RF",xlim=c(-4,13),ylim = c(0,0.3))
lines(Afr.prod.prior.RF$x,Afr.prod.prior.RF$y2, col = "blue", lwd = 1)

plot(density(log(asso.l.mix$BUGSoutput$sims.list$exp.prod_k[,,3,3])),
     col = "red", lwd = 2,type="l",main="n-C33 RF",
     xlim=c(-4,13),ylim = c(0,0.3),axes=F)
lines(Afr.prod.prior.RF$x,Afr.prod.prior.RF$y3, col = "blue", lwd = 1)
#log-scale axis
axis(2,c(0,0.3))
axis(1,log(c(0.01,0.1,1,10,100,1e3,1e4,1e5)))

###CS2: prior vs posterior density of d13C####
#blue is prior, red is posterior
par(mfrow=c(3,3)) #700*800
plot(density(asso.l.mix$BUGSoutput$sims.list$d13C.k[,,1,1]),
     col = "red", lwd = 2,type="l",main="n-C29 GR",xlim=c(-50,-15),ylim = c(0,0.25))
lines(Afr.d13C.prior.GR$x,Afr.d13C.prior.GR$y1, col = "blue", lwd = 1)

plot(density(asso.l.mix$BUGSoutput$sims.list$d13C.k[,,1,2]),
     col = "red", lwd = 2,type="l",main="n-C31 GR",xlim=c(-50,-15),ylim = c(0,0.25))
lines(Afr.d13C.prior.GR$x,Afr.d13C.prior.GR$y2, col = "blue", lwd = 1)

plot(density(asso.l.mix$BUGSoutput$sims.list$d13C.k[,,1,3]),
     col = "red", lwd = 2,type="l",main="n-C33 GR",xlim=c(-50,-15),ylim = c(0,0.25))
lines(Afr.d13C.prior.GR$x,Afr.d13C.prior.GR$y3, col = "blue", lwd = 1)

plot(density(asso.l.mix$BUGSoutput$sims.list$d13C.k[,,2,1]),
     col = "red", lwd = 2,type="l",main="n-C29 SV",xlim=c(-50,-15),ylim = c(0,0.25))
lines(Afr.d13C.prior.SV$x,Afr.d13C.prior.SV$y1, col = "blue", lwd = 1)

plot(density(asso.l.mix$BUGSoutput$sims.list$d13C.k[,,2,2]),
     col = "red", lwd = 2,type="l",main="n-C31 SV",xlim=c(-50,-15),ylim = c(0,0.25))
lines(Afr.d13C.prior.SV$x,Afr.d13C.prior.SV$y2, col = "blue", lwd = 1)

plot(density(asso.l.mix$BUGSoutput$sims.list$d13C.k[,,2,3]),
     col = "red", lwd = 2,type="l",main="n-C33 SV",xlim=c(-50,-15),ylim = c(0,0.25))
lines(Afr.d13C.prior.SV$x,Afr.d13C.prior.SV$y3, col = "blue", lwd = 1)

plot(density(asso.l.mix$BUGSoutput$sims.list$d13C.k[,,3,1]),
     col = "red", lwd = 2,type="l",main="n-C29 RF",xlim=c(-50,-15),ylim = c(0,0.25))
lines(Afr.d13C.prior.RF$x,Afr.d13C.prior.RF$y1, col = "blue", lwd = 1)

plot(density(asso.l.mix$BUGSoutput$sims.list$d13C.k[,,3,2]),
     col = "red", lwd = 2,type="l",main="n-C31 RF",xlim=c(-50,-15),ylim = c(0,0.25))
lines(Afr.d13C.prior.RF$x,Afr.d13C.prior.RF$y2, col = "blue", lwd = 1)

plot(density(asso.l.mix$BUGSoutput$sims.list$d13C.k[,,3,3]),
     col = "red", lwd = 2,type="l",main="n-C33 RF",xlim=c(-50,-15),ylim = c(0,0.25))
lines(Afr.d13C.prior.RF$x,Afr.d13C.prior.RF$y3, col = "blue", lwd = 1)