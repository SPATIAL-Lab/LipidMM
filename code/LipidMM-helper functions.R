library(coda)
library(lattice)
library(rjags)
library(R2jags)
library(msm)
library(ggmcmc)
library(MASS)
library(EnvStats)
library(OpenMx)
library(bayestestR)

####function 1: calculate chain specific mixing ratios from MCMC output f.sum.prod_n_i###
chain.spec.ratio<-function(MCMC.res){
  dim.MCMC <- dim(MCMC.res)
  #typically, the first element is # of interation
  #second element is end member
  #third element is chain
  MCMC.sums<-array(c(0),dim=c(dim.MCMC[1],dim.MCMC[3]))
  for(i in 1:dim.MCMC[2]){
    for(j in 1:dim.MCMC[3]){#extract individual chains from member
      MCMC.sums[,j]<-MCMC.sums[,j] + MCMC.res[,i,j]
    }
    #sums of chains
  }
  MCMC.ratios <- array(c(0),dim = dim.MCMC)
  for(i in 1:dim.MCMC[2]){
    for(j in 1:dim.MCMC[3]){
      #density has x and y 
      MCMC.ratios[,i,j] <- MCMC.res[,i,j]/MCMC.sums[,j]
    }
  }
  return(MCMC.ratios)
}

####function 2: prior multivariate normal density###
#can be used in both d13C and log transformed concentration
#note that this function is optimized for a 3 chain 3 endmember mixing regime

pri.multi.norm.den <- function(X.min,X.max,mu,vcov){
  require(OpenMx)
  X.range <- X.max - X.min
  X.interval<-X.range/511
  X.multinorm <- seq(from = X.min, to = X.max, by =X.interval)
  
  density.1<-rep(0,length(X.multinorm))
  density.2<-rep(0,length(X.multinorm))
  density.3<-rep(0,length(X.multinorm))
  
  for(i in 2:length(X.multinorm)){
    #integrating one interval for density[i]
    density.1[i] <- omxMnor(vcov,mu,c(X.multinorm[i-1],X.min,X.min),
                            c(X.multinorm[i],X.max,X.max))
    density.2[i] <- omxMnor(vcov,mu,c(X.min,X.multinorm[i-1],X.min),
                            c(X.max,X.multinorm[i],X.max))
    density.3[i] <- omxMnor(vcov,mu,c(X.min,X.min,X.multinorm[i-1]),
                            c(X.max,X.max,X.multinorm[i]))
    
  }
  #scaling density values to 1
  density.1.sc<-density.1/X.interval
  density.2.sc<-density.2/X.interval
  density.3.sc<-density.3/X.interval
  #monitoring the density scaling approximation
  #the sum should be close to 1
  print(sum(density.1))
  print(sum(density.2))
  print(sum(density.3))
  
  results<-list(x = X.multinorm, y1 = density.1.sc, y2 = density.2.sc, y3 = density.3.sc)
  return(results)
}

#Function 3: making contour plots
contour.fn<-function(MCMC.f){
  require(MASS)
  f1 <- MCMC.f[,1]
  f2 <- MCMC.f[,2]
  f3 <- MCMC.f[,3]
  contour.1.2 <- kde2d(f1, f2, n = 64)
  contour.1.3 <- kde2d(f1, f3, n = 64)
  contour.2.3 <- kde2d(f2, f3, n = 64)
  results<-list(f1, f2, f3, contour.1.2, contour.1.3, contour.2.3)
  return(results)
}
