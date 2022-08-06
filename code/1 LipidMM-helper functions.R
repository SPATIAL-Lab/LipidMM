library(MASS)
library(OpenMx)

####function 1: calculate FSCn (fractional source contribution to each chain n)
#or chain specific mixing ratios from MCMC output f.sum.conc_n_i###
FSC<-function(MCMC.res){
  dim.MCMC <- dim(MCMC.res)
  #typically, the first element is # of interation
  #second element is source
  #third element is chain
  MCMC.sums<-array(c(0),dim=c(dim.MCMC[1],dim.MCMC[3]))
  for(i in 1:dim.MCMC[2]){
    for(j in 1:dim.MCMC[3]){#extract individual chains from source
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
#can be used in d2H, d13C and log transformed concentration

pri.multi.norm.den <- function(X.min,X.max,mu,vcov){
  require(OpenMx)
  X.range <- X.max - X.min
  X.interval <- X.range/512
  X.multinorm <- seq(from = X.min, to = X.max, by = X.interval)
  dim.vcov <- dim(vcov)
  
  #initiate density matrix
  density <- as.data.frame(matrix(0, nrow=length(X.multinorm), ncol = dim.vcov[1]))
  
  #initiate index matrix
  index <- as.data.frame(matrix(nrow=dim.vcov[1], ncol = 2))
  for (i in 1: dim.vcov[1]){
    index[1,i] <- X.min
    index[2,i] <- X.max
  }
  
  for (i in 1: dim.vcov[1]){
    for(j in 2:length(X.multinorm)){
      #assigning new values
      index[1,i] <- X.multinorm[j-1]
      index[2,i] <- X.multinorm[j]
      #integrating one interval for density[i]
      density[i,j-1] <- omxMnor(vcov,mu,index[1,],index[2,])
    }
    #revert to the default values
    index[1,i] <- X.min
    index[2,i] <- X.max
  }
  
  #scaling density values to 1
  density.sc <- density[,]/X.interval

  #monitoring the density scaling approximation
  #the sum should be close to 1
  
  for (i in 1: dim.vcov[1]){
    print(sum(density[i,]))
  }

  results<-list(x = (X.multinorm[1:512]+ X.interval/2), y = density.sc)
  return(results)
}

#Function 3: making contour plots for FLMCs
#note that this function is optimized for 3 chains
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

#Function 4: calculate MAP and CI% HDI of MCMC result
MAP_HDI <- function(MCMC.res, CI){
  require(bayestestR)
  require(KernSmooth)
  MCMC.MAP <- map_estimate(MCMC.res, method = "KernSmooth")
  MCMC.HDI <- hdi(MCMC.res, ci = CI)
  MCMC.HDIL <- MCMC.HDI[1,2]
  MCMC.HDIH <- MCMC.HDI[1,3]
  return(c(MCMC.MAP, MCMC.HDIL, MCMC.HDIH))
}

#Function 5: plotting points with error bars
PlotPE <- function(x, y, error.y, col){
  if(length(x) != length(y)){
    warning("Error in plot.xy, x and y have different lengths")
  }
  if(length(x) != length(error.y)){
    warning("Error in plot.xy, x and error.y have different lengths")
  }
  if(length(y) != length(error.y)){
    warning("Error in plot.xy, y and error.y have different lengths")
  }
  if(is.null(col)){
    warning("Please specify plotting color")
  }
  
  n <- length(x)

  for(i in 1:n){
    arrows(x[i], y[i], x[i], y[i] + error.y[i], length = 0.03, angle = 90, col = col)
    arrows(x[i], y[i], x[i], y[i] - error.y[i], length = 0.03, angle = 90, col = col)
  }
  points (x, y, pch = 16, col = col)
  
}
