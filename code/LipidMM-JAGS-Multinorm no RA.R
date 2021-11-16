model {
  
  for(n in 1:N){
    
    #d13C.mix[n] ~ dnorm(d13C.mix.m[n], d13C.mean.pre)
    d13C.mix[n] ~ dnorm(d13C.mix.m[n], 1/d13C.mea.sd[n]^2)
    
    #dt (,,1) is a special case of cauchy distribution
    #RA.mix[n] ~ dt(RA.mix.m[n], RA.mix.m.pre, 1)
    
  }
  
  #RA.mix.m.sd <- 0.03
  
  RA.mix.m.pre ~ dgamma(RA.mix.m.pre.shp, RA.mix.m.pre.rate)
  RA.mix.m.pre.shp = 5
  RA.mix.m.pre.rate = 0.05
  
  for(n in 1:N){
    
    #calculate modeled d13C and RA 
    
    RA.mix.m[n] <-  sum(f.sum.prod_n_i[, n]) / sum(f.sum.prod_n_i)
    
    #d13C of chain n
    d13C.mix.m[n] <- sum(f.prod.d13C_n_i[, n]) / sum(f.sum.prod_n_i[, n])
  }
  
  for(i in 1:I){
    
    for(n in 1:N){
      #sum of n-alkane production for each chain n and component i
      sum.prod_n_i[i, n] <- sum(exp.prod_k[, i, n])
      
      #sum of product of production and d13C for each chain n and component i 
      #for d13C of mixing (weighted mean)
      prod.d13C_n_i[i, n] <- sum(prod.d13C_k[, i, n]) 
      
      #leaf area index added here
      
      #fractional sum of production of chain n and component i
      f.sum.prod_n_i[i, n] <- f[i] * sum.prod_n_i[i, n] 
      
      #fractional sum of d13C weighted production of chain n and component i
      f.prod.d13C_n_i[i, n] <- f[i] * prod.d13C_n_i[i, n] 
    }
  }
  
  for(k in 1:K){
    
    for(i in 1:I){
      
      #multivariat normal produces n products for each chain as d13C for the chains
      #for each gram, it is a new draw, but do we want to record this
      
      d13C.k[k, i, 1:N] ~ dmnorm.vcov(d13C.mu.est[i,1:N], d13C.omega[i,1:N,1:N])
      prod_k[k, i, 1:N] ~ dmnorm.vcov(prod.mu.est[i,1:N], prod.omega[i,1:N,1:N])
      
      for(n in 1:N){ #production of each gram of leaf is modeled indvidually
        exp.prod_k[k, i, n] <- exp (prod_k[k, i, n])
        prod.d13C_k[k, i, n] <- exp.prod_k[k, i, n] * d13C.k[k, i, n]
      }
    }
  }
  
  #d13C.omega and prod.omega are the compiled vcov, which provides indexing pattern for omega
  #mu is the dataframe for means, component as rows and chains as columns
  for(i in 1:I){
    
    d13C.omega[i,1:N,1:N] <- d13C.omega.est[c((N*(i-1)+1):(N*(i-1)+N)),1:N]
    prod.omega[i,1:N,1:N] <- prod.omega.est[c((N*(i-1)+1):(N*(i-1)+N)),1:N]
    
  }
  
  #Mixing fractions of component i: fractions f[i] ~ Dirichlet (1,...,1)
  
  f ~ ddirch(alpha)
  
  #Define parameter alpha for Dirichlet distribution
  
  alpha <- rep(1, I)
  
}