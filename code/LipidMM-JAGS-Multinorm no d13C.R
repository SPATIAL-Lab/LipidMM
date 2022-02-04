model {
  #data evaluation, note that evaluation on d13C.mix[n] is commented out
  for(n in 1:N){
    
    #d13C.mix[n] ~ dnorm(d13C.mix.m[n], 1/d13C.mea.sd[n]^2)
    
    #dt (,,1) is a special case of cauchy distribution
    RA.mix[n] ~ dt(RA.mix.m[n], RA.mix.m.pre, 1)
    
  }
  #set precision for RA.mix
  RA.mix.m.pre ~ dgamma(RA.mix.m.pre.shp, RA.mix.m.pre.rate)
  RA.mix.m.pre.shp = 5
  RA.mix.m.pre.rate = 0.05
  
  for(n in 1:N){
    
    #calculate modeled d13C and RA 
    
    RA.mix.m[n] <-  sum(f.sum.conc_n_i[, n]) / sum(f.sum.conc_n_i)
    
    #d13C of chain n
    d13C.mix.m[n] <- sum(f.conc.d13C_n_i[, n]) / sum(f.sum.conc_n_i[, n])
  }
  
  for(i in 1:I){
    
    for(n in 1:N){
      #sum of n-alkane concentration for each chain n and component i
      sum.conc_n_i[i, n] <- sum(exp.conc_k[, i, n])
      
      #sum of product of concentration and d13C for each chain n and component i 
      #for d13C of mixing (weighted mean)
      conc.d13C_n_i[i, n] <- sum(conc.d13C_k[, i, n]) 
      
      #leaf area index added here
      
      #fractional sum of concentration of chain n and component i
      f.sum.conc_n_i[i, n] <- FLMC[i] * sum.conc_n_i[i, n] 
      
      #fractional sum of d13C weighted concentration of chain n and component i
      f.conc.d13C_n_i[i, n] <- FLMC[i] * conc.d13C_n_i[i, n] 
    }
  }
  
  for(k in 1:K){
    
    for(i in 1:I){
      
      #multivariat normal produces n products for each chain as d13C for the chains
      #for each gram, it is a new draw, but do we want to record this
      
      d13C.k[k, i, 1:N] ~ dmnorm.vcov(d13C.mu.est[i,1:N], d13C.omega[i,1:N,1:N])
      conc_k[k, i, 1:N] ~ dmnorm.vcov(conc.mu.est[i,1:N], conc.omega[i,1:N,1:N])
      
      for(n in 1:N){ #concentration of each gram of leaf is modeled indvidually
        exp.conc_k[k, i, n] <- exp (conc_k[k, i, n])
        conc.d13C_k[k, i, n] <- exp.conc_k[k, i, n] * d13C.k[k, i, n]
      }
    }
  }
  
  #d13C.omega and conc.omega are the compiled vcov, which provides indexing pattern for omega
  #mu is the dataframe for means, component as rows and chains as columns
  for(i in 1:I){
    
    d13C.omega[i,1:N,1:N] <- d13C.omega.est[c((N*(i-1)+1):(N*(i-1)+N)),1:N]
    conc.omega[i,1:N,1:N] <- conc.omega.est[c((N*(i-1)+1):(N*(i-1)+N)),1:N]
    
  }
  
  #Mixing fractions of component i: fractions FLMC[i] ~ Dirichlet (1,...,1)
  
  FLMC ~ ddirch(alpha)
  
  #Define parameter alpha for Dirichlet distribution
  
  alpha <- rep(1, I)
  
}