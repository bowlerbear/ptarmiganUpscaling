cat("
  model{
  # JAGS code for SPARTA model plus random walk prior
  # on the year effect of the state model + intercept + halfcauchy hyperpriors
  
  # State model
  for (i in 1:nsite){ 
    for (t in 1:nyear){
      z[i,t] ~ dbern(muZ[i,t]) 
      logit(muZ[i,t])<- int + a[t] + inprod(beta[],occDM[site[i],])
      #
      # year as a fixed factor and site as a random factor and environ variables
    }
  }   
  
  ### Observation Model
  for(j in 1:nvisit) {
    y[j] ~ dbern(Py[j]) #data is Y
    Py[j]<- z[site[j],year[j]]*p[j] #probability to detect = prob of occ * prob of detection

    #detection model:
    logit(p[j]) <-  alpha.p[year[j]] + dtype.p*L[j] 
    #depends on year and list length
    } 
  

  # Derived parameters
  for (t in 1:nyear) {  
    psi.fs[t] <- sum(z[1:nsite, t])/nsite
  } 

  #Priors 

  # State model priors
    int ~ dnorm(0,0.001)
 
    #years
    for(t in 1:nyear){
      a[t] ~ dnorm(0, tau.a)
    }
    tau.a <- 1/(sd.a * sd.a)
    sd.a ~ dt(0, 1, 1)T(0,) 
    
    #site effects
    for (i in 1:nsite) {
      eta[i] ~ dnorm(0, tau2)       
    } 
    tau2 <- 1/(sigma2 * sigma2) 
    sigma2 ~ dt(0, 1, 1)T(0,) 

    #Observation model priors 
    #year effects
    for (t in 1:nyear) {
      alpha.p[t] ~ dnorm(0, tau.lp)            
    }
    tau.lp <- 1 / (sd.lp * sd.lp)                 
    sd.lp ~ dt(0, 1, 1)T(0,)  
    
    for(i in 1:n.covs){
      beta[i] ~ dnorm(0,0.1)
    }

    #observation model covariates
    dtype.p ~ dnorm(0, 0.01)

  }
    ",fill=TRUE,file="BUGS_sparta_variables.txt")