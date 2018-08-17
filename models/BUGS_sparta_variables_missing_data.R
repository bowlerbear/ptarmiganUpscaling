setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/models")
cat("
  model{
  # JAGS code for SPARTA model
  
  # State model
  for (i in 1:nsite){ 
    for (t in 1:nyear){
      z[i,t] ~ dbern(muZ[i,t]) 
      logit(muZ[i,t])<- int.alpha + inprod(beta[],occDM[i,]) 
    }
  }   
  
  ### Observation Model
  for(j in 1:nvisit) {
    y[j] ~ dbern(Py[j]) #data is Y
    Py[j]<- z[site[j],year[j]]*p[j] #probability to detect = prob of occ * prob of detection

    #detection model:
    logit(p[j]) <-  int.d + inprod(beta.det[],occDM[site[j],]) + beta.e * Effort[j]
    } 
  

  #model for missing list length data
  for(j in 1:nvisit){
    L[j] ~ dpois(muL[j])
    L2[j] ~ dpois(muL2[j])
    muL[j] <- intL 
    muL2[j] <- intL2 
    }
    intL ~ dunif(0,100)  
    intL2 ~ dunif(0,100)

  #Priors 

  # State model priors
    ############################
    mean.psi ~ dunif(0, 1)       # Occupancy intercept on prob. scale
    int.alpha <- logit(mean.psi)
 
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

    #adm effects
    for(i in 1:n.adm){
      random.adm[i] ~ dnorm(0,random.adm.tau)
    }
    random.adm.tau <- pow(random.adm.sd,-2)
    random.adm.sd ~ dunif(0,10)

    #Observation model priors
    ##########################
    mean.p ~ dunif(0, 1)         # Detection intercept on prob. scale
    int.d <- logit(mean.p)
    
    beta.sd ~ dnorm(0,0.01)
    beta.e ~ dnorm(0,0.01)

    #year effects
    for (t in 1:nyear) {
      alpha.p[t] ~ dnorm(0, tau.lp)            
    }
    tau.lp <- 1 / (sd.lp * sd.lp)                 
    sd.lp ~ dt(0, 1, 1)T(0,)  
    
    for(i in 1:n.covs){
      beta[i] ~ dnorm(0,0.1)
    }
    
    for(i in 1:n.covs){
      beta.det[i] ~ dnorm(0,0.1)
    }

    #adm effects
    for(i in 1:n.adm){
    random.adm.det[i] ~ dnorm(0,random.adm.det.tau)
    }
    random.adm.det.tau <- pow(random.adm.det.sd,-2)
    random.adm.det.sd ~ dunif(0,10)

  }
    ",fill=TRUE,file="BUGS_sparta_variables_missing.txt")