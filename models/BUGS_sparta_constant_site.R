setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/models")

cat("
  model{
  # JAGS code for SPARTA model 
  # on the year effect of the state model + intercept 
  
  # State model
  for (i in 1:nsite){ 
      z[i] ~ dbern(muZ[i]) 
      logit(muZ[i]) <- int.alpha
    }
  
  ### Observation Model
  for(j in 1:nvisit) {

    y[j] ~ dbern(Py[j]) #data is Y
    
    Py[j]<- z[site[j]]*p[j] #probability to detect = prob of occ * prob of detection

    logit(p[j]) <-  int.d + eta.s.y[siteyear[j]]

  } 

  #Priors 
    mean.psi ~ dunif(0, 1)       # Occupancy intercept on prob. scale
    int.alpha <- logit(mean.psi)
    mean.p ~ dunif(0, 1)         # Detection intercept on prob. scale
    int.d <- logit(mean.p)

    #random year effects
    for(t in 1:nyear){
      a[t] ~ dnorm(0, tau.y)
    }
    tau.y <- 1/(sd.y * sd.y)
    sd.y ~ dunif(0,10) 
    
    #site effects
    for (i in 1:nsite) {
      eta[i] ~ dnorm(0, tau.s)       
    } 
    tau.s <- 1/(sd.s * sd.s) 
    sd.s ~ dunif(0,10)

    #site/year effects
    for (i in 1:nsiteyear) {
      eta.s.y[i] ~ dnorm(0, tau.s.y)       
    } 
    tau.s.y <- 1/(sd.s.y * sd.s.y) 
    sd.s.y ~ dunif(0,10)


  }
    ",fill=TRUE,file="BUGS_sparta_constant_site.txt")