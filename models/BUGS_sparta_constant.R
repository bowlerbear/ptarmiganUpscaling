setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/models")

cat("
  model{
  # JAGS code for SPARTA model 
  # on the year effect of the state model + intercept 
  
  # State model
  for (i in 1:nsite){ 
    for (t in 1:nyear){
      
      z[i,t] ~ dbern(muZ[i,t]) 
      
      logit(muZ[i,t]) <- int.alpha 

    }
  }   
  
  ### Observation Model
  for(j in 1:nvisit) {

    y[j] ~ dbern(Py[j]) #data is Y
    
    Py[j]<- z[site[j],year[j]]*p[j] #probability to detect = prob of occ * prob of detection

    logit(p[j]) <-  int.d 

  } 

  #Priors 
    mean.psi ~ dunif(0, 1)       # Occupancy intercept on prob. scale
    int.alpha <- logit(mean.psi)
    mean.p ~ dunif(0, 1)         # Detection intercept on prob. scale
    int.d <- logit(mean.p)

  }
    ",fill=TRUE,file="BUGS_sparta_constant.txt")