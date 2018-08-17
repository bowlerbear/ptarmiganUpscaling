setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/models")
# Specify model in BUGS language
sink("hurdleModel.txt")
cat("
    model {
    
    #Model to predict site suitability
    ##################################

    # Priors
    mean.p ~ dunif(0, 1)         # Detection intercept on prob. scale
    int.p <- logit(mean.p)      # Detection intercept
  
    mean.psi ~ dunif(0, 1)       # Occupancy intercept on prob. scale
    int.psi <- logit(mean.psi)     # Occupancy intercept


    #random site effect on occupancy
    for(i in 1:nsite){
      random.site[i] ~ dnorm(0,random.site.tau)
    }
    random.site.tau <- pow(random.site.sd,-2)
    random.site.sd ~ dunif(0,10)

    #adm effects on occupancy
    for(i in 1:n.adm){
      random.adm[i] ~ dnorm(0,random.adm.tau)
    }
    random.adm.tau <- pow(random.adm.sd,-2)
    random.adm.sd ~ dunif(0,10)

    for(i in 1:n.adm2){
      random.adm2[i] ~ dnorm(0,random.adm2.tau)
    }
    random.adm2.tau <- pow(random.adm2.sd,-2)
    random.adm2.sd ~ dunif(0,10)

    #random year effect
    for(i in 1:nyear){
      random.o.year[i] ~ dnorm(0,random.o.year.tau)
    }
    random.o.year.tau <- pow(random.o.year.sd,-2)
    random.o.year.sd ~ dunif(0,10)


    #adm effects on detection:

    for(i in 1:n.adm){
      random.det.adm[i] ~ dnorm(0,random.det.adm.tau)
    }
    random.det.adm.tau <- pow(random.det.adm.sd,-2)
    random.det.adm.sd ~ dunif(0,10)
    
    for(i in 1:n.adm2){
    random.det.adm2[i] ~ dnorm(0,random.det.adm2.tau)
    }
    random.det.adm2.tau <- pow(random.det.adm2.sd,-2)
    random.det.adm2.sd ~ dunif(0,10)

    #random year effect
    for(i in 1:nyear){
      random.det.year[i] ~ dnorm(0,random.det.year.tau)
    }
    random.det.year.tau <- pow(random.det.year.sd,-2)
    random.det.year.sd ~ dunif(0,10)

    #other covariates
    beta.e  ~ dnorm(0,0.01)
    beta.e2  ~ dnorm(0,0.01)


    # State model
    for (i in 1:nsite){ 

      for (t in 1:nyear){
    
      z[i,t] ~ dbern(muZ[i,t]) 
    
      logit(muZ[i,t]) <- int.psi + random.site[i] + random.adm[adm[i]] +
                          random.adm2[adm2[i]] + random.o.year[t] 

      }
    }   
    
    ### Observation Model
    for(j in 1:nvisit) {
    
      y[j] ~ dbern(Py[j]) #data is Y
    
      Py[j]<- z[site[j],year[j]]*p[j] #probability to detect = prob of occ * prob of detection
    
      #detection model:
      logit(p[j]) <-  int.p + beta.e * Effort[j] + beta.e2 * Effort2[j] + Density[site[j]]

    }


    #get mean occupancy for each site over time
    for(j in 1:nsite){
      meanZ[j] <- mean(z[j,])
    }

    #Model the abundance at line transect grids 
    ##############################################################

    #intercept
    int.d ~ dnorm(0,0.001)    
    
    #adm effects on abundance
    for(i in 1:n.adm){
    random.adm[i] ~ dnorm(0,random.adm.tau)
    }
    random.adm.tau <- pow(random.adm.sd,-2)
    random.adm.sd ~ dunif(0,10)
    
    for(i in 1:n.adm2){
    random.adm2[i] ~ dnorm(0,random.adm2.tau)
    }
    random.adm2.tau <- pow(random.adm2.sd,-2)
    random.adm2.sd ~ dunif(0,10)
    
    #adm effects on detection
    for(i in 1:n.adm){
    random.det.adm[i] ~ dnorm(0,random.det.adm.tau)
    }
    random.det.adm.tau <- pow(random.det.adm.sd,-2)
    random.det.adm.sd ~ dunif(0,10)
    
    for(i in 1:n.adm2){
    random.det.adm2[i] ~ dnorm(0,random.det.adm2.tau)
    }
    random.det.adm2.tau <- pow(random.det.adm2.sd,-2)
    random.det.adm2.sd ~ dunif(0,10)

    
    #Work out the fraction of the both that was surveyed
    ESW.tauLT <- pow(ESW.sdLT,-2)
    for(j in 1:nsite){
      for(t in 1:nyear){
        predESW[j,t] ~ dnorm(ESW.meanLT[j,t],ESW.tauLT[j,t])
        surveyArea[j,t] <- TransectLengthLT[j,t]/1000 * predESW[j,t]/1000 * 2
        surveyFraction[j,t] <- surveyArea[j,t]/25
      }
    }

    #Observation model:
    for(j in 1:nsite){
      for(t in 1:nyear){
        NuIndivsLT[j,t] ~ dpois(expNuIndivs[j,t])
        expNuIndivs[j,t] <- Density[j,t] * surveyFraction[j,t]
      }
    }
    
    #State model to observed data
    for(j in 1:nsite){
      for(t in 1:nyear){
          log(Density[j,t]) <- int.d + random.site[j] + random.adm[adm[j]] + 
                                random.adm2[adm2[j]] + + random.line.year[t] 
      }
    }

    for(j in 1:nsite){
      meanDensity[j] <- mean(Density[j,])
    }


    #Predict realized abundance, across all sites, i.e. suitability * abundance
    for(j in 1:nsite){
      realAbund[j] ~ dpois(meanDensity[j] * meanZ[j])
    }

    totAbund<-sum(realAbund)

    }
    ",fill = TRUE)
sink()
