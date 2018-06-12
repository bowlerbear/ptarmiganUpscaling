setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/models")
# Specify model in BUGS language
sink("standard_occModel.txt")
cat("
    model {
    
    #Model to predict site suitability

    # Priors
    mean.p ~ dunif(0, 1)         # Detection intercept on prob. scale
    int.p <- logit(mean.p)      # Detection intercept
  
    mean.psi ~ dunif(0, 1)       # Occupancy intercept on prob. scale
    int.psi <- logit(mean.psi)     # Occupancy intercept

    #adm effects
    for(i in 1:n.adm){
      random.adm[i] ~ dnorm(0,random.adm.tau)
    }
    random.adm.tau <- pow(random.adm.sd,-2)
    random.adm.sd ~ dunif(0,10)

    # Likelihood
    for (j in 1:nsite) {
      for( t in 1:nyear){
        # True state model for the partially observed true state
        z[j,t] ~ dbern(psi[j,t])      # True occupancy z at site j in year t
        logit(psi[j,t]) <- int.psi + random.adm[adm[j]] 
    
        for (v in 1:nvisit) {
        # Observation model for the actual observations
        y[j,t,v] ~ dbern(p.eff[j,t,v])    # Detection-nondetection at each v
        p.eff[j,t,v] <- z[j,t] * p[j,t,v]  
        logit(p[j,t,v]) <- int.p 
      }
      }
    }

    #get mean occupancy for each site over time
    for(j in 1:nsite){
      meanZ[j] <- mean(z[j,])
    }

    #Model the abundance, when a site is predicted to be suitable, 

    #intercept
    int.d ~ dnorm(0,0.001) 

    for(i in 1:n.covs){
      beta[i] ~ dnorm(0,0.1)
    }

    #random line effect
    line.d.sd ~ dunif(0,10)
    line.d.tau <- pow(line.d.sd,-2)
    for(j in 1:n.Lines){
      random.d.line[j] ~ dnorm(0,line.d.tau)
    }
    
    #random site effect
    site.d.sd ~ dunif(0,10)
    site.d.tau <- pow(site.d.sd,-2)
    for(j in 1:n.Sites){
      random.d.site[j] ~ dnorm(0,site.d.tau)
    }

    #random site2 effect
    site2.d.sd ~ dunif(0,10)
    site2.d.tau <- pow(site2.d.sd,-2)
    for(j in 1:n.Sites2){
      random.d.site2[j] ~ dnorm(0,site2.d.tau)
    }

    #Observation model:
    for(j in 1:nsite){
      for(t in 1:nyear){
        NuIndivs[j,t] ~ dpois(expNuIndivs[j,t])
        expNuIndivs[j,t] <- (predDensity[j,t] * (TransectLength[j,t]/1000 * predESW[j,t]/1000 * 2))
        predDensity[j,t] ~ dpois(meanZ[j] * Density[j,t])  
      }
    }
    
    #State model
    for(j in 1:nsite){
      for(t in 1:nyear){
        log(Density[j,t]) <- int.d + 
                              random.d.line[j] + 
                              #random.d.site2[site2[j]] +
                              inprod(beta[],occDM[j,]) 
      }
    }

    #get average density per site
    for(j in 1:nsite){
      lineDensity[j] <- mean(Density[j,])
    }
    
    }
    ",fill = TRUE)
sink()
