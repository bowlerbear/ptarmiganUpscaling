# Specify model in BUGS language
sink("models/hurdleModel.txt")
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

    #for(i in 1:n.covs){
    #  beta[i] ~ dnorm(0,0.1)
    #}
    
    #for(i in 1:n.covs){
    #  beta.det[i] ~ dnorm(0,0.1)
    #}


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
      logit(p[j]) <-  int.p + beta.e * Effort[j] + beta.e2 * Effort2[j] + random.det.year[year[j]]
    #random.det.adm[det.adm[j]]
    #random.det.adm2[det.adm2[j]]
    #inprod(beta.det[],occDM[site[j],])

    }


    #get mean occupancy for each site over time
    for(j in 1:nsite){
      meanZ[j] <- mean(z[j,])
    }

    #Model the abundance at line transect grids 
    ##############################################################

    #intercept
    int.d ~ dnorm(0,0.001)    
    
    #for(i in 1:n.covs){
    #  beta.abundance[i] ~ dnorm(0,0.1)
    #} 


    #site effects
    for(i in 1:n.LinesLT){
      random.line.site[i] ~ dnorm(0,random.line.site.tau)
    }
    random.line.site.tau <- pow(random.line.site.sd,-2)
    random.line.site.sd ~ dunif(0,10)

    #adm effects
    for(i in 1:n.SitesLT){
      random.line.adm[i] ~ dnorm(0,random.line.adm.tau)
    }
    random.line.adm.tau <- pow(random.line.adm.sd,-2)
    random.line.adm.sd ~ dunif(0,10)
    
    for(i in 1:n.Sites2LT){
      random.line.adm2[i] ~ dnorm(0,random.line.adm2.tau)
    }
    random.line.adm2.tau <- pow(random.line.adm2.sd,-2)
    random.line.adm2.sd ~ dunif(0,10)

    for(i in 1:n.YearsLT){
      random.line.year[i] ~ dnorm(0,random.line.year.tau)
    }
    random.line.year.tau <- pow(random.line.year.sd,-2)
    random.line.year.sd ~ dunif(0,10)

    #Work out the fraction of the both that was surveyed
    ESW.tauLT <- pow(ESW.sdLT,-2)
    for(j in 1:n.LinesLT){
      for(t in 1:n.YearsLT){
        predESW[j,t] ~ dnorm(ESW.meanLT[j,t],ESW.tauLT[j,t])
        surveyArea[j,t] <- TransectLengthLT[j,t]/1000 * predESW[j,t]/1000 * 2
        surveyFraction[j,t] <- surveyArea[j,t]/25
      }
    }

    #Observation model:
    for(j in 1:n.LinesLT){
      for(t in 1:n.YearsLT){
        NuIndivsLT[j,t] ~ dpois(expNuIndivs[j,t])
        expNuIndivs[j,t] <- Density[j,t] * surveyFraction[j,t]
      }
    }
    
    #State model to observed data
    for(j in 1:n.LinesLT){
      for(t in 1:n.YearsLT){
          log(Density[j,t]) <- int.d + random.line.site[j] + random.line.adm[siteLT[j]] + 
                                random.line.adm2[site2LT[j]] + random.line.year[t]
                              #inprod(beta.abundance[],occDM_abundance[j,])
      }
    }
    
    #Predict density for all sites
    for(j in 1:nsite){
      meanDensity[j] <- exp(int.d) 
    }
    #inprod(beta.occupancy[],occDM[j,])

    #Predict realized abundance, across all sites, i.e. suitability * abundance
    for(j in 1:nsite){
      realAbund[j] ~ dpois(meanDensity[j] * meanZ[j])
    }

    totAbund<-sum(realAbund)

    }
    ",fill = TRUE)
sink()
