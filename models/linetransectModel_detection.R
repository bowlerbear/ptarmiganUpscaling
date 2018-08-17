setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/models")
cat("
    model{

    #################
    #Detection model#
    #################
    pi <- 3.141593
    
    # priors for fixed effect parms for half-normal detection parm sigma
    b.df.0 ~ dunif(0,20)        
    b.group.size ~ dnorm(0,0.001)

    for( i in 1:N){
    #linear predictor
    mu.df[i] <- b.df.0 + b.group.size * detectionGroupSize[i] 
    
    # estimate of sd and var, given coeffs above
    sig.df[i] <- exp(mu.df[i])
    sig2.df[i] <- sig.df[i]*sig.df[i]
    
    # effective strip width
    esw[i] <- sqrt(pi * sig2.df[i] / 2) 
    f0[i] <- 1/esw[i] #assuming all detected on the line
    
    # LIKELIHOOD
    # using zeros trick
    y[i] ~ dunif(0,W) 
    L.f0[i] <- exp(-y[i]*y[i] / (2*sig2.df[i])) * 1/esw[i] #y are the distances
    nlogL.f0[i] <-  -log(L.f0[i])
    zeros.dist[i] ~ dpois(nlogL.f0[i])
    }

    #using this model and predicted group size (below), get predicted ESW for each line and year
    for(t in 1:n.Years){
      for(j in 1:n.Lines){
        pred.sig[j,t] <- exp(b.df.0+ b.group.size * log(predGroupSize[j,t]+1)) 
        pred.sig2[j,t] <- pow(pred.sig[j,t],2)
        predESW[j,t] <- sqrt(pi * pred.sig2[j,t] / 2)
      }
    }


    #Group size model
    #priors
    int.gs ~ dnorm(0,0.001)    
    
    #random line effect
    line.sd ~ dunif(0,10)
    line.tau <- pow(line.sd,-2)
    for(j in 1:n.Lines){
      random.gs.line[j] ~ dnorm(0,line.tau)
    }
    
    #random site2 effect
    site2.sd ~ dunif(0,10)
    site2.tau <- pow(site2.sd,-2)
    for(j in 1:n.Sites2){
      random.gs.site2[j] ~ dnorm(0,site2.tau)
    }
    
    #random year effect
    year.sd ~ dunif(0,10)
    year.tau <- pow( year.sd,-2)
    for(t in 1:n.Years){
      random.gs.year[t] ~ dnorm(0, year.tau)
    }
    
    #random line/year effect
    line.year.sd ~ dunif(0,10)
    line.year.tau <- pow( line.year.sd,-2)
    for(j in 1:n.Lines){
      for(t in 1:n.Years){
        random.gs.line.year[j,t] ~ dnorm(0, line.year.tau)
      }
    }
    
    #random site2/year effect
    site2.year.sd ~ dunif(0,10)
    site2.year.tau <- pow(site2.year.sd,-2)
    for(j in 1:n.Sites2){
      for(t in 1:n.Years){
        random.gs.site2.year[j,t] ~ dnorm(0, site2.year.tau)
      }
    }
    
    #Model
    #for each detection, model group size
    for(i in 1:N){
      GroupSize[i] ~ dpois(expGroupSize[i])
      log(expGroupSize[i]) <- int.gs + random.gs.year[detectionYear[i]] + random.gs.line[detectionLine[i]] + 
                              random.gs.site2[detectionSite[i]]+
                              random.gs.line.year[detectionLine[i],detectionYear[i]] +
                              random.gs.site2.year[detectionSite[i],detectionYear[i]]
    }
    
    #using this model, get predicted group size for each line and year
    for(t in 1:n.Years){
      for(j in 1:n.Lines){
        log(predGroupSize[j,t]) <- int.gs + random.gs.year[t] + random.gs.line[j] + random.gs.site2[site2[j]]+
                                    random.gs.line.year[j,t] +
                                    random.gs.site2.year[site2[j],t]
      }
    }

    }
    ",fill=TRUE,file="linetransectModel_detection.txt")
