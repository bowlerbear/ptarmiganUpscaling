################################################################################

#Generate data for all sites:
library(VGAM)
library(boot)

output<-data.frame()

for(j in 1:1000){

#import siteInfo
lines<-1:100
site2<-rep(1:10,each=10)
site<-rep(1:5,each=20)

#predict occupancy for all sites:
lineOccupancy <- rnorm(length(lines),0,1)
site2Occupancy <- rnorm(length(site2),0,1) 
siteOccupancy <- rnorm(length(site),0,1)

#define average
int.o<-logit(0.1)

logit.occupancy<-as.numeric()
for(i in 1:length(lines)){
  logit.occupancy[i] <- int.o + lineOccupancy[i] + site2Occupancy[site2[i]] + siteOccupancy[site[i]]
}
occupancy <- inv.logit(logit.occupancy)

#predict density for all sites:
lineDensity <- rnorm(length(lines),0,1)
site2Abundance <- rnorm(length(site2),0,1)
siteAbundance <- rnorm(length(site),0,1)
  
#predict abundance for all sites:

#define average
int.a <- log(10)

abundance<-as.numeric()
for(i in 1:length(lines)){
  abundance[i] <- int.a + lineDensity[i] + site2Abundance[site2[i]] + siteAbundance[site[i]]
}  
abundance <- exp(abundance)

#realized abundance
realAbundance<-as.numeric()
for(i in 1:length(lines)){
  realAbundance[i]<- rpois(1,abundance[i] * occupancy[i])  
}
sum(realAbundance)

#combine data
df <- data.frame(realAbundance,lines,site2,site)
df$PA <- ifelse(df$realAbundance>0,1,0)

################################################################################

#sample data
if(type=="basic"){

}else if(type=="line"){
  
  #loose 10% at random of lines
  loose<-ceiling(runif(0.2*nrow(df),0,nrow(df)))
  df$realAbundance[loose]<-NA
  df$PA[loose]<-NA

}else if (type=="site2"){
  
  #loose one site at random
  looseSite2<-ceiling(runif(0.2*length(unique(df$site2)),0,length(unique(df$site2))))
  df$realAbundance[which(df$site2%in%looseSite2)]<-NA
  df$PA[which(df$site2%in%looseSite2)]<-NA
  
  
}else if(type=="site"){
  
  #loose one site at random
  looseSite<-ceiling(runif(0.2*length(unique(df$site)),0,length(unique(df$site))))
  df$realAbundance[which(df$site%in%looseSite)]<-NA
  df$PA[which(df$site%in%looseSite)]<-NA
}

################################################################################
#compile data for model

bugs.data <- list(n.a.lines = nrow(subset(df,!is.na(df$realAbundance))),
                  n.o.lines = nrow(subset(df,!is.na(df$PA))),
                  abund = df$realAbundance,
                  PA = df$PA,
                  lines = df$lines,
                  site2 = df$site2,
                  n.site2 = length(unique(df$site2)),
                  site = df$site,
                  n.site = length(unique(df$site)))

#################################################################################
#Fit model:

setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/models")
# Specify model in BUGS language
sink("simulations_hurdleModel.txt")
cat("
    model {
    
    #Model to predict site suitability
    ##################################

    # Priors
    mean.psi ~ dunif(0, 1)       # Occupancy intercept on prob. scale
    int.psi <- logit(mean.psi)     # Occupancy intercept


    #random line effect on occupancy
    for(i in 1:n.o.lines){
      random.o.line[i] ~ dnorm(0,random.o.line.tau)
    }
    random.o.line.tau <- pow(random.o.line.sd,-2)
    random.o.line.sd ~ dunif(0,10)

    #site effects on occupancy
    for(i in 1:n.site){
      random.o.site[i] ~ dnorm(0,random.o.site.tau)
    }
    random.o.site.tau <- pow(random.o.site.sd,-2)
    random.o.site.sd ~ dunif(0,10)

    for(i in 1:n.site2){
      random.o.site2[i] ~ dnorm(0,random.o.site2.tau)
    }
    random.o.site2.tau <- pow(random.o.site2.sd,-2)
    random.o.site2.sd ~ dunif(0,10)

   
    # State model
    for (i in 1:n.o.lines){ 
      PA[i] ~ dbern(muZ[i]) 
      logit(muZ[i]) <- int.psi + random.o.line[i] + random.o.site[site[i]] + random.o.site2[site2[i]] 
    }   
    
   
    #Model the abundance at line transect grids 
    ##############################################################

    #intercept
    mean.d ~ dunif(0,100)
    int.d <- log(mean.d)    
    
    #random line effect on abundance
    for(i in 1:n.a.lines){
      random.a.line[i] ~ dnorm(0,random.a.line.tau)
    }
    random.a.line.tau <- pow(random.a.line.sd,-2)
    random.a.line.sd ~ dunif(0,10)

    #site effects on occupancy
    for(i in 1:n.site){
      random.a.site[i] ~ dnorm(0,random.a.site.tau)
    }
    random.a.site.tau <- pow(random.a.site.sd,-2)
    random.a.site.sd ~ dunif(0,10)

    for(i in 1:n.site2){
      random.a.site2[i] ~ dnorm(0,random.a.site2.tau)
    }
    random.a.site2.tau <- pow(random.a.site2.sd,-2)
    random.a.site2.sd ~ dunif(0,10)

    #State model to observed data
    for(j in 1:n.a.lines){
          abund[j] ~ dpois(Density[j])
          log(Density[j]) <- int.d +  random.a.line[j] + random.a.site[site[j]] + random.a.site2[site2[j]] 
    }
    

    #Predict realized abundance, across all sites, i.e. suitability * abundance
    for(j in 1:n.a.lines){
      realAbund[j] ~ dpois(Density[j] * muZ[j])
    }

    totAbund<-sum(realAbund)

    }
    ",fill = TRUE)
sink()

#######################################################################################

#run model

library(rjags)
library(jagsUI)

params=c("totAbund","mean.psi","mean.d",
         "random.o.line.sd","random.o.site.sd","random.o.site2.sd",
         "random.a.line.sd","random.a.site.sd","random.a.site2.sd")

#run models
out1 <- jags(bugs.data, inits=NULL, parameters.to.save = params,
             "simulations_hurdleModel.txt", n.thin=5,
             n.chains=3, n.burnin=5000,n.iter=10000,parallel = TRUE)

out1DF<-data.frame(out1$summary)
out1DF$Param<-row.names(out1DF)
out1DF$simNu<-j
output<-rbind(output,out1DF)

}

#############################################################################################
#sim 1 - 100 sites
setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/model-outputs/sim1")
save(output,file="output_Basic.RData")

#loose 10% of lines at random
type="line"
save(output,file="output_Line.RData")

#loose 10% of site2 at random
type="site2"
save(output,file="output_Site2.RData")

#loose 10% of site at random
type="site"
save(output,file="output_Site.RData")

###########################################################################################

setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/model-outputs")

#plot results
files<-c("output_Basic.RData","output_Site.RData","output_Site2.RData","output_Line.RData")
library(plyr)
output<-ldply(files,function(x){
  load(x)
  output$File<-x
  return(output)
})
output$Type <- sapply(output$File,function(x)strsplit(as.character(x),"_")[[1]][2])
output$Type <- sapply(output$Type,function(x)strsplit(as.character(x),"\\.")[[1]][1])

#plot totAbund

#effect on the estimate
library(ggplot2)
ggplot(subset(output,Param=="totAbund"))+
  geom_boxplot(aes(x=Type,y=mean),outlier.shape=NA)+
  ylim(0,2000)
  
#effect on the standard error
ggplot(subset(output,Param=="totAbund"))+
  geom_boxplot(aes(x=Type,y=sd),outlier.shape=NA)+
  ylim(0,200)

#line level sampling most important...

#########################################################################################
#########################################################################################
#assume different data on occupancy and abundance

#Generate data for all sites:
library(VGAM)
library(boot)

output<-data.frame()

for(j in 1:1000){
  
  #import siteInfo
  lines<-1:20
  site2<-rep(1:5,each=4)
  site<-rep(1:2,each=10)
  
  #predict occupancy for all sites:
  lineOccupancy <- rnorm(length(lines),0,1)
  site2Occupancy <- rnorm(length(site2),0,1) 
  siteOccupancy <- rnorm(length(site),0,1)
  
  #define average
  int.o<-cloglog(0.5)
  
  logit.occupancy<-as.numeric()
  for(i in 1:length(lines)){
    logit.occupancy[i] <- int.o + lineOccupancy[i] + site2Occupancy[site2[i]] + siteOccupancy[site[i]]
  }
  occupancy <- cloglog(logit.occupancy,inverse=T)
  
  #predict density for all sites:
  lineDensity <- rnorm(length(lines),0,1)
  site2Abundance <- rnorm(length(site2),0,1)
  siteAbundance <- rnorm(length(site),0,1)
  
  #predict abundance for all sites:
  
  #define average
  int.a <- log(10)
  
  abundance<-as.numeric()
  for(i in 1:length(lines)){
    abundance[i] <- int.a + lineDensity[i] + site2Abundance[site2[i]] + siteAbundance[site[i]]
  }  
  abundance <- exp(abundance)
  
  #realized abundance
  realAbundance<-as.numeric()
  for(i in 1:length(lines)){
    realAbundance[i]<- rpois(1,abundance[i] * occupancy[i])  
  }
  sum(realAbundance)
  
  #combine data
  df <- data.frame(realAbundance,lines,site2,site)
  df$PA <- ifelse(df$realAbundance>0,1,0)
  
  ################################################################################
  
  #sample data
  if(type=="basic"){
    
  }else if(type=="line.a"){
    
    #loose 10% at random of lines
    loose<-ceiling(runif(0.2*nrow(df),0,nrow(df)))
    df$realAbundance[loose]<-NA
    
  }else if (type=="site2.a"){
    
    #loose one site at random
    looseSite2<-ceiling(runif(0.2*length(unique(df$site2)),0,length(unique(df$site2))))
    df$realAbundance[which(df$site2%in%looseSite2)]<-NA

  }else if(type=="site.a"){
    
    #loose one site at random
    looseSite<-ceiling(runif(0.2*length(unique(df$site)),0,length(unique(df$site))))
    df$realAbundance[which(df$site%in%looseSite)]<-NA
  
    }else if(type=="line.o"){
    
    #loose 10% at random of lines
    loose<-ceiling(runif(0.2*nrow(df),0,nrow(df)))
    df$PA[loose]<-NA
    
  }else if (type=="site2.o"){
    
    #loose one site at random
    looseSite2<-ceiling(runif(0.2*length(unique(df$site2)),0,length(unique(df$site2))))
    df$PA[which(df$site2%in%looseSite2)]<-NA
    
  }else if(type=="site.o"){
    
    #loose one site at random
    looseSite<-ceiling(runif(0.2*length(unique(df$site)),0,length(unique(df$site))))
    df$PA[which(df$site%in%looseSite)]<-NA
  }
  
  
  ################################################################################
  #compile data for model
  
  bugs.data <- list(n.a.lines = nrow(df),
                    n.o.lines = nrow(df),
                    abund = df$realAbundance,
                    PA = df$PA,
                    lines = df$lines,
                    site2 = df$site2,
                    n.site2 = length(unique(df$site2)),
                    site = df$site,
                    n.site = length(unique(df$site)))
  
  #################################################################################
  #Fit model:
  
  setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/models")
  # Specify model in BUGS language
  sink("simulations_hurdleModel.txt")
  cat("
      model {
      
      #Model to predict site suitability
      ##################################
      
      # Priors
      mean.psi ~ dunif(0, 1)       # Occupancy intercept on prob. scale
      int.psi <- cloglog(mean.psi)     # Occupancy intercept
      
      
      #random line effect on occupancy
      for(i in 1:n.o.lines){
      random.o.line[i] ~ dnorm(0,random.o.line.tau)
      }
      random.o.line.tau <- pow(random.o.line.sd,-2)
      random.o.line.sd ~ dunif(0,10)
      
      #site effects on occupancy
      for(i in 1:n.site){
      random.o.site[i] ~ dnorm(0,random.o.site.tau)
      }
      random.o.site.tau <- pow(random.o.site.sd,-2)
      random.o.site.sd ~ dunif(0,10)
      
      for(i in 1:n.site2){
      random.o.site2[i] ~ dnorm(0,random.o.site2.tau)
      }
      random.o.site2.tau <- pow(random.o.site2.sd,-2)
      random.o.site2.sd ~ dunif(0,10)
      
      
      # State model
      for (i in 1:n.o.lines){ 
      PA[i] ~ dbern(muZ[i]) 
        cloglog(muZ[i]) <- int.psi + random.o.line[i] + random.o.site[site[i]] + random.o.site2[site2[i]] 
      }   
      
      
      #Model the abundance at line transect grids 
      ##############################################################
      
      #intercept
      mean.d ~ dunif(0,100)
      int.d <- log(mean.d)    
      
      #random line effect on abundance
      for(i in 1:n.a.lines){
      random.a.line[i] ~ dnorm(0,random.a.line.tau)
      }
      random.a.line.tau <- pow(random.a.line.sd,-2)
      random.a.line.sd ~ dunif(0,10)
      
      #site effects on occupancy
      for(i in 1:n.site){
      random.a.site[i] ~ dnorm(0,random.a.site.tau)
      }
      random.a.site.tau <- pow(random.a.site.sd,-2)
      random.a.site.sd ~ dunif(0,10)
      
      for(i in 1:n.site2){
      random.a.site2[i] ~ dnorm(0,random.a.site2.tau)
      }
      random.a.site2.tau <- pow(random.a.site2.sd,-2)
      random.a.site2.sd ~ dunif(0,10)
      
      #State model to observed data
      for(j in 1:n.a.lines){
      abund[j] ~ dpois(Density[j])
      log(Density[j]) <- int.d +  random.a.line[j] + random.a.site[site[j]] + random.a.site2[site2[j]] 
      }
      
      #Predict realized abundance, across all sites, i.e. suitability * abundance
      for(j in 1:n.a.lines){
      realAbund[j] ~ dpois(Density[j] * muZ[j])
      }
      
      totAbund<-sum(realAbund)
      
      }
      ",fill = TRUE)
  sink()
  
  #######################################################################################
  
  #run model
  
  library(rjags)
  library(jagsUI)
  
  params=c("totAbund","mean.psi","mean.d",
           "random.o.line.sd","random.o.site.sd","random.o.site2.sd",
           "random.a.line.sd","random.a.site.sd","random.a.site2.sd")
  
  #run models
  out1 <- jags(bugs.data, inits=NULL, parameters.to.save = params,
               "simulations_hurdleModel.txt", n.thin=5,
               n.chains=3, n.burnin=5000,n.iter=10000,parallel = TRUE)
  
  out1DF<-data.frame(out1$summary)
  out1DF$Param<-row.names(out1DF)
  out1DF$simNu<-j
  output<-rbind(output,out1DF)
  
}

##########################################################################################

#running
setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/model-outputs/sim2")

type="basic"
save(output,file="output_Basic.RData")

type="line.a"
save(output,file="output_Line.a.RData")

type="line.o"
save(output,file="output_Line.o.RData")

type="site2.a"
save(output,file="output_Site2.a.RData")

type="site2.o"
save(output,file="output_Site2.o.RData")

#plot results
files<-list.files()
library(plyr)
output<-ldply(files,function(x){
  load(x)
  output$File<-x
  return(output)
})
output$Type <- sapply(output$File,function(x)strsplit(as.character(x),"_")[[1]][2])
output$Type <- sapply(output$Type,function(x)strsplit(as.character(x),"\\.RData")[[1]][1])

#plot totAbund

#effect on the estimate
library(ggplot2)
ggplot(subset(output,Param=="totAbund"))+
  geom_boxplot(aes(x=Type,y=mean),outlier.shape=NA)+
  ylim(0,2000)

#effect on the standard error
ggplot(subset(output,Param=="totAbund"))+
  geom_boxplot(aes(x=Type,y=sd),outlier.shape=NA)+
  ylim(0,200000)

#loosing abundance data is more important, especially at the level of site2

#######################################################################################