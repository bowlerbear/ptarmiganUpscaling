#load libraries
library(unmarked)
library(rjags)
library(AHMbook)
library(R2WinBUGS)
library(jagsUI)
library(ggplot2)

# JAGS setting b/c otherwise JAGS cannot build a sampler, rec. by M. Plummer
set.factory("bugs::Conjugate", FALSE, type="sampler")

# Default parameters
ni <- 6000   ;   nb <- 2000   ;   nt <- 2   ;   nc <- 3


#model checking
my_ggs_density<-function(D){
  library(ggplot2)
  f <- ggplot(D, aes(x = value))
  f <- f + geom_density(alpha = 0.3,fill="steelblue")+facet_wrap(~Chain)
  return(f)
}

my_ggs_ppmean<-function(D,outcome){
  ppM <- D %>% dplyr::group_by(Iteration) %>% dplyr::summarize(m = mean(value))
  f <- ggplot()+ geom_density(data=ppM,aes(x=m),alpha=0.5,
                              fill="steelblue",
                              colour="steelblue") +  
    geom_vline(xintercept = mean(outcome,na.rm=T))
  return(f)
}


#http://xavier-fim.net/packages/ggmcmc/
#library(ggmcmc)
#bayes.mod.fit.gg <- ggs(out1$samples,family="int.d")
#bayes.mod.fit.gg <- ggs(out1$samples,family="beta.covariateS")
#bayes.mod.fit.gg <- ggs(out1$samples,family="beta.covariateT")
#ggs_density(bayes.mod.fit.gg)
#ggs_histogram(bayes.mod.fit.gg)
#ggs_traceplot(bayes.mod.fit.gg)
#ggs_running(bayes.mod.fit.gg)
#ggs_compare_partial(bayes.mod.fit.gg)
#ggs_autocorrelation(bayes.mod.fit.gg)
#ggmcmc(bayes.mod.fit.gg, file = "bayes_fit_ggmcmc.pdf")


#model fits

getBUGSFits<-function(model,param="Density"){
  modelSummary<-data.frame(model$summary[grepl(param,row.names(model$summary)),])
  #extract bits in brackets 
  modelSummary$ParamNu <- sub(".*\\[([^][]+)].*", "\\1", row.names(modelSummary))
  modelSummary$lineIndex <- sapply(modelSummary$ParamNu,function(x)strsplit(x,",")[[1]][1])
  modelSummary$yearIndex <- sapply(modelSummary$ParamNu,function(x)strsplit(x,",")[[1]][2])
  
  if("LinjeID"%in%names(siteInfo)){
  modelSummary<-merge(modelSummary,siteInfo,by.x="lineIndex",by.y="LinjeID")
  modelSummary$Fylkesnavn <- sapply(as.character(modelSummary$Line),function(x)strsplit(x,"_")[[1]][1])
  }else{
    modelSummary<-merge(modelSummary,siteInfo,by.x="lineIndex",by.y="gridIndex")  
  }
  
  names(modelSummary)[c(4,6)]<-c("lower","upper")
  modelSummary$Year <- as.numeric(as.character(modelSummary$yearIndex)) + 2006
  return(modelSummary)
}

