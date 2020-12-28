setwd("C:/Users/db40fysa/Dropbox/ptarmigan Upscaling")
library(plyr)
library(ggplot2)

#read in dataset
tempDF <- read.delim("siteOccupancy_ptarmigan.txt",as.is=TRUE)
#each row is a grid cell (5 x 5 km) and contains 
#data on observations of ptarmigan and
#ecological covariates of each grid

#model 1 - binary
#response = "species" - 0 or 1 to indicate whether the species was ever seen in the grid
#explanatory variables_ bio1 (mean temp), Open (amount of open habitat) + tree line
glm1<-glm(species ~ scale(bio1) + Open + scale(tree_line_position) + 
            scale(I(tree_line_position^2)),family="binomial",data=tempDF)
summary(glm1)

#plot predicted probabilities of occupancy
tempDF$fits<-predict(glm1,type="response",newdata=tempDF)
#plot fits
library(ggplot2)
ggplot(tempDF)+
  geom_point(aes(x,y,colour=fits),shape=15,size=rel(1))+
  scale_colour_gradient(low="steelblue",high="red")

#model 2 - proportion
#response = nuObs/nuVisits (number of observations of ptarmigan out of the total number of survey visits)
#explanatory variables_ bio1 (mean temp), Open (amount of open habitat) + tree line
glm1 <- glm(cbind(nuObs,nuVisits-nuObs)~
              scale(bio1) + scale(Open) + 
              scale(tree_line_position) + scale(I(tree_line_position^2)),
            data=tempDF,family="binomial")

#plot predicted probabilities
tempDF$fits<-predict(glm1,type="response",newdata=tempDF)
ggplot(tempDF)+
  geom_point(aes(x,y,colour=fits),shape=15,size=rel(1))+
  scale_colour_gradient(low="steelblue",high="red")

#also look at SE of predicted probabilities
tempDF$fits.se<-predict(glm1,type="response",newdata=tempDF,se.fit=TRUE)$se.fit
ggplot(tempDF)+
  geom_point(aes(x,y,colour=fits.se),shape=15,size=rel(1))+
  scale_colour_gradient(low="steelblue",high="red")

#also get predictions on the scale of the link function (i.e.. logits)
tempDF$fitsL <- predict(glm1,type="link",newdata=tempDF)
tempDF$fitsL.se <- predict(glm1,type="link",newdata=tempDF,se.fit=TRUE)$se.fit

#group together all predictiond in one data frame
allPreds <- tempDF[,c("grid","x","y","nuVisits","nuObs",
                      "fits","fits.se","fitsL","fitsL.se")]

#we are interested in the population parameter: 
#total number of occupied grid cells
sum(allPreds$fits)

#do a sensitivity analysis:
#for each grid
#take new samples based on predicted mean and Se of occupancy (logit scale)
#for all other grids, just take fixed mean value
#and get new predicted total sample size for each sample

myGridSurveys <- function(allPreds,mygrid){
  
  #for nonfocal grid
  currentFits <- allPreds$fits[allPreds$grid!=mygrid]
  
  #for focal grid, get samples
  totalNu <- as.numeric()
  require(boot)
  for(i in 1:100){
    newSample <- inv.logit(rnorm(1,allPreds$fitsL[allPreds$grid==mygrid],
                                allPreds$fitsL.se[allPreds$grid==mygrid]))
    
    #new predicted total number of grids
    totalNu[i] <- sum(currentFits,newSample)
  }
  
  #get mean and sd of predicted total occupied grids for each sample
  data.frame(grid=mygrid,mean(totalNu),sd(totalNu))
}

#apply function to each grid
library(plyr)
sensitivityOutput <- ldply(allPreds$grid,function(x){
  myGridSurveys(allPreds,mygrid=x)
})

#merge and plot - this shows where variability in the grid prediction
#has most effect of the variability on the population prediction
allPreds <- merge(allPreds,sensitivityOutput,by="grid")
ggplot(allPreds)+
  geom_point(aes(x,y,colour=sd.totalNu.),shape=15,size=rel(1))+
  scale_colour_gradient2(low="grey",mid="pink",high="purple",
                         midpoint=median(allPreds$sd.totalNu.))

#how is the variability in the total number of occupied grid related to
#number of visits and predicted mean proportion
summary(lm(sqrt(sd.totalNu.) ~ scale(fits) + scale(nuVisits),data=allPreds))

#simulate change in binomial variance
#with different values of n

#variation only partly caused by number of visits
#pull out the predicted effect of another visit?
summary(lm(sqrt(sd.totalNu.) ~ scale(fits)*nuVisits,data=allPreds))

allPreds$fitsB <- ifelse(allPreds$fits < median(allPreds$fits),"L","H")
summary(lm(sqrt(sd.totalNu.) ~ fitsB + fitsB:nuVisits,data=allPreds))
#effect of visit more important in higher density areas

#little simulation to explore properties of binomial
n <- 1:100
p1 <- 0.1
p2 <- 0.5
p3 <- 0.9

varBin <- function(n,p){
  sqrt(p*(1-p)/n)
}

varBin(n=n,p=p1)

df1 <- data.frame(n,varBin(n,p1))
df1$P <- p1
names(df1)<-c("n","var","p")

df2 <- data.frame(n,varBin(n,p2))
df2$P <- p2
names(df2)<-c("n","var","p")

df3 <- data.frame(n,varBin(n,p3))
df3$P <- p3
names(df3)<-c("n","var","p")

allDF <- rbind(df1,df2,df3)
qplot(n,var,data=allDF,facets=~factor(p),geom="line")+
  scale_x_log10()+
  ylab("SE")

###influence measures#######################################################################

#DFBETA measures the difference in each parameter estimate with and without the influential point
temp <- dfbeta(glm1)
#this does it for each parameter in the model

#really want to do this on the population parameter but that is derived,
#so just do it on the model intercept (mean prediction)
tempDF$dfbeta <- temp[,1]
ggplot(tempDF)+
  geom_point(aes(x,y,colour=abs(dfbeta)),shape=15,size=rel(1))+
  scale_colour_gradient2(low="grey",mid="pink",high="purple")


###simulate addding in one extra visit#######################################################

#assume probability of success depends on what was already seen
#obs p
tempDF$obs_p <- tempDF$nuObs/tempDF$nuVisits
summary(tempDF$obs_p)#rare!!!

#for each grid cell, imagine one more visit and get observation on new visit
#refit model and get predicted number of sites and successes

myGridSurveys <- function(tempDF,allPreds,mygrid){
  
  require(boot)
  #get new possible observations - for focal grid
  newObs <- as.numeric()
  for(i in 1:100){
    new_P <- inv.logit(rnorm(1,allPreds$fitsL[allPreds$grid==mygrid],
                                 allPreds$fitsL.se[allPreds$grid==mygrid]))
    
    #new observation based on this p
    newObs[i] <- rbinom(1,1,new_P)
  }
  
  newSum <- as.numeric()
  #update the observations and refit the model
  for(i in 1:100){
    tempDF$nuObs[tempDF$grid==mygrid] <- tempDF$nuObs[tempDF$grid==mygrid]+ newObs[i]
    tempDF$nuVisits[tempDF$grid==mygrid] <- tempDF$nuVisits[tempDF$grid==mygrid] +1
    
    #refit the model
    glm1 <- glm(cbind(nuObs,nuVisits-nuObs)~
                  scale(bio1) + scale(Open) + 
                  scale(tree_line_position) + 
                  scale(I(tree_line_position^2)),
                data=tempDF,family="binomial")
    
    #get new predicted probabilities
    tempDF$fits<-predict(glm1,type="response",newdata=tempDF)
    tempDF$fits.se<-predict(glm1,type="response",newdata=tempDF,se.fit=TRUE)$se.fit
    
    #get sum of new fits
    newSum <- sum(tempDF$fits)
  }
  
  data.frame(meanSum = mean(newSum),sdSum = sd(newSum))
}

#apply function to each grid
library(plyr)
sensitivityOutput <- ldply(allPreds$grid,function(x){
  myGridSurveys(tempDF,allPreds,mygrid=x)
})

#plotting
allPreds <- merge(allPreds,sensitivityOutput,by="grid")
ggplot(allPreds)+
  geom_point(aes(x,y,colour=meanSum),shape=15,size=rel(1))+
  scale_colour_gradient2(low="grey",mid="pink",high="purple",
                         midpoint=median(allPreds$sd.totalNu.))

###work at adm2 levels##########################################################################################

#how many grids in each adms
temp <- ddply(tempDF,.(adm2),summarise,nuGrids = length(unique(grid)))
summary(temp$nuGrids)

tempDF_Agg <- ddply(tempDF,.(adm2),summarise,
                    tree_line_position = mean(tree_line_position,na.rm=T),
                    bio1 = mean(bio1,na.rm=T),
                    Open = mean(Open,na.rm=T),
                    Forest = mean(Forest,na.rm=T),
                    PrefOpen = mean(PrefOpen,na.rm=T),
                    nuObs = sum(nuObs,na.rm=T),
                    nuVisits = sum(nuVisits,na.rm=T))

glm1 <- glm(cbind(nuObs,nuVisits-nuObs)~
              scale(bio1) + scale(Open) + scale(Open^2)+
              scale(Forest) + scale(Forest^2) + scale(PrefOpen)+ scale(PrefOpen^2)+
              scale(tree_line_position) + scale(I(tree_line_position^2)),
            data=tempDF_Agg,family="quasibinomial")

#plot predicted probabilities
tempDF_Agg$fits<-predict(glm1,type="response",newdata=tempDF_Agg)
tempDF$fits <- tempDF_Agg$fits[match(tempDF$adm2,tempDF_Agg$adm2)]
ggplot(tempDF)+
  geom_point(aes(x,y,colour=fits),shape=15,size=rel(1))+
  scale_colour_gradient(low="steelblue",high="red")
sum(tempDF$fits)
#now model is terrible

#how to check fit of an SDM
#where is the most discrepancy between model and observations
tempDF_Agg$obs_p <- tempDF_Agg$nuObs/tempDF_Agg$nuVisits
tempDF$obs_p <- tempDF_Agg$obs_p[match(tempDF$adm2,tempDF_Agg$adm2)]
#check how well predictions match data
myMedian <- median(abs(tempDF$fits-tempDF$obs_p))
ggplot(tempDF)+
  geom_point(aes(x,y,colour=abs(fits-obs_p)),shape=15,size=rel(1))+
  scale_colour_gradient2(low="grey",mid="pink",high="purple",
                         midpoint=0.005)

#also look at SE of predicted probabilities
tempDF_Agg$fits.se<-predict(glm1,type="response",newdata=tempDF_Agg,se.fit=TRUE)$se.fit
tempDF_Agg$fitsL <- predict(glm1,type="link",newdata=tempDF_Agg)
tempDF_Agg$fitsL.se <- predict(glm1,type="link",newdata=tempDF_Agg,se.fit=TRUE)$se.fit
#group together all predictiond in one data frame
allPreds_Agg <- tempDF_Agg[,c("adm2","nuVisits","nuObs",
                      "fits","fits.se","fitsL","fitsL.se")]

myGridSurveys <- function(tempDF_Agg,allPreds_Agg,myadm2){
  
  require(boot)
  #get new possible observations - for focal grid
  newObs <- as.numeric()
  for(i in 1:100){
    new_P <- inv.logit(rnorm(1,allPreds_Agg$fitsL[allPreds_Agg$adm2==myadm2],
                             allPreds_Agg$fitsL.se[allPreds_Agg$adm2==myadm2]))
    
    #new observation based on this p
    newObs[i] <- rbinom(1,1,new_P)
  }
  
  newSum <- as.numeric()
  #update the observations and refit the model
  for(i in 1:100){
    tempDF_Agg$nuObs[tempDF_Agg$adm2==myadm2] <- tempDF_Agg$nuObs[tempDF_Agg$adm2==myadm2]+ newObs[i]
    tempDF_Agg$nuVisits[tempDF_Agg$adm2==myadm2] <- tempDF_Agg$nuVisits[tempDF_Agg$adm2==myadm2] +1
    
    #refit the model
    glm1 <- glm(cbind(nuObs,nuVisits-nuObs)~
                  scale(bio1) + scale(Open) + 
                  scale(tree_line_position) + 
                  scale(I(tree_line_position^2)),
                data=tempDF_Agg,family="binomial")
    
    #get new predicted probabilities
    tempDF_Agg$fits<-predict(glm1,type="response",newdata=tempDF_Agg)
    tempDF_Agg$fits.se<-predict(glm1,type="response",newdata=tempDF_Agg,se.fit=TRUE)$se.fit
    
    #get sum of new fits
    newSum <- sum(tempDF_Agg$fits)
  }
  
  data.frame(adm2 = myadm2, meanSum = mean(newSum),sdSum = sd(newSum))
}

#apply function to each grid
library(plyr)
sensitivityOutput <- ldply(allPreds_Agg$adm2,function(x){
  myGridSurveys(tempDF_Agg,allPreds_Agg,
                myadm2=x)
})
tempDF$fits <- sensitivityOutput$meanSum[match(tempDF$adm2,sensitivityOutput$adm2)]
ggplot(tempDF)+
  geom_point(aes(x,y,colour=fits),shape=15,size=rel(1))+
  scale_colour_gradient(low="steelblue",high="red")


###options

#sample where there is most difference berween model predictions and observations
#sample where more samples are expected to chaneg nationwide parameter
