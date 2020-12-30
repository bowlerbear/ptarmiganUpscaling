tempDF$fitsL <- predict(glm1,type="link",newdata=tempDF)
tempDF$fitsL.se <- predict(glm1,type="link",newdata=tempDF,se.fit=TRUE)$se.fit

#sensitivity analysis
allPreds <- tempDF[,c("grid","x","y","nuVisits",
                      "fits","fits.se","fitsL","fitsL.se")]

#total predicted number of grid cells
sum(allPreds$fits)
sum(allPreds$fits.se)

#for each grid
#all current occupancies for all other grids
#take new samples for focus grid
#get new predicted total sample size

myGridSurveys <- function(allPreds,mygrid){
  
  currentFits <- allPreds$fits[allPreds$grid!=mygrid]
  
  totalNu <- as.numeric()
  require(boot)
  for(i in 1:100){
    newSample <- invlogit(rnorm(1,allPreds$fitsL[allPreds$grid==mygrid],
                                allPreds$fitsL.se[allPreds$grid==mygrid]))
    
    #new predicted total number of grids
    totalNu[i] <- sum(currentFits,newSample)
  }
  
  data.frame(grid=mygrid,mean(totalNu),sd(totalNu))
}

#apply function to each grid
sensitivityOutput <- ldply(allPreds$grid,function(x){
  myGridSurveys(allPreds,mygrid=x)
})

allPreds <- merge(allPreds,sensitivityOutput,by="grid")
ggplot(allPreds)+
  geom_point(aes(x,y,colour=sd.totalNu.),shape=15,size=rel(1))+
  scale_colour_gradient2(low="grey",mid="pink",high="purple",
                         midpoint=0.004)

qplot(fits,sd.totalNu.,data=allPreds)
qplot(sd.totalNu.,nuVisits,data=allPreds)

summary(lm(sqrt(sd.totalNu.) ~ scale(fits) + scale(nuVisits),data=allPreds))

#simulate change in binomial variance
#with different values of n

#variation only partly caused by number of visits
#pull out the predicted effect of another visit?
summary(lm(sqrt(sd.totalNu.) ~ scale(fits)*nuVisits,data=allPreds))

allPreds$fitsB <- ifelse(allPreds$fits < median(allPreds$fits),"L","H")
summary(lm(sqrt(sd.totalNu.) ~ fitsB + fitsB:nuVisits,data=allPreds))
#effect of visit more important in higher density areas


n <- 1:100

p1 <- 0.1
p2 <- 0.5
p3 <- 0.9

varBin <- function(n,p){
  n*p*(1-p)
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
qplot(n,var,data=allDF,facets=~factor(p))

###influence measures#######################################################################

#https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/influence.measures

#Bar Plot of Cook's distance to detect observations that strongly influence fitted values of the mode
#cooks.distance

tempDF$cooks <- cooks.distance(glm1)
ggplot(tempDF)+
  geom_point(aes(x,y,colour=cooks),shape=15,size=rel(1))+
  scale_colour_gradient(low="steelblue",high="red")

#DFBETA measures the difference in each parameter estimate with and without the influential point
#dfbeta - for the overall mean
temp <- dfbeta(glm1)

tempDF$dfbeta <- temp[,1]
ggplot(tempDF)+
  geom_point(aes(x,y,colour=abs(dfbeta)),shape=15,size=rel(1))+
  scale_colour_gradient2(low="grey",mid="pink",high="purple")