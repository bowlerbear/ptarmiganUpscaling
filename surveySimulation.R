#survey simulation

################################################################################################################

#need to specify:

#effects
lineSD<-2
site2SD<-2
siteSD<-2

#sample sizes
n.lineN<-20
n.site2N<-5
n.siteN<-5#much be equal to or less than n.site2N

################################################################################################################

#function to generate abundances

myfun<-function(lineSD,site2SD,siteSD,n.lineN,n.site2N,n.siteN){
#combine into a design matrix
site2N<-sample(1:n.site2N, n.lineN, replace=T)
siteMatch<-data.frame(site2N=1:n.site2N, siteN=sample(1:n.siteN, n.site2N, replace=T))
siteN<-siteMatch$siteN[match(site2N,siteMatch$site2N)]
df<-data.frame(Line=1:n.lineN,siteN,site2N)

#repeat below 1000 times
abund<-as.numeric()
totAbundReps<-as.numeric()

for(n in 1:1000){
  
  #effects
  lineEffect<-rnorm(n.lineN,0,lineSD)
  site2Effect<-rnorm(n.site2N,0,site2SD)
  siteEffect<-rnorm(n.siteN,0,siteSD)
  
  for(i in 1:nrow(df)){
  abund[i]<-lineEffect[df$Line[i]] + site2Effect[df$site2N[i]] + siteEffect[df$siteN[i]]  
  }  
  
#get total abundance
totAbund<-sum(abund)
totAbundReps[n]<-totAbund
}

sd(totAbundReps)#sampling variability

}

################################################################################################################out

#running the function:

#assuming that all have the same variability:

outSameVar_LineN<-sapply(
myfun(lineSD=2,site2SD=2,siteSD=2,n.lineN=20,n.site2N=5,n.siteN=5)

