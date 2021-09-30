require(foreign)
require(ggplot2)
require(MASS)

dat <- read.dta("https://stats.idre.ucla.edu/stat/stata/dae/nb_data.dta")
dat <- within(dat, {
  prog <- factor(prog, levels = 1:3, labels = c("General", "Academic", "Vocational"))
  id <- factor(id)
})

summary(m1 <- glm.nb(daysabs ~ math, data = dat))
dat$preds <- predict(m1,type="response")
qplot(daysabs,preds,data=dat)

#format data for JAGS
bugs.data <- list(
  N = nrow(dat),
  y = dat$daysabs,
  cov1 = dat$math
)


sink("neg_bin_simpleModel.txt")

cat("
model{

  #intercept
  int.d ~ dnorm(0,0.001)    
  
  #covariate effects
  beta1 ~ dnorm(0,0.01)
  beta2 ~ dnorm(0,0.01)
  
  r ~ dunif(0,50)
  
  for(i in 1:N){
      
      y[i] ~ dnegbin(p[i],r)
      
      p[i] <- r/(r+expNuIndivs[i])

      log(expNuIndivs[i]) <- int.d + beta1*cov1[i]

    }
  
  #calculate the Bayesian p-value
  e <- 0.0001
  
  for(i in 1:N){
      
      #Discrepancy
      exp1[i] <- r/p[i]
      exp2[i] <- (p[i] * r)/(1-p[i])
      exp3[i] <- (r*(1-p[i]))/p[i]
      
      #compared to data
      chi2[i] <- pow((y[i] - exp1[i]),2) / (sqrt(exp1[i])+e)
      
      #New samples
      #compared to simulated values
      y.new[i] ~ dnegbin(p[i],r)
      
      chi2.new[i] <- pow((y.new[i] - exp1[i]),2) / (sqrt(exp1[i])+e) 
      
  }
  
  
    fit <- sum(chi2[])                     
    fit.new <- sum(chi2.new[])             
  
}
",fill = TRUE)
sink()


library(rjags)
library(jagsUI)

params <- c("int.d","beta1","fit","fit.new","expNuIndivs","y.new",
            "exp1","exp2","exp3")


#run model
out1 <- jags(bugs.data, 
             inits=NULL, 
             params, 
             "neg_bin_simpleModel.txt", 
             n.thin=3,
             n.chains=3, 
             n.burnin=1000,
             n.iter=2000,
             parallel=T)

mean(out1$sims.list$fit.new > out1$sims.list$fit)
hist(out1$sims.list$fit.new)
abline(v=mean(out1$sims.list$fit))
#good!!!

plot(out1$mean$exp1,out1$mean$exp2)
plot(out1$mean$exp2,out1$mean$exp3)
plot(out1$mean$exp1,out1$mean$exp3)#one and three are the same!!!
plot(dat$preds,out1$mean$exp3)#three matches the predictions of the glm
cor(dat$preds,out1$mean$exp3)

#relationship with link model response
plot(out1$mean$expNuIndivs,out1$mean$exp3)#the same!!!
