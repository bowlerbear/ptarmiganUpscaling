##### Code for binomial simulations #####
##  Simulate my data
# Mydatalenght is set to 100, and  the probabilities for n and x are set to 0.6 and 0.3 respectively, but these can be changed
# specify the length of my vector (i.e. how many pixels)
mydatalength <- 1000
# Approximate average number of visits per square
sim_n <- 10
# Simulate the number of trials
n <- rbinom(mydatalength,sim_n,0.6)
# Simulate the number of successes - for p I use 0.3 because that is somehow close to the Ptarmigan data
x <- rbinom(mydatalength,n,0.3)
# This is to create some covariate data in case we need it (this can also be randomized)
cov <- seq(0.1, 2,length.out=mydatalength)
# Combine simulated data together
mydata1 <- as.data.frame(cbind(n,x,cov))
mydata1

# Base regression - regression model on my base data (the data from which we are going to make assumptions)
glmbase <-glm(cbind(x,n-x) ~ 1, data=mydata1,family="binomial")
summary(glmbase)
#Obtain transformed probabilities
pars_base <-predict(glmbase,type="response",newdata=mydata1)
#Obtain transformed standard errors
se_base <- predict(glmbase,type="response",newdata=mydata1,se.fit=TRUE)$se.fit


#-------------# 
##### Part 1: What happens when I do an extra visit to square 1?
#----#
# Here I add 1 more trial to the first row (first n), but I keep the number of successes the same. (i.e. I visit my first pixel one more time and I see nothing)
mydata1.0 <- as.data.frame(cbind(c(mydata1[1,1]+1,mydata1[2,1],mydata1[3:mydatalength,1]), mydata1[,2],mydata1[,3]))
names(mydata1.0) <- names(mydata1)

# Regression model on my data1.0
glmbase1.0 <- glm(cbind(x,n-x) ~ 1, data=mydata1.0,family="binomial")
#Obtain transformed probabilities
pars_1.0 <- predict(glmbase1.0,type="response",newdata=mydata1.0)
#Obtain transformed standard errors
se_1.0 <-predict(glmbase1.0,type="response",newdata=mydata1.0,se.fit=TRUE)$se.fit
#----#


# Here I add 1 more trial to the first row (first n), and 1 success to the first row (i.e. I visit my first pixel one more time and I see a bird)
mydata1.1 <- as.data.frame(cbind(c(mydata1[1,1]+1,mydata1[2:mydatalength,1]), c(mydata1[1,2]+1,mydata1[2:mydatalength,2]),mydata1[,3]))
names(mydata1.1) <- names(mydata1)
# Regression model on my data1.1
glmbase1.1 <- glm(cbind(x,n-x) ~ 1, data=mydata1.1,family="binomial")
#Obtain transformed probabilities
pars_1.1 <- predict(glmbase1.1,type="response",newdata=mydata1.1)
#Obtain transformed standard errors
se_1.1 <- predict(glmbase1.1,type="response",newdata=mydata1.1,se.fit=TRUE)$se.fit
#----#


## Compare estimates obtained from the 3 data sets using predict to calculating them with the formula
# This lines were just to check that hat_p = sum(x_i)/sum(n_i)
# This is what I obtain from predict() if I do not have covariates
pars_base[1]
#This is the formula for p_hat instead of using a built in function
sum(mydata1[,2])/sum(mydata1[,1])
# This is what I obtain from predict() if I do not have covariates
pars_1.0[1]
#This is the formula for p_hat instead of using a built in function
sum(mydata1.0[,2])/sum(mydata1.0[,1])
# This is what I obtain from predict() if I do not have covariates
pars_1.1[1]
#This is the formula for p_hat instead of using a built in function
sum(mydata1.1[,2])/sum(mydata1.1[,1])
# This lines were just to check se(hat_p) 
se_base[1]
sqrt(pars_base*(1-pars_base)/sum(mydata1[,1]))[1]
se_1.0[1]
sqrt(pars_1.0*(1-pars_1.0)/sum(mydata1.0[,1]))[1]
se_1.1[1]
sqrt(pars_1.1*(1-pars_1.1)/sum(mydata1.1[,1]))[1]


# Comparison of parameter estimates and standard errors between the 3 data sets
# Note here the difference between the parameters will depend a lot on:
# the number of pixels visited (or how large we set mydatalength), 
# on the average number of visits per square (sim_n) and 
#also on the probability of success of x 
# This is because of how p is estimated hat_p = sum(x_i)/sum(n_i) 
myestimates <- cbind(pars_base,pars_1.0,pars_1.1)
mystanderrors <-cbind(se_base,se_1.0,se_1.1)
round(myestimates[1,],4)
round(mystanderrors[1,],4)

#-------------# 
##### Part 2: What happens when I do one more visit to square 1, and I take one visit out from square 2?
# There are 4 outcomes:
# 1. I visit square 1, 1 more time and I see nothing:
#a) The visit I removed from square 2 was a visit with no bird spotted 
mydata1.0.1 <- as.data.frame(cbind(c(mydata1[1,1]+1,mydata1[2,1]-1,mydata1[3:mydatalength,1]), mydata1[,2],mydata1[,3]))
# b) The visit I removed from square 2 was a visit with a bird spotted 
lenghtrows <- dim(mydata1)[1]
fail <-matrix(data=c(rep(0,lenghtrows+1),1, rep(0,(lenghtrows*2)-2)), nrow=dim(mydata1)[1], ncol=dim(mydata1)[2])
mydata1.0.0 <- as.data.frame(mydata1.0.1-fail)
# 2. I visit square 1, 1 more time and I see a bird:
#a) The visit I removed from square 2 was a visit with no bird spotted 
fail1 <-matrix(data=c(0,1,rep(0,(lenghtrows*3-2))), nrow=dim(mydata1)[1], ncol=dim(mydata1)[2])
mydata1.1.1 <- as.data.frame(mydata1.1-fail1)
# b) The visit I removed from square 2 was a visit with a bird spotted 
mydata1.1.0 <- as.data.frame(mydata1.1.1-fail)
###----
# Regresion model on data 1.0.1
glmbase1.0.1 <- glm(cbind(x,n-x) ~ 1, data=mydata1.0.1,family="binomial")
#Obtain transformed probabilities
pars_1.0.1 <- predict(glmbase1.0.1,type="response",newdata=mydata1.0.1)
#Obtain transformed standard errors
se_1.0.1 <-predict(glmbase1.0.1,type="response",newdata=mydata1.0.1,se.fit=TRUE)$se.fit
#----# 
# Regresion model on data 1.0.0
glmbase1.0.0 <- glm(cbind(x,n-x) ~ 1, data=mydata1.0.0,family="binomial")
#Obtain transformed probabilities
pars_1.0.0 <- predict(glmbase1.0.0,type="response",newdata=mydata1.0.0)
#Obtain transformed standard errors
se_1.0.0 <-predict(glmbase1.0.0,type="response",newdata=mydata1.0.0,se.fit=TRUE)$se.fit
# Regresion model on data 1.1.1
glmbase1.1.1 <- glm(cbind(x,n-x) ~ 1, data=mydata1.1.1,family="binomial")
#Obtain transformed probabilities
pars_1.1.1 <- predict(glmbase1.1.1,type="response",newdata=mydata1.1.1)
#Obtain transformed standard errors
se_1.1.1 <-predict(glmbase1.1.1,type="response",newdata=mydata1.1.1,se.fit=TRUE)$se.fit
#----# 
# Regresion model on data 1.1.0
glmbase1.1.0 <- glm(cbind(x,n-x) ~ 1, data=mydata1.1.0,family="binomial")
#Obtain transformed probabilities
pars_1.1.0 <- predict(glmbase1.1.0,type="response",newdata=mydata1.1.0)
#Obtain transformed standard errors
se_1.1.0 <-predict(glmbase1.1.0,type="response",newdata=mydata1.1.0,se.fit=TRUE)$se.fit

# Note here the number of trials never change, only the number of successes change
# having 1 less success in the case of 1.0.0 and 1 more in the case of 1.1.1
myestimates2 <- cbind(pars_base,pars_1.0.1,pars_1.0.0,pars_1.1.1, pars_1.1.0)
mystanderrors2 <-cbind(se_base,se_1.0.1,se_1.0.0,se_1.1.1, se_1.1.0)
round(myestimates2[1,],4)
round(mystanderrors2[1,],4)



