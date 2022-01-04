#choose folder

#myfolder <- "data" #local
myfolder <- "/data/idiv_ess/ptarmiganUpscaling" #HPC

### get data ################################################

#get the occupancy data:
bugs.data_ArtsDaten <- readRDS(paste(myfolder,"bugs.data_ArtsDaten.rds",sep="/"))
zst <- readRDS(paste(myfolder,"zst_ArtsDaten.rds",sep="/"))
inits <- function(){list(z = zst)}


names(bugs.data_ArtsDaten) <- sapply(names(bugs.data_ArtsDaten),
                                       function(x)
                                         paste0(x,"CT"))

### get data ###################################################

#get the abundance data:
bugs.data_LineTransects <- readRDS(paste(myfolder,"bugs.data_linetransects.rds",sep="/")) 

names(bugs.data_LineTransects) <- sapply(names(bugs.data_LineTransects),
                                       function(x)
                                         paste0(x,"LT"))

### combine all ##############################################

bugs.data <- c(bugs.data_ArtsDaten,bugs.data_LineTransects)
names(bugs.data)

### fit model ##################################################

library(rjags)
library(jagsUI)

#n.cores = as.integer(Sys.getenv("NSLOTS", "1")) 
n.cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))

params <- c("total.pop","annual.pop","grid.pop")

modelfile <- paste(myfolder,"fullCombinedModel.txt",sep="/")

n.iterations <- 20000

out1 <- jags(bugs.data, 
             inits=inits, 
             params, 
             modelfile, 
             n.thin = 10,
             n.chains = n.cores, 
             n.burnin = n.iterations/2,
             n.iter = n.iterations,
             parallel = T)

saveRDS(out1$summary,file="outSummary_fullCombinedModel.rds")

#######################################################################################
# #https://stats.stackexchange.com/questions/234887/predicting-future-values-with-hurdle-poisson-model
# #https://stats.stackexchange.com/questions/15967/how-can-i-set-up-a-zero-inflated-poisson-in-jags
# #https://stats.stackexchange.com/questions/320924/is-a-hurdle-model-really-one-model-or-just-two-separate-sequential-models
# 
# 
# ### hurdle model formula
# 
# #The basic idea is 'success probability of binomial distribution' * 'lambda (= an expected value) of #poisson distribution'. But you have to consider that the poisson model in count part never returns 0.
# #phi_zero <- 1 / ( 1 + exp(-(model$coef$zero[[1]] + model$coef$zero[[2]] * c(23, 26, 29))))
# 
# 
# #Second, you calculate a param (= expected value) of a poisson model, mu, and the > 0 probability, #phi_count
# #mu <- exp(model$coef$count[[1]] + model$coef$count[[2]] * c(23, 26, 29))
# #phi_count <- ppois(0, lambda = mu, lower.tail = F)  # not 0 probability
# 
# #integrate both values
# phi <- phi_zero / phi_count  # because there isn't 0 coming from poisson.
# rval <- phi * mu
# 
# #zeros only comes from first dataset
# 
# #file:///C:/Users/db40fysa/AppData/Local/Temp/v27i08.pdf
# #p6
# #Hurdle models combine a left-truncated count component with a right-censored hurdle com-ponent
# #They are two-component models
# #A truncated count component, such as Poisson, geometric or negative binomial, is employed for positive #counts, and a hurdle component models zero vs. larger counts.  For thelatter, either a binomial model or #a censored count distribution can be employed.
# 
# #play with data
# set.seed(1839)
# # simulate poisson with many zeros
# x <- rnorm(100)
# e <- rnorm(100)
# y <- rpois(100, exp(-1.5 + x + e))
# 
# # how many zeroes?
# table(y == 0)
# 
# #hurdle
# library(pscl)
# hurdle1 <- pscl::hurdle(y ~ x)
# summary(hurdle1)
# predict(hurdle1)
# 
# #separate models
# library(VGAM)
# summary(VGAM::vglm(y[y > 0] ~ x[y > 0], family = pospoisson()))
# summary(glm(I(y == 0) ~ x, family = binomial))
# 
# data("bioChemists", package = "pscl")
# 
# ## logit-poisson
# ## "art ~ ." is the same as "art ~ . | .", i.e.
# ## "art ~ fem + mar + kid5 + phd + ment | fem + mar + kid5 + phd + ment"
# 
# fm_hp1 <- hurdle(art ~ fem + mar + kid5, data = bioChemists)
# summary(fm_hp1)
# bioChemists$hurdle_preds <- predict(fm_hp1)
# 
# #binomial preds
# fm_hp1 <- glm(I(art==0) ~ fem + mar + kid5, data = bioChemists,family=binomial)
# bioChemists$binomial_preds <- predict(fm_hp1,type="response")
# 
# #pospoisson model
# fm_hp1 <- vglm(art ~ fem + mar + kid5, data = subset(bioChemists,art>0), family = pospoisson())
# bioChemists$pospoison_preds <- NA
# bioChemists$pospoison_preds[bioChemists$art>0] <- predict(fm_hp1,type="response")
# 
# 
# 
# #other dataset
# library(pscl)
# data <- data.frame(y = c(8, 0, 3, 7), width = c(34.40, 22.50, 28.34, 32.22))
# data$success <- ifelse(data$y>0,1,0)
# model <- hurdle(y ~ width, data = data)
# model
# 
# #binomial model - succes probability
# phi_zero <- 1 / ( 1 + exp(-(model$coef$zero[[1]] + model$coef$zero[[2]] * c(23, 26, 29))))
# phi_zero_coef <- model$coef$zero[[1]] + model$coef$zero[[2]] * c(23, 26, 29)
# # p0_zero <- log(phi_zero)
# 
# #count model
# mu <- exp(model$coef$count[[1]] + model$coef$count[[2]] * c(23, 26, 29))
# phi_count <- ppois(0, lambda = mu, lower.tail = F)  # not 0 probability
# #1-exp(-lambda)
# 
# # p0_count <- ppois(0, lambda = mu, lower.tail = F, log.p = T)
# 
# #integration
# phi <- phi_zero / phi_count  # because there isn't 0 coming from poisson.
# rval <- phi * mu
# # logphi <- p0_zero - p0_count
# # rval2 <- exp(logphi + log(mu))
# rval
# # [1] 8.051324e-09 2.429582e+00 3.674317e+00
# 
# #compare
# predict(model, data.frame(width = c(23, 26, 29)), type = "response")
# 
# #breaking it down
# #binomial model
# binModel <- glm(success ~ width, data=data,family="binomial")
# summary(binModel)
# successProb_link <- predict(binModel,data.frame(width = c(23, 26, 29)),type="link")
# successProb <- predict(binModel,data.frame(width = c(23, 26, 29)),type="response")
# #same as phi_zero
# 
# #poisson model
# poisModel <- vglm(y ~ width, data = subset(data,y>0), family = pospoisson())
# summary(poisModel)#matches with "model" above
# countProb <- predict(poisModel,data.frame(width = c(23, 26, 29)),type="response")
# #Coefficients: 
# #            Estimate Std. Error z value Pr(>|z|)
# #(Intercept)  -3.4520     3.6332      NA       NA
# #width         0.1629     0.1110   1.467    0.142
# 
# #successProb * countProb - already accounts for the zeros....
# #8.051323e-09
# #2.429581e+00
# #3.674316e+00
# 
# #zero truncated poisson in JAGS
# #y[i] ~ dpois(mu[i]) T(1, )
