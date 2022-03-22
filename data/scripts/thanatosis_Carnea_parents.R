###############################
### Carnea thanatosis analysis
### Katie Taylor and Tim Farkas
### Fall 2016
###############################

#### Preliminaries ####
# load libraries
library(R2jags)

# load and manipulate data
td <- read.csv("~/Dropbox/Projects/Carnea_Thanatosis/data/thanatosis_Carnea_parents.csv")
td$c.sex <- ifelse(is.na(td$sex), NA, ifelse(td$sex == 0, -.5, .5))
td$m.sex <- ifelse(is.na(td$sex), NA, ifelse(td$sex == 0, 1, 0))
td$f.sex <- td$sex
td$starve.tmnt <- ifelse(td$tmnt=="fed", .5, -.25)
td$junk.tmnt <- ifelse(td$tmnt=="one", -.25, ifelse(td$tmnt=="fed", 0, .25))
td$tmnt.one <- ifelse(td$tmnt=="fed", 0, ifelse(td$tmnt=="two", 0, 1))
td$tmnt.two <- ifelse(td$tmnt=="fed", 0, ifelse(td$tmnt=="one", 0, 1))
td_feign <- td[!is.na(td$feign) & !is.na(td$sex),]
td_prod <- td[!is.na(td$prod) & !is.na(td$sex),]

contrasts(td$tmnt) <- cbind(one=c(0,1,0), two=c(0,0,1))
contrasts(td$tmnt) <- cbind(starve=c(.5, -.25, -.25), junk=c(0, -.25, .25))

# code GLMM

an1 <- glm(feign ~ tmnt * c.sex, family=quasibinomial, data=td)
summary(an1)
Anova(an1, type=3)

an2 <- glm(prod ~ tmnt * sex, family=quasipoisson, data=td_prod)
summary(an2)
Anova(an2, type=3)

#### JAGS Model: Feigning Yes/NO ####

sink(file="~/Dropbox/Projects/Carnea_Thanatosis/data/feign_JAGS_model.txt")
cat(
  "model{
    
    # priors on beta
      
      b0 ~ dnorm(0, 0.0001)
      #b.tmnt_starve ~ dnorm(0, 0.0001)
      #b.tmnt_junk ~ dnorm(0, 0.0001)
      b.tmnt_one ~ dnorm(0, 0.0001)
      b.tmnt_two ~ dnorm(0, 0.0001)
      #b.sex ~ dnorm(0, 0.0001)

    # likelihood

      for (i in 1:enn) {
        
        y[i] ~ dbern(p[i])
        logit(p[i]) <- eta[i]
        #eta[i] <- b0 + b.tmnt_starve * tmnt_starve[i] + b.tmnt_junk * tmnt_junk[i] + b.sex * c.sex[i]
        #eta[i] <- b0 + b.tmnt_starve * tmnt_starve[i] + b.tmnt_junk * tmnt_junk[i] 
        eta[i] <- b0 + b.tmnt_one * tmnt_one[i] + b.tmnt_two * tmnt_two[i]

      }

    # discrepancy measures
      
      for (i in 1:enn) {  
        exp[i] <- p[i]
        Var[i] <- p[i] * (1 - p[i])
        e[i] <- (y[i] - exp[i]) / sqrt(Var[i]) # pearson residuals
        d[i] <- pow(e[i], 2)
      
      }

      # simulate new data
      
      for (i in 1:enn) {

        ynew[i] ~ dbern(p[i])
        enew[i] <- (ynew[i] - exp[i]) / sqrt(Var[i])
        dnew[i] <- pow(enew[i], 2)

      }

    fit <- sum(d[1:enn])
    fitnew <- sum(dnew[1:enn])

  }",
fill=TRUE)
sink()

#### Run JAGS Model: Feigning Yes/NO ####

#parms <- c("b0", "b.tmnt_starve", "b.tmnt_junk", "b.sex", "fit", "fitnew", "e", "exp")
parms <- c("b0", "b.tmnt_one", "b.tmnt_two", "fit", "fitnew", "e", "exp")

inits <- function() {
  list(b0=rnorm(1, 0, 0.01),
       b.tmnt_one=rnorm(1, 0, 0.01),
       b.tmnt_two=rnorm(1, 0, 0.01))
}

# data <- list(y = td_feign$feign, 
#              tmnt_starve=td_feign$starve.tmnt,
#              tmnt_junk=td_feign$junk.tmnt,
#              c.sex = td_feign$c.sex, 
#              enn = nrow(td_feign))

data <- list(y = td_feign$feign, 
             tmnt_one=td_feign$tmnt.one,
             tmnt_two=td_feign$tmnt.two,
             enn = nrow(td_feign))

mod1 <- jags(model.file = "~/Dropbox/Projects/Carnea_Thanatosis/data/feign_JAGS_model.txt",
             parameters.to.save = parms,
             data = data,
             n.iter = 5000, 
             n.burnin = 1000, 
             n.chains = 3, 
             n.thin = 2)

#### Validation: Feigning Yes/No ####

out1 <- mod1$BUGSoutput
parms2print <- parms[1:3]
print(out1$summary[parms2print,], digits=3)

# MCMC validation

traceplot(mod1, varname=parms2print)

# overdispersion

ep <- mod1$BUGSoutput$mean$e
enn <- nrow(td_feign)
pm <- 4
disp <- sum(ep^2)/(enn - 4) # 1.09

# model fit?
mean(mod1$BUGSoutput$sims.list$fit > mod1$BUGSoutput$sims.list$fitnew) # 0.668

#### Plot Results: Feigning Yes/No ####

# plot posterior

b0_post <- mod1$BUGSoutput$sims.list$b0
b.tmnt_starve_post <- mod1$BUGSoutput$sims.list$b.tmnt_starve
b.sex_post <- mod1$BUGSoutput$sims.list$b.sex

b.tmnt_one_post <- out1$sims.list$b.tmnt_one
b.tmnt_two_post <- out1$sims.list$b.tmnt_two

hist(b0_post)
hist(b.tmnt_starve_post)
hist(b.sex_post)

# sketch fitted values

mydata <- data.frame(one=c(0, 1, 0), two=c(0, 0, 1))
mat <- model.matrix( ~ one + two, data=mydata)
betas <- c(out1$mean$b0, out1$mean$b.tmnt_one, out1$mean$b.tmnt_two)
mu <- inv.logit(mat %*% betas)

# another way with 95% CIs

mcmc.betas <- cbind(b0_post, b.tmnt_one_post, b.tmnt_two_post)
mcmc.mus <- inv.logit(mat %*% t(mcmc.betas))

plotsum <- matrix(nrow=nrow(mcmc.mus), ncol=4)

for (i in 1:nrow(mcmc.mus)) {
  
  xi <- mcmc.mus[i,]
  plotsum[i, 3:4] <- quantile(xi, probs=c(0.025, 0.975))
  plotsum[i, 1] <- mean(xi)
  plotsum[i, 2] <- sd(xi)
  
}
  
colnames(plotsum) <- c("mean", "sd", "2.5%", "97.5%")
rownames(plotsum) <- c("0", "1", "2")

bp1 <- barplot(height=plotsum[,"mean"], ylim=c(0, 1), 
        xlab = "days starved", ylab = "proportion feigning death", las=1, 
        cex.axis=1.3, cex=1.3)
arrows(x0=bp1, 
       y0=plotsum[,"2.5%"], 
       x1=bp1,
       y1=plotsum[,"97.5%"], 
       angle=90, code=3, lwd=1.3)
text(x=bp1, y=0.05, labels=c("A", "B", "B"), cex=2)

#### JAGS Model: Prods Poisson ####

sink(file="~/Dropbox/Projects/Carnea_Thanatosis/data/prod_poisson_JAGS_model.txt")
cat(
  "model{
  
  # priors on beta
  
  b0 ~ dnorm(0, 0.0001)
  b.tmnt_starve ~ dnorm(0, 0.0001)
  b.tmnt_junk ~ dnorm(0, 0.0001)
  #b.tmnt_one ~ dnorm(0, 0.0001)
  #b.tmnt_two ~ dnorm(0, 0.0001)
  b.sex ~ dnorm(0, 0.0001)
  b.int_starve ~ dnorm(0, 0.0001)
  b.int_junk ~ dnorm(0, 0.0001)
  # b.int_one ~ dnorm(0, 0.0001)
  # b.int_two ~ dnorm(0, 0.0001)

  # likelihood
  
  for (i in 1:enn) {
  
  y[i] ~ dpois(mu[i])

  #eta[i] <- b0 + b.tmnt_starve * tmnt_starve[i] + b.tmnt_junk * tmnt_junk[i] + b.sex * c.sex[i]
  #eta[i] <- b0 + b.tmnt_starve * tmnt_starve[i] + b.tmnt_junk * tmnt_junk[i] 
  #eta[i] <- b0 + b.tmnt_one * tmnt_one[i] + b.tmnt_two * tmnt_two[i]

  log(mu[i]) <- b0 + b.tmnt_starve * tmnt_starve[i] + b.tmnt_junk * tmnt_junk[i] + b.sex * c.sex[i] + b.int_starve * tmnt_starve[i] * c.sex[i] + b.int_junk * tmnt_junk[i] * c.sex[i]

  }
  
  # discrepancy measures
  
  for (i in 1:enn) {  
  exp[i] <- mu[i]
  Var[i] <- mu[i] 
  e[i] <- (y[i] - exp[i]) / sqrt(Var[i]) # pearson residuals
  d[i] <- pow(e[i], 2)
  
  }
  
  # simulate new data
  
  for (i in 1:enn) {
  
  ynew[i] ~ dpois(mu[i])
  enew[i] <- (ynew[i] - exp[i]) / sqrt(Var[i])
  dnew[i] <- pow(enew[i], 2)
  
  }
  
  fit <- sum(d[1:enn])
  fitnew <- sum(dnew[1:enn])
  
  }",
fill=TRUE)
sink()

#### JAGS Model: Prods Negative Binomial ####

sink(file="~/Dropbox/Projects/Carnea_Thanatosis/data/prod_negbin_JAGS_model.txt")
cat(
  "model{
  
  # priors on beta
  
    b0 ~ dnorm(0, 0.0001)
    b.tmnt_starve ~ dnorm(0, 0.0001)
    b.tmnt_junk ~ dnorm(0, 0.0001)
    #b.tmnt_one ~ dnorm(0, 0.0001)
    #b.tmnt_two ~ dnorm(0, 0.0001)
    b.sex ~ dnorm(0, 0.0001)
    b.int_starve ~ dnorm(0, 0.0001)
    b.int_junk ~ dnorm(0, 0.0001)
    # b.int_one ~ dnorm(0, 0.0001)
    # b.int_two ~ dnorm(0, 0.0001)

  # prior on size parameter for neg bin
  
    size ~ dunif(0.0001, 5)

  # likelihood
  
    for (i in 1:enn) {
  
      y[i] ~ dnegbin(p[i], size)
      p[i]  <- size / (size + mu[i])
      log(mu[i]) <- eta[i]
  
      #eta[i] <- b0 + b.tmnt_starve * tmnt_starve[i] + b.tmnt_junk * tmnt_junk[i] + b.sex * c.sex[i]
      #eta[i] <- b0 + b.tmnt_starve * tmnt_starve[i] + b.tmnt_junk * tmnt_junk[i] 
      #eta[i] <- b0 + b.tmnt_one * tmnt_one[i] + b.tmnt_two * tmnt_two[i]
  
      eta[i] <- b0 +
                b.tmnt_starve * tmnt_starve[i] +
                b.tmnt_junk * tmnt_junk[i] +
                b.sex * c.sex[i] +
                b.int_starve * tmnt_starve[i] * c.sex[i] +
                b.int_junk * tmnt_junk[i] * c.sex[i]

  # discrepancy measures

      exp[i] <- mu[i]
      Var[i] <- mu[i] + mu[i] * mu[i] / size
      e[i] <- (y[i] - exp[i]) / sqrt(Var[i]) # pearson residuals
      d[i] <- pow(e[i], 2)
  
  # simulate new data
  
      ynew[i] ~ dnegbin(p[i], size)
      enew[i] <- (ynew[i] - exp[i]) / sqrt(Var[i])
      dnew[i] <- pow(enew[i], 2)
  
    }
  
  fit <- sum(d[1:enn])
  fitnew <- sum(dnew[1:enn])
  
  }",
fill=TRUE)
sink()
#### Run JAGS Model: Prods -- NegBin ####

#parms <- c("b0", "b.tmnt_starve", "b.tmnt_junk", "b.sex", "fit", "fitnew", "e", "exp")
parms <- c("b0", 
           "b.tmnt_starve", 
           "b.tmnt_junk", 
           "b.sex",
           "b.int_starve", 
           "b.int_junk",
           "fit", "fitnew", "e", "exp")

# inits <- function() {
#   list(b0=rnorm(1, 0, 0.01),
#        b.tmnt_one=rnorm(1, 0, 0.01),
#        b.tmnt_two=rnorm(1, 0, 0.01))
# }

# data <- list(y = td_feign$feign, 
#              tmnt_starve=td_feign$starve.tmnt,
#              tmnt_junk=td_feign$junk.tmnt,
#              c.sex = td_feign$c.sex, 
#              enn = nrow(td_feign))

data <- list(y = td_prod$prod, 
             tmnt_starve=td_prod$starve.tmnt,
             tmnt_junk=td_prod$junk.tmnt,
             c.sex=td_prod$f.sex,
             enn = nrow(td_prod))

mod3 <- jags(model.file = "~/Dropbox/Projects/Carnea_Thanatosis/data/prod_negbin_JAGS_model.txt",
             parameters.to.save = parms,
             data = data,
             n.iter = 4000, 
             n.burnin = 1000, 
             n.chains = 3, 
             n.thin = 2)
mod3 <- update(mod2, n.iter = 20000, n.thin = 10)

#### Validation: Prods ####
mod3 <- mod2
out3 <- mod3$BUGSoutput
parms2print <- parms[1:6]
print(out3$summary[parms2print,], digits=3)

# MCMC validation

traceplot(mod3, varname=parms2print)

# overdispersion

ep <- mod3$BUGSoutput$mean$e
enn <- nrow(td_prod)
pm <- 6
disp <- sum(ep^2)/(enn - pm) # 1.04! problem solved

# model fit?
mean(mod3$BUGSoutput$sims.list$fit > mod3$BUGSoutput$sims.list$fitnew) # 0.617

#### Plot Results: Prods ####

# plot posterior

b0_post <- mod1$BUGSoutput$sims.list$b0
b.tmnt_starve_post <- mod1$BUGSoutput$sims.list$b.tmnt_starve
b.sex_post <- mod1$BUGSoutput$sims.list$b.sex

b.tmnt_one_post <- out1$sims.list$b.tmnt_one
b.tmnt_two_post <- out1$sims.list$b.tmnt_two

hist(b0_post)
hist(b.tmnt_starve_post)
hist(b.sex_post)

# sketch fitted values

mydata <- data.frame(one=c(0, 1, 0), two=c(0, 0, 1))
mat <- model.matrix( ~ one + two, data=mydata)
betas <- c(out1$mean$b0, out1$mean$b.tmnt_one, out1$mean$b.tmnt_two)
mu <- inv.logit(mat %*% betas)

# another way with 95% CIs

mcmc.betas <- cbind(b0_post, b.tmnt_one_post, b.tmnt_two_post)
mcmc.mus <- inv.logit(mat %*% t(mcmc.betas))

plotsum <- matrix(nrow=nrow(mcmc.mus), ncol=4)

for (i in 1:nrow(mcmc.mus)) {
  
  xi <- mcmc.mus[i,]
  plotsum[i, 3:4] <- quantile(xi, probs=c(0.025, 0.975))
  plotsum[i, 1] <- mean(xi)
  plotsum[i, 2] <- sd(xi)
  
}

colnames(plotsum) <- c("mean", "sd", "2.5%", "97.5%")
rownames(plotsum) <- c("0", "1", "2")

bp1 <- barplot(height=plotsum[,"mean"], ylim=c(0, 1), 
               xlab = "days starved", ylab = "proportion feigning death", las=1, 
               cex.axis=1.3, cex=1.3)
arrows(x0=bp1, 
       y0=plotsum[,"2.5%"], 
       x1=bp1,
       y1=plotsum[,"97.5%"], 
       angle=90, code=3, lwd=1.3)
text(x=bp1, y=0.05, labels=c("A", "B", "B"), cex=2)

#### JAGS Model: Duration -- Beta ####

td_dur <- td[!is.na(td$time) & !is.na(td$sex),]
hist(td_dur$time/max(td_dur$time))
td_dur$b_time <- td_dur$time/max(td_dur$time)
td_dur$b_time <- ifelse(td_dur$b_time==1, .99999, td_dur$b_time)
td_dur$jsex.m <- td_dur$junk.tmnt * td_dur$m.sex
td_dur$jsex.f <- td_dur$junk.tmnt * td_dur$f.sex
td_dur$ssex.m <- td_dur$starve.tmnt * td_dur$m.sex
td_dur$ssex.f <- td_dur$starve.tmnt * td_dur$f.sex

sink(file="~/Dropbox/Projects/Carnea_Thanatosis/data/dur_beta_JAGS_model.txt")
cat(
  "model{
  
  # priors on beta
  
  b0 ~ dnorm(0, 0.0001)
  b.tmnt_starve ~ dnorm(0, 0.0001)
  b.tmnt_junk ~ dnorm(0, 0.0001)
  b.sex ~ dnorm(0, 0.0001)
  b.int_junk ~ dnorm(0, 0.0001)
  b.int_starve ~ dnorm(0, 0.0001)
  b.tmnt_one ~ dnorm(0, 0.0001)
  b.tmnt_two ~ dnorm(0, 0.0001)
  b.int_one ~ dnorm(0, 0.0001)
  b.int_two ~ dnorm(0, 0.0001)
  
  # prior on shape parameters for beta
  
  phi ~ dunif(0,10)
  
  # likelihood
  
  for (i in 1:enn) {
  
  y[i] ~ dbeta(shp1[i], shp2[i])
  shp1[i] <- mu[i] * phi
  shp2[i] <- (1 - mu[i]) * phi
  logit(mu[i]) <- eta[i]
  
  # eta[i] <- b0 +
  #           b.tmnt_starve * tmnt_starve[i] +
  #           b.tmnt_junk * tmnt_junk[i] +
  #           b.sex * sex[i] +
  #           b.int_starve * tmnt_starve[i] * sex[i] + 
  #           b.int_junk * tmnt_junk[i] * sex[i] 

  eta[i] <- b0 +
            b.tmnt_one * tmnt_one[i] +
            b.tmnt_two * tmnt_two[i] +
            b.sex * sex[i] +
            b.int_one * tmnt_one[i] * sex[i] + 
            b.int_two * tmnt_two[i] * sex[i] 

  }
 
}",
fill=TRUE)
sink()
#### Run JAGS Model: Duration ####

#parms <- c("b0", "b.tmnt_starve", "b.tmnt_junk", "b.sex", "fit", "fitnew", "e", "exp")
parms <- c("b0", 
           "b.tmnt_starve", 
           "b.tmnt_junk", 
           "b.sex", 
           "fit", "fitnew", "e", "exp")
parms <- c("b0",
           "b.tmnt_one", 
           "b.tmnt_two", 
           "b.sex",
           "b.int_one", 
           "b.int_two")

# inits <- function() {
#   list(b0=rnorm(1, 0, 0.01),
#        b.tmnt_one=rnorm(1, 0, 0.01),
#        b.tmnt_two=rnorm(1, 0, 0.01))
# }

inits <- function() {list("phi"=runif(1, 0, 10))}


# data <- list(y = td_feign$feign, 
#              tmnt_starve=td_feign$starve.tmnt,
#              tmnt_junk=td_feign$junk.tmnt,
#              c.sex = td_feign$c.sex, 
#              enn = nrow(td_feign))

data <- list(y = td_dur$b_time, 
             tmnt_one=td_dur$tmnt.one,
             tmnt_two=td_dur$tmnt.two,
             sex=td_dur$f.sex,
             enn = nrow(td_dur))
#data <- list(y = td_dur$b_time, 
#             enn = nrow(td_dur))

data$y <- ifelse(data$y==1, 0.999999, data$y)
mod3 <- jags(model.file = "~/Dropbox/Projects/Carnea_Thanatosis/data/dur_beta_JAGS_model.txt",
             parameters.to.save = parms,
             data = data,
             n.iter = 4000, 
             n.burnin = 1000, 
             n.thin = 2,
             n.chains=3)
             
mod3 <- update(mod3, n.iter = 20000, n.thin = 10)

#### Validation: Duration ####
mod3 <- mod2
out3 <- mod3$BUGSoutput
parms2print <- parms[1:6]
print(out3$summary[parms2print,], digits=3)

## GLM beta frequentist

an1 <- betareg(b_time ~ f.sex * tmnt, data=td_dur)
summary(an1)
# MCMC validation

traceplot(mod3, varname=parms2print)

# overdispersion

ep <- mod3$BUGSoutput$mean$e
enn <- nrow(td_prod)
pm <- 6
disp <- sum(ep^2)/(enn - pm) # 1.04! problem solved

# model fit?
mean(mod3$BUGSoutput$sims.list$fit > mod3$BUGSoutput$sims.list$fitnew) # 0.617

#### Plot Results: Duration ####

# plot posterior

b0_post <- mod1$BUGSoutput$sims.list$b0
b.tmnt_starve_post <- mod1$BUGSoutput$sims.list$b.tmnt_starve
b.sex_post <- mod1$BUGSoutput$sims.list$b.sex

b.tmnt_one_post <- out1$sims.list$b.tmnt_one
b.tmnt_two_post <- out1$sims.list$b.tmnt_two

b.int_on_post <- out3$sims.list$b.int_one

hist(b0_post)
hist(b.tmnt_starve_post)
hist(b.sex_post)
hist(b.int_on_post)

# sketch fitted values

mydata <- data.frame(one=c(0, 1, 0), two=c(0, 0, 1))
mat <- model.matrix( ~ one + two, data=mydata)
betas <- c(out1$mean$b0, out1$mean$b.tmnt_one, out1$mean$b.tmnt_two)
mu <- inv.logit(mat %*% betas)

# another way with 95% CIs

mcmc.betas <- cbind(b0_post, b.tmnt_one_post, b.tmnt_two_post)
mcmc.mus <- inv.logit(mat %*% t(mcmc.betas))

plotsum <- matrix(nrow=nrow(mcmc.mus), ncol=4)

for (i in 1:nrow(mcmc.mus)) {
  
  xi <- mcmc.mus[i,]
  plotsum[i, 3:4] <- quantile(xi, probs=c(0.025, 0.975))
  plotsum[i, 1] <- mean(xi)
  plotsum[i, 2] <- sd(xi)
  
}

colnames(plotsum) <- c("mean", "sd", "2.5%", "97.5%")
rownames(plotsum) <- c("0", "1", "2")

bp1 <- barplot(height=plotsum[,"mean"], ylim=c(0, 1), 
               xlab = "days starved", ylab = "proportion feigning death", las=1, 
               cex.axis=1.3, cex=1.3)
arrows(x0=bp1, 
       y0=plotsum[,"2.5%"], 
       x1=bp1,
       y1=plotsum[,"97.5%"], 
       angle=90, code=3, lwd=1.3)
text(x=bp1, y=0.05, labels=c("A", "B", "B"), cex=2)
