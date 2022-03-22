############################
#### Test of MVN vs.looping in JAGS
#### Dec 2016, T Farkas
############################

sma <- matrix(.5, nrow=2, ncol=2)
diag(sma) <- c(2, .5)

dd <- rmvnorm(1, sigma=sma)


library(R2jags)
#### Simulate data

y <- rnorm(100, 10)

#### Looping JAGS model

sink(file="~/Dropbox/Projects/Carnea_Thanatosis/mvn_test_model_loop.txt")
cat(
  "model{

    # prior on the variance

      sigma2 ~ dgamma(0.0001, 0.0001)
      tau <- 1 / sigma2

    # prior on the mean
      
      #ybar <- mean(y[])
      b0 ~ dnorm(ybar, 0.0001)

    # likelihood

      for (i in 1:enn) {
  
      y[i] ~ dnorm(mu[i], tau)
      mu[i] <- b0

      }

  }",
  fill=TRUE
  )
sink()

data <- list(y=y, 
             ybar=mean(y), 
             enn=length(y)
             )

parms <- c("b0")

time.l <- system.time(
mod.l <- jags(model.file="~/Dropbox/Projects/Carnea_Thanatosis/mvn_test_model_loop.txt",
              data=data, 
              parameters.to.save=parms, 
              n.iter=10000, 
              n.burnin=1000,
              n.thin = 10,
              n.chains=3)
)

#### MVN JAGS model: Super duper slow (67x slower!), but low deviance

sink(file="~/Dropbox/Projects/Carnea_Thanatosis/mvn_test_model_mvn.txt")
cat(
  "model{
  
  # prior on the variance
  
  sigma2 ~ dgamma(0.0001, 0.0001)
  tau <- 1 / sigma2
  tau.mat <- tau * eye
  
  # prior on the mean
  
  #ybar <- mean(y[])
  b0 ~ dnorm(ybar, 0.0001)
  
  # likelihood
  
  y[] ~ dmnorm(mu[], tau.mat[,] )
  mu <- b0 * ones
  
  }",
  fill=TRUE
)
sink()

data <- list(y=y, 
             ybar=mean(y), 
             #enn=length(y),
             eye=diag(length(y)),
             ones=matrix(1, nrow=length(y))
)

parms <- c("b0")

time.mvn <- system.time(
  mod.mvn <- jags(model.file="~/Dropbox/Projects/Carnea_Thanatosis/mvn_test_model_mvn.txt",
                data=data, 
                parameters.to.save=parms, 
                n.iter=10000, 
                n.burnin=1000,
                n.thin = 10,
                n.chains=3)
)

#### Error Out JAGS model: Not quite ...

sink(file="~/Dropbox/Projects/Carnea_Thanatosis/mvn_test_model_errorout.txt")
cat(
  "model{
  
  # prior on the variance
  
  sigma2 ~ dgamma(0.0001, 0.0001)
  tau <- 1 / sigma2
  
  # prior on the mean
  
  #ybar <- mean(y[])
  b0 ~ dnorm(ybar, 0.0001)

  # prior on the error

  for (i in 1:enn) {

  e[i] ~ dnorm(0, tau)
  err[i] <- y[i] - b0
  }
  
  }",
  fill=TRUE
)
sink()

data <- list(y=y, 
             ybar=mean(y), 
             enn=length(y)
)

parms <- c("b0")

time.eout <- system.time(
  mod.eout <- jags(model.file="~/Dropbox/Projects/Carnea_Thanatosis/mvn_test_model_errorout.txt",
                data=data, 
                parameters.to.save=parms, 
                n.iter=10000, 
                n.burnin=1000,
                n.thin = 10,
                n.chains=3)
)

