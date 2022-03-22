#######################################
#### Test script to understand how to 
#### estimate heritability with binary 
#### data. T. Farkas Dec 2016.
#######################################

# h2.sim function is attempt to generate phenotypic data for parents and their 
# offspring given a known heritability for a single trait. use can specify the 
# number of parents (enn), average breeding value for parents (mu.bv), standard
# deviation for the breeding values (sd.bv), narrow-sense heritability (h2), 
# and number of offspring per breeding pair (noff). parents are randomly paired
# in a single mating event. function returns the phenotypic data for all 
# individuals, a pedigree for all individuals, and the empirically calculated 
# h2. when noff is low, sampling error can lead to h2 > 1.

h2.sim <- function(enn = 100, mu.bv = 100, sd.bv=mu.bv/100, h2=.5, noff=5) {
  
  bvs.p <- rnorm(enn, mu.bv, sd.bv) # breeding values
  va.p <- var(bvs.p) # additive genetic variance
  vp1.p <- va.p / h2 # phenotypic variance
  ve.p <- vp1.p - va.p # enivronmental variance
  y.p <- bvs.p + rnorm(enn, sd=sqrt(ve.p)) # parent phenotypes
  
  vp.p <- var(y.p)
  e.h2 <- va.p / vp.p
  e.h2
  y.p.m <- bvs.p
  m <- 1:length(y.p.m)
  p.match <- matrix(ncol=2, nrow=length(bvs.p)/2)
  
  for (i in 1:(length(bvs.p)/2)) {
    
    p.m <- sample(m, size=2)
    m <- m[-c(p.m)]
    p.match[i,] <- p.m
    
  }
  
  # matches <- as.data.frame(matrix(bvs.p[p.match], ncol=2, byrow=FALSE))
  # matches$f1.hat <- rowMeans(matches)
  # matches <- data.frame(matches, as.data.frame(matrix(y.p[p.match], ncol=2, byrow=FALSE)))
  # matches <- data.frame(pair.id=1:nrow(matches), matches)
  # colnames(matches) <- c("pair.id", "bv.p1", "bv.p2", "f1.hat", "y.p1", "y.p2")
  # 
  # f1.hat <- matrix(matches$f1.hat, ncol=1)
  # f1.hats <- as.vector(t(do.call(cbind, replicate(noff, f1.hat, simplify=FALSE))))
  # #f1.ps <- f1.hats + ifelse(h2==1, 0, rnorm(noff, mean=0, sd=sqrt(ve.p)))
  # 
  # all.p <- c(y.p, f1.ps)
  # all.p <- data.frame(animal=1:length(all.p), y = all.p)
  
  sires <- as.vector(t(do.call(cbind, replicate(noff, p.match[,1], simplify=FALSE))))
  dams <- as.vector(t(do.call(cbind, replicate(noff, p.match[,2], simplify=FALSE))))
  
  all.ped <- cbind(sires, dams)
  na.p <- matrix(NA, nrow=enn, ncol=2)
  all.ped <- rbind(na.p, all.ped)
  all.ped <- cbind(1:nrow(all.ped), all.ped)
  
  #return(list(ys=all.p, ped=all.ped, h2=e.h2))
  return(all.ped)
  
}


size_bias <- lapply(c(10, 20, 40, 80), function(q) {
  
  h2bias <- sapply(seq(0, 1, by=.1), function(z) {
    
    h2bars <- sapply(1:10, function(why) {
      
      h2s <- sapply(1:100, function(x) {
        
        ped <- h2.sim(enn=q, noff=10, h2=z)
        K <- 2 * kinship2::kinship(ped[, 1], ped[, 2], ped[, 3])
        Va <- z
        A <- K * Va
        a <- rmvnorm(1, sigma= A)
        e  <- rmvnorm(1, sigma= (1-Va) * diag(length(a)))
        mu <- 0
        y <- mu + a + e 
        
        e.h2 <- var(t(a))/var(t(y))
        return(e.h2)
      })
      return(mean(h2s))
    })
    return(mean(h2bars))
  })
  return(h2bias)
})

names(size_bias) <- colnames(size)

save(size_bias, file="size_bias.RData")
  
plot(seq(0, 1, by=.1), y=h2bias)
abline(a=0, b=1, lty=2)








# h2.sim function is attempt to generate phenotypic data for parents and their 
# offspring given a known heritability for a single trait. use can specify the 
# number of parents (enn), average breeding value for parents (mu.bv), standard
# deviation for the breeding values (sd.bv), narrow-sense heritability (h2), 
# and number of offspring per breeding pair (noff). parents are randomly paired
# in a single mating event. function returns the phenotypic data for all 
# individuals, a pedigree for all individuals, and the empirically calculated 
# h2. when noff is low, sampling error can lead to h2 > 1.

h2.sim <- function(enn = 100, mu.bv = 100, sd.bv=mu.bv/100, h2=.5, noff=5) {
  
  bvs.p <- rnorm(enn, mu.bv, sd.bv) # breeding values
  va.p <- var(bvs.p) # additive genetic variance
  vp1.p <- va.p / h2 # phenotypic variance
  ve.p <- vp1.p - va.p # enivronmental variance
  y.p <- bvs.p + rnorm(enn, sd=sqrt(ve.p)) # parent phenotypes
  
  vp.p <- var(y.p)
  e.h2 <- va.p / vp.p
  e.h2
  y.p.m <- bvs.p
  m <- 1:length(y.p.m)
  p.match <- matrix(ncol=2, nrow=length(bvs.p)/2)
  
  for (i in 1:(length(bvs.p)/2)) {
    
    p.m <- sample(m, size=2)
    m <- m[-c(p.m)]
    p.match[i,] <- p.m
    
  }
  
  # matches <- as.data.frame(matrix(bvs.p[p.match], ncol=2, byrow=FALSE))
  # matches$f1.hat <- rowMeans(matches)
  # matches <- data.frame(matches, as.data.frame(matrix(y.p[p.match], ncol=2, byrow=FALSE)))
  # matches <- data.frame(pair.id=1:nrow(matches), matches)
  # colnames(matches) <- c("pair.id", "bv.p1", "bv.p2", "f1.hat", "y.p1", "y.p2")
  # 
  # f1.hat <- matrix(matches$f1.hat, ncol=1)
  # f1.hats <- as.vector(t(do.call(cbind, replicate(noff, f1.hat, simplify=FALSE))))
  # #f1.ps <- f1.hats + ifelse(h2==1, 0, rnorm(noff, mean=0, sd=sqrt(ve.p)))
  # 
  # all.p <- c(y.p, f1.ps)
  # all.p <- data.frame(animal=1:length(all.p), y = all.p)
  
  sires <- as.vector(t(do.call(cbind, replicate(noff, p.match[,1], simplify=FALSE))))
  dams <- as.vector(t(do.call(cbind, replicate(noff, p.match[,2], simplify=FALSE))))
  
  all.ped <- cbind(sires, dams)
  na.p <- matrix(NA, nrow=enn, ncol=2)
  all.ped <- rbind(na.p, all.ped)
  all.ped <- cbind(1:nrow(all.ped), all.ped)
  
  #return(list(ys=all.p, ped=all.ped, h2=e.h2))
  return(all.ped)
  
}

matches <- h2.sim(h2=.6)

h2.sim(h2=1)
#### MCMCglmm with BTdata ####

bd <- BTdata
bp <- BTped

prior <- list(R = list(V=1, nu=0.002), G = list(G1 = list(V=1, nu=0.002)))
mod <- MCMCglmm(Trait ~ 1, random=~ID, family="gaussian", prior=prior, 
                pedigree=ped[,1:3], data=ped, nitt=10000, burnin=1000, thin=10)

plot(mod$Sol)
plot(mod$VCV)
autocorr.diag(mod$Sol)
autocorr.diag(mod$VCV)
effectiveSize(mod$Sol)
effectiveSize(mod$VCV)
heidel.diag(mod$VCV)

summary(mod)

herit <- mod$VCV[,"animal"] / (mod$VCV[,"animal"] + mod$VCV[,"units"])

effectiveSize(herit)
HPDinterval(herit)
mean(herit)

#### MCMCglmm with simulated data ####
sims <- h2.sim(enn=200, h2=.5, noff=20)
sims$h2

bd <- sims$ys
bp <- sims$ped

prior <- list(R = list(V=1, nu=0.002), G = list(G1 = list(V=1, nu=0.002)))
mod <- MCMCglmm(y ~ 1, random=~animal, family="gaussian", prior=prior, 
                pedigree=bp, data=bd, nitt=10000, burnin=1000, thin=10)

plot(mod$Sol)
plot(mod$VCV)
autocorr.diag(mod$Sol)
autocorr.diag(mod$VCV)
effectiveSize(mod$Sol)
effectiveSize(mod$VCV)
heidel.diag(mod$VCV)

summary(mod)

herit <- mod$VCV[,"animal"] / (mod$VCV[,"animal"] + mod$VCV[,"units"])

effectiveSize(herit)
HPDinterval(herit)
mean(herit)
sims$h2


#### Try it with JAGS ####

sink("~/Dropbox/Projects/Carnea_Thanatosis/jags_animal_model.txt")
cat(
  "model{

  # half_Cauchy hyperprior for additive genetic variance (Va) and environmental 
  # variance
    
    Va ~ dunif(0,p.var)
    #gam ~ dgamma(0.0001, 0.0001)
    #Va <- 1/gam
    #num_Va ~ dnorm(0, 0.0016)
    #den_Va ~ dnorm(0, 1)
    #sig_Va <- abs(num_Va/den_Va)
    #Va <- pow(sig_Va, 2)
    
    var2_a <- A * Va
    tau_a <- inverse(var2_a)

    a[1:enn] ~ dmnorm(mu[], tau_a[,])

    Ve ~ dunif(0, p.var)
    #gam2 ~ dgamma(0.0001, 0.0001)
    #Ve <- 1 / gam2
    #num_Ve ~ dnorm(0, 0.0016)
    #den_Ve ~ dnorm(0, 1)
    #sig_Ve <- abs(num_Ve/den_Ve)
    #Ve <- pow(sig_Ve, 2)
    #var2_Ve <- eye * var_Ve
    #tau_e <- inverse(var_e)
    tau_e <- 1 / Ve
    
  # prior on intercept (mean)
    
    b0 ~ dnorm(p.mu, 0.001)

  # likelihood 

    for (i in 1:enn) {
      
      y[i] ~ dnorm(yhat[i], tau_e)
      yhat[i] <- b0 + a[i]  # animal model 

    }

  # derived quantities
      
    Vp <- Va + Ve # total phenotypic variance
    h2 <- Va/Vp # heritability

  }", 
  fill=TRUE)
sink()

#### Run JAGS animal model on simulated data

sims <- h2.sim(enn=50, h2=.5, noff=10)
sims$h2

bd <- sims$ys
bp <- sims$ped

A <- 2 * kinship2::kinship(bp[, 1], bp[, 2], bp[, 3])

eye <- diag(nrow(bp))
mu <- rep(0, nrow(A))

data <- list(y=bd$y, 
             A=A,
             #eye=eye,
             mu= mu,
             p.mu=mean(bd$y),
             p.var=var(bd$y),
             enn=nrow(bp))

parms <- c("b0", "Va", "Ve", "h2")

inits <- function() {list(b0     = rnorm(1, mean(bd$y), 0.01),
                          Vp = rnorm(1, var(bd$y)/2))}

time <- system.time(
  mod <- jags(model.file="~/Dropbox/Projects/Carnea_Thanatosis/jags_animal_model.txt",
              data=data, 
              parameters.to.save = parms,
              n.iter=10000, 
              n.burnin=1000,
              n.thin=10, 
              n.chains=2, 
              #inits=inits
  )
)

traceplot(mod)

# try it with MCMCglmm
prior <- list(R = list(V=1, nu=0.002), G = list(G1 = list(V=1, nu=0.002)))
mod2 <- MCMCglmm(y ~ 1, random=~animal, family="gaussian", prior=prior, 
                 pedigree=bp, data=bd, nitt=10000, burnin=1000, thin=10)

plot(mod$Sol)
plot(mod$VCV)
autocorr.diag(mod$Sol)
autocorr.diag(mod$VCV)
effectiveSize(mod$Sol)
effectiveSize(mod$VCV)
heidel.diag(mod$VCV)

summary(mod2)

herit <- mod2$VCV[,"animal"] / (mod2$VCV[,"animal"] + mod2$VCV[,"units"])

effectiveSize(herit)
HPDinterval(herit)
mean(herit)

