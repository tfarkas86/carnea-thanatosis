nd <- datafedonly %>% 
  select(animal, feign, mom) 

reps <- 100 # number of simulations

priorm <- list(R = list(V = 1, fix = 1), # priors
               G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1), 
                        G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

# a matrix to collect model results
out <- matrix(NA, nrow=reps, ncol=4) 

# loop to simulate
for(i in 1:reps) {
  
  # permute the motherhood relationships
  nd$mom <- sample(nd$mom, size=nrow(nd), replace=FALSE)
  
  # run model
  modelm2<-MCMCglmm(feign ~ 1, random= ~animal+mom,
                    family="categorical", prior=priorm, pedigree=ped, data=nd,
                    nitt = 11000, burnin = 1000, thin = 10, verbose=FALSE)
  
  
  an_mean <- mean(modelm2$VCV[, "animal"]) # mean for Va
  an_low <- quantile(modelm2$VCV[, "animal"], prob=0.025) # lower 0.025%
  an_high <- quantile(modelm2$VCV[, "animal"], prob=0.975) # upper 0.975%
  DIC <- modelm2$DIC # DIC
  
  out[i, ] <- cbind(an_mean, an_low, an_high, DIC)
  
}

save(sim_maternal, file="~/Dropbox/Projects/Carnea_Thanatosis/data/sim_maternal.RData")
