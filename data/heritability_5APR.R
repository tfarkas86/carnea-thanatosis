##########################################################
# Data analysis for heritabilty of death feigning 
# Katie Taylor and Tim Farkas April 2017
##########################################################

rd <- read.csv("~/Dropbox/Projects/Carnea_Thanatosis/data/DeathFeignMaster_9APR.csv")
rd$pdiet <- ifelse(is.na(rd$crossf), rd$diet, ifelse(rd$crossf < 200, 1, ifelse(rd$crossf < 300, 2, 3)))
rd <- rd[rd$generation!="N",]
rd$feign01 <- ifelse(rd$feign=="N", 0, 1)

# remove non-familial cases ("") and half-sib families
dd <- rd[!rd$family %in% c("", "E", "J", "L", "M", "N", "Q", "KK"),]

##### fed group only #####
fd <- dd[dd$pdiet==1, ]

# get # feigned or not for each family offspring

f1 <- aggregate(fd$feign01[fd$generation=="O"], by=list(fd$family[fd$generation=="O"]), FUN=sum)
f0 <- aggregate(fd$feign01[fd$generation=="O"], by=list(fd$family[fd$generation=="O"]), 
                FUN=function(x) length(x) - sum(x))
sires <- aggregate(fd$crossm[fd$generation=="O"], by=list(fd$family[fd$generation=="O"]), FUN=mean)
dams <- aggregate(fd$crossf[fd$generation=="O"], by=list(fd$family[fd$generation=="O"]), FUN=mean)

fd2 <- data.frame(fam=f1[,1], yes=f1[,2], no=f0[,2], sires=sires[,2], dams=dams[,2])

sires01 <- fd$feign01[match(fd2$sires, fd$ID)]
dams01 <- fd$feign01[match(fd2$dams, fd$ID)]

fd2$sires01 <- sires01
fd2$dams01 <- dams01
fd2$mid <- rowMeans(cbind(fd2$sires01, fd2$dams01))
#fd2 <- fd2[fd2$mid!=.5,]

# logistic regresion with single parent

# sire (one-parent regression)
an1 <- glm(cbind(fd2$yes, fd2$no) ~ fd2$sires01, family=binomial)
summary(an1)

h2(an1) * 2 # h2 = 0.487

# dam (one-parent regression)
an1 <- glm(cbind(fd2$yes, fd2$no) ~ fd2$dams01, family=binomial)
summary(an1)

h2(an1) * 2 # h2 = 0.0179

# mid-parent regression

an3  <- glm(cbind(fd2$yes, fd2$no) ~ fd2$mid, family=binomial)
summary(an3)

h2(an3) # h2 = 0.15
fd2

##### 1-day starved group only ####
fd <- dd[dd$pdiet==2, ]

# get # feigned or not for each family offspring

f1 <- aggregate(fd$feign01[fd$generation=="O"], by=list(fd$family[fd$generation=="O"]), FUN=sum)
f0 <- aggregate(fd$feign01[fd$generation=="O"], by=list(fd$family[fd$generation=="O"]), 
                FUN=function(x) length(x) - sum(x))
sires <- aggregate(fd$crossm[fd$generation=="O"], by=list(fd$family[fd$generation=="O"]), FUN=mean)
dams <- aggregate(fd$crossf[fd$generation=="O"], by=list(fd$family[fd$generation=="O"]), FUN=mean)

fd2 <- data.frame(fam=f1[,1], yes=f1[,2], no=f0[,2], sires=sires[,2], dams=dams[,2])

sires01 <- fd$feign01[match(fd2$sires, fd$ID)]
dams01 <- fd$feign01[match(fd2$dams, fd$ID)]

fd2$sires01 <- sires01
fd2$dams01 <- dams01
fd2$mid <- rowMeans(cbind(fd2$sires01, fd2$dams01))
#fd2 <- fd2[fd2$mid!=.5,]

# logistic regresion with single parent

# sire (one-parent regression)
an1 <- glm(cbind(fd2$yes, fd2$no) ~ fd2$sires01, family=binomial)
summary(an1)

h2(an1) * 2 # h2 = 0.487

# dam (one-parent regression)
an1 <- glm(cbind(fd2$yes, fd2$no) ~ fd2$dams01, family=binomial)
summary(an1)

h2(an1) * 2 # h2 = 0.0179

# mid-parent regression

an3  <- glm(cbind(fd2$yes, fd2$no) ~ fd2$mid, family=binomial)
summary(an3)

h2(an3) # h2 = 0.15
fd2


##### 2-day starved group only #####
fd <- dd[dd$pdiet==3, ]

# get # feigned or not for each family offspring

f1 <- aggregate(fd$feign01[fd$generation=="O"], by=list(fd$family[fd$generation=="O"]), FUN=sum)
f0 <- aggregate(fd$feign01[fd$generation=="O"], by=list(fd$family[fd$generation=="O"]), 
                FUN=function(x) length(x) - sum(x))
sires <- aggregate(fd$crossm[fd$generation=="O"], by=list(fd$family[fd$generation=="O"]), FUN=mean)
dams <- aggregate(fd$crossf[fd$generation=="O"], by=list(fd$family[fd$generation=="O"]), FUN=mean)

fd2 <- data.frame(fam=f1[,1], yes=f1[,2], no=f0[,2], sires=sires[,2], dams=dams[,2])

sires01 <- fd$feign01[match(fd2$sires, fd$ID)]
dams01 <- fd$feign01[match(fd2$dams, fd$ID)]

fd2$sires01 <- sires01
fd2$dams01 <- dams01
fd2$mid <- rowMeans(cbind(fd2$sires01, fd2$dams01))
#fd2 <- fd2[fd2$mid!=.5,]

# logistic regresion with single parent

# sire (one-parent regression)
an1 <- glm(cbind(fd2$yes, fd2$no) ~ fd2$sires01, family=binomial)
summary(an1)

h2(an1) * 2 # h2 = 0.487

# dam (one-parent regression)
an1 <- glm(cbind(fd2$yes, fd2$no) ~ fd2$dams01, family=binomial)
summary(an1)

h2(an1) * 2 # h2 = 0.0179

# mid-parent regression

an3  <- glm(cbind(fd2$yes, fd2$no) ~ fd2$mid, family=binomial)
summary(an3)

h2(an3) # h2 = 0.15
fd2

##### Animal Model #####

library(MCMCglmm)



