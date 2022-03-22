#####################################################
##### Carnea ploribunda thanatosis pilot study ######
#####################################################
##### performed by K. Taylor and T. Farkas ##########
#####################################################
##### September 7, 2016 #############################
#####################################################

##### input and data manipulation #####

pd <- read.csv("~/Dropbox/Projects/Carnea_Thanatosis/ploribunda_pilot.csv")

pd <- pd[-c(1,2,4),] # remove eraser method trials
pd[1,1] <- "002" # rename 002b
pd <- pd[,-2] # remove method column
pd$prods_1 <- ifelse(pd$prods_1 >= 20, 20, pd$prods_1)

# make long form data for analysis
pd_long <- reshape(pd[pd$prods_1 < 20,], varying=list(c(2, 4, 6), c(3, 5, 7)), 
                   direction="long",
                   idvar="id")

##### summarize data #####

# histogram of prods needed to elicit death feigning
pdf(file="~/Dropbox/Projects/Carnea_Thanatosis/Figures/Cp_pilot_prods_hist.pdf")
h_prods <- hist(pd$prods_1, breaks=-1:20+.5, main="prods at first trial", 
                ylab="number of individuals", xlab="number of prods")
dev.off()

# histogram of duration death feigning
pdf(file="~/Dropbox/Projects/Carnea_Thanatosis/Figures/Cp_pilot_dur_hist.pdf")
h_dur <- hist(pd$dur_1[pd$dur_1 > -1], breaks=seq(-20, 220, by=20), main="duration", 
              xlab="duration immobile (seconds)", ylab="number of individuals")
dev.off()

# trends through time for number of prods and duration
avgs <- colMeans(pd[pd$prods_1<20,2:7])
ses <- apply(pd[pd$prods_1 < 20, 2:7], MARGIN=2, 
             FUN=function(x) (sd(x)/sqrt(length(x)-1)))

p_avgs <- avgs[c(1,3,5)]
p_ses <- ses[c(1,3,5)]
d_avgs <- avgs[c(2,4,6)]
d_ses <- ses[c(2,4,6)]

# barplots

pdf(file="~/Dropbox/Projects/Carnea_Thanatosis/Figures/Cp_pilot_prods_repeat.pdf")
bp1 <- barplot(p_avgs, ylim=c(0, 8), xlab="trial", ylab="number of prods")
arrows(bp1, p_avgs + p_ses, bp1, p_avgs - p_ses, angle=90, code=3)
dev.off()

pdf(file="~/Dropbox/Projects/Carnea_Thanatosis/Figures/Cp_pilot_dur_repeat.pdf")
bp2 <- barplot(d_avgs, ylim=c(0,70), xlab="trial", ylab="duration of immobility (seconds)")
arrows(bp2, d_avgs + d_ses, bp2, d_avgs - d_ses, angle=90, code=3)
dev.off()

##### Analysis #####

# habituation of prods
an1 <- lm(prods_1 ~ factor(time), data=pd_long)
summary(an1)
Anova(an1) # p = 0.535

# habituation of duration
an1 <- lm(dur_1 ~ factor(time), data=pd_long)
summary(an1)
Anova(an1) # p = 0.876

# repeatability for prods

an1 <- lm(prods_1 ~ time + id, data=pd_long)
summary(an1)
aov1 <- anova(an1)

btw <- aov1$`Mean Sq`[1]
wthn <- (aov1$`Mean Sq`[1] - aov1$`Mean Sq`[2]) / 3
rep <- wthn/(wthn + btw) # 0.0136 = not repeatable

# repeatability for duration

an1 <- lm(dur_1 ~ id, data=pd_long)
summary(an1)
aov1 <- anova(an1)

btw <- aov1$`Mean Sq`[1]
wthn <- (aov1$`Mean Sq`[1] - aov1$`Mean Sq`[2]) / 3
rep <- wthn/(wthn + btw) # 0.075 = more repeatable



