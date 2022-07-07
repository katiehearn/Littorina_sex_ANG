library(bbmle)
library(hzar)
library(tidyverse)

#### fitting inversion clines

## data already loaded from previous scripts or will need to be read in 
# set working directory
region1_females_invfreqscol <- read.csv("region1_females_invfreqs.csv")
region2a_females_invfreqscol <- read.csv("region2a_females_invfreqs.csv")
region2b_females_invfreqscol <- read.csv("region2b_females_invfreqs.csv")
region3_females_invfreqscol <- read.csv("region3_females_invfreqs.csv")
region1_females <- read.csv("region1_females_PCandcluster.csv")
region2a_females <- read.csv("region2a_females_PCandcluster.csv")
region2b_females <- read.csv("region2b_females_PCandcluster.csv")
region3_females <- read.csv("region3_females_PCandcluster.csv")
region1_males_invfreqscol <- read.csv("region1_males_invfreqs.csv")
region2a_males_invfreqscol <- read.csv("region2a_males_invfreqs.csv")
region2b_males_invfreqscol <- read.csv("region2b_males_invfreqs.csv")
region3_males_invfreqscol <- read.csv("region3_males_invfreqs.csv")
region1_males <- read.csv("region1_males_PCandcluster.csv")
region2a_males <- read.csv("region2a_males_PCandcluster.csv")
region2b_males <- read.csv("region2b_males_PCandcluster.csv")
region3_males <- read.csv("region3_males_PCandcluster.csv")



#### clean copy of cline script

cline <- function(hap,pos,c,lw,lp_c,lp_w){  # using logit frequencies and log(width), to avoid boundary issues
  d <- pos-c
  
  p_x <- 1/(1+exp(0-4*(d/exp(lw))))  
  
  crab <- exp(lp_c)/(1+exp(lp_c)) # converts logit(freq) back to frequency
  wave <- exp(lp_w)/(1+exp(lp_w))
  z_x <- crab + (wave-crab)*p_x  
  # z_x is expected frequency
  
  minusll <- -sum(dbinom(hap,2,z_x,log=T))
  
  return(minusll)
}

theta.init <- list(c=left,lw=log(10),lp_c=start_c,lp_w=start_w) # start_c and start_w need to be logit scale, ‘left’ is just the environmental transition or the expected cline centre

mle.cline <- mle2(cline, theta.init,  
                  control=list(parscale=abs(unlist(theta.init))),
                  data=list(pos=dist,hap=inv))

summary(mle.cline)
confint(mle.cline)
aic_cline <- AIC(mle.cline)


#################################################################

#### next copy to edit and use

## region1 female
ggplot(region1_females_invfreqscol, aes(x=mean_dist, y=prop_x)) + geom_point() + scale_y_continuous(limits=c(0,1))

theta.init <- list(c=77.8,lw=log(10),lp_c=qlogis(0.001), lp_w=qlogis(0.32)) # start_c and start_w need to be logit scale, ‘left’ is just the environmental transition or the expected cline centre
#theta.init <- list(c=133.2,lw=-0.81,lp_c=0.33,lp_w=2.2) # alternative if needed for a second fit – used if confint profiling finds a better fit than the original fit

### need to change coding to 0,1,2 and also remove any NAs 
region1_females <- region1_females %>% mutate(cluster2=ifelse(cluster==1, 0,
                                                              ifelse(cluster==2, 1,
                                                                     ifelse(cluster==3, 2, NA))))
region1_females_nona <- region1_females %>% filter(!is.na(DistAlongPath))


mle.cline <- mle2(cline, theta.init,  
                  control=list(parscale=abs(unlist(theta.init))),
                  data=list(pos=region1_females_nona$DistAlongPath, hap=as.numeric(region1_females_nona$cluster2)))

summary(mle.cline)
confint(mle.cline)
aic_cline <- AIC(mle.cline)

# find fitted cline to plot and to get residuals
pars.m <- coef(mle.cline) # get coefficients for fitted model
# back transform the parameters
#pars.m[2] <- exp(pars.m[2])
#pars.m[3] <- exp(pars.m[3])/(1+exp(pars.m[3]))
#pars.m[4] <- exp(pars.m[4])/(1+exp(pars.m[4]))
# save paramaters
r1_f_pars <- pars.m
r1_f_summary <- summary(mle.cline)
r1_f_confint <- confint(mle.cline)
r1_f_aiccline <- AIC(mle.cline)

pos <- c(0:160) # ordered distance vector (or you can use actual snail distances if you want to find residuals)
dist.m <- pos-pars.m[1]
freq.m <- 1/(1+exp(0-4*(dist.m)/(exp(pars.m[2]))))
mean.m <- (exp(pars.m[3])/(1+exp(pars.m[3])))+((exp(pars.m[4])/(1+exp(pars.m[4])))-(exp(pars.m[3])/(1+exp(pars.m[3]))))*freq.m # fitted phenotype mean
#mean.m <- pars.m[3]+(pars.m[4]-pars.m[3])*freq.m
region1_female_fitted <- data.frame(dist.m, mean.m)
ggplot(region1_female_fitted, aes(x=dist.m+87.5, y=1-mean.m)) + geom_line() + 
  scale_x_continuous(limits=c(0,160)) +
  scale_y_continuous(limits=c(0,1)) + 
  geom_point(data=region1_females_invfreqscol, aes(x=mean_dist, y=prop_x), alpha=0.5, size=3, colour="deepskyblue4") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank()) +
  xlab("Distance along path") + ylab("Inversion frequency")

r1_f_pars_back <- r1_f_pars[1]
r1_f_pars_back[2] <- exp(r1_f_pars[2])
r1_f_pars_back[3] <- exp(r1_f_pars[3])/(1+exp(r1_f_pars[3]))
r1_f_pars_back[4] <- exp(r1_f_pars[4])/(1+exp(r1_f_pars[4]))

r1_f_confint_back <- data.frame(4,2)
r1_f_confint_back[1,] <- r1_f_confint[1,]
r1_f_confint_back[2,] <- exp(r1_f_confint[2,])
r1_f_confint_back[3,] <- exp(r1_f_confint[3,])/(1+exp(r1_f_confint[3,]))
r1_f_confint_back[4,] <- exp(r1_f_confint[4,])/(1+exp(r1_f_confint[4,]))
###########################################################################

###########################################################################
### region 2a female

ggplot(region2a_females_invfreqscol, aes(x=mean_dist, y=prop_x)) + geom_point() + scale_y_continuous(limits=c(0,1))

theta.init <- list(c=100,lw=log(10),lp_c=qlogis(0.499), lp_w=qlogis(0.55))
region2a_females <- region2a_females %>% mutate(cluster2=ifelse(cluster==1, 0,
                                                              ifelse(cluster==2, 1,
                                                                     ifelse(cluster==3, 2, NA))))
region2a_females_nona <- region2a_females %>% filter(!is.na(DistAlongPath) & !is.na(cluster2))


mle.cline <- mle2(cline, theta.init,  
                  control=list(parscale=abs(unlist(theta.init))),
                  data=list(pos=region2a_females_nona$DistAlongPath, hap=as.numeric(region2a_females_nona$cluster2)))

summary(mle.cline)
confint(mle.cline)
aic_cline <- AIC(mle.cline)

# save new fit from confint
#r2a_f_newfit <- confint(mle.cline)
#summary(r2a_f_newfit)

# find fitted cline to plot and to get residuals
#pars.m <- coef(r2a_f_newfit) # get coefficients for fitted model
pars.m <- coef(mle.cline)
# back transform the parameters
pars.m[2] <- exp(pars.m[2])
pars.m[3] <- exp(pars.m[3])/(1+exp(pars.m[3]))
pars.m[4] <- exp(pars.m[4])/(1+exp(pars.m[4]))
# save paramaters
r2a_f_pars <- pars.m
r2a_f_summary <- summary(mle.cline)
r2a_f_confint <- confint(mle.cline)
r2a_f_aiccline <- AIC(mle.cline)


pos <- c(0:160) # ordered distance vector (or you can use actual snail distances if you want to find residuals)
dist.m <- pos-pars.m[1]
freq.m <- 1/(1+exp(0-4*(dist.m)/(exp(pars.m[2]))))
#mean.m <- (exp(pars.m[3])/(1+exp(pars.m[3])))+((exp(pars.m[4])/(1+exp(pars.m[4])))-(exp(pars.m[3])/(1+exp(pars.m[3]))))*freq.m # fitted phenotype mean
mean.m <- pars.m[3]+(pars.m[4]-pars.m[3])*freq.m
region2a_female_fitted <- data.frame(dist.m, mean.m)
ggplot(region2a_female_fitted, aes(x=dist.m+86.32, y=1-mean.m)) + geom_line() + 
  scale_x_continuous(limits=c(0,160)) +
  scale_y_continuous(limits=c(0,1)) + 
  geom_point(data=region2a_females_invfreqscol, aes(x=mean_dist, y=prop_x), alpha=0.5, size=3, colour="deepskyblue4") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank()) +
  xlab("Distance along path") + ylab("Inversion frequency") +
  geom_vline(xintercept=77.82, colour="red4", linetype=2, alpha=0.5)

r2a_f_pars_back <- r2a_f_pars[1]
#r2a_f_pars_back[2] <- exp(r2a_f_pars[2])
#r2a_f_pars_back[3] <- exp(r2a_f_pars[3])/(1+exp(r2a_f_pars[3]))
#r2a_f_pars_back[4] <- exp(r2a_f_pars[4])/(1+exp(r2a_f_pars[4]))

# not working as cline is so flat?
r2a_f_confint_back <- data.frame(4,2)
r2a_f_confint_back[1,] <- r2a_f_confint[1,]
r2a_f_confint_back[2,] <- exp(r2a_f_confint[2,])
r2a_f_confint_back[3,] <- exp(r2a_f_confint[3,])/(1+exp(r2a_f_confint[3,]))
r2a_f_confint_back[4,] <- exp(r2a_f_confint[4,])/(1+exp(r2a_f_confint[4,]))
###########################################################################

###########################################################################
### region 2b female

ggplot(region2b_females_invfreqscol, aes(x=mean_dist, y=1-prop_x)) + geom_point() + scale_y_continuous(limits=c(0,1))

theta.init <- list(c=85,lw=log(10),lp_c=qlogis(0.499), lp_w=qlogis(0.9))
region2b_females <- region2b_females %>% mutate(cluster2=ifelse(cluster==1, 0,
                                                                ifelse(cluster==2, 1,
                                                                       ifelse(cluster==3, 2, NA))))
region2b_females_nona <- region2b_females %>% filter(!is.na(DistAlongPath) & !is.na(cluster2))


mle.cline <- mle2(cline, theta.init,  
                  control=list(parscale=abs(unlist(theta.init))),
                  data=list(pos=region2b_females_nona$DistAlongPath, hap=as.numeric(region2b_females_nona$cluster2)))

summary(mle.cline)
confint(mle.cline)
aic_cline <- AIC(mle.cline)

# find fitted cline to plot and to get residuals
pars.m <- coef(mle.cline) # get coefficients for fitted model
# back transform the parameters
#pars.m[2] <- exp(pars.m[2])
#pars.m[3] <- exp(pars.m[3])/(1+exp(pars.m[3]))
#pars.m[4] <- exp(pars.m[4])/(1+exp(pars.m[4]))
# save paramaters
r2b_f_pars <- pars.m
r2b_f_summary <- summary(mle.cline)
r2b_f_confint <- confint(mle.cline)
r2b_f_aiccline <- AIC(mle.cline)

pos <- c(0:160) # ordered distance vector (or you can use actual snail distances if you want to find residuals)
dist.m <- pos-pars.m[1]
freq.m <- 1/(1+exp(0-4*(dist.m)/(exp(pars.m[2]))))
mean.m <- (exp(pars.m[3])/(1+exp(pars.m[3])))+((exp(pars.m[4])/(1+exp(pars.m[4])))-(exp(pars.m[3])/(1+exp(pars.m[3]))))*freq.m # fitted phenotype mean
#mean.m <- pars.m[3]+(pars.m[4]-pars.m[3])*freq.m
region2b_female_fitted <- data.frame(dist.m, mean.m)
ggplot(region2b_female_fitted, aes(x=dist.m+90.23, y=1-mean.m)) + geom_line() + 
  scale_x_continuous(limits=c(0,160)) +
  scale_y_continuous(limits=c(0,1)) + 
  geom_point(data=region2b_females_invfreqscol, aes(x=mean_dist, y=prop_x), alpha=0.5, size=3, colour="deepskyblue4") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank()) +
  xlab("Distance along path") + ylab("Inversion frequency") +
  geom_vline(xintercept=77.82, colour="red4", linetype=2, alpha=0.5)

r2b_f_pars_back <- r2b_f_pars[1]
r2b_f_pars_back[2] <- exp(r2b_f_pars[2])
r2b_f_pars_back[3] <- exp(r2b_f_pars[3])/(1+exp(r2b_f_pars[3]))
r2b_f_pars_back[4] <- exp(r2b_f_pars[4])/(1+exp(r2b_f_pars[4]))

r2b_f_confint_back <- data.frame(4,2)
r2b_f_confint_back[1,] <- r2b_f_confint[1,]
r2b_f_confint_back[2,] <- exp(r2b_f_confint[2,])
r2b_f_confint_back[3,] <- exp(r2b_f_confint[3,])/(1+exp(r2b_f_confint[3,]))
r2b_f_confint_back[4,] <- exp(r2b_f_confint[4,])/(1+exp(r2b_f_confint[4,]))
###########################################################################

###########################################################################
### region 3 female

ggplot(region3_females_invfreqscol, aes(x=mean_dist, y=1-prop_x)) + geom_point() + scale_y_continuous(limits=c(0,1))

theta.init <- list(c=85,lw=log(10),lp_c=qlogis(0.6), lp_w=qlogis(0.1))
region3_females <- region3_females %>% mutate(cluster2=ifelse(cluster==1, 0,
                                                                ifelse(cluster==2, 1,
                                                                       ifelse(cluster==3, 2, NA))))
region3_females_nona <- region3_females %>% filter(!is.na(DistAlongPath) & !is.na(cluster2))


mle.cline <- mle2(cline, theta.init,  
                  control=list(parscale=abs(unlist(theta.init))),
                  data=list(pos=region3_females_nona$DistAlongPath, hap=as.numeric(region3_females_nona$cluster2)))

summary(mle.cline)
confint(mle.cline)
aic_cline <- AIC(mle.cline)

# find fitted cline to plot and to get residuals
pars.m <- coef(mle.cline) # get coefficients for fitted model
# back transform the parameters
#pars.m[2] <- exp(pars.m[2])
#pars.m[3] <- exp(pars.m[3])/(1+exp(pars.m[3]))
#pars.m[4] <- exp(pars.m[4])/(1+exp(pars.m[4]))
# save paramaters
r3_f_pars <- pars.m
r3_f_summary <- summary(mle.cline)
r3_f_confint <- confint(mle.cline)
r3_f_aiccline <- AIC(mle.cline)

pos <- c(0:160) # ordered distance vector (or you can use actual snail distances if you want to find residuals)
dist.m <- pos-pars.m[1]
freq.m <- 1/(1+exp(0-4*(dist.m)/(exp(pars.m[2]))))
mean.m <- (exp(pars.m[3])/(1+exp(pars.m[3])))+((exp(pars.m[4])/(1+exp(pars.m[4])))-(exp(pars.m[3])/(1+exp(pars.m[3]))))*freq.m # fitted phenotype mean
region3_female_fitted <- data.frame(dist.m, mean.m)
ggplot(region3_female_fitted, aes(x=dist.m+89.98, y=mean.m)) + geom_line() + 
  scale_x_continuous(limits=c(0,160)) +
  scale_y_continuous(limits=c(0,1)) + 
  geom_point(data=region3_females_invfreqscol, aes(x=mean_dist, y=1-prop_x), alpha=0.5, size=3, colour="deepskyblue4") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank()) +
  xlab("Distance along path") + ylab("Inversion frequency") +
  geom_vline(xintercept=77.82, colour="red4", linetype=2, alpha=0.5)

r3_f_pars_back <- r3_f_pars[1]
r3_f_pars_back[2] <- exp(r3_f_pars[2])
r3_f_pars_back[3] <- exp(r3_f_pars[3])/(1+exp(r3_f_pars[3]))
r3_f_pars_back[4] <- exp(r3_f_pars[4])/(1+exp(r3_f_pars[4]))

r3_f_confint_back <- data.frame(4,2)
r3_f_confint_back[1,] <- r3_f_confint[1,]
r3_f_confint_back[2,] <- exp(r3_f_confint[2,])
r3_f_confint_back[3,] <- exp(r3_f_confint[3,])/(1+exp(r3_f_confint[3,]))
r3_f_confint_back[4,] <- exp(r3_f_confint[4,])/(1+exp(r3_f_confint[4,]))
###########################################################################


###########################################################################
### region 1 male

ggplot(region1_males_invfreqscol, aes(x=mean_dist, y=1-prop_x)) + geom_point() + scale_y_continuous(limits=c(0,1))

theta.init <- list(c=81,lw=log(10),lp_c=qlogis(0.001), lp_w=qlogis(0.2))
region1_males <- region1_males %>% mutate(cluster2=ifelse(cluster==1, 0,
                                                              ifelse(cluster==2, 1,
                                                                     ifelse(cluster==3, 2, NA))))
region1_males_nona <- region1_males %>% filter(!is.na(DistAlongPath) & !is.na(cluster2))


mle.cline <- mle2(cline, theta.init,  
                  control=list(parscale=abs(unlist(theta.init))),
                  data=list(pos=region1_males_nona$DistAlongPath, hap=as.numeric(region1_males_nona$cluster2)))

summary(mle.cline)
confint(mle.cline)
aic_cline <- AIC(mle.cline)

# find fitted cline to plot and to get residuals
pars.m <- coef(mle.cline) # get coefficients for fitted model
# back transform the parameters
#pars.m[2] <- exp(pars.m[2])
#pars.m[3] <- exp(pars.m[3])/(1+exp(pars.m[3]))
#pars.m[4] <- exp(pars.m[4])/(1+exp(pars.m[4]))
# save paramaters
r1_m_pars <- pars.m
r1_m_summary <- summary(mle.cline)
r1_m_confint <- confint(mle.cline)
r1_m_aiccline <- AIC(mle.cline)

pos <- c(0:160) # ordered distance vector (or you can use actual snail distances if you want to find residuals)
dist.m <- pos-pars.m[1]
freq.m <- 1/(1+exp(0-4*(dist.m)/(exp(pars.m[2]))))
mean.m <- (exp(pars.m[3])/(1+exp(pars.m[3])))+((exp(pars.m[4])/(1+exp(pars.m[4])))-(exp(pars.m[3])/(1+exp(pars.m[3]))))*freq.m # fitted phenotype mean
region1_male_fitted <- data.frame(dist.m, mean.m)
ggplot(region1_male_fitted, aes(x=dist.m+85.12, y=1-mean.m)) + geom_line() + 
  scale_x_continuous(limits=c(0,160)) +
  scale_y_continuous(limits=c(0,1)) + 
  geom_point(data=region1_males_invfreqscol, aes(x=mean_dist, y=prop_x), alpha=0.5, size=3, colour="deepskyblue4") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank()) +
  xlab("Distance along path") + ylab("Inversion frequency") +
  geom_vline(xintercept=77.82, colour="red4", linetype=2, alpha=0.5)

r1_m_pars_back <- r1_m_pars[1]
r1_m_pars_back[2] <- exp(r1_m_pars[2])
r1_m_pars_back[3] <- exp(r1_m_pars[3])/(1+exp(r1_m_pars[3]))
r1_m_pars_back[4] <- exp(r1_m_pars[4])/(1+exp(r1_m_pars[4]))

r1_m_confint_back <- data.frame(4,2)
r1_m_confint_back[1,] <- r1_m_confint[1,]
r1_m_confint_back[2,] <- exp(r1_m_confint[2,])
r1_m_confint_back[3,] <- exp(r1_m_confint[3,])/(1+exp(r1_m_confint[3,]))
r1_m_confint_back[4,] <- exp(r1_m_confint[4,])/(1+exp(r1_m_confint[4,]))
###########################################################################

###########################################################################
### region 2a male

ggplot(region2a_males_invfreqscol, aes(x=mean_dist, y=1-prop_x)) + geom_point() + scale_y_continuous(limits=c(0,1))

theta.init <- list(c=85,lw=log(10),lp_c=qlogis(0.1), lp_w=qlogis(0.3))
region2a_males <- region2a_males %>% mutate(cluster2=ifelse(cluster==1, 0,
                                                          ifelse(cluster==2, 1,
                                                                 ifelse(cluster==3, 2, NA))))
region2a_males_nona <- region2a_males %>% filter(!is.na(DistAlongPath) & !is.na(cluster2))


mle.cline <- mle2(cline, theta.init,  
                  control=list(parscale=abs(unlist(theta.init))),
                  data=list(pos=region2a_males_nona$DistAlongPath, hap=as.numeric(region2a_males_nona$cluster2)))

summary(mle.cline)
confint(mle.cline)
aic_cline <- AIC(mle.cline)

# find fitted cline to plot and to get residuals
pars.m <- coef(mle.cline) # get coefficients for fitted model
# back transform the parameters
#pars.m[2] <- exp(pars.m[2])
#pars.m[3] <- exp(pars.m[3])/(1+exp(pars.m[3]))
#pars.m[4] <- exp(pars.m[4])/(1+exp(pars.m[4]))
# save paramaters
r2a_m_pars <- pars.m
r2a_m_summary <- summary(mle.cline)
r2a_m_confint <- confint(mle.cline)
r2a_m_aiccline <- AIC(mle.cline)

pos <- c(0:160) # ordered distance vector (or you can use actual snail distances if you want to find residuals)
dist.m <- pos-pars.m[1]
freq.m <- 1/(1+exp(0-4*(dist.m)/(exp(pars.m[2]))))
mean.m <- (exp(pars.m[3])/(1+exp(pars.m[3])))+((exp(pars.m[4])/(1+exp(pars.m[4])))-(exp(pars.m[3])/(1+exp(pars.m[3]))))*freq.m # fitted phenotype mean
region2a_male_fitted <- data.frame(dist.m, mean.m)
ggplot(region2a_male_fitted, aes(x=dist.m+83.76, y=1-mean.m)) + geom_line() + 
  scale_x_continuous(limits=c(0,160)) +
  scale_y_continuous(limits=c(0,1)) + 
  geom_point(data=region2a_males_invfreqscol, aes(x=mean_dist, y=prop_x), alpha=0.5, size=3, colour="deepskyblue4") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank()) +
  xlab("Distance along path") + ylab("Inversion frequency") +
  geom_vline(xintercept=77.82, colour="red4", linetype=2, alpha=0.5)

r2a_m_pars_back <- r2a_m_pars[1]
r2a_m_pars_back[2] <- exp(r2a_m_pars[2])
r2a_m_pars_back[3] <- exp(r2a_m_pars[3])/(1+exp(r2a_m_pars[3]))
r2a_m_pars_back[4] <- exp(r2a_m_pars[4])/(1+exp(r2a_m_pars[4]))

# not working as confit finds better fit
r2a_m_confint_back <- data.frame(4,2)
r2a_m_confint_back[1,] <- r2a_m_confint[1,]
r2a_m_confint_back[2,] <- exp(r2a_m_confint[2,])
r2a_m_confint_back[3,] <- exp(r2a_m_confint[3,])/(1+exp(r2a_m_confint[3,]))
r2a_m_confint_back[4,] <- exp(r2a_m_confint[4,])/(1+exp(r2a_m_confint[4,]))
###########################################################################

###########################################################################
### region 2b male

ggplot(region2b_males_invfreqscol, aes(x=mean_dist, y=1-prop_x)) + geom_point() + scale_y_continuous(limits=c(0,1))

theta.init <- list(c=80,lw=log(10),lp_c=qlogis(0.99), lp_w=qlogis(0.9))
region2b_males <- region2b_males %>% mutate(cluster2=ifelse(cluster==1, 0,
                                                            ifelse(cluster==2, 1,
                                                                   ifelse(cluster==3, 2, NA))))
region2b_males_nona <- region2b_males %>% filter(!is.na(DistAlongPath) & !is.na(cluster2))


mle.cline <- mle2(cline, theta.init,  
                  control=list(parscale=abs(unlist(theta.init))),
                  data=list(pos=region2b_males_nona$DistAlongPath, hap=as.numeric(region2b_males_nona$cluster2)))

summary(mle.cline)
confint(mle.cline)
aic_cline <- AIC(mle.cline)

# find fitted cline to plot and to get residuals
pars.m <- coef(mle.cline) # get coefficients for fitted model
# back transform the parameters
#pars.m[2] <- exp(pars.m[2])
#pars.m[3] <- exp(pars.m[3])/(1+exp(pars.m[3]))
#pars.m[4] <- exp(pars.m[4])/(1+exp(pars.m[4]))
# save paramaters
r2b_m_pars <- pars.m
r2b_m_summary <- summary(mle.cline)
r2b_m_confint <- confint(mle.cline)
r2b_m_aiccline <- AIC(mle.cline)

pos <- c(0:160) # ordered distance vector (or you can use actual snail distances if you want to find residuals)
dist.m <- pos-pars.m[1]
freq.m <- 1/(1+exp(0-4*(dist.m)/(exp(pars.m[2]))))
mean.m <- (exp(pars.m[3])/(1+exp(pars.m[3])))+((exp(pars.m[4])/(1+exp(pars.m[4])))-(exp(pars.m[3])/(1+exp(pars.m[3]))))*freq.m # fitted phenotype mean
region2b_male_fitted <- data.frame(dist.m, mean.m)
ggplot(region2b_male_fitted, aes(x=dist.m+73.24, y=1-mean.m)) + geom_line() + 
  scale_x_continuous(limits=c(0,160)) +
  scale_y_continuous(limits=c(0,1)) + 
  geom_point(data=region2b_males_invfreqscol, aes(x=mean_dist, y=prop_x), alpha=0.5, size=3, colour="deepskyblue4") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank()) +
  xlab("Distance along path") + ylab("Inversion frequency") +
  geom_vline(xintercept=77.82, colour="red4", linetype=2, alpha=0.5)

r2b_m_pars_back <- r2b_m_pars[1]
r2b_m_pars_back[2] <- exp(r2b_m_pars[2])
r2b_m_pars_back[3] <- exp(r2b_m_pars[3])/(1+exp(r2b_m_pars[3]))
r2b_m_pars_back[4] <- exp(r2b_m_pars[4])/(1+exp(r2b_m_pars[4]))

r2b_m_confint_back <- data.frame(4,2)
r2b_m_confint_back[1,] <- r2b_m_confint[1,]
r2b_m_confint_back[2,] <- exp(r2b_m_confint[2,])
r2b_m_confint_back[3,] <- exp(r2b_m_confint[3,])/(1+exp(r2b_m_confint[3,]))
r2b_m_confint_back[4,] <- exp(r2b_m_confint[4,])/(1+exp(r2b_m_confint[4,]))
###########################################################################

###########################################################################
### region 3 male

ggplot(region3_males_invfreqscol, aes(x=mean_dist, y=1-prop_x)) + geom_point() + scale_y_continuous(limits=c(0,1))

theta.init <- list(c=85,lw=log(20),lp_c=qlogis(0.3), lp_w=qlogis(0.15))
region3_males <- region3_males %>% mutate(cluster2=ifelse(cluster==1, 0,
                                                            ifelse(cluster==2, 1,
                                                                   ifelse(cluster==3, 2, NA))))
region3_males_nona <- region3_males %>% filter(!is.na(DistAlongPath) & !is.na(cluster2))


mle.cline <- mle2(cline, theta.init,  
                  control=list(parscale=abs(unlist(theta.init))),
                  data=list(pos=region3_males_nona$DistAlongPath, hap=as.numeric(region3_males_nona$cluster2)))

summary(mle.cline)
confint(mle.cline)
aic_cline <- AIC(mle.cline)

# find fitted cline to plot and to get residuals
pars.m <- coef(mle.cline) # get coefficients for fitted model
# back transform the parameters
#pars.m[2] <- exp(pars.m[2])
#pars.m[3] <- exp(pars.m[3])/(1+exp(pars.m[3]))
#pars.m[4] <- exp(pars.m[4])/(1+exp(pars.m[4]))
# save paramaters
r3_m_pars <- pars.m
r3_m_summary <- summary(mle.cline)
r3_m_confint <- confint(mle.cline)
r3_m_aiccline <- AIC(mle.cline)

pos <- c(0:160) # ordered distance vector (or you can use actual snail distances if you want to find residuals)
dist.m <- pos-pars.m[1]
freq.m <- 1/(1+exp(0-4*(dist.m)/(exp(pars.m[2]))))
mean.m <- (exp(pars.m[3])/(1+exp(pars.m[3])))+((exp(pars.m[4])/(1+exp(pars.m[4])))-(exp(pars.m[3])/(1+exp(pars.m[3]))))*freq.m # fitted phenotype mean
region3_male_fitted <- data.frame(dist.m, mean.m)
ggplot(region3_male_fitted, aes(x=dist.m+93.5, y=mean.m)) + geom_line() + 
  scale_x_continuous(limits=c(0,160)) +
  scale_y_continuous(limits=c(0,1)) + 
  geom_point(data=region3_males_invfreqscol, aes(x=mean_dist, y=1-prop_x), alpha=0.5, size=3, colour="deepskyblue4") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank()) +
  xlab("Distance along path") + ylab("Inversion frequency") +
  geom_vline(xintercept=77.82, colour="red4", linetype=2, alpha=0.5)

r3_m_pars_back <- r3_m_pars[1]
r3_m_pars_back[2] <- exp(r3_m_pars[2])
r3_m_pars_back[3] <- exp(r3_m_pars[3])/(1+exp(r3_m_pars[3]))
r3_m_pars_back[4] <- exp(r3_m_pars[4])/(1+exp(r3_m_pars[4]))

r3_m_confint_back <- data.frame(4,2)
r3_m_confint_back[1,] <- r3_m_confint[1,]
r3_m_confint_back[2,] <- exp(r3_m_confint[2,])
r3_m_confint_back[3,] <- exp(r3_m_confint[3,])/(1+exp(r3_m_confint[3,]))
r3_m_confint_back[4,] <- exp(r3_m_confint[4,])/(1+exp(r3_m_confint[4,]))
###########################################################################


###########################################################################
### plots with male and female clines for each region together

## region 1
ggplot(region1_male_fitted, aes(x=dist.m+85.12, y=1-mean.m)) + geom_line(colour="dodgerblue4") + 
  geom_line(data=region1_female_fitted, aes(x=dist.m+87.59, y=1-mean.m), colour="brown4") +
  scale_x_continuous(limits=c(0,160)) +
  scale_y_continuous(limits=c(0,1)) + 
  geom_point(data=region1_females_invfreqscol, aes(x=mean_dist, y=prop_x), alpha=0.5, size=3, colour="#F8766D") +
  geom_point(data=region1_males_invfreqscol, aes(x=mean_dist, y=prop_x), alpha=0.5, size=3, colour="#619CFF") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank()) +
  xlab("Distance along path") + ylab("Inversion frequency") +
  geom_vline(xintercept=67, colour="orange2", linetype=2, alpha=0.8) +
  geom_vline(xintercept=85, colour="orange2", linetype=2, alpha=0.8) +
  geom_vline(xintercept=91.8, colour="green4", linetype=2, alpha=0.8) +
  theme(aspect.ratio = 1)

## region 2a
ggplot(region2a_male_fitted, aes(x=dist.m+83.76, y=1-mean.m)) + geom_line(colour="dodgerblue4") + 
  geom_line(data=region2a_female_fitted, aes(x=dist.m+86.32, y=1-mean.m), colour="brown4") +
  scale_x_continuous(limits=c(0,160)) +
  scale_y_continuous(limits=c(0,1)) + 
  geom_point(data=region2a_females_invfreqscol, aes(x=mean_dist, y=prop_x), alpha=0.5, size=3, colour="#F8766D") +
  geom_point(data=region2a_males_invfreqscol, aes(x=mean_dist, y=prop_x), alpha=0.5, size=3, colour="#619CFF") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank()) +
  xlab("Distance along path") + ylab("Inversion frequency") +
  geom_vline(xintercept=77.82, colour="red4", linetype=2, alpha=0.5)

## region 2b
ggplot(region2b_male_fitted, aes(x=dist.m+73.24, y=1-mean.m)) + geom_line(colour="dodgerblue4") + 
  geom_line(data=region2b_female_fitted, aes(x=dist.m+90.23, y=1-mean.m), colour="brown4") +
  scale_x_continuous(limits=c(0,160)) +
  scale_y_continuous(limits=c(0,1)) + 
  geom_point(data=region2b_females_invfreqscol, aes(x=mean_dist, y=prop_x), alpha=0.5, size=3, colour="#F8766D") +
  geom_point(data=region2b_males_invfreqscol, aes(x=mean_dist, y=prop_x), alpha=0.5, size=3, colour="#619CFF") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank()) +
  xlab("Distance along path") + ylab("Inversion frequency") +
  geom_vline(xintercept=77.82, colour="red4", linetype=2, alpha=0.5)

## region 3
ggplot(region3_male_fitted, aes(x=dist.m+93.46, y=mean.m)) + geom_line(colour="dodgerblue4") + 
  geom_line(data=region3_female_fitted, aes(x=dist.m+89.98, y=mean.m), colour="brown4") +
  scale_x_continuous(limits=c(0,160)) +
  scale_y_continuous(limits=c(0,1)) + 
  geom_point(data=region3_females_invfreqscol, aes(x=mean_dist, y=1-prop_x), alpha=0.5, size=3, colour="#F8766D") +
  geom_point(data=region3_males_invfreqscol, aes(x=mean_dist, y=1-prop_x), alpha=0.5, size=3, colour="#619CFF") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank()) +
  xlab("Distance along path") + ylab("Inversion frequency") +
  geom_vline(xintercept=77.82, colour="red4", linetype=2, alpha=0.5)


