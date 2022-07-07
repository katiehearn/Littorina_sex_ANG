#####################
#####################
##
## follow-on script to inversion_cline_fits.R
## 
## control clines- no cline model, comparing male and female clines (separate vs combined, constrained frequencies and cline centres/widths)
##
#####################
#####################



##### null cline - one parameter- logit(f) - where f is overall allele frequency across the cline

# example format
cline_control <- function(hap, pos, f){
  z_x <- exp(f)/(1+exp(f))
  minusll <- -sum(dbinom(hap, 2, z_x,log=T))
  return(minusll)
}

theta.init <- list(f=qlogis(freq_est))

mle.cline <- mle2(cline, theta.init,  
                  control=list(parscale=abs(unlist(theta.init))),
                  data=list(pos=dist,hap=inv))

summary(mle.cline)
confint(mle.cline)
AIC(mle.cline)


####################################################################################
## region1 females

counts <- region1_females %>% group_by(cluster) %>% count()
freq <- ((counts[1,2])+(counts[2,2]/2))/(sum(counts$n))

theta.init <- list(f=qlogis(0.99))

mle.cline <- mle2(cline_control, theta.init,  
                  control=list(parscale=abs(unlist(theta.init))),
                  data=list(pos=region1_females_nona$DistAlongPath, hap=as.numeric(region1_females_nona$cluster2)))

summary(mle.cline)
confint(mle.cline)
AIC(mle.cline)

# save the mean freq from mle and other things
r1_f_control_pars <- coef(mle.cline)
r1_f_control_summary <- summary(mle.cline)
r1_f_control_confint <- confint(mle.cline)
r1_f_control_aiccline <- AIC(mle.cline)

# back transform parameter and confidence intervals
r1_f_control_pars_back <- exp(r1_f_control_pars)/(1+exp(r1_f_control_pars))
r1_f_control_confint_back <- exp(r1_f_control_confint)/(1+exp(r1_f_control_confint))

# add control cline to plot (copy plot script over from other file and just add the line)
ggplot(region1_female_fitted, aes(x=dist.m+87.5, y=1-mean.m)) + geom_line() + 
  scale_x_continuous(limits=c(0,160)) +
  geom_line(aes(y=1-r1_f_control_pars_back), alpha=0.7, colour="orange") +
  scale_y_continuous(limits=c(0,1)) + 
  geom_point(data=region1_females_invfreqscol, aes(x=mean_dist, y=prop_x), alpha=0.5, size=3, colour="deepskyblue4") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank()) +
  xlab("Distance along path") + ylab("Inversion frequency") +
  geom_vline(xintercept=77.82, colour="red4", linetype="longdash", alpha=0.7)

####################################################################################

####################################################################################
## region 2a females

counts <- region2a_females %>% group_by(cluster) %>% count()
freq <- ((counts[1,2])+(counts[2,2]/2))/(sum(counts$n))
freq

theta.init <- list(f=qlogis(0.49))

mle.cline <- mle2(cline_control, theta.init,  
                  control=list(parscale=abs(unlist(theta.init))),
                  data=list(pos=region2a_females_nona$DistAlongPath, hap=as.numeric(region2a_females_nona$cluster2)))

summary(mle.cline)
confint(mle.cline)
AIC(mle.cline)

# save the mean freq from mle and other things
r2a_f_control_pars <- coef(mle.cline)
r2a_f_control_summary <- summary(mle.cline)
r2a_f_control_confint <- confint(mle.cline)
r2a_f_control_aiccline <- AIC(mle.cline)

# back transform parameter and confidence intervals
r2a_f_control_pars_back <- exp(r2a_f_control_pars)/(1+exp(r2a_f_control_pars))
r2a_f_control_confint_back <- exp(r2a_f_control_confint)/(1+exp(r2a_f_control_confint))

# add control cline to plot (copy plot script over from other file and just add the line)
ggplot(region2a_female_fitted, aes(x=dist.m+86.32, y=1-mean.m)) + geom_line() + 
  scale_x_continuous(limits=c(0,160)) +
  scale_y_continuous(limits=c(0,1)) + 
  geom_line(aes(y=1-r2a_f_control_pars_back), alpha=0.7, colour="orange") +
  geom_point(data=region2a_females_invfreqscol, aes(x=mean_dist, y=prop_x), alpha=0.5, size=3, colour="deepskyblue4") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank()) +
  xlab("Distance along path") + ylab("Inversion frequency") +
  geom_vline(xintercept=77.82, colour="red4", linetype="longdash", alpha=0.7)

####################################################################################

####################################################################################
## region 2b females

counts <- region2b_females %>% group_by(cluster) %>% count()
freq <- ((counts[1,2])+(counts[2,2]/2))/(sum(counts$n))
freq

theta.init <- list(f=qlogis(0.35))

mle.cline <- mle2(cline_control, theta.init,  
                  control=list(parscale=abs(unlist(theta.init))),
                  data=list(pos=region2b_females_nona$DistAlongPath, hap=as.numeric(region2b_females_nona$cluster2)))

summary(mle.cline)
confint(mle.cline)
AIC(mle.cline)

# save the mean freq from mle and other things
r2b_f_control_pars <- coef(mle.cline)
r2b_f_control_summary <- summary(mle.cline)
r2b_f_control_confint <- confint(mle.cline)
r2b_f_control_aiccline <- AIC(mle.cline)

# back transform parameter and confidence intervals
r2b_f_control_pars_back <- exp(r2b_f_control_pars)/(1+exp(r2b_f_control_pars))
r2b_f_control_confint_back <- exp(r2b_f_control_confint)/(1+exp(r2b_f_control_confint))

# add control cline to plot (copy plot script over from other file and just add the line)
ggplot(region2b_female_fitted, aes(x=dist.m+90.23, y=1-mean.m)) + geom_line() + 
  scale_x_continuous(limits=c(0,160)) +
  scale_y_continuous(limits=c(0,1)) + 
  geom_line(aes(y=1-r2b_f_control_pars_back), alpha=0.7, colour="orange") +
  geom_point(data=region2b_females_invfreqscol, aes(x=mean_dist, y=prop_x), alpha=0.5, size=3, colour="deepskyblue4") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank()) +
  xlab("Distance along path") + ylab("Inversion frequency") +
  geom_vline(xintercept=77.82, colour="red4", linetype="longdash", alpha=0.7)

####################################################################################

####################################################################################
## region 3 females

counts <- region3_females %>% group_by(cluster) %>% count()
freq <- ((counts[1,2])+(counts[2,2]/2))/(sum(counts$n))
freq

theta.init <- list(f=qlogis(0.55))

mle.cline <- mle2(cline_control, theta.init,  
                  control=list(parscale=abs(unlist(theta.init))),
                  data=list(pos=region3_females_nona$DistAlongPath, hap=as.numeric(region3_females_nona$cluster2)))

summary(mle.cline)
confint(mle.cline)
AIC(mle.cline)

# save the mean freq from mle and other things
r3_f_control_pars <- coef(mle.cline)
r3_f_control_summary <- summary(mle.cline)
r3_f_control_confint <- confint(mle.cline)
r3_f_control_aiccline <- AIC(mle.cline)

# back transform parameter and confidence intervals
r3_f_control_pars_back <- exp(r3_f_control_pars)/(1+exp(r3_f_control_pars))
r3_f_control_confint_back <- exp(r3_f_control_confint)/(1+exp(r3_f_control_confint))

# add control cline to plot (copy plot script over from other file and just add the line)
ggplot(region3_female_fitted, aes(x=dist.m+89.98, y=mean.m)) + geom_line() + 
  scale_x_continuous(limits=c(0,160)) +
  scale_y_continuous(limits=c(0,1)) + 
  geom_line(aes(y=r3_f_control_pars_back), alpha=0.7, colour="orange") +
  geom_point(data=region3_females_invfreqscol, aes(x=mean_dist, y=1-prop_x), alpha=0.5, size=3, colour="deepskyblue4") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank()) +
  xlab("Distance along path") + ylab("Inversion frequency") +
  geom_vline(xintercept=77.82, colour="red4", linetype="longdash", alpha=0.7)

####################################################################################

####################################################################################
## region1 males

counts <- region1_males %>% group_by(cluster) %>% count()
freq <- ((counts[1,2])+(counts[2,2]/2))/(sum(counts$n))
freq

theta.init <- list(f=qlogis(0.89))

mle.cline <- mle2(cline_control, theta.init,  
                  control=list(parscale=abs(unlist(theta.init))),
                  data=list(pos=region1_males_nona$DistAlongPath, hap=as.numeric(region1_males_nona$cluster2)))

summary(mle.cline)
confint(mle.cline)
AIC(mle.cline)

# save the mean freq from mle and other things
r1_m_control_pars <- coef(mle.cline)
r1_m_control_summary <- summary(mle.cline)
r1_m_control_confint <- confint(mle.cline)
r1_m_control_aiccline <- AIC(mle.cline)

# back transform parameter and confidence intervals
r1_m_control_pars_back <- exp(r1_m_control_pars)/(1+exp(r1_m_control_pars))
r1_m_control_confint_back <- exp(r1_m_control_confint)/(1+exp(r1_m_control_confint))

# add control cline to plot (copy plot script over from other file and just add the line)
ggplot(region1_male_fitted, aes(x=dist.m+85.12, y=1-mean.m)) + geom_line() + 
  scale_x_continuous(limits=c(0,160)) +
  geom_line(aes(y=1-r1_m_control_pars_back), alpha=0.7, colour="orange") +
  scale_y_continuous(limits=c(0,1)) + 
  geom_point(data=region1_males_invfreqscol, aes(x=mean_dist, y=prop_x), alpha=0.5, size=3, colour="deepskyblue4") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank()) +
  xlab("Distance along path") + ylab("Inversion frequency") +
  geom_vline(xintercept=77.82, colour="red4", linetype="longdash", alpha=0.7)

####################################################################################

####################################################################################
## region 2a males

counts <- region2a_males %>% group_by(cluster) %>% count()
freq <- ((counts[1,2])+(counts[2,2]/2))/(sum(counts$n))
freq

theta.init <- list(f=qlogis(0.85))

mle.cline <- mle2(cline_control, theta.init,  
                  control=list(parscale=abs(unlist(theta.init))),
                  data=list(pos=region2a_males_nona$DistAlongPath, hap=as.numeric(region2a_males_nona$cluster2)))

summary(mle.cline)
confint(mle.cline)
AIC(mle.cline)

# save the mean freq from mle and other things
r2a_m_control_pars <- coef(mle.cline)
r2a_m_control_summary <- summary(mle.cline)
r2a_m_control_confint <- confint(mle.cline)
r2a_m_control_aiccline <- AIC(mle.cline)

# back transform parameter and confidence intervals
r2a_m_control_pars_back <- exp(r2a_m_control_pars)/(1+exp(r2a_m_control_pars))
r2a_m_control_confint_back <- exp(r2a_m_control_confint)/(1+exp(r2a_m_control_confint))

# add control cline to plot (copy plot script over from other file and just add the line)
ggplot(region2a_male_fitted, aes(x=dist.m+83.76, y=1-mean.m)) + geom_line() + 
  scale_x_continuous(limits=c(0,160)) +
  geom_line(aes(y=1-r2a_m_control_pars_back), alpha=0.7, colour="orange") +
  scale_y_continuous(limits=c(0,1)) + 
  geom_point(data=region2a_males_invfreqscol, aes(x=mean_dist, y=prop_x), alpha=0.5, size=3, colour="deepskyblue4") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank()) +
  xlab("Distance along path") + ylab("Inversion frequency") +
  geom_vline(xintercept=77.82, colour="red4", linetype="longdash", alpha=0.7)

####################################################################################

####################################################################################
## region 2b males

counts <- region2b_males %>% group_by(cluster) %>% count()
freq <- ((counts[1,2])+(counts[2,2]/2))/(sum(counts$n))
freq

theta.init <- list(f=qlogis(0.04))

mle.cline <- mle2(cline_control, theta.init,  
                  control=list(parscale=abs(unlist(theta.init))),
                  data=list(pos=region2b_males_nona$DistAlongPath, hap=as.numeric(region2b_males_nona$cluster2)))

summary(mle.cline)
confint(mle.cline)
AIC(mle.cline)

# save the mean freq from mle and other things
r2b_m_control_pars <- coef(mle.cline)
r2b_m_control_summary <- summary(mle.cline)
r2b_m_control_confint <- confint(mle.cline)
r2b_m_control_aiccline <- AIC(mle.cline)

# back transform parameter and confidence intervals
r2b_m_control_pars_back <- exp(r2b_m_control_pars)/(1+exp(r2b_m_control_pars))
r2b_m_control_confint_back <- exp(r2b_m_control_confint)/(1+exp(r2b_m_control_confint))

# add control cline to plot (copy plot script over from other file and just add the line)
ggplot(region2b_male_fitted, aes(x=dist.m+73.24, y=1-mean.m)) + geom_line() + 
  scale_x_continuous(limits=c(0,160)) +
  geom_line(aes(y=1-r2b_m_control_pars_back), alpha=0.7, colour="orange") +
  scale_y_continuous(limits=c(0,1)) + 
  geom_point(data=region2b_males_invfreqscol, aes(x=mean_dist, y=prop_x), alpha=0.5, size=3, colour="deepskyblue4") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank()) +
  xlab("Distance along path") + ylab("Inversion frequency") +
  geom_vline(xintercept=77.82, colour="red4", linetype="longdash", alpha=0.7)

####################################################################################

####################################################################################
## region 3 males

counts <- region3_males %>% group_by(cluster) %>% count()
freq <- ((counts[1,2])+(counts[2,2]/2))/(sum(counts$n))
freq

theta.init <- list(f=qlogis(0.75))

mle.cline <- mle2(cline_control, theta.init,  
                  control=list(parscale=abs(unlist(theta.init))),
                  data=list(pos=region3_males_nona$DistAlongPath, hap=as.numeric(region3_males_nona$cluster2)))

summary(mle.cline)
confint(mle.cline)
AIC(mle.cline)

# save the mean freq from mle and other things
r3_m_control_pars <- coef(mle.cline)
r3_m_control_summary <- summary(mle.cline)
r3_m_control_confint <- confint(mle.cline)
r3_m_control_aiccline <- AIC(mle.cline)

# back transform parameter and confidence intervals
r3_m_control_pars_back <- exp(r3_m_control_pars)/(1+exp(r3_m_control_pars))
r3_m_control_confint_back <- exp(r3_m_control_confint)/(1+exp(r3_m_control_confint))

# add control cline to plot (copy plot script over from other file and just add the line)
ggplot(region3_male_fitted, aes(x=dist.m+93.5, y=mean.m)) + geom_line() + 
  scale_x_continuous(limits=c(0,160)) +
  geom_line(aes(y=r3_m_control_pars_back), alpha=0.7, colour="orange") +
  scale_y_continuous(limits=c(0,1)) + 
  geom_point(data=region3_males_invfreqscol, aes(x=mean_dist, y=1-prop_x), alpha=0.5, size=3, colour="deepskyblue4") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank()) +
  xlab("Distance along path") + ylab("Inversion frequency") +
  geom_vline(xintercept=77.82, colour="red4", linetype="longdash", alpha=0.7)

####################################################################################


########
########
######### cline fitted to combined male and female data - to test whether cline parameters differ between males and females
########
########


### combine the two datasets then run as normal
### cline function should already be ready from previous script


####################################################################################

## region 1 

region1_combined <- region1_males %>% select(-Axis3, -Axis4) %>% bind_rows(., region1_females)

ggplot(region1_combined, aes(x=DistAlongPath, y=cluster, colour=sex_dissection)) + 
  geom_jitter(alpha=0.5, size=4, width=0, height=0.15) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())

ggplot(region1_females_invfreqscol, aes(x=mean_dist, y=prop_x)) + geom_point(alpha=0.5, colour="#F8766D") +
  geom_point(data=region1_males_invfreqscol, aes(x=mean_dist, y=prop_x), alpha=0.5, colour="#619CFF") +
  scale_y_continuous(limits=c(0,1)) + scale_x_continuous(limits=c(0,130))

theta.init <- list(c=80,lw=log(10),lp_c=qlogis(0.01), lp_w=qlogis(0.25))

region1_combined_nona <- region1_combined %>% filter(!is.na(DistAlongPath) & !is.na(cluster2))

mle.cline <- mle2(cline, theta.init,  
                  control=list(parscale=abs(unlist(theta.init))),
                  data=list(pos=region1_combined_nona$DistAlongPath, hap=as.numeric(region1_combined_nona$cluster2)))

summary(mle.cline)
confint(mle.cline)
AIC(mle.cline)

pars.m <- coef(mle.cline)

pos <- c(0:160) 
dist.m <- pos-pars.m[1]
freq.m <- 1/(1+exp(0-4*(dist.m)/(exp(pars.m[2]))))
mean.m <- (exp(pars.m[3])/(1+exp(pars.m[3])))+((exp(pars.m[4])/(1+exp(pars.m[4])))-(exp(pars.m[3])/(1+exp(pars.m[3]))))*freq.m

region1_combined_fitted <- data.frame(dist.m, mean.m)

ggplot(region1_combined_fitted, aes(x=dist.m+85.27, y=1-mean.m)) + geom_line() + 
  scale_x_continuous(limits=c(0,160)) +
  scale_y_continuous(limits=c(0,1)) + 
  geom_point(data=region1_females_invfreqscol, aes(x=mean_dist, y=prop_x), alpha=0.5, size=3, colour="#F8766D") +
  geom_point(data=region1_males_invfreqscol, aes(x=mean_dist, y=prop_x), alpha=0.5, size=3, colour="#619CFF") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank()) +
  xlab("Distance along path") + ylab("Inversion frequency") +
  geom_vline(xintercept=77.82, colour="red4", linetype=2, alpha=0.5)

r1_comb_pars <- pars.m
r1_comb_summary <- summary(mle.cline)
r1_comb_confint <- confint(mle.cline)
r1_comb_aiccline <- AIC(mle.cline)

r1_comb_pars_back <- r1_comb_pars[1]
r1_comb_pars_back[2] <- exp(r1_comb_pars[2])
r1_comb_pars_back[3] <- exp(r1_comb_pars[3])/(1+exp(r1_comb_pars[3]))
r1_comb_pars_back[4] <- exp(r1_comb_pars[4])/(1+exp(r1_comb_pars[4]))

#r1_comb_confint_back <- data.frame(4,2)
#r1_comb_confint_back[1,] <- r1_comb_confint[1,]
#r1_comb_confint_back[2,] <- exp(r1_comb_confint[2,])
#r1_comb_confint_back[3,] <- exp(r1_comb_confint[3,])/(1+exp(r1_comb_confint[3,]))
#r1_comb_confint_back[4,] <- exp(r1_comb_confint[4,])/(1+exp(r1_comb_confint[4,]))

####################################################################################

####################################################################################

## region 2a 

region2a_combined <- region2a_males %>% select(-Axis3) %>% bind_rows(., region2a_females)

ggplot(region2a_combined, aes(x=DistAlongPath, y=cluster, colour=sex_dissection)) + 
  geom_jitter(alpha=0.5, size=4, width=0, height=0.15) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())

ggplot(region2a_females_invfreqscol, aes(x=mean_dist, y=prop_x)) + geom_point(alpha=0.5, colour="#F8766D") +
  geom_point(data=region2a_males_invfreqscol, aes(x=mean_dist, y=prop_x), alpha=0.5, colour="#619CFF") +
  scale_y_continuous(limits=c(0,1)) + scale_x_continuous(limits=c(0,130))

theta.init <- list(c=85,lw=log(10),lp_c=qlogis(0.25), lp_w=qlogis(0.4))

region2a_combined_nona <- region2a_combined %>% filter(!is.na(DistAlongPath) & !is.na(cluster2))

mle.cline <- mle2(cline, theta.init,  
                  control=list(parscale=abs(unlist(theta.init))),
                  data=list(pos=region2a_combined_nona$DistAlongPath, hap=as.numeric(region2a_combined_nona$cluster2)))

summary(mle.cline)
confint(mle.cline)
AIC(mle.cline)

pars.m <- coef(mle.cline)

pos <- c(0:160) 
dist.m <- pos-pars.m[1]
freq.m <- 1/(1+exp(0-4*(dist.m)/(exp(pars.m[2]))))
mean.m <- (exp(pars.m[3])/(1+exp(pars.m[3])))+((exp(pars.m[4])/(1+exp(pars.m[4])))-(exp(pars.m[3])/(1+exp(pars.m[3]))))*freq.m

region2a_combined_fitted <- data.frame(dist.m, mean.m)

ggplot(region2a_combined_fitted, aes(x=dist.m+83.59, y=1-mean.m)) + geom_line() + 
  scale_x_continuous(limits=c(0,160)) +
  scale_y_continuous(limits=c(0,1)) + 
  geom_point(data=region2a_females_invfreqscol, aes(x=mean_dist, y=prop_x), alpha=0.5, size=3, colour="#F8766D") +
  geom_point(data=region2a_males_invfreqscol, aes(x=mean_dist, y=prop_x), alpha=0.5, size=3, colour="#619CFF") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank()) +
  xlab("Distance along path") + ylab("Inversion frequency") +
  geom_vline(xintercept=77.82, colour="red4", linetype=2, alpha=0.5)

r2a_comb_pars <- pars.m
r2a_comb_summary <- summary(mle.cline)
r2a_comb_confint <- confint(mle.cline)
r2a_comb_aiccline <- AIC(mle.cline)

r2a_comb_pars_back <- r2a_comb_pars[1]
r2a_comb_pars_back[2] <- exp(r2a_comb_pars[2])
r2a_comb_pars_back[3] <- exp(r2a_comb_pars[3])/(1+exp(r2a_comb_pars[3]))
r2a_comb_pars_back[4] <- exp(r2a_comb_pars[4])/(1+exp(r2a_comb_pars[4]))

#r2a_comb_confint_back <- data.frame(4,2)
#r2a_comb_confint_back[1,] <- r2a_comb_confint[1,]
#r2a_comb_confint_back[2,] <- exp(r2a_comb_confint[2,])
#r2a_comb_confint_back[3,] <- exp(r2a_comb_confint[3,])/(1+exp(r2a_comb_confint[3,]))
#r2a_comb_confint_back[4,] <- exp(r2a_comb_confint[4,])/(1+exp(r2a_comb_confint[4,]))

####################################################################################

####################################################################################

## region 2b 

region2b_combined <- region2b_males %>% bind_rows(., region2b_females)

ggplot(region2b_combined, aes(x=DistAlongPath, y=cluster, colour=sex_dissection)) + 
  geom_jitter(alpha=0.5, size=4, width=0, height=0.15) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())

ggplot(region2b_females_invfreqscol, aes(x=mean_dist, y=prop_x)) + geom_point(alpha=0.5, colour="#F8766D") +
  geom_point(data=region2b_males_invfreqscol, aes(x=mean_dist, y=prop_x), alpha=0.5, colour="#619CFF") +
  scale_y_continuous(limits=c(0,1)) + scale_x_continuous(limits=c(0,130))

theta.init <- list(c=85,lw=log(10),lp_c=qlogis(0.7), lp_w=qlogis(0.9))

region2b_combined_nona <- region2b_combined %>% filter(!is.na(DistAlongPath) & !is.na(cluster2))

mle.cline <- mle2(cline, theta.init,  
                  control=list(parscale=abs(unlist(theta.init))),
                  data=list(pos=region2b_combined_nona$DistAlongPath, hap=as.numeric(region2b_combined_nona$cluster2)))

summary(mle.cline)
confint(mle.cline)
AIC(mle.cline)

pars.m <- coef(mle.cline)

pos <- c(0:160) 
dist.m <- pos-pars.m[1]
freq.m <- 1/(1+exp(0-4*(dist.m)/(exp(pars.m[2]))))
mean.m <- (exp(pars.m[3])/(1+exp(pars.m[3])))+((exp(pars.m[4])/(1+exp(pars.m[4])))-(exp(pars.m[3])/(1+exp(pars.m[3]))))*freq.m

region2b_combined_fitted <- data.frame(dist.m, mean.m)

ggplot(region2b_combined_fitted, aes(x=dist.m+93.75, y=1-mean.m)) + geom_line() + 
  scale_x_continuous(limits=c(0,160)) +
  scale_y_continuous(limits=c(0,1)) + 
  geom_point(data=region2b_females_invfreqscol, aes(x=mean_dist, y=prop_x), alpha=0.5, size=3, colour="#F8766D") +
  geom_point(data=region2b_males_invfreqscol, aes(x=mean_dist, y=prop_x), alpha=0.5, size=3, colour="#619CFF") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank()) +
  xlab("Distance along path") + ylab("Inversion frequency") +
  geom_vline(xintercept=77.82, colour="red4", linetype=2, alpha=0.5)

r2b_comb_pars <- pars.m
r2b_comb_summary <- summary(mle.cline)
r2b_comb_confint <- confint(mle.cline)
r2b_comb_aiccline <- AIC(mle.cline)

r2b_comb_pars_back <- r2b_comb_pars[1]
r2b_comb_pars_back[2] <- exp(r2b_comb_pars[2])
r2b_comb_pars_back[3] <- exp(r2b_comb_pars[3])/(1+exp(r2b_comb_pars[3]))
r2b_comb_pars_back[4] <- exp(r2b_comb_pars[4])/(1+exp(r2b_comb_pars[4]))

#r2b_comb_confint_back <- data.frame(4,2)
#r2b_comb_confint_back[1,] <- r2b_comb_confint[1,]
#r2b_comb_confint_back[2,] <- exp(r2b_comb_confint[2,])
#r2b_comb_confint_back[3,] <- exp(r2b_comb_confint[3,])/(1+exp(r2b_comb_confint[3,]))
#r2b_comb_confint_back[4,] <- exp(r2b_comb_confint[4,])/(1+exp(r2b_comb_confint[4,]))

####################################################################################

####################################################################################

## region 3 

region3_combined <- region3_females %>% select(-Axis3) %>% bind_rows(., region3_males)

ggplot(region3_combined, aes(x=DistAlongPath, y=cluster, colour=sex_dissection)) + 
  geom_jitter(alpha=0.5, size=4, width=0, height=0.15) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())

ggplot(region3_females_invfreqscol, aes(x=mean_dist, y=prop_x)) + geom_point(alpha=0.5, colour="#F8766D") +
  geom_point(data=region3_males_invfreqscol, aes(x=mean_dist, y=prop_x), alpha=0.5, colour="#619CFF") +
  scale_y_continuous(limits=c(0,1)) + scale_x_continuous(limits=c(0,130))

theta.init <- list(c=80,lw=log(10),lp_c=qlogis(0.45), lp_w=qlogis(0.15))

region3_combined_nona <- region3_combined %>% filter(!is.na(DistAlongPath) & !is.na(cluster2))

mle.cline <- mle2(cline, theta.init,  
                  control=list(parscale=abs(unlist(theta.init))),
                  data=list(pos=region3_combined_nona$DistAlongPath, hap=as.numeric(region3_combined_nona$cluster2)))

summary(mle.cline)
confint(mle.cline)
AIC(mle.cline)

pars.m <- coef(mle.cline)

pos <- c(0:160) 
dist.m <- pos-pars.m[1]
freq.m <- 1/(1+exp(0-4*(dist.m)/(exp(pars.m[2]))))
mean.m <- (exp(pars.m[3])/(1+exp(pars.m[3])))+((exp(pars.m[4])/(1+exp(pars.m[4])))-(exp(pars.m[3])/(1+exp(pars.m[3]))))*freq.m

region3_combined_fitted <- data.frame(dist.m, mean.m)

ggplot(region3_combined_fitted, aes(x=dist.m+90.6, y=mean.m)) + geom_line() + 
  scale_x_continuous(limits=c(0,160)) +
  scale_y_continuous(limits=c(0,1)) + 
  geom_point(data=region3_females_invfreqscol, aes(x=mean_dist, y=1-prop_x), alpha=0.5, size=3, colour="#F8766D") +
  geom_point(data=region3_males_invfreqscol, aes(x=mean_dist, y=1-prop_x), alpha=0.5, size=3, colour="#619CFF") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank()) +
  xlab("Distance along path") + ylab("Inversion frequency") +
  geom_vline(xintercept=77.82, colour="red4", linetype=2, alpha=0.5)

r3_comb_pars <- pars.m
r3_comb_summary <- summary(mle.cline)
r3_comb_confint <- confint(mle.cline)
r3_comb_aiccline <- AIC(mle.cline)

r3_comb_pars_back <- r3_comb_pars[1]
r3_comb_pars_back[2] <- exp(r3_comb_pars[2])
r3_comb_pars_back[3] <- exp(r3_comb_pars[3])/(1+exp(r3_comb_pars[3]))
r3_comb_pars_back[4] <- exp(r3_comb_pars[4])/(1+exp(r3_comb_pars[4]))

#r3_comb_confint_back <- data.frame(4,2)
#r3_comb_confint_back[1,] <- r3_comb_confint[1,]
#r3_comb_confint_back[2,] <- exp(r3_comb_confint[2,])
#r3_comb_confint_back[3,] <- exp(r3_comb_confint[3,])/(1+exp(r3_comb_confint[3,]))
#r1_comb_confint_back[4,] <- exp(r3_comb_confint[4,])/(1+exp(r3_comb_confint[4,]))

####################################################################################




########
########
######### constrained cline - do end frequencies differ (m vs f) when centre and width are constrained to be the same
########
########


####################################################################################

# new cline function with separate male and female frequency parameters

cline_mf <- function(hap,pos,sex,c,lw,lp_cf,lp_wf,lp_cm,lp_wm){  # now separate male and female frequency parameters, need to pass a variable with sex information
  d <- pos-c
  
  p_x <- 1/(1+exp(0-4*(d/exp(lw))))  
  
  crab_f <- exp(lp_cf)/(1+exp(lp_cf)) # separate male and female again
  wave_f <- exp(lp_wf)/(1+exp(lp_wf))
  crab_m <- exp(lp_cm)/(1+exp(lp_cm)) 
  wave_m <- exp(lp_wm)/(1+exp(lp_wm))
  
  z_x_f <- crab_f + (wave_f-crab_f)*p_x    # z_x_f is expected frequency if female 
  z_x_m <- crab_m + (wave_m-crab_m)*p_x    # z_x_m is expected frequency if male
  
  minusll <- -sum(dbinom(hap[sex == "f"],2,z_x_f[sex == "f"],log=T)) - sum(dbinom(hap[sex == "m"],2,z_x_m[sex == "m"],log=T)) # -sums of probabilities for females and males separately
  
  return(minusll)
}

####################################################################################

## region 1

region1_combined_nona <- region1_combined_nona %>% mutate(sex=ifelse(sex_dissection=="male", "m",
                                                            ifelse(sex_dissection=="female", "f", NA)))

ggplot(region1_combined, aes(x=DistAlongPath, y=cluster, colour=sex_dissection)) + 
  geom_jitter(alpha=0.5, size=4, width=0, height=0.15) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())


ggplot(region1_females_invfreqscol, aes(x=mean_dist, y=prop_x)) + geom_point(alpha=0.5, colour="#F8766D") +
  geom_point(data=region1_males_invfreqscol, aes(x=mean_dist, y=prop_x), alpha=0.5, colour="#619CFF") +
  scale_y_continuous(limits=c(0,1)) + scale_x_continuous(limits=c(0,130))


theta.init <- list(c=80,lw=log(10), lp_cf=qlogis(0.01), lp_wf=qlogis(0.25), lp_cm=qlogis(0.01), lp_wm=qlogis(0.25))

mle.cline <- mle2(cline_mf, theta.init,  
                  control=list(parscale=abs(unlist(theta.init))),
                  data=list(pos=region1_combined_nona$DistAlongPath, hap=as.numeric(region1_combined_nona$cluster2),
                            sex=region1_combined_nona$sex))

summary(mle.cline)
confint(mle.cline)
AIC(mle.cline)

pars.m <- coef(mle.cline)

pos <- c(0:160) 
dist.m <- pos-pars.m[1]
freq.m <- 1/(1+exp(0-4*(dist.m)/(exp(pars.m[2]))))
mean.m_f <- (exp(pars.m[3])/(1+exp(pars.m[3])))+((exp(pars.m[4])/(1+exp(pars.m[4])))-(exp(pars.m[3])/(1+exp(pars.m[3]))))*freq.m
mean.m_m <- (exp(pars.m[5])/(1+exp(pars.m[5])))+((exp(pars.m[6])/(1+exp(pars.m[6])))-(exp(pars.m[5])/(1+exp(pars.m[5]))))*freq.m

region1_mf_fitted <- data.frame(dist.m, mean.m_f, mean.m_m)

ggplot(region1_mf_fitted, aes(x=dist.m+85.41, y=1-mean.m_m)) + geom_line(colour="dodgerblue4") + 
  geom_line(aes(y=1-mean.m_f), colour="brown4") +
  scale_x_continuous(limits=c(0,160)) +
  scale_y_continuous(limits=c(0,1)) + 
  geom_point(data=region1_females_invfreqscol, aes(x=mean_dist, y=prop_x), alpha=0.5, size=3, colour="#F8766D") +
  geom_point(data=region1_males_invfreqscol, aes(x=mean_dist, y=prop_x), alpha=0.5, size=3, colour="#619CFF") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank()) +
  xlab("Distance along path") + ylab("Inversion frequency") +
  geom_vline(xintercept=77.82, colour="red4", linetype=2, alpha=0.5)

r1_mf_pars <- pars.m
r1_mf_summary <- summary(mle.cline)
r1_mf_confint <- confint(mle.cline)
r1_mf_aiccline <- AIC(mle.cline)

r1_mf_pars_back <- r1_mf_pars[1]
r1_mf_pars_back[2] <- exp(r1_mf_pars[2])
r1_mf_pars_back[3] <- exp(r1_mf_pars[3])/(1+exp(r1_mf_pars[3]))
r1_mf_pars_back[4] <- exp(r1_mf_pars[4])/(1+exp(r1_mf_pars[4]))
r1_mf_pars_back[5] <- exp(r1_mf_pars[5])/(1+exp(r1_mf_pars[5]))
r1_mf_pars_back[6] <- exp(r1_mf_pars[6])/(1+exp(r1_mf_pars[6]))

####################################################################################

####################################################################################

## region 2a

region2a_combined_nona <- region2a_combined_nona %>% mutate(sex=ifelse(sex_dissection=="male", "m",
                                                                     ifelse(sex_dissection=="female", "f", NA)))

ggplot(region2a_combined, aes(x=DistAlongPath, y=cluster, colour=sex_dissection)) + 
  geom_jitter(alpha=0.5, size=4, width=0, height=0.15) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())


ggplot(region2a_females_invfreqscol, aes(x=mean_dist, y=prop_x)) + geom_point(alpha=0.5, colour="#F8766D") +
  geom_point(data=region2a_males_invfreqscol, aes(x=mean_dist, y=prop_x), alpha=0.5, colour="#619CFF") +
  scale_y_continuous(limits=c(0,1)) + scale_x_continuous(limits=c(0,130))


theta.init <- list(c=80,lw=log(10), lp_cf=qlogis(0.49), lp_wf=qlogis(0.45), lp_cm=qlogis(0.01), lp_wm=qlogis(0.25))

mle.cline <- mle2(cline_mf, theta.init,  
                  control=list(parscale=abs(unlist(theta.init))),
                  data=list(pos=region2a_combined_nona$DistAlongPath, hap=as.numeric(region2a_combined_nona$cluster2),
                            sex=region2a_combined_nona$sex))

summary(mle.cline)
confint(mle.cline)
AIC(mle.cline)

pars.m <- coef(mle.cline)

pos <- c(0:160) 
dist.m <- pos-pars.m[1]
freq.m <- 1/(1+exp(0-4*(dist.m)/(exp(pars.m[2]))))
mean.m_f <- (exp(pars.m[3])/(1+exp(pars.m[3])))+((exp(pars.m[4])/(1+exp(pars.m[4])))-(exp(pars.m[3])/(1+exp(pars.m[3]))))*freq.m
mean.m_m <- (exp(pars.m[5])/(1+exp(pars.m[5])))+((exp(pars.m[6])/(1+exp(pars.m[6])))-(exp(pars.m[5])/(1+exp(pars.m[5]))))*freq.m

region2a_mf_fitted <- data.frame(dist.m, mean.m_f, mean.m_m)

ggplot(region2a_mf_fitted, aes(x=dist.m+85.06, y=1-mean.m_m)) + geom_line(colour="dodgerblue4") + 
  geom_line(aes(y=1-mean.m_f), colour="brown4") +
  scale_x_continuous(limits=c(0,160)) +
  scale_y_continuous(limits=c(0,1)) + 
  geom_point(data=region2a_females_invfreqscol, aes(x=mean_dist, y=prop_x), alpha=0.5, size=3, colour="#F8766D") +
  geom_point(data=region2a_males_invfreqscol, aes(x=mean_dist, y=prop_x), alpha=0.5, size=3, colour="#619CFF") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank()) +
  xlab("Distance along path") + ylab("Inversion frequency") +
  geom_vline(xintercept=67, colour="orange2", linetype=2, alpha=0.8) +
  geom_vline(xintercept=85, colour="orange2", linetype=2, alpha=0.8) +
  geom_vline(xintercept=91.8, colour="green4", linetype=2, alpha=0.8) +
  theme(aspect.ratio = 1)

r2a_mf_pars <- pars.m
r2a_mf_summary <- summary(mle.cline)
r2a_mf_confint <- confint(mle.cline)
r2a_mf_aiccline <- AIC(mle.cline)

r2a_mf_pars_back <- r2a_mf_pars[1]
r2a_mf_pars_back[2] <- exp(r2a_mf_pars[2])
r2a_mf_pars_back[3] <- exp(r2a_mf_pars[3])/(1+exp(r2a_mf_pars[3]))
r2a_mf_pars_back[4] <- exp(r2a_mf_pars[4])/(1+exp(r2a_mf_pars[4]))
r2a_mf_pars_back[5] <- exp(r2a_mf_pars[5])/(1+exp(r2a_mf_pars[5]))
r2a_mf_pars_back[6] <- exp(r2a_mf_pars[6])/(1+exp(r2a_mf_pars[6]))

####################################################################################

####################################################################################

## region 2b

region2b_combined_nona <- region2b_combined_nona %>% mutate(sex=ifelse(sex_dissection=="male", "m",
                                                                       ifelse(sex_dissection=="female", "f", NA)))

ggplot(region2b_combined, aes(x=DistAlongPath, y=cluster, colour=sex_dissection)) + 
  geom_jitter(alpha=0.5, size=4, width=0, height=0.15) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())


ggplot(region2b_females_invfreqscol, aes(x=mean_dist, y=prop_x)) + geom_point(alpha=0.5, colour="#F8766D") +
  geom_point(data=region2b_males_invfreqscol, aes(x=mean_dist, y=prop_x), alpha=0.5, colour="#619CFF") +
  scale_y_continuous(limits=c(0,1)) + scale_x_continuous(limits=c(0,130))


theta.init <- list(c=80,lw=log(10), lp_cf=qlogis(0.49), lp_wf=qlogis(0.9), lp_cm=qlogis(0.99), lp_wm=qlogis(0.95))

mle.cline <- mle2(cline_mf, theta.init,  
                  control=list(parscale=abs(unlist(theta.init))),
                  data=list(pos=region2b_combined_nona$DistAlongPath, hap=as.numeric(region2b_combined_nona$cluster2),
                            sex=region2b_combined_nona$sex))

summary(mle.cline)
confint(mle.cline)
AIC(mle.cline)

pars.m <- coef(mle.cline)

pos <- c(0:160) 
dist.m <- pos-pars.m[1]
freq.m <- 1/(1+exp(0-4*(dist.m)/(exp(pars.m[2]))))
mean.m_f <- (exp(pars.m[3])/(1+exp(pars.m[3])))+((exp(pars.m[4])/(1+exp(pars.m[4])))-(exp(pars.m[3])/(1+exp(pars.m[3]))))*freq.m
mean.m_m <- (exp(pars.m[5])/(1+exp(pars.m[5])))+((exp(pars.m[6])/(1+exp(pars.m[6])))-(exp(pars.m[5])/(1+exp(pars.m[5]))))*freq.m

region2b_mf_fitted <- data.frame(dist.m, mean.m_f, mean.m_m)

ggplot(region2b_mf_fitted, aes(x=dist.m+88.75, y=1-mean.m_m)) + geom_line(colour="dodgerblue4") + 
  geom_line(aes(y=1-mean.m_f), colour="brown4") +
  scale_x_continuous(limits=c(0,160)) +
  scale_y_continuous(limits=c(0,1)) + 
  geom_point(data=region2b_females_invfreqscol, aes(x=mean_dist, y=prop_x), alpha=0.5, size=3, colour="#F8766D") +
  geom_point(data=region2b_males_invfreqscol, aes(x=mean_dist, y=prop_x), alpha=0.5, size=3, colour="#619CFF") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank()) +
  xlab("Distance along path") + ylab("Inversion frequency")

r2b_mf_pars <- pars.m
r2b_mf_summary <- summary(mle.cline)
r2b_mf_confint <- confint(mle.cline)
r2b_mf_aiccline <- AIC(mle.cline)

r2b_mf_pars_back <- r2b_mf_pars[1]
r2b_mf_pars_back[2] <- exp(r2b_mf_pars[2])
r2b_mf_pars_back[3] <- exp(r2b_mf_pars[3])/(1+exp(r2b_mf_pars[3]))
r2b_mf_pars_back[4] <- exp(r2b_mf_pars[4])/(1+exp(r2b_mf_pars[4]))
r2b_mf_pars_back[5] <- exp(r2b_mf_pars[5])/(1+exp(r2b_mf_pars[5]))
r2b_mf_pars_back[6] <- exp(r2b_mf_pars[6])/(1+exp(r2b_mf_pars[6]))


r2b_mf_confint_back <- data.frame(6,2)
r2b_mf_confint_back[1,] <- r2b_mf_confint[1,]
r2b_mf_confint_back[2,] <- exp(r2b_mf_confint[2,])
r2b_mf_confint_back[3,] <- exp(r2b_mf_confint[3,])/(1+exp(r2b_mf_confint[3,]))
r2b_mf_confint_back[4,] <- exp(r2b_mf_confint[4,])/(1+exp(r2b_mf_confint[4,]))
r2b_mf_confint_back[5,] <- exp(r2b_mf_confint[5,])/(1+exp(r2b_mf_confint[5,]))
r2b_mf_confint_back[6,] <- exp(r2b_mf_confint[6,])/(1+exp(r2b_mf_confint[6,]))

####################################################################################

####################################################################################

## region 3

region3_combined_nona <- region3_combined_nona %>% mutate(sex=ifelse(sex_dissection=="male", "m",
                                                                       ifelse(sex_dissection=="female", "f", NA)))

ggplot(region3_combined, aes(x=DistAlongPath, y=cluster, colour=sex_dissection)) + 
  geom_jitter(alpha=0.5, size=4, width=0, height=0.15) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())


ggplot(region3_females_invfreqscol, aes(x=mean_dist, y=1-prop_x)) + geom_point(alpha=0.5, colour="#F8766D") +
  geom_point(data=region3_males_invfreqscol, aes(x=mean_dist, y=1-prop_x), alpha=0.5, colour="#619CFF") +
  scale_y_continuous(limits=c(0,1)) + scale_x_continuous(limits=c(0,130))


theta.init <- list(c=80,lw=log(10), lp_cf=qlogis(0.55), lp_wf=qlogis(0.1), lp_cm=qlogis(0.3), lp_wm=qlogis(0.1))

mle.cline <- mle2(cline_mf, theta.init,  
                  control=list(parscale=abs(unlist(theta.init))),
                  data=list(pos=region3_combined_nona$DistAlongPath, hap=as.numeric(region3_combined_nona$cluster2),
                            sex=region3_combined_nona$sex))

summary(mle.cline)
confint(mle.cline)
AIC(mle.cline)

pars.m <- coef(mle.cline)

pos <- c(0:160) 
dist.m <- pos-pars.m[1]
freq.m <- 1/(1+exp(0-4*(dist.m)/(exp(pars.m[2]))))
mean.m_f <- (exp(pars.m[3])/(1+exp(pars.m[3])))+((exp(pars.m[4])/(1+exp(pars.m[4])))-(exp(pars.m[3])/(1+exp(pars.m[3]))))*freq.m
mean.m_m <- (exp(pars.m[5])/(1+exp(pars.m[5])))+((exp(pars.m[6])/(1+exp(pars.m[6])))-(exp(pars.m[5])/(1+exp(pars.m[5]))))*freq.m

region3_mf_fitted <- data.frame(dist.m, mean.m_f, mean.m_m)

ggplot(region3_mf_fitted, aes(x=dist.m+90.37, y=mean.m_m)) + geom_line(colour="dodgerblue4") + 
  geom_line(aes(y=mean.m_f), colour="brown4") +
  scale_x_continuous(limits=c(0,160)) +
  scale_y_continuous(limits=c(0,1)) + 
  geom_point(data=region3_females_invfreqscol, aes(x=mean_dist, y=1-prop_x), alpha=0.5, size=3, colour="#F8766D") +
  geom_point(data=region3_males_invfreqscol, aes(x=mean_dist, y=1-prop_x), alpha=0.5, size=3, colour="#619CFF") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank()) +
  xlab("Distance along path") + ylab("Inversion frequency") +
  geom_vline(xintercept=77.82, colour="red4", linetype=2, alpha=0.5)

r3_mf_pars <- pars.m
r3_mf_summary <- summary(mle.cline)
r3_mf_confint <- confint(mle.cline)
r3_mf_aiccline <- AIC(mle.cline)

r3_mf_pars_back <- r3_mf_pars[1]
r3_mf_pars_back[2] <- exp(r3_mf_pars[2])
r3_mf_pars_back[3] <- exp(r3_mf_pars[3])/(1+exp(r3_mf_pars[3]))
r3_mf_pars_back[4] <- exp(r3_mf_pars[4])/(1+exp(r3_mf_pars[4]))
r3_mf_pars_back[5] <- exp(r3_mf_pars[5])/(1+exp(r3_mf_pars[5]))
r3_mf_pars_back[6] <- exp(r3_mf_pars[6])/(1+exp(r3_mf_pars[6]))


r3_mf_confint_back <- data.frame(6,2)
r3_mf_confint_back[1,] <- r3_mf_confint[1,]
r3_mf_confint_back[2,] <- exp(r3_mf_confint[2,])
r3_mf_confint_back[3,] <- exp(r3_mf_confint[3,])/(1+exp(r3_mf_confint[3,]))
r3_mf_confint_back[4,] <- exp(r3_mf_confint[4,])/(1+exp(r3_mf_confint[4,]))
r3_mf_confint_back[5,] <- exp(r3_mf_confint[5,])/(1+exp(r3_mf_confint[5,]))
r3_mf_confint_back[6,] <- exp(r3_mf_confint[6,])/(1+exp(r3_mf_confint[6,]))

####################################################################################



########
########
######### Wave-constrained - constrain M and F wave freqs to be the same and see if the model is any better/worse than the current best model
######### - tests whether male and female Wave freqs differ
######### combined cline is best fit for 12.1 so no need to test wave-constrained
######### for inv12.2 and 12.3 one sex best model is no cline- incorporate into wave-constrained model for these (Fs no cline in 12.2 and Ms no cline in 12.3)
########
########

####################################################################################

# new cline function with no cline (one freq) for females, males have a cline but Wave freq is same as female freq

cline_wave <- function(hap, pos, sex, f_fwm, c_m, lw_m, lp_cm){     #f_fwm is the frequency in females and wave males; c_m is cline centre for males; lw_m is the cline width for males; lp_cm is the frequency for crab males
  d <- pos-c_m
  p_x <- 1/(1+exp(0-4*(d/exp(lw_m)))) 
  crab_m <- exp(lp_cm)/(1+exp(lp_cm))
  z_x <- exp(f_fwm)/(1+exp(f_fwm))
  
  z_x_f <- z_x
  z_x_m <- crab_m + (z_x-crab_m)*p_x
  
  minusll <- -sum(dbinom(hap[sex=="f"],2, z_x_f,log=T)) - sum(dbinom(hap[sex=="m"],2,z_x_m[sex=="m"],log=T))
  return(minusll)
}

####################################################################################

####################################################################################

## region 2a

#region2a_combined_nona <- region2a_combined_nona %>% mutate(sex=ifelse(sex_dissection=="male", "m",
#                                                                       ifelse(sex_dissection=="female", "f", NA)))

ggplot(region2a_combined, aes(x=DistAlongPath, y=cluster, colour=sex_dissection)) + 
  geom_jitter(alpha=0.5, size=4, width=0, height=0.15) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())


ggplot(region2a_females_invfreqscol, aes(x=mean_dist, y=prop_x)) + geom_point(alpha=0.5, colour="#F8766D") +
  geom_point(data=region2a_males_invfreqscol, aes(x=mean_dist, y=prop_x), alpha=0.5, colour="#619CFF") +
  scale_y_continuous(limits=c(0,1)) + scale_x_continuous(limits=c(0,130))


theta.init <- list(f_fwm=qlogis(0.49), c_m=80, lw_m=log(10), lp_cm=qlogis(0.01))

mle.cline <- mle2(cline_wave, theta.init,  
                  control=list(parscale=abs(unlist(theta.init))),
                  data=list(pos=region2a_combined_nona$DistAlongPath, hap=as.numeric(region2a_combined_nona$cluster2),
                            sex=region2a_combined_nona$sex))

summary(mle.cline)
confint(mle.cline)
AIC(mle.cline)

pars.m <- coef(mle.cline)

pos <- c(0:160) 
dist.m <- pos-pars.m[2]
freq.m <- 1/(1+exp(0-4*(dist.m)/(exp(pars.m[3]))))
mean.m_f <- (exp(pars.m[1])/(1+exp(pars.m[1])))+((exp(pars.m[1])/(1+exp(pars.m[1])))-(exp(pars.m[1])/(1+exp(pars.m[1]))))*freq.m
mean.m_m <- (exp(pars.m[4])/(1+exp(pars.m[4])))+((exp(pars.m[1])/(1+exp(pars.m[1])))-(exp(pars.m[4])/(1+exp(pars.m[4]))))*freq.m

region2a_wave_fitted <- data.frame(dist.m, mean.m_f, mean.m_m)

ggplot(region2a_wave_fitted, aes(x=dist.m+85.48, y=1-mean.m_m)) + geom_line(colour="dodgerblue4") + 
  geom_line(aes(y=1-mean.m_f), colour="brown4") +
  scale_x_continuous(limits=c(0,160)) +
  scale_y_continuous(limits=c(0,1)) + 
  geom_point(data=region2a_females_invfreqscol, aes(x=mean_dist, y=prop_x), alpha=0.5, size=3, colour="#F8766D") +
  geom_point(data=region2a_males_invfreqscol, aes(x=mean_dist, y=prop_x), alpha=0.5, size=3, colour="#619CFF") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank()) +
  xlab("Distance along path") + ylab("Inversion frequency") +
  geom_vline(xintercept=77.82, colour="red4", linetype=2, alpha=0.5)

r2a_wave_pars <- pars.m
r2a_wave_summary <- summary(mle.cline)
r2a_wave_confint <- confint(mle.cline)
r2a_wave_aiccline <- AIC(mle.cline)

r2a_wave_pars_back <- exp(r2a_wave_pars[1])/(1+exp(r2a_wave_pars[1]))
r2a_wave_pars_back[2] <- r2a_wave_pars[2]
r2a_wave_pars_back[3] <- exp(r2a_wave_pars[3])
r2a_wave_pars_back[4] <- exp(r2a_wave_pars[4])/(1+exp(r2a_wave_pars[4]))

####################################################################################

####################################################################################

# new cline function with no cline (one freq) for males, females have a cline but Wave freq is same as male freq

cline_wave2 <- function(hap, pos, sex, f_mwf, c_f, lw_f, lp_cf){     #f_fwm is the frequency in females and wave males; c_m is cline centre for males; lw_m is the cline width for males; lp_cm is the frequency for crab males
  d <- pos-c_f
  p_x <- 1/(1+exp(0-4*(d/exp(lw_f)))) 
  crab_f <- exp(lp_cf)/(1+exp(lp_cf))
  z_x <- exp(f_mwf)/(1+exp(f_mwf))
  
  z_x_m <- z_x
  z_x_f <- crab_f + (z_x-crab_f)*p_x
  
  minusll <- -sum(dbinom(hap[sex=="m"],2, z_x_m,log=T)) - sum(dbinom(hap[sex=="f"],2,z_x_f[sex=="f"],log=T))
  return(minusll)
}

####################################################################################

####################################################################################

## region 2b

#region2b_combined_nona <- region2b_combined_nona %>% mutate(sex=ifelse(sex_dissection=="male", "m",
#                                                                       ifelse(sex_dissection=="female", "f", NA)))

ggplot(region2b_combined, aes(x=DistAlongPath, y=cluster, colour=sex_dissection)) + 
  geom_jitter(alpha=0.5, size=4, width=0, height=0.15) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())


ggplot(region2b_females_invfreqscol, aes(x=mean_dist, y=prop_x)) + geom_point(alpha=0.5, colour="#F8766D") +
  geom_point(data=region2b_males_invfreqscol, aes(x=mean_dist, y=prop_x), alpha=0.5, colour="#619CFF") +
  scale_y_continuous(limits=c(0,1)) + scale_x_continuous(limits=c(0,130))


theta.init <- list(f_mwf=qlogis(0.95), c_f=80, lw_f=log(10), lp_cf=qlogis(0.95))

mle.cline <- mle2(cline_wave2, theta.init,  
                  control=list(parscale=abs(unlist(theta.init))),
                  data=list(pos=region2b_combined_nona$DistAlongPath, hap=as.numeric(region2b_combined_nona$cluster2),
                            sex=region2b_combined_nona$sex))

summary(mle.cline)
confint(mle.cline)
AIC(mle.cline)

pars.m <- coef(mle.cline)

pos <- c(0:160) 
dist.m <- pos-pars.m[2]
freq.m <- 1/(1+exp(0-4*(dist.m)/(exp(pars.m[3]))))
mean.m_m <- (exp(pars.m[1])/(1+exp(pars.m[1])))+((exp(pars.m[1])/(1+exp(pars.m[1])))-(exp(pars.m[1])/(1+exp(pars.m[1]))))*freq.m
mean.m_f <- (exp(pars.m[4])/(1+exp(pars.m[4])))+((exp(pars.m[1])/(1+exp(pars.m[1])))-(exp(pars.m[4])/(1+exp(pars.m[4]))))*freq.m

region2b_wave_fitted <- data.frame(dist.m, mean.m_f, mean.m_m)

ggplot(region2b_wave_fitted, aes(x=dist.m+90.33, y=1-mean.m_m)) + geom_line(colour="dodgerblue4") + 
  geom_line(aes(y=1-mean.m_f), colour="brown4") +
  scale_x_continuous(limits=c(0,160)) +
  scale_y_continuous(limits=c(0,1)) + 
  geom_point(data=region2b_females_invfreqscol, aes(x=mean_dist, y=prop_x), alpha=0.5, size=3, colour="#F8766D") +
  geom_point(data=region2b_males_invfreqscol, aes(x=mean_dist, y=prop_x), alpha=0.5, size=3, colour="#619CFF") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank()) +
  xlab("Distance along path") + ylab("Inversion frequency") +
  geom_vline(xintercept=67, colour="orange2", linetype=2, alpha=0.8) +
  geom_vline(xintercept=85, colour="orange2", linetype=2, alpha=0.8) +
  geom_vline(xintercept=91.8, colour="green4", linetype=2, alpha=0.8) +
  theme(aspect.ratio = 1)

r2b_wave_pars <- pars.m
r2b_wave_summary <- summary(mle.cline)
r2b_wave_confint <- confint(mle.cline)
r2b_wave_aiccline <- AIC(mle.cline)

r2b_wave_pars_back <- exp(r2b_wave_pars[1])/(1+exp(r2b_wave_pars[1]))
r2b_wave_pars_back[2] <- r2b_wave_pars[2]
r2b_wave_pars_back[3] <- exp(r2b_wave_pars[3])
r2b_wave_pars_back[4] <- exp(r2b_wave_pars[4])/(1+exp(r2b_wave_pars[4]))

r2b_wave_confint_back <- data.frame(4,2)
r2b_wave_confint_back[1,] <- exp(r2b_wave_confint[1,])/(1+exp(r2b_wave_confint[1,]))
r2b_wave_confint_back[2,] <- r2b_wave_confint[2,]
r2b_wave_confint_back[3,] <- exp(r2b_wave_confint[3,])
r2b_wave_confint_back[4,] <- exp(r2b_wave_confint[4,])/(1+exp(r2b_wave_confint[4,]))

####################################################################################


####################################################################################

# new cline function where cline centres, widths and wave freq are the same for males and females (only crab freq differs between males and females)

cline_wave3 <- function(hap, pos, sex, c, lw, lp_cf, lp_cm, lp_w){
  d <- pos-c
  p_x <- 1/(1+exp(0-4*(d/exp(lw)))) 
  
  crab_f <- exp(lp_cf)/(1+exp(lp_cf))
  crab_m <- exp(lp_cm)/(1+exp(lp_cm))
  wave <- exp(lp_w)/(1+exp(lp_w))
  
  z_x_f <- crab_f + (wave-crab_f)*p_x 
  z_x_m <- crab_m + (wave-crab_m)*p_x
  
  minusll <- -sum(dbinom(hap[sex == "f"],2,z_x_f[sex == "f"],log=T)) - sum(dbinom(hap[sex == "m"],2,z_x_m[sex == "m"],log=T)) # -sums of probabilities for females and males separately
  
  return(minusll)
}

####################################################################################

####################################################################################

## region 3

#region3_combined_nona <- region3_combined_nona %>% mutate(sex=ifelse(sex_dissection=="male", "m",
#                                                                       ifelse(sex_dissection=="female", "f", NA)))

ggplot(region3_combined, aes(x=DistAlongPath, y=cluster, colour=sex_dissection)) + 
  geom_jitter(alpha=0.5, size=4, width=0, height=0.15) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())


ggplot(region3_females_invfreqscol, aes(x=mean_dist, y=prop_x)) + geom_point(alpha=0.5, colour="#F8766D") +
  geom_point(data=region3_males_invfreqscol, aes(x=mean_dist, y=prop_x), alpha=0.5, colour="#619CFF") +
  scale_y_continuous(limits=c(0,1)) + scale_x_continuous(limits=c(0,130))


theta.init <- list(c=80, lw=log(10), lp_cf=qlogis(0.55), lp_cm=qlogis(0.3), lp_w=qlogis(0.1))

mle.cline <- mle2(cline_wave3, theta.init,  
                  control=list(parscale=abs(unlist(theta.init))),
                  data=list(pos=region3_combined_nona$DistAlongPath, hap=as.numeric(region3_combined_nona$cluster2),
                            sex=region3_combined_nona$sex))

summary(mle.cline)
confint(mle.cline)
AIC(mle.cline)

pars.m <- coef(mle.cline)

pos <- c(0:160) 
dist.m <- pos-pars.m[1]
freq.m <- 1/(1+exp(0-4*(dist.m)/(exp(pars.m[2]))))
mean.m_f <- (exp(pars.m[3])/(1+exp(pars.m[3])))+((exp(pars.m[5])/(1+exp(pars.m[5])))-(exp(pars.m[3])/(1+exp(pars.m[3]))))*freq.m
mean.m_m <- (exp(pars.m[4])/(1+exp(pars.m[4])))+((exp(pars.m[5])/(1+exp(pars.m[5])))-(exp(pars.m[4])/(1+exp(pars.m[4]))))*freq.m

region3_wave_fitted <- data.frame(dist.m, mean.m_f, mean.m_m)

ggplot(region3_wave_fitted, aes(x=dist.m+90.43, y=mean.m_m)) + geom_line(colour="dodgerblue4") + 
  geom_line(aes(y=mean.m_f), colour="brown4") +
  scale_x_continuous(limits=c(0,160)) +
  scale_y_continuous(limits=c(0,1)) + 
  geom_point(data=region3_females_invfreqscol, aes(x=mean_dist, y=1-prop_x), alpha=0.5, size=3, colour="#F8766D") +
  geom_point(data=region3_males_invfreqscol, aes(x=mean_dist, y=1-prop_x), alpha=0.5, size=3, colour="#619CFF") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank()) +
  xlab("Distance along path") + ylab("Inversion frequency") +
  geom_vline(xintercept=67, colour="orange2", linetype=2, alpha=0.8) +
  geom_vline(xintercept=85, colour="orange2", linetype=2, alpha=0.8) +
  geom_vline(xintercept=91.8, colour="green4", linetype=2, alpha=0.8) +
  theme(aspect.ratio = 1)

r3_wave_pars <- pars.m
r3_wave_summary <- summary(mle.cline)
r3_wave_confint <- confint(mle.cline)
r3_wave_aiccline <- AIC(mle.cline)


r3_wave_pars_back <- r3_wave_pars[1]
r3_wave_pars_back[2] <- exp(r3_wave_pars[2])
r3_wave_pars_back[3] <- exp(r3_wave_pars[3])/(1+exp(r3_wave_pars[3]))
r3_wave_pars_back[4] <- exp(r3_wave_pars[4])/(1+exp(r3_wave_pars[4]))
r3_wave_pars_back[5] <- exp(r3_wave_pars[5])/(1+exp(r3_wave_pars[5]))

r3_wave_confint_back <- data.frame(4,2)
r3_wave_confint_back[1,] <- r3_wave_confint[1,]
r3_wave_confint_back[2,] <- exp(r3_wave_confint[2,])
r3_wave_confint_back[3,] <- exp(r3_wave_confint[3,])/(1+exp(r3_wave_confint[3,]))
r3_wave_confint_back[4,] <- exp(r3_wave_confint[4,])/(1+exp(r3_wave_confint[4,]))
r3_wave_confint_back[5,] <- exp(r3_wave_confint[5,])/(1+exp(r3_wave_confint[5,]))

####################################################################################


####################
## plot for inversion 12.2 that is female no cline and male full cline (the best fitting combination)

ggplot(region2a_male_fitted, aes(x=dist.m+83.76, y=1-mean.m)) + geom_line(colour="dodgerblue4") + 
  geom_line(aes(y=1-r2a_f_control_pars_back), alpha=0.7, colour="brown4") +
  scale_x_continuous(limits=c(0,160)) +
  scale_y_continuous(limits=c(0,1)) + 
  geom_point(data=region2a_females_invfreqscol, aes(x=mean_dist, y=prop_x), alpha=0.5, size=3, colour="#F8766D") +
  geom_point(data=region2a_males_invfreqscol, aes(x=mean_dist, y=prop_x), alpha=0.5, size=3, colour="#619CFF") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank()) +
  xlab("Distance along path") + ylab("Inversion frequency") +
  geom_vline(xintercept=67, colour="orange2", linetype=2, alpha=0.8) +
  geom_vline(xintercept=85, colour="orange2", linetype=2, alpha=0.8) +
  geom_vline(xintercept=91.8, colour="green4", linetype=2, alpha=0.8) +
  theme(aspect.ratio = 1)



