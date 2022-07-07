rm(list=ls())
# set working directory
library(tidyverse)

######################################################
##                                                  ##
##             Pi and dxy analysis                  ##
##          Data generated on Sharc                 ##
##	(Using scripts from Simon Martin)		 ##
##                                                  ##
######################################################

## read in the stats data files
inv1_stats <- read.csv("inv1_stats_50000_25.csv")
inv2a_stats <- read.csv("inv2a_stats_50000_25.csv")
inv2b_stats <- read.csv("inv2b_stats_50000_25.csv")
inv3_stats <- read.csv("inv3_stats_50000_25.csv")

## read in the CW split data files
inv1_stats_CW <- read.csv("inv1_stats_CW.csv")
inv2a_stats_CW <- read.csv("inv2a_stats_CW.csv")
inv2b_stats_CW <- read.csv("inv2b_stats_CW.csv")
inv3_stats_CW <- read.csv("inv3_stats_CW.csv")

## read in the geno split only data files
inv1_stats_genoonly <- read.csv("inv1_stats_genoonly.csv") %>% slice(-1)
inv2a_stats_genoonly <- read.csv("inv2a_stats_genoonly.csv") %>% slice(-1)
inv2b_stats_genoonly <- read.csv("inv2b_stats_genoonly.csv") %>% slice(-1)
inv3_stats_genoonly <- read.csv("inv3_stats_genoonly.csv") %>% slice(-1)

## read in other data files needed (map positions etc)
map <- read.csv("LG12_mappos_bycontig.csv") %>% separate(contig_ID, c("contig", "pos"), sep="_") %>%
  select(-pos) %>% distinct()


## make data tidy
## first gather into group(s) and stat value
## then split group(s) col into stat type and (group(s)
inv1_stats_tidy <- inv1_stats %>% gather(group, stat_value, 7:42) %>% 
  separate(group, c("stat_type", "group"), sep="_", extra="merge")
inv2a_stats_tidy <- inv2a_stats %>% gather(group, stat_value, 7:42) %>% 
  separate(group, c("stat_type", "group"), sep="_", extra="merge")
inv2b_stats_tidy <- inv2b_stats %>% gather(group, stat_value, 7:42) %>% 
  separate(group, c("stat_type", "group"), sep="_", extra="merge")
inv3_stats_tidy <- inv3_stats %>% gather(group, stat_value, 7:42) %>% 
  separate(group, c("stat_type", "group"), sep="_", extra="merge")

inv1_stats_tidy_CW <- inv1_stats_CW %>% gather(group, stat_value, 7:70) %>% 
  separate(group, c("stat_type", "group"), sep="_", extra="merge")
inv2a_stats_tidy_CW <- inv2a_stats_CW %>% gather(group, stat_value, 7:127) %>% 
  separate(group, c("stat_type", "group"), sep="_", extra="merge")
inv2b_stats_tidy_CW <- inv2b_stats_CW %>% gather(group, stat_value, 7:106) %>% 
  separate(group, c("stat_type", "group"), sep="_", extra="merge")
inv3_stats_tidy_CW <- inv3_stats_CW %>% gather(group, stat_value, 7:127) %>% 
  separate(group, c("stat_type", "group"), sep="_", extra="merge")

## add map info
inv1_stats_tidy <- inv1_stats_tidy %>% left_join(., map, by=c("scaffold"="contig"))
inv2a_stats_tidy <- inv2a_stats_tidy %>% left_join(., map, by=c("scaffold"="contig"))
inv2b_stats_tidy <- inv2b_stats_tidy %>% left_join(., map, by=c("scaffold"="contig"))
inv3_stats_tidy <- inv3_stats_tidy %>% left_join(., map, by=c("scaffold"="contig"))

inv1_stats_tidy_CW <- inv1_stats_tidy_CW %>% left_join(., map, by=c("scaffold"="contig"))
inv2a_stats_tidy_CW <- inv2a_stats_tidy_CW %>% left_join(., map, by=c("scaffold"="contig"))
inv2b_stats_tidy_CW <- inv2b_stats_tidy_CW %>% left_join(., map, by=c("scaffold"="contig"))
inv3_stats_tidy_CW <- inv3_stats_tidy_CW %>% left_join(., map, by=c("scaffold"="contig"))

#####
# recode to the RR/RA/AA genotypes for inversions
inv1_stats_tidy$group <- str_replace_all(inv1_stats_tidy$group, c("het"="RA", "hom1"="RR", "hom2"="AA"))
inv2a_stats_tidy$group <- str_replace_all(inv2a_stats_tidy$group, c("het"="RA", "hom1"="RR", "hom2"="AA"))
inv2b_stats_tidy$group <- str_replace_all(inv2b_stats_tidy$group, c("het"="RA", "hom1"="RR", "hom2"="AA"))
inv3_stats_tidy$group <- str_replace_all(inv3_stats_tidy$group, c("het"="RA", "hom1"="AA", "hom2"="RR"))

inv1_stats_tidy_CW$group <- str_replace_all(inv1_stats_tidy_CW$group, c("het"="RA", "hom1"="RR", "hom2"="AA"))
inv2a_stats_tidy_CW$group <- str_replace_all(inv2a_stats_tidy_CW$group, c("het"="RA", "hom1"="RR", "hom2"="AA"))
inv2b_stats_tidy_CW$group <- str_replace_all(inv2b_stats_tidy_CW$group, c("het"="RA", "hom1"="RR", "hom2"="AA"))
inv3_stats_tidy_CW$group <- str_replace_all(inv3_stats_tidy_CW$group, c("het"="RA", "hom1"="AA", "hom2"="RR"))

inv3_stats_tidy$group <- str_replace_all(inv3_stats_tidy$group, c("F_AA_F_RA"="F_RA_F_AA", 
                                                                  "F_AA_F_RR"="F_RR_F_AA", 
                                                                  "F_RA_F_RR"="F_RR_F_RA",
                                                                  "M_AA_M_RA"="M_RA_M_AA",
                                                                  "M_AA_M_RR"="M_RR_M_AA",
                                                                  "M_RA_M_RR"="M_RR_M_RA"))
inv3_stats_tidy_CW$group <- str_replace_all(inv3_stats_tidy_CW$group, c("F_AA_C_F_RR_C"="F_RR_C_F_AA_C",
                                                                        "F_AA_C_F_RR_W"="F_RR_W_F_AA_C",
                                                                        "F_AA_W_F_RR_C"="F_RR_C_F_AA_W",
                                                                        "F_AA_W_F_RR_W"="F_RR_W_F_AA_C",
                                                                        "M_AA_C_M_RR_C"="M_RR_C_M_AA_C",
                                                                        "M_AA_C_M_RR_W"="M_RR_W_M_AA_C",
                                                                        "M_AA_W_M_RR_C"="M_RR_C_M_AA_W",
                                                                        "M_AA_W_M_RR_W"="M_RR_W_M_AA_C"))




## combined data frame 
inv1_stats_tidy <- inv1_stats_tidy %>% mutate(inv="1")
inv2a_stats_tidy <- inv2a_stats_tidy %>% mutate(inv="2")
inv2b_stats_tidy <- inv2b_stats_tidy %>% mutate(inv="3")
inv3_stats_tidy <- inv3_stats_tidy %>% mutate(inv="4")
inv_stats_tidy <- bind_rows(inv1_stats_tidy, inv2a_stats_tidy) %>% bind_rows(., inv2b_stats_tidy) %>% 
  bind_rows(., inv3_stats_tidy)

inv1_stats_tidy_CW <- inv1_stats_tidy_CW %>% mutate(inv="1")
inv2a_stats_tidy_CW <- inv2a_stats_tidy_CW %>% mutate(inv="2")
inv2b_stats_tidy_CW <- inv2b_stats_tidy_CW %>% mutate(inv="3")
inv3_stats_tidy_CW <- inv3_stats_tidy_CW %>% mutate(inv="4")
inv_stats_tidy_CW <- bind_rows(inv1_stats_tidy_CW, inv2a_stats_tidy_CW) %>% bind_rows(., inv2b_stats_tidy_CW) %>% 
  bind_rows(., inv3_stats_tidy_CW)



## remove comparisons/groups i don't need
reduced_set <- inv_stats_tidy %>% filter(!grepl("RA", group))
reduced_set_CW <- inv_stats_tidy_CW %>% filter(!grepl("RA",group))


#############
## use the set with hets included
#############

## add grouping info (male/female)
## add colour-coding for comparison type (between-sex, between-ecotype, between-arrangement)
## remove unnecessary dxy comparisons (where more than one above factor changes)


## same but not using the means per map position
reduced_set %>% filter(stat_type=="pi") %>% ggplot(., aes(x=group, y=stat_value, colour=group)) +
  geom_boxplot() + facet_wrap(~inv, nrow=1)

reduced_set %>% filter(stat_type=="dxy") %>% ggplot(., aes(x=group, y=stat_value, colour=group)) +
  geom_boxplot() + facet_wrap(~inv, nrow=1) + theme(axis.text.x = element_text(angle = 90))

reduced_set_CW %>% filter(stat_type=="pi") %>% ggplot(., aes(x=group, y=stat_value, colour=group)) +
  geom_boxplot() + facet_wrap(~inv, nrow=1) + theme(axis.text.x = element_text(angle = 90))

reduced_set_CW %>% filter(stat_type=="dxy") %>% ggplot(., aes(x=group, y=stat_value, colour=group)) +
  geom_boxplot() + facet_wrap(~inv) + theme(axis.text.x = element_text(angle = 90))


########
# further reduce the dataset to only the most important dxy comparisons
extra_reduced_set <- reduced_set %>% filter(group!="F_AA_M_RR" & group!="F_RR_M_AA")
extra_reduced_set_CW <- reduced_set_CW %>% filter(group!="F_AA_C_M_AA_W" & 
                                                                group!="F_AA_C_M_RR_C" &
                                                                group!="F_AA_C_M_RR_W" &
                                                                group!="F_AA_W_M_AA_C" &
                                                                group!="F_AA_W_M_RR_C" &
                                                                group!="F_AA_W_M_RR_W" &
                                                                group!="F_RR_C_F_AA_W" &
                                                                group!="F_RR_C_M_AA_C" &
                                                                group!="F_RR_C_M_AA_W" &
                                                                group!="F_RR_C_M_RR_W" &
                                                                group!="F_RR_W_F_AA_C" &
                                                                group!="F_RR_W_M_AA_C" &
                                                                group!="F_RR_W_M_AA_W" &
                                                                group!="F_RR_W_M_RR_C" &
                                                                group!="M_RR_C_M_AA_W")


extra_reduced_set %>% filter(stat_type=="dxy") %>% ggplot(., aes(x=group, y=stat_value, colour=group)) +
  geom_boxplot() + facet_wrap(~inv, nrow=1) + theme(text = element_text(size=20), axis.text.x = element_text(angle = 90))

extra_reduced_set_CW %>% filter(stat_type=="dxy") %>% ggplot(., aes(x=group, y=stat_value, colour=group)) +
  geom_boxplot() + facet_wrap(~inv) + theme(text = element_text(size=15), axis.text.x = element_text(angle = 90), legend.position = "none")

## add colour-coding for between-sex, between-ecotype, between-arrangement comparisons

extra_reduced_set <- extra_reduced_set %>% 
  mutate(comp_type=ifelse(stat_type=="dxy" & !grepl("F", group), "between_arrangement",
                          ifelse(stat_type=="dxy" & !grepl("M", group), "between_arrangement", "between_sex")))

extra_reduced_set_CW <- extra_reduced_set_CW %>% 
  mutate(comp_type=ifelse(stat_type=="dxy" & !grepl("F", group) & !grepl("C", group), "between_arrangement",
                          ifelse(stat_type=="dxy" & !grepl("F", group) & !grepl("W", group), "between_arrangement", 
                          ifelse(stat_type=="dxy" & !grepl("M", group) & !grepl("C", group), "between_arrangement",
                          ifelse(stat_type=="dxy" & !grepl("M", group) & !grepl("W", group), "between_arrangement",
                          ifelse(stat_type=="dxy" & !grepl("F", group) & !grepl("RR", group), "between_ecotype",
                          ifelse(stat_type=="dxy" & !grepl("F", group) & !grepl("AA", group), "between_ecotype",
                          ifelse(stat_type=="dxy" & !grepl("M", group) & !grepl("RR", group), "between_ecotype",
                          ifelse(stat_type=="dxy" & !grepl("M", group) & !grepl("AA", group), "between_ecotype", "between_sex")))))))))


extra_reduced_set %>% filter(stat_type=="dxy") %>% ggplot(., aes(x=group, y=stat_value, colour=comp_type)) +
  geom_boxplot() + facet_wrap(~inv, nrow=1) + theme(text = element_text(size=20), axis.text.x = element_text(angle = 90), legend.text=element_text(size=10))

extra_reduced_set_CW %>% filter(stat_type=="dxy") %>% ggplot(., aes(x=group, y=stat_value, colour=comp_type)) +
  geom_boxplot() + facet_wrap(~inv) + theme(text = element_text(size=15), axis.text.x = element_text(angle = 90), legend.text=element_text(size=10))




# other plots 
inv_stats_tidy %>% filter(stat_type=="pi") %>% 
  ggplot(., aes(x=group, y=stat_value, colour=inv)) + geom_boxplot() + facet_wrap(~av) +
  theme(axis.text.x = element_text(angle = 90))

inv_stats_tidy_CW %>% filter(stat_type=="pi") %>% 
  ggplot(., aes(x=group, y=stat_value, colour=inv)) + geom_boxplot() + facet_wrap(~av) +
  theme(axis.text.x = element_text(angle = 90))

inv_stats_tidy %>% ggplot(., aes(x=group, y=stat_value)) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90)) + facet_wrap(~stat_type)


reduced_set_CW %>% filter(stat_type=="pi") %>% ggplot(., aes(x=av, y=stat_value, colour=inv)) + geom_point(alpha=0.5) + facet_wrap(~group)





###############################################################

###############################################################

#### 
## using the pi for genotype groups (no sex/ecotype split) and calculating pi/dxy for R and A
## (see notes)

## first get genoonly datasets into tidy form

# gather and add map info
inv1_stats_genoonly_tidy <- inv1_stats_genoonly %>% gather(group, stat_value, 7:15) %>% 
  separate(group, c("stat_type", "group"), sep="_", extra="merge")%>%
  left_join(., map, by=c("scaffold"="contig"))
inv2a_stats_genoonly_tidy <- inv2a_stats_genoonly %>% gather(group, stat_value, 7:15) %>% 
  separate(group, c("stat_type", "group"), sep="_", extra="merge")%>%
  left_join(., map, by=c("scaffold"="contig"))
inv2b_stats_genoonly_tidy <- inv2b_stats_genoonly %>% gather(group, stat_value, 7:15) %>% 
  separate(group, c("stat_type", "group"), sep="_", extra="merge")%>%
  left_join(., map, by=c("scaffold"="contig"))
inv3_stats_genoonly_tidy <- inv3_stats_genoonly %>% gather(group, stat_value, 7:15) %>% 
  separate(group, c("stat_type", "group"), sep="_", extra="merge")%>%
  left_join(., map, by=c("scaffold"="contig"))

# recode to R/A
inv1_stats_genoonly_tidy$group <- str_replace_all(inv1_stats_genoonly_tidy$group, c("het"="RA", "hom1"="RR", "hom2"="AA"))
inv2a_stats_genoonly_tidy$group <- str_replace_all(inv2a_stats_genoonly_tidy$group, c("het"="RA", "hom1"="RR", "hom2"="AA"))
inv2b_stats_genoonly_tidy$group <- str_replace_all(inv2b_stats_genoonly_tidy$group, c("het"="RA", "hom1"="RR", "hom2"="AA"))
inv3_stats_genoonly_tidy$group <- str_replace_all(inv3_stats_genoonly_tidy$group, c("het"="RA", "hom1"="AA", "hom2"="RR"))

# make into one data frame
inv1_stats_genoonly_tidy <- inv1_stats_genoonly_tidy %>% mutate(inv="1")
inv2a_stats_genoonly_tidy <- inv2a_stats_genoonly_tidy %>% mutate(inv="2")
inv2b_stats_genoonly_tidy <- inv2b_stats_genoonly_tidy %>% mutate(inv="3")
inv3_stats_genoonly_tidy <- inv3_stats_genoonly_tidy %>% mutate(inv="4")
inv_stats_genoonly_tidy <- bind_rows(inv1_stats_genoonly_tidy, inv2a_stats_genoonly_tidy) %>% 
  bind_rows(., inv2b_stats_genoonly_tidy) %>% bind_rows(., inv3_stats_genoonly_tidy)


# quick look at boxplot of pi per genotype
inv_stats_genoonly_tidy %>% filter(stat_type=="pi") %>% ggplot(aes(x=group, y=stat_value)) + 
  geom_boxplot() + facet_wrap(~inv)


##########
## calculating dxy between arrangements 
## dxy = 2*(pi_het - pi_AA/4 - pi_RR/4)

inv1_genoonly_dxyarr <- inv1_stats_genoonly %>% mutate(dxy_A_R=2*(pi_het-(pi_hom1/4)-(pi_hom2/4))) 
inv2a_genoonly_dxyarr <- inv2a_stats_genoonly %>% mutate(dxy_A_R=2*(pi_het-(pi_hom1/4)-(pi_hom2/4))) 
inv2b_genoonly_dxyarr <- inv2b_stats_genoonly %>% mutate(dxy_A_R=2*(pi_het-(pi_hom1/4)-(pi_hom2/4))) 
inv3_genoonly_dxyarr <- inv3_stats_genoonly %>% mutate(dxy_A_R=2*(pi_het-(pi_hom1/4)-(pi_hom2/4))) 

# use previous code to add map info, get data tidy etc
inv1_genoonly_dxyarr_tidy <- inv1_genoonly_dxyarr %>% gather(group, stat_value, 7:16) %>% 
  separate(group, c("stat_type", "group"), sep="_", extra="merge")%>%
  left_join(., map, by=c("scaffold"="contig"))
inv2a_genoonly_dxyarr_tidy <- inv2a_genoonly_dxyarr %>% gather(group, stat_value, 7:16) %>% 
  separate(group, c("stat_type", "group"), sep="_", extra="merge")%>%
  left_join(., map, by=c("scaffold"="contig"))
inv2b_genoonly_dxyarr_tidy <- inv2b_genoonly_dxyarr %>% gather(group, stat_value, 7:16) %>% 
  separate(group, c("stat_type", "group"), sep="_", extra="merge")%>%
  left_join(., map, by=c("scaffold"="contig"))
inv3_genoonly_dxyarr_tidy <- inv3_genoonly_dxyarr %>% gather(group, stat_value, 7:16) %>% 
  separate(group, c("stat_type", "group"), sep="_", extra="merge")%>%
  left_join(., map, by=c("scaffold"="contig"))

inv1_genoonly_dxyarr_tidy$group <- str_replace_all(inv1_genoonly_dxyarr_tidy$group, c("het"="RA", "hom1"="RR", "hom2"="AA"))
inv2a_genoonly_dxyarr_tidy$group <- str_replace_all(inv2a_genoonly_dxyarr_tidy$group, c("het"="RA", "hom1"="RR", "hom2"="AA"))
inv2b_genoonly_dxyarr_tidy$group <- str_replace_all(inv2b_genoonly_dxyarr_tidy$group, c("het"="RA", "hom1"="RR", "hom2"="AA"))
inv3_genoonly_dxyarr_tidy$group <- str_replace_all(inv3_genoonly_dxyarr_tidy$group, c("het"="RA", "hom1"="AA", "hom2"="RR"))

inv1_genoonly_dxyarr_tidy <- inv1_genoonly_dxyarr_tidy %>% mutate(inv="1")
inv2a_genoonly_dxyarr_tidy <- inv2a_genoonly_dxyarr_tidy %>% mutate(inv="2")
inv2b_genoonly_dxyarr_tidy <- inv2b_genoonly_dxyarr_tidy %>% mutate(inv="3")
inv3_genoonly_dxyarr_tidy <- inv3_genoonly_dxyarr_tidy %>% mutate(inv="4")
inv_genoonly_dxyarr_tidy <- bind_rows(inv1_genoonly_dxyarr_tidy, inv2a_genoonly_dxyarr_tidy) %>% 
  bind_rows(., inv2b_genoonly_dxyarr_tidy) %>% bind_rows(., inv3_genoonly_dxyarr_tidy)

means_genoonly_dxyarr <- inv_genoonly_dxyarr_tidy %>% group_by(inv, av, stat_type, group) %>% 
  mutate(mean_stat_value=mean(stat_value, na.rm=TRUE))

# plots 
inv_genoonly_dxyarr_tidy %>% filter(stat_type=="dxy") %>% ggplot(aes(x=group, y=stat_value)) + 
  geom_boxplot() + facet_wrap(~inv)
means_genoonly_dxyarr %>% filter(stat_type=="dxy") %>% ggplot(aes(x=group, y=mean_stat_value)) + 
  geom_boxplot() + facet_wrap(~inv, nrow=1)







######################
######################


#####
## linear models testing
#####

#library(lme4)


#######################
#######################

### try the MuMIn package for automated model checking
library(MuMIn)

###
##############################
###
###
### Remove groups with <2 snails and reorder groups so intercept is a group w snails
### On all of them make F_RA_W the intercept (order F_RA_W, M_RA_W, F_RA_C, M_RA_C, F_AA_W, M_AA_W, F_AA_C, M_AA_C, F_RR_W, M_RR_W, F_RR_C, M_RR_C)
###


##### inversion 1 - model fitting 
## remove F_AA_C, M_AA_C, M_AA_W, F_RA_C, M_RA_C
inv1_filtered <- inv1_stats_tidy_CW %>% filter(stat_type=="pi") %>%
  mutate(log_stat=log(stat_value)) %>% 
  filter(group!="F_AA_C" & group!="M_AA_C" & group!="M_AA_W" & group!="F_RA_C" & group!="M_RA_C") %>%
  separate(group, c("sex", "genotype", "ecotype"), sep="_") %>%
  filter(!is.nan(log_stat) & !is.infinite(log_stat))
inv1_filtered$sex <- factor(inv1_filtered$sex, levels=c("F", "M"))
inv1_filtered$ecotype <- factor(inv1_filtered$ecotype, levels=c("W", "C"))
inv1_filtered$genotype <- factor(inv1_filtered$genotype, levels=c("RA", "AA", "RR"))
inv1_filtered[order(inv1_filtered$genotype, inv1_filtered$sex, inv1_filtered$ecotype),]
options(na.action="na.omit")
inv1_filtered_model <- lmer(log_stat ~ genotype*sex*ecotype + (1 + genotype*sex*ecotype | av), data=inv1_filtered)
options(na.action="na.fail")
inv1_filtered_modelsel <- dredge(inv1_filtered_model, trace=TRUE, rank="AICc", REML=FALSE)
inv1_fmList_filtered <- get.models(inv1_filtered_modelsel, subset=delta<2)
inv1_top <- subset(inv1_filtered_modelsel, delta<2)
#inv1_filtered_modelavg <- model.avg(inv1_fmList_filtered, revised.var=TRUE)
#inv1_filtered_modelavg_summary <- summary(inv1_filtered_modelavg)
#inv1_df <- as.data.frame(inv1_filtered_modelavg_summary$coefmat.full)
inv1_filtered_estimates <- as.data.frame(coef(summary(inv1_fmList_filtered[[1]])))
inv1_filtered_estimates <- rownames_to_column(inv1_filtered_estimates, "Coefficients")
names(inv1_filtered_estimates) <- gsub(" ", "", names(inv1_filtered_estimates))

# convert estimates from model to estimates for groups
inv1_filtered_groupestimates <- inv1_filtered_estimates %>% select(Coefficients, Estimate) %>% 
  spread(Coefficients, Estimate) %>%
  rename(Intercept=`(Intercept)`) %>%
  mutate(RA=Intercept,
         AA=Intercept+genotypeAA,
         RR=Intercept+genotypeRR) %>%
  gather(Group, Estimate, 4:6) %>%
  mutate(Estimate_exp=exp(Estimate))


##### inversion 2 - model fitting 
## remove F_AA_C, M_AA_C, F_RR_C
inv2_filtered <- inv2a_stats_tidy_CW %>% filter(stat_type=="pi") %>%
  mutate(log_stat=log(stat_value)) %>% 
  filter(group!="F_AA_C" & group!="M_AA_C" & group!="F_RR_C") %>%
  separate(group, c("sex", "genotype", "ecotype"), sep="_") %>%
  filter(!is.nan(log_stat) & !is.infinite(log_stat))
inv2_filtered$sex <- factor(inv2_filtered$sex, levels=c("F", "M"))
inv2_filtered$ecotype <- factor(inv2_filtered$ecotype, levels=c("W", "C"))
inv2_filtered$genotype <- factor(inv2_filtered$genotype, levels=c("RA", "AA", "RR"))
inv2_filtered[order(inv2_filtered$genotype, inv2_filtered$sex, inv2_filtered$ecotype),]
options(na.action="na.omit")
inv2_filtered_model <- lmer(log_stat ~ genotype*sex*ecotype + (1 + genotype*sex*ecotype | av), data=inv2_filtered)
options(na.action="na.fail")
inv2_filtered_modelsel <- dredge(inv2_filtered_model, trace=TRUE, rank="AICc", REML=FALSE)
inv2_fmList_filtered <- get.models(inv2_filtered_modelsel, subset=delta<2)
inv2_top <- subset(inv2_filtered_modelsel, delta<2)
#inv2_filtered_modelavg <- model.avg(inv2_fmList_filtered, revised.var=TRUE)
#inv2_filtered_modelavg_summary <- summary(inv2_filtered_modelavg)
#inv2_df <- as.data.frame(inv2_filtered_modelavg_summary$coefmat.full)
inv2_filtered_estimates <- as.data.frame(coef(summary(inv2_fmList_filtered[[1]])))
inv2_filtered_estimates <- rownames_to_column(inv2_filtered_estimates, "Coefficients")
names(inv2_filtered_estimates) <- gsub(" ", "", names(inv2_filtered_estimates))

# convert estimates from model to estimates for groups
inv2_filtered_groupestimates <- inv2_filtered_estimates %>% select(Coefficients, Estimate) %>% 
  spread(Coefficients, Estimate) %>%
  rename(Intercept=`(Intercept)`) %>%
  mutate(#F_AA_C,
         F_AA_W=Intercept+genotypeAA,
         F_RA_C=Intercept+ecotypeC,
         F_RA_W=Intercept,
         #F_RR_C,
         F_RR_W=Intercept+genotypeRR,
         #M_AA_C,
         M_AA_W=Intercept+genotypeAA+sexM,
         M_RA_C=Intercept+sexM+ecotypeC+`ecotypeC:sexM`,
         M_RA_W=Intercept+sexM,
         M_RR_C=Intercept+genotypeRR+ecotypeC+sexM+`ecotypeC:genotypeRR`+`ecotypeC:sexM`,
         M_RR_W=Intercept+genotypeRR+sexM) %>%
  gather(Group, Estimate, 8:16) %>%
  mutate(Estimate_exp=exp(Estimate))


##### inversion 3 - model fitting 
## remove F_AA_C, F_RR_C, F_RR_W, M_RR_C, M_RR_W
inv3_filtered <- inv2b_stats_tidy_CW %>% filter(stat_type=="pi") %>%
  mutate(log_stat=log(stat_value)) %>% 
  filter(group!="F_AA_C" & group!="F_RR_C" & group!="F_RR_W" & group!="M_RR_C" & group!="M_RR_W") %>%
  separate(group, c("sex", "genotype", "ecotype"), sep="_") %>%
  filter(!is.nan(log_stat) & !is.infinite(log_stat))
inv3_filtered$sex <- factor(inv3_filtered$sex, levels=c("F", "M"))
inv3_filtered$ecotype <- factor(inv3_filtered$ecotype, levels=c("W", "C"))
inv3_filtered$genotype <- factor(inv3_filtered$genotype, levels=c("RA", "AA", "RR"))
inv3_filtered[order(inv3_filtered$genotype, inv3_filtered$sex, inv3_filtered$ecotype),]
options(na.action="na.omit")
inv3_filtered_model <- lmer(log_stat ~ genotype*sex*ecotype + (1 + genotype*sex*ecotype | av), data=inv3_filtered)
options(na.action="na.fail")
inv3_filtered_modelsel <- dredge(inv3_filtered_model, trace=TRUE, rank="AICc", REML=FALSE)
inv3_fmList_filtered <- get.models(inv3_filtered_modelsel, subset=delta<2)
inv3_top <- subset(inv3_filtered_modelsel, delta<2)
inv3_filtered_modelavg <- model.avg(inv3_fmList_filtered, revised.var=TRUE)
inv3_filtered_modelavg_summary <- summary(inv3_filtered_modelavg)
inv3_filtered_estimates <- as.data.frame(inv3_filtered_modelavg_summary$coefmat.full)
#inv3_filtered_estimates <- as.data.frame(coef(summary(inv3_fmList_filtered[[1]])))
inv3_filtered_estimates <- rownames_to_column(inv3_filtered_estimates, "Coefficients")
names(inv3_filtered_estimates) <- gsub(" ", "", names(inv3_filtered_estimates))

# convert estimates from model to estimates for groups
inv3_filtered_groupestimates <- inv3_filtered_estimates %>% select(Coefficients, Estimate) %>% 
  spread(Coefficients, Estimate) %>%
  rename(Intercept=`(Intercept)`) %>%
  mutate(#F_AA_C,
    F_AA_W=Intercept+genotypeAA,
    F_RA_C=Intercept+ecotypeC,
    F_RA_W=Intercept,
    #F_RR_C,
    #F_RR_W=Intercept+genotypeRR,
    M_AA_C=Intercept+sexM+genotypeAA+ecotypeC+`genotypeAA:sexM`+`ecotypeC:sexM`,
    M_AA_W=Intercept+genotypeAA+sexM+`genotypeAA:sexM`,
    M_RA_C=Intercept+sexM+ecotypeC+`ecotypeC:sexM`,
    M_RA_W=Intercept+sexM,
    #M_RR_C=Intercept+genotypeRR+ecotypeC+sexM+`ecotypeC:genotypeRR`+`ecotypeC:sexM`,
    #M_RR_W=Intercept+genotypeRR+sexM
    ) %>%
  gather(Group, Estimate, 7:13) %>%
  mutate(Estimate_exp=exp(Estimate))


##### inversion 4 - model fitting 
## remove F_AA_C, F_RR_W, M_RR_W
inv4_filtered <- inv3_stats_tidy_CW %>% filter(stat_type=="pi") %>%
  mutate(log_stat=log(stat_value)) %>% 
  filter(group!="F_AA_C" & group!="F_RR_W" & group!="M_RR_W") %>%
  separate(group, c("sex", "genotype", "ecotype"), sep="_") %>%
  filter(!is.nan(log_stat) & !is.infinite(log_stat))
inv4_filtered$sex <- factor(inv4_filtered$sex, levels=c("F", "M"))
inv4_filtered$ecotype <- factor(inv4_filtered$ecotype, levels=c("W", "C"))
inv4_filtered$genotype <- factor(inv4_filtered$genotype, levels=c("RA", "AA", "RR"))
inv4_filtered[order(inv4_filtered$genotype, inv4_filtered$sex, inv4_filtered$ecotype),]
options(na.action="na.omit")
inv4_filtered_model <- lmer(log_stat ~ genotype*sex*ecotype + (1 + genotype*sex*ecotype | av), data=inv4_filtered)
options(na.action="na.fail")
inv4_filtered_modelsel <- dredge(inv4_filtered_model, trace=TRUE, rank="AICc", REML=FALSE)
inv4_fmList_filtered <- get.models(inv4_filtered_modelsel, subset=delta<2)
inv4_top <- subset(inv4_filtered_modelsel, delta<2)
inv4_filtered_modelavg <- model.avg(inv4_fmList_filtered, revised.var=TRUE)
inv4_filtered_modelavg_summary <- summary(inv4_filtered_modelavg)
inv4_filtered_estimates <- as.data.frame(inv4_filtered_modelavg_summary$coefmat.full)
#inv4_filtered_estimates <- as.data.frame(coef(summary(inv4_fmList_filtered[[1]])))
inv4_filtered_estimates <- rownames_to_column(inv4_filtered_estimates, "Coefficients")
names(inv4_filtered_estimates) <- gsub(" ", "", names(inv4_filtered_estimates))

# convert estimates from model to estimates for groups
inv4_filtered_groupestimates <- inv4_filtered_estimates %>% select(Coefficients, Estimate) %>% 
  spread(Coefficients, Estimate) %>%
  rename(Intercept=`(Intercept)`) %>%
  mutate(#F_AA_C,
    F_AA_W=Intercept+genotypeAA,
    F_RA_C=Intercept+ecotypeC,
    F_RA_W=Intercept,
    F_RR_C=Intercept+genotypeRR+ecotypeC,
    #F_RR_W,
    M_AA_C=Intercept+sexM+genotypeAA+ecotypeC+`ecotypeC:genotypeAA`,
    M_AA_W=Intercept+genotypeAA+sexM,
    M_RA_C=Intercept+sexM+ecotypeC,
    M_RA_W=Intercept+sexM,
    M_RR_C=Intercept+genotypeRR+ecotypeC+sexM+`ecotypeC:genotypeAA`
    #M_RR_W
  ) %>%
  gather(Group, Estimate, 7:15) %>%
  mutate(Estimate_exp=exp(Estimate))



###########
## export tables of best models and estimates
write.table(inv1_top, file="inv1_models.csv", row.names=F, sep=",")
write.table(inv1_filtered_estimates, file="inv1_estimates.csv", row.names=F, sep=",")
write.table(inv2_top, file="inv2_models.csv", row.names=F, sep=",")
write.table(inv2_filtered_estimates, file="inv2_estimates.csv", row.names=F, sep=",")
write.table(inv3_top, file="inv3_models.csv", row.names=F, sep=",")
write.table(inv3_filtered_estimates, file="inv3_estimates.csv", row.names=F, sep=",")
write.table(inv4_top, file="inv4_models.csv", row.names=F, sep=",")
write.table(inv4_filtered_estimates, file="inv4_estimates.csv", row.names=F, sep=",")


###########
### make the plots for SI - pi and dxy with sex and eco groups

# pi - add column about arr type
supp_pi_plot_data <- inv_stats_tidy_CW %>% filter(stat_type=="pi") %>%
  mutate(arrangement=ifelse(grepl("RR", group), "RR",
                            ifelse(grepl("RA", group), "RA", "AA")))
ggplot(supp_pi_plot_data, aes(x=group, y=stat_value, colour=arrangement)) +
  geom_boxplot() +
  facet_wrap(~inv, nrow=1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) 

# dxy - add column about comp type
supp_dxy_plot_data <- inv_stats_tidy_CW %>% filter(stat_type=="dxy") %>%
  mutate(comp_type=ifelse(stat_type=="dxy" & !grepl("F", group) & !grepl("C", group), "between_arrangement",
                   ifelse(stat_type=="dxy" & !grepl("F", group) & !grepl("W", group), "between_arrangement", 
                   ifelse(stat_type=="dxy" & !grepl("M", group) & !grepl("C", group), "between_arrangement",
                   ifelse(stat_type=="dxy" & !grepl("M", group) & !grepl("W", group), "between_arrangement",
                   ifelse(stat_type=="dxy" & !grepl("F", group) & !grepl("RR", group) & !grepl("RA", group), "between_ecotype",
                   ifelse(stat_type=="dxy" & !grepl("F", group) & !grepl("AA", group) & !grepl("RR", group), "between_ecotype",
                   ifelse(stat_type=="dxy" & !grepl("F", group) & !grepl("AA", group) & !grepl("RA", group), "between_ecotype",
                   ifelse(stat_type=="dxy" & !grepl("M", group) & !grepl("RR", group) & !grepl("RA", group), "between_ecotype",
                   ifelse(stat_type=="dxy" & !grepl("M", group) & !grepl("AA", group) & !grepl("RR", group), "between_ecotype",
                   ifelse(stat_type=="dxy" & !grepl("M", group) & !grepl("AA", group) & !grepl("RA", group), "between_ecotype",
                   ifelse(stat_type=="dxy" & !grepl("C", group) & !grepl("RR", group) & !grepl("RA", group), "between_sex",
                   ifelse(stat_type=="dxy" & !grepl("C", group) & !grepl("AA", group) & !grepl("RR", group), "between_sex",
                   ifelse(stat_type=="dxy" & !grepl("C", group) & !grepl("AA", group) & !grepl("RA", group), "between_sex",
                   ifelse(stat_type=="dxy" & !grepl("W", group) & !grepl("RR", group) & !grepl("RA", group), "between_sex",
                   ifelse(stat_type=="dxy" & !grepl("W", group) & !grepl("AA", group) & !grepl("RR", group), "between_sex",
                   ifelse(stat_type=="dxy" & !grepl("W", group) & !grepl("AA", group) & !grepl("RA", group), "between_sex",
                   ifelse(stat_type=="dxy","other", NA)))))))))))))))))) %>%
  filter(comp_type!="other")
ggplot(supp_dxy_plot_data, aes(x=group, y=stat_value, colour=comp_type)) +
  geom_boxplot() +
  geom_abline(aes(intercept=dxymean$mean, slope=0)) +
  facet_wrap(~inv, nrow=2) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) 
#added line for overall mean dxy per contig across all groups, just out of interest
dxymean <- supp_dxy_plot_data %>% summarise(mean=mean(stat_value, na.rm=T))

View(inv4_filtered_groupestimates)
