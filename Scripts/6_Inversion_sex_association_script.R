rm(list=ls())
# set working directory
library(tidyverse)
library(zoo)
library(tidyquant)
library(reshape)

########
###
### Associations between inversions (and sex)
### 09/03/21
###
########


## need table of inversion genotypes for each inversion, plus sex and ecotype info
## use inversion genotypes from PC clusters and convert to Rs and As and export new inversion genotype files

# read in data
inv1_female <- read.csv("region1_females_PCandcluster.csv")
inv1_male <- read.csv("region1_males_PCandcluster.csv")
inv2_female <- read.csv("region2a_females_PCandcluster.csv")
inv2_male <- read.csv("region2a_males_PCandcluster.csv")
inv3_female <- read.csv("region2b_females_PCandcluster.csv")
inv3_male <- read.csv("region2b_males_PCandcluster.csv")
inv4_female <- read.csv("region3_females_PCandcluster.csv")
inv4_male <- read.csv("region3_males_PCandcluster.csv")
ANG_snail_info <- read.csv("ANG_snail_info.csv")
map <- read.csv("LG12_mappos_bycontig.csv")

# select relevent columns and recode clusters to inversion genotypes (AA/RA/RR)
inv1_female <- inv1_female %>% select(snail_ID, DistAlongPath, cluster) %>% dplyr::rename("genotype"="cluster")
inv1_female$genotype <- recode(inv1_female$genotype, `1`="RR", `2`="RA", `3`="AA")
inv1_male <- inv1_male %>% select(snail_ID, DistAlongPath, cluster) %>% dplyr::rename("genotype"="cluster")
inv1_male$genotype <- recode(inv1_male$genotype, `1`="RR", `2`="RA", `3`="AA")

inv2_female <- inv2_female %>% select(snail_ID, DistAlongPath, cluster) %>% dplyr::rename("genotype"="cluster")
inv2_female$genotype <- recode(inv2_female$genotype, `1`="RR", `2`="RA", `3`="AA")
inv2_male <- inv2_male %>% select(snail_ID, DistAlongPath, cluster) %>% dplyr::rename("genotype"="cluster")
inv2_male$genotype <- recode(inv2_male$genotype, `1`="RR", `2`="RA", `3`="AA")

inv3_female <- inv3_female %>% select(snail_ID, DistAlongPath, cluster) %>% dplyr::rename("genotype"="cluster")
inv3_female$genotype <- recode(inv3_female$genotype, `1`="RR", `2`="RA", `3`="AA")
inv3_male <- inv3_male %>% select(snail_ID, DistAlongPath, cluster) %>% dplyr::rename("genotype"="cluster")
inv3_male$genotype <- recode(inv3_male$genotype, `1`="RR", `2`="RA", `3`="AA")

inv4_female <- inv4_female %>% select(snail_ID, DistAlongPath, cluster) %>% dplyr::rename("genotype"="cluster")
inv4_female$genotype <- recode(inv4_female$genotype, `3`="RR", `2`="RA", `1`="AA")
inv4_male <- inv4_male %>% select(snail_ID, DistAlongPath, cluster) %>% dplyr::rename("genotype"="cluster")
inv4_male$genotype <- recode(inv4_male$genotype, `3`="RR", `2`="RA", `1`="AA")

# make and export data frame of inversion genotypes
inv1_female_genos <- inv1_female %>% select(-DistAlongPath) %>% dplyr::rename(inv1=genotype)
inv1_male_genos <- inv1_male %>% select(-DistAlongPath) %>% dplyr::rename(inv1=genotype)
inv2_female_genos <- inv2_female %>% select(-DistAlongPath) %>% dplyr::rename(inv2=genotype)
inv2_male_genos <- inv2_male %>% select(-DistAlongPath) %>% dplyr::rename(inv2=genotype)
inv3_female_genos <- inv3_female %>% select(-DistAlongPath) %>% dplyr::rename(inv3=genotype)
inv3_male_genos <- inv3_male %>% select(-DistAlongPath) %>% dplyr::rename(inv3=genotype)
inv4_female_genos <- inv4_female %>% select(-DistAlongPath) %>% dplyr::rename(inv4=genotype)
inv4_male_genos <- inv4_male %>% select(-DistAlongPath) %>% dplyr::rename(inv4=genotype)

female_genos <- left_join(inv1_female_genos, inv2_female_genos) %>% left_join(., inv3_female_genos) %>%
  left_join(., inv4_female_genos)
male_genos <- left_join(inv1_male_genos, inv2_male_genos) %>% left_join(., inv3_male_genos) %>%
  left_join(., inv4_male_genos)

#write.csv(female_genos, file="ANG_females_inversion_genotypes.csv", row.names = F)
#write.csv(male_genos, file="ANG_males_inversion_genotypes.csv", row.names = F)

## add snail info to the merged data frames
females_info <- read.csv("region1_females_PCandcluster.csv") %>% 
  select(snail_ID, sex_dissection, ecotype.y) %>%
  dplyr::rename(sex=sex_dissection, ecotype=ecotype.y)
males_info <- read.csv("region1_males_PCandcluster.csv") %>% 
  select(snail_ID, sex_dissection, ecotype.y) %>%
  dplyr::rename(sex=sex_dissection, ecotype=ecotype.y)

female_genos <- female_genos %>% left_join(., females_info, by="snail_ID")
male_genos <- male_genos %>% left_join(., males_info, by="snail_ID")

## split into ecotypes
crab_female_genos <- female_genos %>% filter(ecotype=="Crab")
wave_female_genos <- female_genos %>% filter(ecotype=="Wave")
crab_male_genos <- male_genos %>% filter(ecotype=="Crab")
wave_male_genos <- male_genos %>% filter(ecotype=="Wave")

## do associations between inversions analysis separately for these four groups
## i want to know whether the genotypes at each inversion are linked, so...
##    contingency table should be just inversions vs genotypes, with counts as values in the table (??)
## it should be a 3x3 table- rows are genos for one inv and columns are genos for another inv

## so, calculate counts of each genotype for each inversion
crab_female_counts <- crab_female_genos %>% select(-ecotype) %>% gather(inversion, genotype, 2:5) 
crab_female_counts$inversion <- as.factor(crab_female_counts$inversion)
crab_female_counts$genotype <- as.factor(crab_female_counts$genotype)
crab_female_counts <- crab_female_counts %>% group_by(inversion, genotype, .drop=FALSE) %>% count()

crab_male_counts <- crab_male_genos %>% select(-ecotype) %>% gather(inversion, genotype, 2:5) 
crab_male_counts$inversion <- as.factor(crab_male_counts$inversion)
crab_male_counts$genotype <- as.factor(crab_male_counts$genotype)
crab_male_counts <- crab_male_counts %>% group_by(inversion, genotype, .drop=FALSE) %>% count()


wave_female_counts <- wave_female_genos %>% select(-ecotype) %>% gather(inversion, genotype, 2:5) %>%
  group_by(inversion, genotype) %>% count()

wave_male_counts <- wave_male_genos %>% select(-ecotype) %>% gather(inversion, genotype, 2:5) 
wave_male_counts$inversion <- as.factor(wave_male_counts$inversion)
wave_male_counts$genotype <- as.factor(wave_male_counts$genotype)
wave_male_counts <- wave_male_counts %>% group_by(inversion, genotype, .drop=FALSE) %>% count()

## spread into contingency table 
crab_female_table <- crab_female_counts %>% filter(genotype=="RA" | genotype=="RR" | genotype=="AA") %>%
  spread(genotype, n) 
crab_male_table <- crab_male_counts %>% filter(genotype=="RA" | genotype=="RR" | genotype=="AA") %>%
  spread(genotype, n) 
wave_female_table <- wave_female_counts %>% filter(genotype=="RA" | genotype=="RR" | genotype=="AA") %>%
  spread(genotype, n) 
wave_male_table <- wave_male_counts %>% filter(genotype=="RA" | genotype=="RR" | genotype=="AA") %>%
  spread(genotype, n) 

## add 1 to all values (to avoid the '0' issue in the contingency tables)
crab_female_table <- crab_female_table %>% mutate(AA=AA+1, RA=RA+1, RR=RR+1)
crab_male_table <- crab_male_table %>% mutate(AA=AA+1, RA=RA+1, RR=RR+1)
wave_female_table <- wave_female_table %>% mutate(AA=AA+1, RA=RA+1, RR=RR+1)
wave_male_table <- wave_male_table %>% mutate(AA=AA+1, RA=RA+1, RR=RR+1)

## do chisquare test
# this is overall association of genotype across inversions - one p val
crab_female_matrix <- crab_female_table %>% column_to_rownames("inversion") %>% as.matrix(.)
crab_female_chi <- chisq.test(crab_female_matrix)

## do pairwise comparisons to get a p val for each pair of inversions

#  crab female
crab_female_1vs2 <- crab_female_genos %>% select(-sex, -ecotype, -inv3, -inv4) %>%
  group_by(inv1, inv2, .drop=FALSE) %>% count() %>% spread(inv2, n) %>%
  ungroup %>%
  add_row(inv1="RA", AA=0, RA=0, RR=0) %>%
  add_row(inv1="AA", AA=0, RA=0, RR=0) %>% column_to_rownames("inv1") 
crab_female_1vs2 <- crab_female_1vs2 %>% mutate(AA=AA+1, RA=RA+1, RR=RR+1) %>%
  as.matrix(.)
crab_female_1vs2_chi <- chisq.test(crab_female_1vs2)
crab_female_1vs2_pval <- crab_female_1vs2_chi$p.value

crab_female_1vs3 <- crab_female_genos %>% select(-sex, -ecotype, -inv2, -inv4) %>%
  group_by(inv1, inv3, .drop=FALSE) %>% count() %>% spread(inv3, n) %>%
  ungroup %>%
  add_row(inv1="RA", AA=0, RA=0, RR=0) %>%
  add_row(inv1="AA", AA=0, RA=0, RR=0) %>% column_to_rownames("inv1") 
crab_female_1vs3 <- crab_female_1vs3 %>% mutate(AA=AA+1, RA=RA+1, RR=RR+1) %>%
  as.matrix(.)
crab_female_1vs3_chi <- chisq.test(crab_female_1vs3)
crab_female_1vs3_pval <- crab_female_1vs3_chi$p.value

crab_female_1vs4 <- crab_female_genos %>% select(-sex, -ecotype, -inv2, -inv3) %>%
  group_by(inv1, inv4, .drop=FALSE) %>% count() %>% spread(inv4, n) %>%
  ungroup %>%
  add_row(inv1="RA", AA=0, RA=0, RR=0) %>%
  add_row(inv1="AA", AA=0, RA=0, RR=0) %>% column_to_rownames("inv1") %>% select(-`<NA>`)
crab_female_1vs4 <- crab_female_1vs4 %>% mutate(AA=AA+1, RA=RA+1, RR=RR+1) %>%
  as.matrix(.)
crab_female_1vs4_chi <- chisq.test(crab_female_1vs4)
crab_female_1vs4_pval <- crab_female_1vs4_chi$p.value

crab_female_2vs3 <- crab_female_genos %>% select(-sex, -ecotype, -inv1, -inv4) %>%
  group_by(inv2, inv3, .drop=FALSE) %>% count() %>% spread(inv3, n) %>%
  ungroup %>% column_to_rownames("inv2") 
crab_female_2vs3[is.na(crab_female_2vs3)] <- 0 
crab_female_2vs3 <- crab_female_2vs3 %>% mutate(AA=AA+1, RA=RA+1, RR=RR+1) %>%
  as.matrix(.)
crab_female_2vs3_chi <- chisq.test(crab_female_2vs3)
crab_female_2vs3_pval <- crab_female_2vs3_chi$p.value

crab_female_2vs4 <- crab_female_genos %>% select(-sex, -ecotype, -inv1, -inv3) %>%
  group_by(inv2, inv4, .drop=FALSE) %>% count() %>% spread(inv4, n) %>%
  select(-`<NA>`) %>%
  ungroup %>% column_to_rownames("inv2") 
crab_female_2vs4[is.na(crab_female_2vs4)] <- 0 
crab_female_2vs4 <- crab_female_2vs4 %>% mutate(AA=AA+1, RA=RA+1, RR=RR+1) %>%
  as.matrix(.)
crab_female_2vs4_chi <- chisq.test(crab_female_2vs4)
crab_female_2vs4_pval <- crab_female_2vs4_chi$p.value

crab_female_3vs4 <- crab_female_genos %>% select(-sex, -ecotype, -inv1, -inv2) %>%
  group_by(inv3, inv4, .drop=FALSE) %>% count() %>% spread(inv4, n) %>%
  select(-`<NA>`) %>%
  ungroup %>% column_to_rownames("inv3") 
crab_female_3vs4[is.na(crab_female_3vs4)] <- 0 
crab_female_3vs4 <- crab_female_3vs4 %>% mutate(AA=AA+1, RA=RA+1, RR=RR+1) %>%
  as.matrix(.)
crab_female_3vs4_chi <- chisq.test(crab_female_3vs4)
crab_female_3vs4_pval <- crab_female_3vs4_chi$p.value

crab_female_1vs1 <- crab_female_genos %>% select(-sex, -ecotype, -inv2, -inv3, -inv4) %>%
  mutate(inv1a=inv1) %>%
  group_by(inv1, inv1a, .drop=FALSE) %>% count() %>% spread(inv1a, n) %>%
  ungroup %>%
  add_row(inv1="RA", RR=0) %>%
  add_row(inv1="AA", RR=0) %>% 
  mutate(RA=0, AA=0) %>% column_to_rownames("inv1") 
crab_female_1vs1 <- crab_female_1vs1 %>% mutate(AA=AA+1, RA=RA+1, RR=RR+1) %>%
  as.matrix(.)
crab_female_1vs1_chi <- chisq.test(crab_female_1vs1)
crab_female_1vs1_pval <- crab_female_1vs1_chi$p.value

crab_female_2vs2 <- crab_female_genos %>% select(-sex, -ecotype, -inv1, -inv3, -inv4) %>%
  mutate(inv2a=inv2) %>%
  group_by(inv2, inv2a, .drop=FALSE) %>% count() %>% spread(inv2a, n) %>%
  ungroup %>% column_to_rownames("inv2") 
crab_female_2vs2[is.na(crab_female_2vs2)] <- 0 
crab_female_2vs2 <- crab_female_2vs2 %>% mutate(AA=AA+1, RA=RA+1, RR=RR+1) %>%
  as.matrix(.)
crab_female_2vs2_chi <- chisq.test(crab_female_2vs2)
crab_female_2vs2_pval <- crab_female_2vs2_chi$p.value

crab_female_3vs3 <- crab_female_genos %>% select(-sex, -ecotype, -inv1, -inv2, -inv4) %>%
  mutate(inv3a=inv3) %>%
  group_by(inv3, inv3a, .drop=FALSE) %>% count() %>% spread(inv3a, n) %>%
  ungroup %>% column_to_rownames("inv3") 
crab_female_3vs3[is.na(crab_female_3vs3)] <- 0 
crab_female_3vs3 <- crab_female_3vs3 %>% mutate(AA=AA+1, RA=RA+1, RR=RR+1) %>%
  as.matrix(.)
crab_female_3vs3_chi <- chisq.test(crab_female_3vs3)
crab_female_3vs3_pval <- crab_female_3vs3_chi$p.value

crab_female_4vs4 <- crab_female_genos %>% select(-sex, -ecotype, -inv1, -inv2, -inv3) %>%
  filter(!is.na(inv4)) %>%
  mutate(inv4a=inv4) %>%
  group_by(inv4, inv4a, .drop=FALSE) %>% count() %>% spread(inv4a, n) %>%
  ungroup %>% column_to_rownames("inv4") 
crab_female_4vs4[is.na(crab_female_4vs4)] <- 0 
crab_female_4vs4 <- crab_female_4vs4 %>% mutate(AA=AA+1, RA=RA+1, RR=RR+1) %>%
  as.matrix(.)
crab_female_4vs4_chi <- chisq.test(crab_female_4vs4)
crab_female_4vs4_pval <- crab_female_4vs4_chi$p.value



#  crab male
crab_male_1vs2 <- crab_male_genos %>% select(-sex, -ecotype, -inv3, -inv4) %>%
  group_by(inv1, inv2, .drop=FALSE) %>% count() %>% spread(inv2, n) %>%
  ungroup %>%
  add_row(inv1="RA", RA=0, RR=0) %>%
  add_row(inv1="AA", RA=0, RR=0) %>% 
  mutate(AA=0) %>% column_to_rownames("inv1") 
crab_male_1vs2 <- crab_male_1vs2 %>% mutate(AA=AA+1, RA=RA+1, RR=RR+1) %>%
  as.matrix(.)
crab_male_1vs2_chi <- chisq.test(crab_male_1vs2)
crab_male_1vs2_pval <- crab_male_1vs2_chi$p.value

crab_male_1vs3 <- crab_male_genos %>% select(-sex, -ecotype, -inv2, -inv4) %>%
  group_by(inv1, inv3, .drop=FALSE) %>% count() %>% spread(inv3, n) %>%
  ungroup %>%
  add_row(inv1="RA", AA=0, RA=0) %>%
  add_row(inv1="AA", AA=0, RA=0) %>% 
  mutate(RR=0) %>%column_to_rownames("inv1") 
crab_male_1vs3 <- crab_male_1vs3 %>% mutate(AA=AA+1, RA=RA+1, RR=RR+1) %>%
  as.matrix(.)
crab_male_1vs3_chi <- chisq.test(crab_male_1vs3)
crab_male_1vs3_pval <- crab_male_1vs3_chi$p.value

crab_male_1vs4 <- crab_male_genos %>% select(-sex, -ecotype, -inv2, -inv3) %>%
  group_by(inv1, inv4, .drop=FALSE) %>% count() %>% spread(inv4, n) %>%
  ungroup %>%
  add_row(inv1="RA", AA=0, RA=0, RR=0) %>%
  add_row(inv1="AA", AA=0, RA=0, RR=0) %>% column_to_rownames("inv1")
crab_male_1vs4 <- crab_male_1vs4 %>% mutate(AA=AA+1, RA=RA+1, RR=RR+1) %>%
  as.matrix(.)
crab_male_1vs4_chi <- chisq.test(crab_male_1vs4)
crab_male_1vs4_pval <- crab_male_1vs4_chi$p.value

crab_male_2vs3 <- crab_male_genos %>% select(-sex, -ecotype, -inv1, -inv4) %>%
  group_by(inv2, inv3, .drop=FALSE) %>% count() %>% spread(inv3, n) %>%
  ungroup %>% 
  add_row(inv2="AA", AA=0, RA=0) %>%
  mutate(RR=0) %>%
  column_to_rownames("inv2") 
crab_male_2vs3[is.na(crab_male_2vs3)] <- 0 
crab_male_2vs3 <- crab_male_2vs3 %>% mutate(AA=AA+1, RA=RA+1, RR=RR+1) %>%
  as.matrix(.)
crab_male_2vs3_chi <- chisq.test(crab_male_2vs3)
crab_male_2vs3_pval <- crab_male_2vs3_chi$p.value

crab_male_2vs4 <- crab_male_genos %>% select(-sex, -ecotype, -inv1, -inv3) %>%
  group_by(inv2, inv4, .drop=FALSE) %>% count() %>% spread(inv4, n) %>%
  ungroup %>%
  add_row(inv2="AA", AA=0, RA=0, RR=0) %>%
  column_to_rownames("inv2") 
crab_male_2vs4[is.na(crab_male_2vs4)] <- 0 
crab_male_2vs4 <- crab_male_2vs4 %>% mutate(AA=AA+1, RA=RA+1, RR=RR+1) %>%
  as.matrix(.)
crab_male_2vs4_chi <- chisq.test(crab_male_2vs4)
crab_male_2vs4_pval <- crab_male_2vs4_chi$p.value

crab_male_3vs4 <- crab_male_genos %>% select(-sex, -ecotype, -inv1, -inv2) %>%
  group_by(inv3, inv4, .drop=FALSE) %>% count() %>% spread(inv4, n) %>%
  ungroup %>%
  add_row(inv3="RR", AA=0, RA=0, RR=0) %>% column_to_rownames("inv3") 
crab_male_3vs4[is.na(crab_male_3vs4)] <- 0 
crab_male_3vs4 <- crab_male_3vs4 %>% mutate(AA=AA+1, RA=RA+1, RR=RR+1) %>%
  as.matrix(.)
crab_male_3vs4_chi <- chisq.test(crab_male_3vs4)
crab_male_3vs4_pval <- crab_male_3vs4_chi$p.value

crab_male_1vs1 <- crab_male_genos %>% select(-sex, -ecotype, -inv2, -inv3, -inv4) %>%
  mutate(inv1a=inv1) %>%
  group_by(inv1, inv1a, .drop=FALSE) %>% count() %>% spread(inv1a, n) %>%
  ungroup %>%
  add_row(inv1="RA", RR=0) %>%
  add_row(inv1="AA", RR=0) %>% 
  mutate(RA=0, AA=0) %>% column_to_rownames("inv1") 
crab_male_1vs1 <- crab_male_1vs1 %>% mutate(AA=AA+1, RA=RA+1, RR=RR+1) %>%
  as.matrix(.)
crab_male_1vs1_chi <- chisq.test(crab_male_1vs1)
crab_male_1vs1_pval <- crab_male_1vs1_chi$p.value

crab_male_2vs2 <- crab_male_genos %>% select(-sex, -ecotype, -inv1, -inv3, -inv4) %>%
  mutate(inv2a=inv2) %>%
  group_by(inv2, inv2a, .drop=FALSE) %>% count() %>% spread(inv2a, n) %>%
  ungroup %>% 
  add_row(inv2="AA", RA=0, RR=0) %>% 
  mutate(AA=0) %>% column_to_rownames("inv2") 
crab_male_2vs2[is.na(crab_male_2vs2)] <- 0 
crab_male_2vs2 <- crab_male_2vs2 %>% mutate(AA=AA+1, RA=RA+1, RR=RR+1) %>%
  as.matrix(.)
crab_male_2vs2_chi <- chisq.test(crab_male_2vs2)
crab_male_2vs2_pval <- crab_male_2vs2_chi$p.value

crab_male_3vs3 <- crab_male_genos %>% select(-sex, -ecotype, -inv1, -inv2, -inv4) %>%
  mutate(inv3a=inv3) %>%
  group_by(inv3, inv3a, .drop=FALSE) %>% count() %>% spread(inv3a, n) %>%
  ungroup %>%
  add_row(inv3="RR", AA=0, RA=0) %>% 
  mutate(RR=0) %>%
  column_to_rownames("inv3") 
crab_male_3vs3[is.na(crab_male_3vs3)] <- 0 
crab_male_3vs3 <- crab_male_3vs3 %>% mutate(AA=AA+1, RA=RA+1, RR=RR+1) %>%
  as.matrix(.)
crab_male_3vs3_chi <- chisq.test(crab_male_3vs3)
crab_male_3vs3_pval <- crab_male_3vs3_chi$p.value

crab_male_4vs4 <- crab_male_genos %>% select(-sex, -ecotype, -inv1, -inv2, -inv3) %>%
  filter(!is.na(inv4)) %>%
  mutate(inv4a=inv4) %>%
  group_by(inv4, inv4a, .drop=FALSE) %>% count() %>% spread(inv4a, n) %>%
  ungroup %>% column_to_rownames("inv4") 
crab_male_4vs4[is.na(crab_male_4vs4)] <- 0 
crab_male_4vs4 <- crab_male_4vs4 %>% mutate(AA=AA+1, RA=RA+1, RR=RR+1) %>%
  as.matrix(.)
crab_male_4vs4_chi <- chisq.test(crab_male_4vs4)
crab_male_4vs4_pval <- crab_male_4vs4_chi$p.value


#  wave female
wave_female_1vs2 <- wave_female_genos %>% select(-sex, -ecotype, -inv3, -inv4) %>%
  group_by(inv1, inv2, .drop=FALSE) %>% count() %>% spread(inv2, n) %>%
  ungroup %>% column_to_rownames("inv1") 
wave_female_1vs2[is.na(wave_female_1vs2)] <- 0
wave_female_1vs2 <- wave_female_1vs2 %>% mutate(AA=AA+1, RA=RA+1, RR=RR+1) %>%
  as.matrix(.)
wave_female_1vs2_chi <- chisq.test(wave_female_1vs2)
wave_female_1vs2_pval <- wave_female_1vs2_chi$p.value

wave_female_1vs3 <- wave_female_genos %>% select(-sex, -ecotype, -inv2, -inv4) %>%
  group_by(inv1, inv3, .drop=FALSE) %>% count() %>% spread(inv3, n) %>%
  ungroup %>%
  column_to_rownames("inv1") 
wave_female_1vs3[is.na(wave_female_1vs3)] <- 0
wave_female_1vs3 <- wave_female_1vs3 %>% mutate(AA=AA+1, RA=RA+1, RR=RR+1) %>%
  as.matrix(.)
wave_female_1vs3_chi <- chisq.test(wave_female_1vs3)
wave_female_1vs3_pval <- wave_female_1vs3_chi$p.value

wave_female_1vs4 <- wave_female_genos %>% select(-sex, -ecotype, -inv2, -inv3) %>%
  group_by(inv1, inv4, .drop=FALSE) %>% count() %>% spread(inv4, n) %>%
  ungroup %>%
  column_to_rownames("inv1")
wave_female_1vs4[is.na(wave_female_1vs4)] <- 0
wave_female_1vs4 <- wave_female_1vs4 %>% mutate(AA=AA+1, RA=RA+1, RR=RR+1) %>%
  as.matrix(.)
wave_female_1vs4_chi <- chisq.test(wave_female_1vs4)
wave_female_1vs4_pval <- wave_female_1vs4_chi$p.value

wave_female_2vs3 <- wave_female_genos %>% select(-sex, -ecotype, -inv1, -inv4) %>%
  group_by(inv2, inv3, .drop=FALSE) %>% count() %>% spread(inv3, n) %>%
  ungroup %>% column_to_rownames("inv2") 
wave_female_2vs3[is.na(wave_female_2vs3)] <- 0 
wave_female_2vs3 <- wave_female_2vs3 %>% mutate(AA=AA+1, RA=RA+1, RR=RR+1) %>%
  as.matrix(.)
wave_female_2vs3_chi <- chisq.test(wave_female_2vs3)
wave_female_2vs3_pval <- wave_female_2vs3_chi$p.value

wave_female_2vs4 <- wave_female_genos %>% select(-sex, -ecotype, -inv1, -inv3) %>%
  group_by(inv2, inv4, .drop=FALSE) %>% count() %>% spread(inv4, n) %>%
  ungroup %>% column_to_rownames("inv2") 
wave_female_2vs4[is.na(wave_female_2vs4)] <- 0 
wave_female_2vs4 <- wave_female_2vs4 %>% mutate(AA=AA+1, RA=RA+1, RR=RR+1) %>%
  as.matrix(.)
wave_female_2vs4_chi <- chisq.test(wave_female_2vs4)
wave_female_2vs4_pval <- wave_female_2vs4_chi$p.value

wave_female_3vs4 <- wave_female_genos %>% select(-sex, -ecotype, -inv1, -inv2) %>%
  group_by(inv3, inv4, .drop=FALSE) %>% count() %>% spread(inv4, n) %>%
  ungroup %>% column_to_rownames("inv3") 
wave_female_3vs4[is.na(wave_female_3vs4)] <- 0 
wave_female_3vs4 <- wave_female_3vs4 %>% mutate(AA=AA+1, RA=RA+1, RR=RR+1) %>%
  as.matrix(.)
wave_female_3vs4_chi <- chisq.test(wave_female_3vs4)
wave_female_3vs4_pval <- wave_female_3vs4_chi$p.value

wave_female_1vs1 <- wave_female_genos %>% select(-sex, -ecotype, -inv2, -inv3, -inv4) %>%
  mutate(inv1a=inv1) %>%
  group_by(inv1, inv1a, .drop=FALSE) %>% count() %>% spread(inv1a, n) %>%
  ungroup %>% column_to_rownames("inv1") 
wave_female_1vs1[is.na(wave_female_1vs1)] <- 0 
wave_female_1vs1 <- wave_female_1vs1 %>% mutate(AA=AA+1, RA=RA+1, RR=RR+1) %>%
  as.matrix(.)
wave_female_1vs1_chi <- chisq.test(wave_female_1vs1)
wave_female_1vs1_pval <- wave_female_1vs1_chi$p.value

wave_female_2vs2 <- wave_female_genos %>% select(-sex, -ecotype, -inv1, -inv3, -inv4) %>%
  mutate(inv2a=inv2) %>%
  group_by(inv2, inv2a, .drop=FALSE) %>% count() %>% spread(inv2a, n) %>%
  ungroup %>% column_to_rownames("inv2") 
wave_female_2vs2[is.na(wave_female_2vs2)] <- 0 
wave_female_2vs2 <- wave_female_2vs2 %>% mutate(AA=AA+1, RA=RA+1, RR=RR+1) %>%
  as.matrix(.)
wave_female_2vs2_chi <- chisq.test(wave_female_2vs2)
wave_female_2vs2_pval <- wave_female_2vs2_chi$p.value

wave_female_3vs3 <- wave_female_genos %>% select(-sex, -ecotype, -inv1, -inv2, -inv4) %>%
  mutate(inv3a=inv3) %>%
  group_by(inv3, inv3a, .drop=FALSE) %>% count() %>% spread(inv3a, n) %>%
  ungroup %>% column_to_rownames("inv3") 
wave_female_3vs3[is.na(wave_female_3vs3)] <- 0 
wave_female_3vs3 <- wave_female_3vs3 %>% mutate(AA=AA+1, RA=RA+1, RR=RR+1) %>%
  as.matrix(.)
wave_female_3vs3_chi <- chisq.test(wave_female_3vs3)
wave_female_3vs3_pval <- wave_female_3vs3_chi$p.value

wave_female_4vs4 <- wave_female_genos %>% select(-sex, -ecotype, -inv1, -inv2, -inv3) %>%
  filter(!is.na(inv4)) %>%
  mutate(inv4a=inv4) %>%
  group_by(inv4, inv4a, .drop=FALSE) %>% count() %>% spread(inv4a, n) %>%
  ungroup %>% column_to_rownames("inv4") 
wave_female_4vs4[is.na(wave_female_4vs4)] <- 0 
wave_female_4vs4 <- wave_female_4vs4 %>% mutate(AA=AA+1, RA=RA+1, RR=RR+1) %>%
  as.matrix(.)
wave_female_4vs4_chi <- chisq.test(wave_female_4vs4)
wave_female_4vs4_pval <- wave_female_4vs4_chi$p.value


#  wave male
wave_male_1vs2 <- wave_male_genos %>% select(-sex, -ecotype, -inv3, -inv4) %>%
  group_by(inv1, inv2, .drop=FALSE) %>% count() %>% spread(inv2, n) %>%
  ungroup %>% column_to_rownames("inv1") 
wave_male_1vs2[is.na(wave_male_1vs2)] <- 0
wave_male_1vs2 <- wave_male_1vs2 %>% mutate(AA=AA+1, RA=RA+1, RR=RR+1) %>%
  as.matrix(.)
wave_male_1vs2_chi <- chisq.test(wave_male_1vs2)
wave_male_1vs2_pval <- wave_male_1vs2_chi$p.value

wave_male_1vs3 <- wave_male_genos %>% select(-sex, -ecotype, -inv2, -inv4) %>%
  group_by(inv1, inv3, .drop=FALSE) %>% count() %>% spread(inv3, n) %>%
  ungroup %>%
  mutate(RR=0) %>%
  column_to_rownames("inv1") 
wave_male_1vs3[is.na(wave_male_1vs3)] <- 0
wave_male_1vs3 <- wave_male_1vs3 %>% mutate(AA=AA+1, RA=RA+1, RR=RR+1) %>%
  as.matrix(.)
wave_male_1vs3_chi <- chisq.test(wave_male_1vs3)
wave_male_1vs3_pval <- wave_male_1vs3_chi$p.value

wave_male_1vs4 <- wave_male_genos %>% select(-sex, -ecotype, -inv2, -inv3) %>%
  group_by(inv1, inv4, .drop=FALSE) %>% count() %>% spread(inv4, n) %>%
  ungroup %>%
  mutate(RR=0) %>%
  column_to_rownames("inv1")
wave_male_1vs4[is.na(wave_male_1vs4)] <- 0
wave_male_1vs4 <- wave_male_1vs4 %>% mutate(AA=AA+1, RA=RA+1, RR=RR+1) %>%
  as.matrix(.)
wave_male_1vs4_chi <- chisq.test(wave_male_1vs4)
wave_male_1vs4_pval <- wave_male_1vs4_chi$p.value

wave_male_2vs3 <- wave_male_genos %>% select(-sex, -ecotype, -inv1, -inv4) %>%
  group_by(inv2, inv3, .drop=FALSE) %>% count() %>% spread(inv3, n) %>%
  ungroup %>%
  mutate(RR=0) %>% column_to_rownames("inv2") 
wave_male_2vs3[is.na(wave_male_2vs3)] <- 0 
wave_male_2vs3 <- wave_male_2vs3 %>% mutate(AA=AA+1, RA=RA+1, RR=RR+1) %>%
  as.matrix(.)
wave_male_2vs3_chi <- chisq.test(wave_male_2vs3)
wave_male_2vs3_pval <- wave_male_2vs3_chi$p.value

wave_male_2vs4 <- wave_male_genos %>% select(-sex, -ecotype, -inv1, -inv3) %>%
  group_by(inv2, inv4, .drop=FALSE) %>% count() %>% spread(inv4, n) %>%
  ungroup %>%
  mutate(RR=0) %>% column_to_rownames("inv2") 
wave_male_2vs4[is.na(wave_male_2vs4)] <- 0 
wave_male_2vs4 <- wave_male_2vs4 %>% mutate(AA=AA+1, RA=RA+1, RR=RR+1) %>%
  as.matrix(.)
wave_male_2vs4_chi <- chisq.test(wave_male_2vs4)
wave_male_2vs4_pval <- wave_male_2vs4_chi$p.value

wave_male_3vs4 <- wave_male_genos %>% select(-sex, -ecotype, -inv1, -inv2) %>%
  group_by(inv3, inv4, .drop=FALSE) %>% count() %>% spread(inv4, n) %>%
  ungroup %>% 
  add_row(inv3="RR", AA=0, RA=0) %>%
  mutate(RR=0) %>% column_to_rownames("inv3") 
wave_male_3vs4[is.na(wave_male_3vs4)] <- 0 
wave_male_3vs4 <- wave_male_3vs4 %>% mutate(AA=AA+1, RA=RA+1, RR=RR+1) %>%
  as.matrix(.)
wave_male_3vs4_chi <- chisq.test(wave_male_3vs4)
wave_male_3vs4_pval <- wave_male_3vs4_chi$p.value

wave_male_1vs1 <- wave_male_genos %>% select(-sex, -ecotype, -inv2, -inv3, -inv4) %>%
  mutate(inv1a=inv1) %>%
  group_by(inv1, inv1a, .drop=FALSE) %>% count() %>% spread(inv1a, n) %>%
  ungroup %>% column_to_rownames("inv1") 
wave_male_1vs1[is.na(wave_male_1vs1)] <- 0 
wave_male_1vs1 <- wave_male_1vs1 %>% mutate(AA=AA+1, RA=RA+1, RR=RR+1) %>%
  as.matrix(.)
wave_male_1vs1_chi <- chisq.test(wave_male_1vs1)
wave_male_1vs1_pval <- wave_male_1vs1_chi$p.value

wave_male_2vs2 <- wave_male_genos %>% select(-sex, -ecotype, -inv1, -inv3, -inv4) %>%
  mutate(inv2a=inv2) %>%
  group_by(inv2, inv2a, .drop=FALSE) %>% count() %>% spread(inv2a, n) %>%
  ungroup %>% column_to_rownames("inv2") 
wave_male_2vs2[is.na(wave_male_2vs2)] <- 0 
wave_male_2vs2 <- wave_male_2vs2 %>% mutate(AA=AA+1, RA=RA+1, RR=RR+1) %>%
  as.matrix(.)
wave_male_2vs2_chi <- chisq.test(wave_male_2vs2)
wave_male_2vs2_pval <- wave_male_2vs2_chi$p.value

wave_male_3vs3 <- wave_male_genos %>% select(-sex, -ecotype, -inv1, -inv2, -inv4) %>%
  mutate(inv3a=inv3) %>%
  group_by(inv3, inv3a, .drop=FALSE) %>% count() %>% spread(inv3a, n) %>%
  ungroup %>% 
  add_row(inv3="RR", AA=0, RA=0) %>%
  mutate(RR=0) %>% column_to_rownames("inv3") 
wave_male_3vs3[is.na(wave_male_3vs3)] <- 0 
wave_male_3vs3 <- wave_male_3vs3 %>% mutate(AA=AA+1, RA=RA+1, RR=RR+1) %>%
  as.matrix(.)
wave_male_3vs3_chi <- chisq.test(wave_male_3vs3)
wave_male_3vs3_pval <- wave_male_3vs3_chi$p.value

wave_male_4vs4 <- wave_male_genos %>% select(-sex, -ecotype, -inv1, -inv2, -inv3) %>%
  filter(!is.na(inv4)) %>%
  mutate(inv4a=inv4) %>%
  group_by(inv4, inv4a, .drop=FALSE) %>% count() %>% spread(inv4a, n) %>%
  ungroup %>% 
  add_row(inv4="RR", AA=0, RA=0) %>%
  mutate(RR=0) %>% column_to_rownames("inv4") 
wave_male_4vs4[is.na(wave_male_4vs4)] <- 0 
wave_male_4vs4 <- wave_male_4vs4 %>% mutate(AA=AA+1, RA=RA+1, RR=RR+1) %>%
  as.matrix(.)
wave_male_4vs4_chi <- chisq.test(wave_male_4vs4)
wave_male_4vs4_pval <- wave_male_4vs4_chi$p.value



## make the p vals into a table - turn each into a data frame and merge
crab_female_1vs2_pval <- as.data.frame(crab_female_1vs2_pval) %>% dplyr::rename(inv1_inv2=crab_female_1vs2_pval)
crab_female_1vs3_pval <- as.data.frame(crab_female_1vs3_pval) %>% dplyr::rename(inv1_inv3=crab_female_1vs3_pval)
crab_female_1vs4_pval <- as.data.frame(crab_female_1vs4_pval) %>% dplyr::rename(inv1_inv4=crab_female_1vs4_pval)
crab_female_2vs3_pval <- as.data.frame(crab_female_2vs3_pval) %>% dplyr::rename(inv2_inv3=crab_female_2vs3_pval)
crab_female_2vs4_pval <- as.data.frame(crab_female_2vs4_pval) %>% dplyr::rename(inv2_inv4=crab_female_2vs4_pval)
crab_female_3vs4_pval <- as.data.frame(crab_female_3vs4_pval) %>% dplyr::rename(inv3_inv4=crab_female_3vs4_pval)
crab_female_1vs1_pval <- as.data.frame(crab_female_1vs1_pval) %>% dplyr::rename(inv1_inv1=crab_female_1vs1_pval)
crab_female_2vs2_pval <- as.data.frame(crab_female_2vs2_pval) %>% dplyr::rename(inv2_inv2=crab_female_2vs2_pval)
crab_female_3vs3_pval <- as.data.frame(crab_female_3vs3_pval) %>% dplyr::rename(inv3_inv3=crab_female_3vs3_pval)
crab_female_4vs4_pval <- as.data.frame(crab_female_4vs4_pval) %>% dplyr::rename(inv4_inv4=crab_female_4vs4_pval)

crab_male_1vs2_pval <- as.data.frame(crab_male_1vs2_pval) %>% dplyr::rename(inv1_inv2=crab_male_1vs2_pval)
crab_male_1vs3_pval <- as.data.frame(crab_male_1vs3_pval) %>% dplyr::rename(inv1_inv3=crab_male_1vs3_pval)
crab_male_1vs4_pval <- as.data.frame(crab_male_1vs4_pval) %>% dplyr::rename(inv1_inv4=crab_male_1vs4_pval)
crab_male_2vs3_pval <- as.data.frame(crab_male_2vs3_pval) %>% dplyr::rename(inv2_inv3=crab_male_2vs3_pval)
crab_male_2vs4_pval <- as.data.frame(crab_male_2vs4_pval) %>% dplyr::rename(inv2_inv4=crab_male_2vs4_pval)
crab_male_3vs4_pval <- as.data.frame(crab_male_3vs4_pval) %>% dplyr::rename(inv3_inv4=crab_male_3vs4_pval)
crab_male_1vs1_pval <- as.data.frame(crab_male_1vs1_pval) %>% dplyr::rename(inv1_inv1=crab_male_1vs1_pval)
crab_male_2vs2_pval <- as.data.frame(crab_male_2vs2_pval) %>% dplyr::rename(inv2_inv2=crab_male_2vs2_pval)
crab_male_3vs3_pval <- as.data.frame(crab_male_3vs3_pval) %>% dplyr::rename(inv3_inv3=crab_male_3vs3_pval)
crab_male_4vs4_pval <- as.data.frame(crab_male_4vs4_pval) %>% dplyr::rename(inv4_inv4=crab_male_4vs4_pval)

wave_female_1vs2_pval <- as.data.frame(wave_female_1vs2_pval) %>% dplyr::rename(inv1_inv2=wave_female_1vs2_pval)
wave_female_1vs3_pval <- as.data.frame(wave_female_1vs3_pval) %>% dplyr::rename(inv1_inv3=wave_female_1vs3_pval)
wave_female_1vs4_pval <- as.data.frame(wave_female_1vs4_pval) %>% dplyr::rename(inv1_inv4=wave_female_1vs4_pval)
wave_female_2vs3_pval <- as.data.frame(wave_female_2vs3_pval) %>% dplyr::rename(inv2_inv3=wave_female_2vs3_pval)
wave_female_2vs4_pval <- as.data.frame(wave_female_2vs4_pval) %>% dplyr::rename(inv2_inv4=wave_female_2vs4_pval)
wave_female_3vs4_pval <- as.data.frame(wave_female_3vs4_pval) %>% dplyr::rename(inv3_inv4=wave_female_3vs4_pval)
wave_female_1vs1_pval <- as.data.frame(wave_female_1vs1_pval) %>% dplyr::rename(inv1_inv1=wave_female_1vs1_pval)
wave_female_2vs2_pval <- as.data.frame(wave_female_2vs2_pval) %>% dplyr::rename(inv2_inv2=wave_female_2vs2_pval)
wave_female_3vs3_pval <- as.data.frame(wave_female_3vs3_pval) %>% dplyr::rename(inv3_inv3=wave_female_3vs3_pval)
wave_female_4vs4_pval <- as.data.frame(wave_female_4vs4_pval) %>% dplyr::rename(inv4_inv4=wave_female_4vs4_pval)

wave_male_1vs2_pval <- as.data.frame(wave_male_1vs2_pval) %>% dplyr::rename(inv1_inv2=wave_male_1vs2_pval)
wave_male_1vs3_pval <- as.data.frame(wave_male_1vs3_pval) %>% dplyr::rename(inv1_inv3=wave_male_1vs3_pval)
wave_male_1vs4_pval <- as.data.frame(wave_male_1vs4_pval) %>% dplyr::rename(inv1_inv4=wave_male_1vs4_pval)
wave_male_2vs3_pval <- as.data.frame(wave_male_2vs3_pval) %>% dplyr::rename(inv2_inv3=wave_male_2vs3_pval)
wave_male_2vs4_pval <- as.data.frame(wave_male_2vs4_pval) %>% dplyr::rename(inv2_inv4=wave_male_2vs4_pval)
wave_male_3vs4_pval <- as.data.frame(wave_male_3vs4_pval) %>% dplyr::rename(inv3_inv4=wave_male_3vs4_pval)
wave_male_1vs1_pval <- as.data.frame(wave_male_1vs1_pval) %>% dplyr::rename(inv1_inv1=wave_male_1vs1_pval)
wave_male_2vs2_pval <- as.data.frame(wave_male_2vs2_pval) %>% dplyr::rename(inv2_inv2=wave_male_2vs2_pval)
wave_male_3vs3_pval <- as.data.frame(wave_male_3vs3_pval) %>% dplyr::rename(inv3_inv3=wave_male_3vs3_pval)
wave_male_4vs4_pval <- as.data.frame(wave_male_4vs4_pval) %>% dplyr::rename(inv4_inv4=wave_male_4vs4_pval)

crab_female_pvals <- bind_cols(crab_female_1vs2_pval, crab_female_1vs3_pval) %>% 
  bind_cols(., crab_female_1vs4_pval) %>% bind_cols(., crab_female_2vs3_pval) %>% 
  bind_cols(., crab_female_2vs4_pval) %>% bind_cols(., crab_female_3vs4_pval) %>%
  bind_cols(., crab_female_1vs1_pval) %>% bind_cols(., crab_female_2vs2_pval) %>%
  bind_cols(., crab_female_3vs3_pval) %>% bind_cols(., crab_female_4vs4_pval)
crab_male_pvals <- bind_cols(crab_male_1vs2_pval, crab_male_1vs3_pval) %>% 
  bind_cols(., crab_male_1vs4_pval) %>% bind_cols(., crab_male_2vs3_pval) %>% 
  bind_cols(., crab_male_2vs4_pval) %>% bind_cols(., crab_male_3vs4_pval)%>%
  bind_cols(., crab_male_1vs1_pval) %>% bind_cols(., crab_male_2vs2_pval) %>%
  bind_cols(., crab_male_3vs3_pval) %>% bind_cols(., crab_male_4vs4_pval)
wave_female_pvals <- bind_cols(wave_female_1vs2_pval, wave_female_1vs3_pval) %>% 
  bind_cols(., wave_female_1vs4_pval) %>% bind_cols(., wave_female_2vs3_pval) %>% 
  bind_cols(., wave_female_2vs4_pval) %>% bind_cols(., wave_female_3vs4_pval) %>%
  bind_cols(., wave_female_1vs1_pval) %>% bind_cols(., wave_female_2vs2_pval) %>%
  bind_cols(., wave_female_3vs3_pval) %>% bind_cols(., wave_female_4vs4_pval)
wave_male_pvals <- bind_cols(wave_male_1vs2_pval, wave_male_1vs3_pval) %>% 
  bind_cols(., wave_male_1vs4_pval) %>% bind_cols(., wave_male_2vs3_pval) %>% 
  bind_cols(., wave_male_2vs4_pval) %>% bind_cols(., wave_male_3vs4_pval) %>%
  bind_cols(., wave_male_1vs1_pval) %>% bind_cols(., wave_male_2vs2_pval) %>%
  bind_cols(., wave_male_3vs3_pval) %>% bind_cols(., wave_male_4vs4_pval)

row.names(crab_female_pvals) <- "crab_female"
row.names(crab_male_pvals) <- "crab_male"
row.names(wave_female_pvals) <- "wave_female"
row.names(wave_male_pvals) <- "wave_male"

pvals <- union(crab_female_pvals, crab_male_pvals) %>% union(., wave_female_pvals) %>% union(., wave_male_pvals)

## make data tidy
pvals_tidy <- pvals %>% rownames_to_column("group") %>% gather(comp, pval, 2:11) %>% 
  separate(group, c("ecotype", "sex"), sep="_")



### do the inversion-sex associations
# for each inversion, do sex vs genotype

all_genos <- bind_rows(female_genos, male_genos)

crab_inv1_genos <- all_genos %>% filter(ecotype=="Crab") %>% select(snail_ID, inv1, sex)
crab_inv2_genos <- all_genos %>% filter(ecotype=="Crab") %>% select(snail_ID, inv2, sex)
crab_inv3_genos <- all_genos %>% filter(ecotype=="Crab") %>% select(snail_ID, inv3, sex)
crab_inv4_genos <- all_genos %>% filter(ecotype=="Crab") %>% select(snail_ID, inv4, sex)

wave_inv1_genos <- all_genos %>% filter(ecotype=="Wave") %>% select(snail_ID, inv1, sex)
wave_inv2_genos <- all_genos %>% filter(ecotype=="Wave") %>% select(snail_ID, inv2, sex)
wave_inv3_genos <- all_genos %>% filter(ecotype=="Wave") %>% select(snail_ID, inv3, sex)
wave_inv4_genos <- all_genos %>% filter(ecotype=="Wave") %>% select(snail_ID, inv4, sex)

# do counts and make into contingency table then into matrix
# and sort out NAs/missing values
# add 1 to everything to sort out the 0s issue
crab_inv1_counts <- crab_inv1_genos %>% group_by(sex, inv1) %>% count() %>% spread(inv1, n) %>%
  mutate(RA=0, AA=0) %>% mutate(AA=AA+1, RA=RA+1, RR=RR+1) %>% column_to_rownames("sex") %>% as.matrix(.)
crab_inv2_counts <- crab_inv2_genos %>% group_by(sex, inv2) %>% count() %>% spread(inv2, n) %>%
   mutate(AA=AA+1, RA=RA+1, RR=RR+1) %>% column_to_rownames("sex") %>% as.matrix(.)
crab_inv2_counts[is.na(crab_inv2_counts)] <- 1
crab_inv3_counts <- crab_inv3_genos %>% group_by(sex, inv3) %>% count() %>% spread(inv3, n) %>%
  mutate(AA=AA+1, RA=RA+1, RR=RR+1) %>% column_to_rownames("sex") %>% as.matrix(.)
crab_inv3_counts[is.na(crab_inv3_counts)] <- 1
crab_inv4_counts <- crab_inv4_genos %>% group_by(sex, inv4) %>% count() %>% spread(inv4, n) %>%
  mutate(AA=AA+1, RA=RA+1, RR=RR+1) %>% column_to_rownames("sex") %>% select(-`<NA>`) %>% as.matrix(.)

wave_inv1_counts <- wave_inv1_genos %>% group_by(sex, inv1) %>% count() %>% spread(inv1, n) %>%
  column_to_rownames("sex") %>% as.matrix(.)
wave_inv2_counts <- wave_inv2_genos %>% group_by(sex, inv2) %>% count() %>% spread(inv2, n) %>%
  column_to_rownames("sex") %>% as.matrix(.)
wave_inv3_counts <- wave_inv3_genos %>% group_by(sex, inv3) %>% count() %>% spread(inv3, n) %>%
  column_to_rownames("sex") %>% as.matrix(.)
wave_inv3_counts[is.na(wave_inv3_counts)] <- 0
wave_inv4_counts <- wave_inv4_genos %>% group_by(sex, inv4) %>% count() %>% spread(inv4, n) %>%
  column_to_rownames("sex") %>% as.matrix(.)
wave_inv4_counts[is.na(wave_inv4_counts)] <- 0

# do the chisq tests and save the output/pval
crab_inv1_chi <- chisq.test(crab_inv1_counts)
crab_inv1_pval <- crab_inv1_chi$p.value %>% as.data.frame(.) %>% dplyr::rename(inv1=1)
crab_inv2_chi <- chisq.test(crab_inv2_counts)
crab_inv2_pval <- crab_inv2_chi$p.value %>% as.data.frame(.) %>% dplyr::rename(inv2=1)
crab_inv3_chi <- chisq.test(crab_inv3_counts)
crab_inv3_pval <- crab_inv3_chi$p.value %>% as.data.frame(.) %>% dplyr::rename(inv3=1)
crab_inv4_chi <- chisq.test(crab_inv4_counts)
crab_inv4_pval <- crab_inv4_chi$p.value %>% as.data.frame(.) %>% dplyr::rename(inv4=1)

wave_inv1_chi <- chisq.test(wave_inv1_counts)
wave_inv1_pval <- wave_inv1_chi$p.value %>% as.data.frame(.) %>% dplyr::rename(inv1=1)
wave_inv2_chi <- chisq.test(wave_inv2_counts)
wave_inv2_pval <- wave_inv2_chi$p.value %>% as.data.frame(.) %>% dplyr::rename(inv2=1)
wave_inv3_chi <- chisq.test(wave_inv3_counts)
wave_inv3_pval <- wave_inv3_chi$p.value %>% as.data.frame(.) %>% dplyr::rename(inv3=1)
wave_inv4_chi <- chisq.test(wave_inv4_counts)
wave_inv4_pval <- wave_inv4_chi$p.value %>% as.data.frame(.) %>% dplyr::rename(inv4=1)

crab_pvals <- bind_cols(crab_inv1_pval, crab_inv2_pval) %>% bind_cols(., crab_inv3_pval) %>%
  bind_cols(., crab_inv4_pval)
wave_pvals <- bind_cols(wave_inv1_pval, wave_inv2_pval) %>% bind_cols(., wave_inv3_pval) %>%
  bind_cols(wave_inv4_pval)

row.names(crab_pvals) <- "crab"
row.names(wave_pvals) <- "wave"

sex_comp_pvals <- union(crab_pvals, wave_pvals)
sex_comp_pvals_tidy <- sex_comp_pvals %>% rownames_to_column("ecotype") %>% gather(inversion, pval, 2:5) %>%
  mutate(comp="M_F")
sex_comp_pvals_tidy <- sex_comp_pvals_tidy %>% mutate(log_p=-(log(pval)))

                     


##########################################################


## doing r^2 correlation calculations/plots for measure of association rather than just 
##    significance of association

## have genotypes as 0,1,2 and sex as 0,1,2
## same as before, separate table for sex/ecotype and for inversion pairs
## then do the correlation test

# (1 is RR for inv1,2,3; 3 is RR for inv4)

#### # inv-inv
crab_female_corr <- crab_female_genos %>% select(snail_ID, inv1, inv2, inv3, inv4)
crab_female_corr$inv1 <- recode(crab_female_corr$inv1, "RR"=1, "RA"=2, "AA"=3)
crab_female_corr$inv2 <- recode(crab_female_corr$inv2, "RR"=1, "RA"=2, "AA"=3)
crab_female_corr$inv3 <- recode(crab_female_corr$inv3, "RR"=1, "RA"=2, "AA"=3)
crab_female_corr$inv4 <- recode(crab_female_corr$inv4, "RR"=1, "RA"=2, "AA"=3)
crab_female_corr <- crab_female_corr %>% column_to_rownames("snail_ID")
crab_female_cormat <- round(cor(crab_female_corr),2)

crab_male_corr <- crab_male_genos %>% select(snail_ID, inv1, inv2, inv3, inv4)
crab_male_corr$inv1 <- recode(crab_male_corr$inv1, "RR"=1, "RA"=2, "AA"=3)
crab_male_corr$inv2 <- recode(crab_male_corr$inv2, "RR"=1, "RA"=2, "AA"=3)
crab_male_corr$inv3 <- recode(crab_male_corr$inv3, "RR"=1, "RA"=2, "AA"=3)
crab_male_corr$inv4 <- recode(crab_male_corr$inv4, "RR"=1, "RA"=2, "AA"=3)
crab_male_corr <- crab_male_corr %>% column_to_rownames("snail_ID")
crab_male_cormat <- round(cor(crab_male_corr),2)

wave_female_corr <- wave_female_genos %>% select(snail_ID, inv1, inv2, inv3, inv4)
wave_female_corr$inv1 <- recode(wave_female_corr$inv1, "RR"=1, "RA"=2, "AA"=3)
wave_female_corr$inv2 <- recode(wave_female_corr$inv2, "RR"=1, "RA"=2, "AA"=3)
wave_female_corr$inv3 <- recode(wave_female_corr$inv3, "RR"=1, "RA"=2, "AA"=3)
wave_female_corr$inv4 <- recode(wave_female_corr$inv4, "RR"=1, "RA"=2, "AA"=3)
wave_female_corr <- wave_female_corr %>% column_to_rownames("snail_ID")
wave_female_cormat <- round(cor(wave_female_corr),2)

wave_male_corr <- wave_male_genos %>% select(snail_ID, inv1, inv2, inv3, inv4)
wave_male_corr$inv1 <- recode(wave_male_corr$inv1, "RR"=1, "RA"=2, "AA"=3)
wave_male_corr$inv2 <- recode(wave_male_corr$inv2, "RR"=1, "RA"=2, "AA"=3)
wave_male_corr$inv3 <- recode(wave_male_corr$inv3, "RR"=1, "RA"=2, "AA"=3)
wave_male_corr$inv4 <- recode(wave_male_corr$inv4, "RR"=1, "RA"=2, "AA"=3)
wave_male_corr <- wave_male_corr %>% column_to_rownames("snail_ID")
wave_male_cormat <- round(cor(wave_male_corr),2)

get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

crab_female_cormat_upper <- get_upper_tri(crab_female_cormat)
crab_male_cormat_upper <- get_upper_tri(crab_male_cormat)
wave_female_cormat_upper <- get_upper_tri(wave_female_cormat)
wave_male_cormat_upper <- get_upper_tri(wave_male_cormat)

crab_female_melted <- melt(crab_female_cormat_upper, na.rm=TRUE)
crab_male_melted <- melt(crab_male_cormat_upper, na.rm=TRUE)
wave_female_melted <- melt(wave_female_cormat_upper, na.rm=TRUE)
wave_male_melted <- melt(wave_male_cormat_upper, na.rm=TRUE)

View(wave_male_cormat_upper)
# NAs for a lot of the crab ones are bc the correlation couldn't be calculated, not bc i removed them


### inv-sex associations

genos <- bind_rows(female_genos, male_genos) %>% filter(ecotype=="Crab" | ecotype=="Wave")
crab_genos <- genos %>% filter(ecotype=="Crab") %>% select(-ecotype)
wave_genos <- genos %>% filter(ecotype=="Wave") %>% select(-ecotype)

crab_genos$inv1 <- recode(crab_genos$inv1, "RR"=1, "RA"=2, "AA"=3)
crab_genos$inv2 <- recode(crab_genos$inv2, "RR"=1, "RA"=2, "AA"=3)
crab_genos$inv3 <- recode(crab_genos$inv3, "RR"=1, "RA"=2, "AA"=3)
crab_genos$inv4 <- recode(crab_genos$inv4, "RR"=1, "RA"=2, "AA"=3)
crab_genos$sex <- recode(crab_genos$sex, "female"=0, "male"=1)

wave_genos$inv1 <- recode(wave_genos$inv1, "RR"=1, "RA"=2, "AA"=3)
wave_genos$inv2 <- recode(wave_genos$inv2, "RR"=1, "RA"=2, "AA"=3)
wave_genos$inv3 <- recode(wave_genos$inv3, "RR"=1, "RA"=2, "AA"=3)
wave_genos$inv4 <- recode(wave_genos$inv4, "RR"=1, "RA"=2, "AA"=3)
wave_genos$sex <- recode(wave_genos$sex, "female"=0, "male"=1)


crab_corr <- crab_genos %>% column_to_rownames("snail_ID") %>% filter(!is.na(inv4)) 
crab_cormat <- round(cor(crab_corr),2)
crab_cormat_upper <- crab_cormat %>% get_upper_tri(.)
crab_melted <- melt(crab_cormat_upper, na.rm=TRUE)
crab_melted <- crab_melted %>% mutate(ecotype="Crab")

wave_corr <- wave_genos %>% column_to_rownames("snail_ID") %>% filter(!is.na(inv4))
wave_cormat <- round(cor(wave_corr),2)
wave_cormat_upper <- wave_cormat %>% get_upper_tri(.)
wave_melted <- melt(wave_cormat_upper, na.rm=TRUE)
wave_melted <- wave_melted %>% mutate(ecotype="Wave")

sexcomp_melted <- bind_rows(crab_melted, wave_melted)
sexcomp_melted <- sexcomp_melted %>% mutate(r2=value*value)

sexcomp_melted_vals <- sexcomp_melted %>% filter(X2=="sex" & X1!="sex")
