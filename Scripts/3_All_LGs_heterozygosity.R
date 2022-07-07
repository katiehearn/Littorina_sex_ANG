rm(list=ls())
library(tidyverse)
library(Hmisc)
library(adegenet)
library(hexbin)
# set working directory

##############################################
##
##
## Testing all linkage groups for sex-linkage in Wave
##    Heterozygosity plots
##
##
##############################################

# read in map
map <- read.table("map_v11.txt", header=T)

# read in data for ANG
ANG_snail_info <- read.csv("ANG_snail_info.csv")
ANG_sex <- read.csv("ANG_dissections_20150723.csv") %>% select(snail_ID, sex_dissection)
ANG_path <- read.csv("ANG_SnailDataLitPath_20160702.csv") %>% select(snail_ID, DistAlongPath)
ANG_allLGs <- read.table("CZ013_ANG.GT.FORMAT", header=T)

# make data tidy and add sex, ecotype, map info
ANG_allLGs_tidy <- ANG_allLGs %>% rename("chrom"="CHROM", "pos"="POS") %>%
  gather(snail_ID, genotype, 3:382) %>% 
  mutate(geno2=ifelse(genotype=="0/0", 0,
                      ifelse(genotype=="0/1", 1,
                             ifelse(genotype=="1/1", 2, NA)))) %>%
  select(-genotype) %>% rename("genotype"="geno2") %>%
  left_join(., ANG_sex, by="snail_ID") %>%
  left_join(., ANG_path, by="snail_ID") %>% 
  mutate(Ecotype=ifelse(DistAlongPath<=68, "Crab",
                        ifelse(DistAlongPath>=108, "Wave", "Hybrid")))
#ANG_allLGs_tidy <- ANG_allLGs_tidy %>% separate(contig_ID, c("contig", "pos"), sep="_")
map_contigs <- map %>% select(-pos) %>% rename("chrom"="contig")
ANG_allLGs_tidy <- ANG_allLGs_tidy %>% left_join(., map_contigs, by="chrom") %>% filter(!is.na(LG))

ANG_allLGs %>% rename("chrom"="CHROM") %>% left_join(., map_contigs, by="chrom") %>% filter(!is.na(LG)) %>%
  select(chrom) %>% unique() %>% View()

# calculate genotype props for each ecotype and sex
ANG_allLGs_hets <- ANG_allLGs_tidy %>% unite(contig_ID, c("chrom", "pos"), sep="_") %>%
  filter(Ecotype!="Hybrid") %>% filter(!is.na(genotype)) %>%
  group_by(LG, as.factor(Ecotype), as.factor(sex_dissection), as.factor(contig_ID), as.factor(genotype), .drop=FALSE) %>%
  count() %>%
  group_by(LG, `as.factor(Ecotype)`, `as.factor(sex_dissection)`, `as.factor(contig_ID)`) %>%
  mutate(total=sum(n), prop=n/total) %>%
  rename("Ecotype"="as.factor(Ecotype)", "sex_dissection"="as.factor(sex_dissection)", 
         "contig_ID"="as.factor(contig_ID)", "genotype"="as.factor(genotype)") 

# keep hets only
ANG_allLGs_hets_only <- ANG_allLGs_hets %>% filter(genotype==1) %>% select(-n, -total) %>% spread(sex_dissection, prop)
ANG_allLGs_hets_only <- ANG_allLGs_hets_only %>% filter(female>=0)

# plot female vs male hets

# all LGs
#ggplot(ANG_allLGs_hets_only, aes(x=female, y=male)) + geom_point(shape=21, alpha=0.5, size=0.5) + 
#  facet_grid(cols=vars(LG), rows=vars(Ecotype)) + theme_bw() + 
#  xlab("Proportion of female heterozygotes") +
#  ylab("Proportion of male heterozygotes") +
#  theme(aspect.ratio = 1,
#        panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank(),
#        strip.text = element_text(size=8, face="bold"),
#        axis.text = element_text(size=5),
#        axis.title = element_text(size=6, face="bold")) +
#  scale_x_continuous(limits=c(0,1)) +
#  scale_y_continuous(limits=c(0,1)) +
#  geom_abline(slope=1, intercept=0, colour="darkred", linetype=2, size=0.5)

## hexplot to show density without using alpha transparency
ANG_allLGs_hets_only %>%
  ggplot(., aes(x=female, y=male)) + geom_hex() + 
  facet_wrap(~LG+Ecotype) + theme_bw() + 
  scale_fill_gradient(high="grey88", low="black", trans="log10") +
  xlab("Proportion of female heterozygotes") +
  ylab("Proportion of male heterozygotes") +
  theme(aspect.ratio = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=8, face="bold"),
        axis.text = element_text(size=5),
        axis.title = element_text(size=6, face="bold")) +
  scale_x_continuous(limits=c(0,1)) +
  scale_y_continuous(limits=c(0,1)) +
  geom_abline(slope=1, intercept=0, colour="darkred", linetype=2, size=0.5)
  

#### calculate top 1% of residuals for SNPs- in Wave- then plot those on the hex plot in a different colour

# calculate residuals
ANG_allLGs_hets_only <- ANG_allLGs_hets_only %>% mutate(resid=male-female)
# mark out top 1% in Wave (really bottom 1% as we want the most negative residuals)
# 113156 SNPs so top 1% is 1132 SNPs 
# get list of top SNPs then mutate the full dataframe to label them
top1pc_list <- ANG_allLGs_hets_only %>% filter(Ecotype=="Wave") %>% dplyr::arrange(resid) %>%
  slice(1:1132) %>% pull(., contig_ID)
top1pc <- ANG_allLGs_hets_only %>% filter(Ecotype=="Wave") %>% dplyr::arrange(resid) %>%
  slice(1:1132)
  
ANG_allLGs_hets_only %>%
  ggplot(., aes(x=female, y=male)) + geom_hex() + geom_point(data=top1pc, shape=21, colour="orange") +
  facet_wrap(~LG+Ecotype) + theme_bw() + 
  scale_fill_gradient(high="grey88", low="black", trans="log10") +
  xlab("Proportion of female heterozygotes") +
  ylab("Proportion of male heterozygotes") +
  theme(aspect.ratio = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=8, face="bold"),
        axis.text = element_text(size=5),
        axis.title = element_text(size=6, face="bold")) +
  scale_x_continuous(limits=c(0,1)) +
  scale_y_continuous(limits=c(0,1)) +
  geom_abline(slope=1, intercept=0, colour="darkred", linetype=2, size=0.5)


# barplot of percentage of top 1% SNPs on each linkage group

top1pc <- top1pc %>% group_by(LG) %>% mutate(count=count(LG), pc=count/1132*100)

ggplot(top1pc, aes(x=LG, y=pc)) + geom_col(colour="orange") +
 theme_bw() +
 xlab("Linkage group") +
 ylab("Percentage of SNPs") +
 theme(aspect.ratio = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=8, face="bold"),
        axis.text = element_text(size=5),
        axis.title = element_text(size=6, face="bold"))

