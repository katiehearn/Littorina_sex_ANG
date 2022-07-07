###########
### Plots to examine halpotypes / patterns of genotypes along LG12 and along the transect
###########

rm(list=ls())
# set working directory
library(tidyverse)
library(groupdata2)


female <- read.csv("ANGfemale.csv")
male <- read.csv("ANGmale.csv")
snail_info <- read.csv("ANG_snail_info.csv")
mappos <- read.csv("LG12_mappos_bycontig.csv")
dist <- read.csv("ANG_SnailDataLitPath_20160702.csv")


## make data tidy
female_tidy <- female %>% unite(contig_ID, CHROM, POS, sep="_") %>% gather(snail_ID, genotype, -1)

## put snps and snails in order
# first add info
dist <- dist %>% select(snail_ID, DistAlongPath)
female_tidy <- female_tidy %>% left_join(., dist, by="snail_ID") %>% 
  left_join(., mappos, by="contig_ID") 
female_tidy <- female_tidy %>% separate(contig_ID, c("contig", "pos"), sep="_") %>%
  arrange(av, contig, pos) %>% group(n=209, method="greedy", col_name="snporder")

female_tidy <- female_tidy %>% unite(contig_ID, contig, pos, sep="_") %>% 
  arrange(DistAlongPath) %>% group(n=8657, method="greedy", col_name="snailorder")

## filter to retain SNPs in focal clusters only

female_tidy$genotype <- as.factor(female_tidy$genotype)

cluster2317_0.89_snps <- read.csv("cluster2317_0.89.csv")
cluster2318_0.89_snps <- read.csv("cluster2318_0.89.csv")
cluster2398_0.88_snps <- read.csv("cluster2398_0.88.csv")
cluster3414_0.69_snps <- read.csv("cluster3414_0.69.csv")
cluster3508_0.67_snps <- read.csv("cluster3508_0.67.csv")

clustersnps_genos <- female_tidy %>% filter(contig_ID%in%cluster2317_0.89_snps$contig_ID |
                                          contig_ID%in%cluster2318_0.89_snps$contig_ID |
                                          contig_ID%in%cluster2398_0.88_snps$contig_ID |
                                          contig_ID%in%cluster3414_0.69_snps$contig_ID |
                                          contig_ID%in%cluster3508_0.67_snps$contig_ID)

clustersnps_genos <- clustersnps_genos %>% arrange(snporder) %>%
  group(n=209, method="greedy", col_name="snporder2")

# plot
cols <- c("0" = "palegreen", "1" = "green", "2" = "darkgreen")
ggplot(clustersnps_genos, aes(x=snailorder, y=as.numeric(snporder2))) + 
  geom_tile(aes(fill=genotype)) + 
  scale_fill_manual(values=cols) +
  theme_bw()



#### repeat for males

## make data tidy
male_tidy <- male %>% unite(contig_ID, CHROM, POS, sep="_") %>% gather(snail_ID, genotype, -1)

## put snps and snails in order
# first add info
male_tidy <- male_tidy %>% left_join(., dist, by="snail_ID") %>% 
  left_join(., mappos, by="contig_ID") 
male_tidy <- male_tidy %>% separate(contig_ID, c("contig", "pos"), sep="_") %>%
  arrange(av, contig, pos) %>% group(n=171, method="greedy", col_name="snporder")
male_tidy <- male_tidy %>% unite(contig_ID, contig, pos, sep="_") %>% 
  arrange(DistAlongPath) %>% group(n=8657, method="greedy", col_name="snailorder")
male_tidy$genotype <- as.factor(male_tidy$genotype)

clustersnps_m_genos <- male_tidy %>% filter(contig_ID%in%cluster2317_0.89_snps$contig_ID |
                                              contig_ID%in%cluster2318_0.89_snps$contig_ID |
                                              contig_ID%in%cluster2398_0.88_snps$contig_ID |
                                              contig_ID%in%cluster3414_0.69_snps$contig_ID |
                                              contig_ID%in%cluster3508_0.67_snps$contig_ID)
clustersnps_m_genos <- clustersnps_m_genos %>% arrange(snporder) %>%
  group(n=171, method="greedy", col_name="snporder2")
ggplot(clustersnps_m_genos, aes(x=snailorder, y=as.numeric(snporder2))) + 
  geom_tile(aes(fill=genotype)) + 
  scale_fill_manual(values=cols) +
  theme_bw()

