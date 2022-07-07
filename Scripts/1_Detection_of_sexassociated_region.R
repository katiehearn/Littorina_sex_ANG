rm(list=ls())
library(tidyverse)
library(zoo)
# set working directory


#####################
#####################
### Detection of a sex-associated region
### Sex differences in heterozygosity on LG12
####################
####################


#########
## ANG
#########

ANG_genotypes <- read.csv("SL.filt2.LG12.GT.csv", header=TRUE)
ANG_dissections <- read.csv("ANG_dissections_20150723.csv", header=TRUE)
ANG_pathdistance <- read.csv("ANG_SnailDataLitPath_20160702.csv", header=TRUE)
map1 <- read.table("map_v11.txt", header=TRUE)

## extract sex from dissections and exclude juveniles (and "u"s)
ANG_sex_info <- select(ANG_dissections, snail_ID, sex_dissection)
ANG_sex_info_filtered <- ANG_sex_info %>% filter(!is.na(sex_dissection), sex_dissection!="juvenile")

## extract distance along path from distance files
ANG_distances <- ANG_pathdistance %>% 
  select(snail_ID, DistAlongPath) %>%
  filter(DistAlongPath!="NA")

## merge sex with distance along path
ANG_merged_data <- inner_join(ANG_sex_info_filtered, ANG_distances, by="snail_ID")

##### remove 'hybrid' individuals- 10m either side of the transition
ANG_merged_data <- ANG_merged_data %>% filter(DistAlongPath>=88 | DistAlongPath<=68)

## classify individuals as Wave or Crab
ANG_merged_data <- ANG_merged_data %>%
  mutate(ecotype=ifelse(DistAlongPath<78, "Crab", "Wave"))


## tidy genotype info ready for merging - merge columns, transpose, extract column names,
#   merge names to data frame, make names column heading
ANG_merged_genotypes <- ANG_genotypes %>% unite(snail_ID, CHROM, POS, sep="_")
ANG_genotypes_transposed <- as.tibble(t(ANG_merged_genotypes))
names <- colnames(ANG_merged_genotypes)
ANG_genotypes_tidy <- data.frame(names, ANG_genotypes_transposed)
names(ANG_genotypes_tidy)
colnames(ANG_genotypes_tidy)=ANG_genotypes_tidy[1,]
ANG_genotypes_tidy <- rename(ANG_genotypes_tidy, "snail_ID"="383")
ANG_genotypes_tidy <- slice(ANG_genotypes_tidy, -1)

## merge genotype data to other merged data
ANG_merged_full <- inner_join(ANG_merged_data, ANG_genotypes_tidy, by="snail_ID")

## make the data tidy, and arrange by snail ID
ANG_merged_fulltidy <- gather(ANG_merged_full, key="contig_ID", value="genotype", -c(1:4))
ANG_merged_fulltidy <- arrange(ANG_merged_fulltidy, snail_ID)
ANG_merged_fulltidy <- rename(ANG_merged_fulltidy, sex=sex_dissection)


## filter out SNPs on other LGs using contigs from the map
ANGtester <- ANG_merged_fulltidy %>% separate(contig_ID, c("contig","pos"), sep="_")
testmap <- map1 %>% filter(LG=="12")
ANGtestfilter <- ANGtester %>% filter(contig %in% testmap$contig)
ANGtestfilter <- ANGtestfilter %>% unite(contig_ID, contig, pos, sep="_")

## now group by sex and ecotype and do proportion of heterozygotes
#   count of each genotype
ANG_genotype_counts <- ANGtestfilter %>% group_by(sex, ecotype, contig_ID) %>%
  count(genotype) %>%
  rename(count=n)
ANG_genotype_counts <- ANG_genotype_counts %>% filter(genotype!="juvenile")
# impute misssing genotypes with count of 0
ANG_genotype_counts <- ANG_genotype_counts %>% 
  ungroup() %>%  
  complete(sex, ecotype, contig_ID, genotype, fill=list(count=0))
# remove the rows with gender 'u' that got added for some reason
ANG_genotype_counts <- ANG_genotype_counts %>% filter(sex!="juvenile")
#   remove rows with missing genotypes (ie ./.)
ANG_genotype_counts_filtered <- filter(ANG_genotype_counts, !is.na(genotype))
#   calculate proportions of genotypes
ANG_genotype_proportions <- ANG_genotype_counts_filtered %>% 
  group_by(sex, ecotype, contig_ID) %>%
  mutate(total=sum(count),
         proportion=count/total) %>%
  arrange(contig_ID)
#   calculate proportion heterozygotes for each SNP for each sex
#ANG_genotype_proportions <- ANG_genotype_proportions %>% rename("0"="0",
#                                                                "1"="1",
#                                                                "2"="2")
ANG_heterozygote_proportions <- ANG_genotype_proportions %>%
  filter(genotype==" 1") %>%
  select(-count, -total,-genotype) %>%
  spread(sex,proportion) %>%
  arrange(contig_ID)
# change NAs to 0s - they're NA as all SNPs are homozygous for one allele or another
#ANG_heterozygote_proportions[is.na(ANG_heterozygote_proportions)] <- 0

## plot proportion of heterozygotes
ANG_heterozygotes_plot <- ggplot(ANG_heterozygote_proportions, aes(x=female, y=male)) +
  geom_point(shape=21) +
  facet_wrap(~ecotype, scales="free") +
  xlab("proportion of heterozygotes (females)") +
  ylab("proportion of heterozygotes (males)") +
  theme_linedraw() +
  theme(aspect.ratio=1, 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


## add map info
ANG_heterozygote_proportions <- ANG_heterozygote_proportions %>% rename(f_het_prop=female, m_het_prop=male) %>% 
separate(contig_ID, c("contig", "pos"), sep="_") %>% ungroup()


## add map positions from hernan's list (for contigs with multiple map positions just uses most common av)
mapone <- read.table("mapposhernan.txt", header=TRUE)
mapone <- mapone %>% select(LG, contig, pos, av, Mpos_majority) %>% 
  rename(av_new=Mpos_majority) %>% filter(LG==12) %>% select(-av, -LG, -pos)
mapone$contig <- as.character(mapone$contig)
ANGmap_SNPs <- left_join(ANG_heterozygote_proportions, mapone, by="contig")
ANGmap_SNPs <- ANGmap_SNPs %>% rename(av=av_new)
## correct number of rows now (no duplicates) and all have a map pos (no NAs)


## calculate residuals

# add 1:1 line to original graph
ggplot(ANGmap_SNPs, aes(x=f_het_prop, y=m_het_prop)) +
  geom_point(shape=21) +
  facet_wrap(~ecotype, scales="free") +
  scale_x_continuous(limits=c(0,1)) +
  scale_y_continuous(limits=c(0,1)) +
  theme_linedraw() +
  theme(aspect.ratio=1) +
  geom_abline(slope=1, intercept=0, colour="blue", linetype=2, size=1.5) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

## calculate residuals and plot against map position
##  r = y - ax - b

ANGresiduals_SNPs <- ANGmap_SNPs %>% mutate(residual=m_het_prop-f_het_prop)

ggplot(ANGresiduals_SNPs, aes(x=as.numeric(av), y=residual)) +
  geom_point(shape=21) +
  xlab("map position") +
  ylab("residual from 1:1 line") +
  scale_x_continuous(limits=c(0,60)) +
  theme_linedraw() +
  facet_wrap(~ecotype) +
  theme(aspect.ratio=0.5) +
  scale_y_continuous(limits=c(-1,1)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
