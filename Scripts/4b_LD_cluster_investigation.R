##########
## inspecting clusters
## to be run alongside/after the genetics_LDna script- requires that script to be loaded run

#### look at clusters
## first look at the four main clusters 

# extract list of SNPs in each cluster
cluster2317_0.89 <- clusters1$clusters$`2317_0.89` %>% data.frame(.) %>% rename("contig_ID"=".")
cluster2318_0.89 <- clusters1$clusters$`2318_0.89` %>% data.frame(.) %>% rename("contig_ID"=".")
cluster2398_0.88 <- clusters1$clusters$`2398_0.88` %>% data.frame(.) %>% rename("contig_ID"=".")
cluster3414_0.69 <- clusters1$clusters$`3414_0.69` %>% data.frame(.) %>% rename("contig_ID"=".")
cluster3508_0.67 <- clusters1$clusters$`3508_0.67` %>% data.frame(.) %>% rename("contig_ID"=".")

# add map positions to these SNPs
SNPslist <- read.csv("ANGfemale.csv") %>% select(CHROM, POS) %>% rename("contig"="CHROM", "pos"="POS")
SNPslist <- SNPslist %>% left_join(., map, by="contig")
SNPslist <- SNPslist %>% unite(contig_ID, c("contig", "pos"), sep="_")

cluster2317_0.89 <- cluster2317_0.89 %>% left_join(., SNPslist, by="contig_ID")
cluster2318_0.89 <- cluster2318_0.89 %>% left_join(., SNPslist, by="contig_ID")
cluster2398_0.88 <- cluster2398_0.88 %>% left_join(., SNPslist, by="contig_ID")
cluster3414_0.69 <- cluster3414_0.69 %>% left_join(., SNPslist, by="contig_ID")

# look at distribution of SNPs on LG12
ggplot(cluster2317_0.89, aes(x=av)) + geom_bar()
ggplot(cluster2318_0.89, aes(x=av)) + geom_bar()
ggplot(cluster2398_0.88, aes(x=av)) + geom_bar()
ggplot(cluster3414_0.69, aes(x=av)) + geom_bar()
ggplot(cluster3508_0.67, aes(x=av)) + geom_bar()

# make one dataset so can plot all on one plot
c1 <- cluster2317_0.89 %>% mutate(cluster="cluster2317_0.89")
c2 <- cluster2318_0.89 %>% mutate(cluster="cluster2318_0.89") 
c3 <- cluster2398_0.88 %>% mutate(cluster="cluster2398_0.88")
c4 <- cluster3414_0.69 %>% mutate(cluster="cluster3414_0.69")
c5 <- cluster3508_0.67 %>% mutate(cluster="cluster3508_0.67")
clusters4 <- bind_rows(c1, c2, c3, c4, c5)
ggplot(clusters4, aes(x=av, fill=cluster)) + geom_bar(width=0.6) + 
  facet_wrap(~cluster, ncol=1) +
  theme_bw()

# save list of SNPs for each cluster
#write.csv(cluster2317_0.89, file="cluster2317_0.89.csv", row.names = F)
#write.csv(cluster2318_0.89, file="cluster2318_0.89.csv", row.names = F)
#write.csv(cluster2398_0.88, file="cluster2398_0.88.csv", row.names = F)
#write.csv(cluster3414_0.69, file="cluster3414_0.69.csv", row.names = F)


# 2 more clusters - from edges400 and phi2-6
cluster3508_0.67 <- clusters1$clusters$`3508_0.67` %>% data.frame(.) %>% rename("contig_ID"=".")
cluster4583_0.44 <- clusters1$clusters$`4583_0.44` %>% data.frame(.) %>% rename("contig_ID"=".")

cluster3508_0.67 <- cluster3508_0.67 %>% left_join(., SNPslist, by="contig_ID")
cluster4583_0.44 <- cluster4583_0.44 %>% left_join(., SNPslist, by="contig_ID")

ggplot(cluster3508_0.67, aes(x=av)) + geom_bar()
ggplot(cluster4583_0.44, aes(x=av)) + geom_bar()

#write.csv(cluster3508_0.67, file="cluster3508_0.67.csv", row.names = F)
#write.csv(cluster4583_0.44, file="cluster4583_0.44.csv", row.names = F)


# cluster from edges200-400, phi1-2
cluster2043_0.92 <- clusters1$clusters$`2043_0.92` %>% data.frame(.) %>% rename("contig_ID"=".")
cluster2043_0.92 <- cluster2043_0.92 %>% left_join(., SNPslist, by="contig_ID")
ggplot(cluster2043_0.92, aes(x=av)) + geom_bar()
#write.csv(cluster2043_0.92, file="cluster2043_0.92.csv", row.names = F)

# cluster from edges 400-600, phi1-3
cluster2667_0.84 <- clusters1$clusters$`2667_0.84` %>% data.frame(.) %>% rename("contig_ID"=".")
cluster2667_0.84 <- cluster2667_0.84 %>% left_join(., SNPslist, by="contig_ID")
ggplot(cluster2667_0.84, aes(x=av)) + geom_bar()
#write.csv(cluster2667_0.84, file="cluster2667_0.84.csv", row.names = F)

# cluster from edges 200, phi1
cluster2401_0.88 <- clusters1$clusters$`2401_0.88` %>% data.frame(.) %>% rename("contig_ID"=".")
cluster2401_0.88 <- cluster2401_0.88 %>% left_join(., SNPslist, by="contig_ID")
ggplot(cluster2401_0.88, aes(x=av)) + geom_bar()
#write.csv(cluster2401_0.88, file="cluster2401_0.88.csv", row.names = F)

# cluster from edges 200, phi1-8
cluster2640_0.85 <- clusters1$clusters$`2640_0.85` %>% data.frame(.) %>% rename("contig_ID"=".")
cluster2640_0.85 <- cluster2640_0.85 %>% left_join(., SNPslist, by="contig_ID")
ggplot(cluster2640_0.85, aes(x=av)) + geom_bar()
#write.csv(cluster2640_0.85, file="cluster2640_0.85.csv", row.names = F)

# cluster from edges 200, phi1-3
cluster3415_0.69 <- clusters1$clusters$`3415_0.69` %>% data.frame(.) %>% rename("contig_ID"=".")
cluster3415_0.69 <- cluster3415_0.69 %>% left_join(., SNPslist, by="contig_ID")
ggplot(cluster3415_0.69, aes(x=av)) + geom_bar()
#write.csv(cluster3415_0.69, file="cluster3415_0.69.csv", row.names = F)

# cluster from edges 200, phi1
cluster4278_0.5 <- clusters1$clusters$`4278_0.5` %>% data.frame(.) %>% rename("contig_ID"=".")
cluster4278_0.5 <- cluster4278_0.5 %>% left_join(., SNPslist, by="contig_ID")
ggplot(cluster4278_0.5, aes(x=av)) + geom_bar()
#write.csv(cluster4278_0.5, file="cluster4278_0.5.csv", row.names = F)






#
#
#
#
#
#
###
### investigate main five clusters using PCAs
###
# set working directory
library(tidyverse)
library(Hmisc)
library(adegenet)

# read in genotype table
ANGfemale <- read.csv("ANGfemale.csv")

# transpose dataframe as adegenet needs SNPs in columns
ANGfemale <- ANGfemale %>% unite(contig_ID, c("CHROM", "POS"), sep="_")
test <- ANGfemale %>% select(-contig_ID)
testtranspose <- as_tibble(t(test))
names <- colnames(test)
testtranspose <- data.frame(names, testtranspose)
names(testtranspose)

# get list of SNPs to be colnames
SNPslist_names <- SNPslist[,1]
snails <- testtranspose %>% select(names)
testtranspose <- testtranspose %>% select(-names)
colnames(testtranspose) <- SNPslist_names
ANG_female <- cbind(snails, testtranspose)
ANG_female <- ANG_female %>% rename("snail_ID"="names")
# make snail names the rownames - so genind can take them
rownames(ANG_female) <- ANG_female$snail_ID
ANG_female <- ANG_female %>% select(-snail_ID)

# change genotype coding to adegenet format
ANG_female[ANG_female==0] <- "0/0"
ANG_female[ANG_female==1] <- "0/1"
ANG_female[ANG_female==2] <- "1/1"

# make a dataframe of ecotype info for the snails to use to colour the PCAs
ANG_female_info <- read.csv("ANG_snail_info.csv") %>% filter(sex_dissection=="female") %>% 
  select(-sex_dissection, -DistAlongPath) %>% arrange(snail_ID)

# filter the genotypes to only the ones in the focal cluster 
#cluster2317_0.89_SNPs <- cluster2317_0.89 %>% pull(., contig_ID)
cluster2317_0.89_SNPs <- SNPslist %>% filter(av>=32.98692 & av<=43.8010) %>% pull(., contig_ID)
cluster2317_0.89_genos <- ANG_female %>% select(all_of(cluster2317_0.89_SNPs))

# make genind object
cluster2317_0.89_genind <- df2genind(cluster2317_0.89_genos, ploidy=2, sep="/")
summary(cluster2317_0.89_genind)

# scale
cluster2317_0.89_scaled <- scaleGen(cluster2317_0.89_genind, NA.method="mean")
class(cluster2317_0.89_scaled)
dim(cluster2317_0.89_scaled)
cluster2317_0.89_scaled[1:5, 1:5]

# PCA
cluster2317_0.89_PCA <- dudi.pca(cluster2317_0.89_scaled, cent=FALSE, scale=FALSE, scannf=TRUE)
barplot(cluster2317_0.89_PCA$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))

# save PCs 
cluster2317_0.89_PCs <- cluster2317_0.89_PCA$li

# add grouping info to add colour to plot - first make rownames a col for the PCs
cluster2317_0.89_PCs <- rownames_to_column(cluster2317_0.89_PCs, var="snail_ID")
cluster2317_0.89_PCs <- cluster2317_0.89_PCs %>% left_join(., ANG_female_info, by="snail_ID")

# plot with colour 
ggplot(cluster2317_0.89_PCs, aes(x=Axis1, y=Axis2, colour=ecotype)) + 
  geom_point(alpha=0.4, size=5) + 
  theme_bw() +
  xlab("PC1") + ylab("PC2") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(text=element_text(size=14))

# also a zoom in of the clustered section
ggplot(cluster2317_0.89_PCs, aes(x=Axis1, y=Axis2, colour=ecotype)) + 
  geom_point(alpha=0.4, size=5) + 
  theme_bw() +
  xlab("PC1") + ylab("PC2") +
  #guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(text=element_text(size=14)) +
  scale_x_continuous(limits=c(-20, 20)) + scale_y_continuous(limits=c(-20, 5))



################
### next cluster
### 2318_0.89

#cluster2318_0.89_SNPs <- cluster2318_0.89 %>% pull(., contig_ID)
cluster2318_0.89_SNPs <- SNPslist %>% filter(av>=33.624 & av<=60.235) %>% pull(., contig_ID)
cluster2318_0.89_genos <- ANG_female %>% select(all_of(cluster2318_0.89_SNPs))

cluster2318_0.89_genind <- df2genind(cluster2318_0.89_genos, ploidy=2, sep="/")

cluster2318_0.89_scaled <- scaleGen(cluster2318_0.89_genind, NA.method="mean")

cluster2318_0.89_PCA <- dudi.pca(cluster2318_0.89_scaled, cent=FALSE, scale=FALSE, scannf=TRUE)

cluster2318_0.89_PCs <- cluster2318_0.89_PCA$li

cluster2318_0.89_PCs <- rownames_to_column(cluster2318_0.89_PCs, var="snail_ID")
cluster2318_0.89_PCs <- cluster2318_0.89_PCs %>% left_join(., ANG_female_info, by="snail_ID")

ggplot(cluster2318_0.89_PCs, aes(x=Axis1, y=Axis2, colour=ecotype)) + 
  geom_point(alpha=0.4, size=5) + 
  theme_bw() +
  xlab("PC1") + ylab("PC2") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(text=element_text(size=14))

ggplot(cluster2318_0.89_PCs, aes(x=Axis1, y=Axis2, colour=ecotype)) + 
  geom_point(alpha=0.4, size=5) + 
  theme_bw() +
  xlab("PC1") + ylab("PC2") +
  #guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(text=element_text(size=14)) +
  scale_x_continuous(limits=c(-75, 15)) + scale_y_continuous(limits=c(-20, 20))



################
### next cluster
### 2398_0.88

#cluster2398_0.88_SNPs <- cluster2398_0.88 %>% pull(., contig_ID)
cluster2398_0.88_SNPs <- SNPslist %>% filter(av>=48.712 & av<=60.235) %>% pull(., contig_ID)
cluster2398_0.88_genos <- ANG_female %>% select(all_of(cluster2398_0.88_SNPs))

cluster2398_0.88_genind <- df2genind(cluster2398_0.88_genos, ploidy=2, sep="/")

cluster2398_0.88_scaled <- scaleGen(cluster2398_0.88_genind, NA.method="mean")

cluster2398_0.88_PCA <- dudi.pca(cluster2398_0.88_scaled, cent=FALSE, scale=FALSE, scannf=TRUE)

cluster2398_0.88_PCs <- cluster2398_0.88_PCA$li


cluster2398_0.88_PCs <- rownames_to_column(cluster2398_0.88_PCs, var="snail_ID")
cluster2398_0.88_PCs <- cluster2398_0.88_PCs %>% left_join(., ANG_female_info, by="snail_ID") %>%
  rename(Ecotype=ecotype)
cluster2398_0.88_PCs[is.na(cluster2398_0.88_PCs)] <- "Hybrid"
cluster2398_0.88_PCs$Ecotype <- factor(cluster2398_0.88_PCs$Ecotype, levels=c("Crab", "Wave", "Hybrid"))
cluster2398_0.88_var <- round((cluster2398_0.88_PCA$eig/(sum(cluster2398_0.88_PCA$eig)))*100, 2)
cluster2398_0.88_var[1:2]

ggplot(cluster2398_0.88_PCs, aes(x=Axis1, y=Axis2, colour=Ecotype)) + 
  geom_point(alpha=0.4, size=5) + 
  theme_bw() +
  xlab("PC1 (20.9%)") + ylab("PC2 (2.9%)") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(text=element_text(size=14)) +
  scale_colour_manual(labels=c("Crab", "Wave", "Hybrid"), values=c("#F8766D", "#00BFC4", "grey70"))

ggplot(cluster2398_0.88_PCs, aes(x=Axis1, y=Axis2, colour=ecotype)) + 
  geom_point(alpha=0.4, size=5) + 
  theme_bw() +
  xlab("PC1") + ylab("PC2") +
  #guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(text=element_text(size=14)) +
  scale_x_continuous(limits=c(40, 65)) + scale_y_continuous(limits=c(0, 20))



################
### next cluster
### 3414_0.69

#cluster3414_0.69_SNPs <- cluster3414_0.69 %>% pull(., contig_ID)
cluster3414_0.69_SNPs <- SNPslist %>% filter(av>=0 & av<=32.785) %>% pull(., contig_ID)
cluster3414_0.69_genos <- ANG_female %>% select(all_of(cluster3414_0.69_SNPs))

cluster3414_0.69_genind <- df2genind(cluster3414_0.69_genos, ploidy=2, sep="/")

cluster3414_0.69_scaled <- scaleGen(cluster3414_0.69_genind, NA.method="mean")

cluster3414_0.69_PCA <- dudi.pca(cluster3414_0.69_scaled, cent=FALSE, scale=FALSE, scannf=TRUE)

cluster3414_0.69_PCs <- cluster3414_0.69_PCA$li

cluster3414_0.69_PCs <- rownames_to_column(cluster3414_0.69_PCs, var="snail_ID")
cluster3414_0.69_PCs <- cluster3414_0.69_PCs %>% left_join(., ANG_female_info, by="snail_ID") %>%
  rename(Ecotype=ecotype)
cluster3414_0.69_PCs[is.na(cluster3414_0.69_PCs)] <- "Hybrid"
cluster3414_0.69_PCs$Ecotype <- factor(cluster3414_0.69_PCs$Ecotype, levels=c("Crab", "Wave", "Hybrid"))
cluster3414_0.69_var <- round((cluster3414_0.69_PCA$eig/(sum(cluster3414_0.69_PCA$eig)))*100, 2)
cluster3414_0.69_var[1:2]

ggplot(cluster3414_0.69_PCs, aes(x=Axis1, y=Axis2, colour=Ecotype)) + 
  geom_point(alpha=0.4, size=5) + 
  theme_bw() +
  xlab("PC1 (18.0%)") + ylab("PC2 (3.6%)") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(text=element_text(size=14)) +
  scale_colour_manual(labels=c("Crab", "Wave", "Hybrid"), values=c("#F8766D", "#00BFC4", "grey70"))

ggplot(cluster3414_0.69_PCs, aes(x=Axis1, y=Axis2, colour=ecotype)) + 
  geom_point(alpha=0.4, size=5) + 
  theme_bw() +
  xlab("PC1") + ylab("PC2") +
  #guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(text=element_text(size=14)) +
  scale_x_continuous(limits=c(-25, 0)) + scale_y_continuous(limits=c(-20, 20))



################
### next cluster
### 3508_0.67

#cluster3508_0.67_SNPs <- cluster3508_0.67 %>% pull(., contig_ID)
cluster3508_0.67_SNPs <- SNPslist %>% filter(av>=33.904 & av<=43.801) %>% pull(., contig_ID)
cluster3508_0.67_genos <- ANG_female %>% select(all_of(cluster3508_0.67_SNPs))

cluster3508_0.67_genind <- df2genind(cluster3508_0.67_genos, ploidy=2, sep="/")

cluster3508_0.67_scaled <- scaleGen(cluster3508_0.67_genind, NA.method="mean")

cluster3508_0.67_PCA <- dudi.pca(cluster3508_0.67_scaled, cent=FALSE, scale=FALSE, scannf=TRUE)

cluster3508_0.67_PCs <- cluster3508_0.67_PCA$li

cluster3508_0.67_PCs <- rownames_to_column(cluster3508_0.67_PCs, var="snail_ID")
cluster3508_0.67_PCs <- cluster3508_0.67_PCs %>% left_join(., ANG_female_info, by="snail_ID")

ggplot(cluster3508_0.67_PCs, aes(x=Axis1, y=Axis2, colour=ecotype)) + 
  geom_point(alpha=0.4, size=5) + 
  theme_bw() +
  xlab("PC1") + ylab("PC2") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(text=element_text(size=14))

ggplot(cluster3508_0.67_PCs, aes(x=Axis1, y=Axis2, colour=ecotype)) + 
  geom_point(alpha=0.4, size=5) + 
  theme_bw() +
  xlab("PC1") + ylab("PC2") +
  #guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(text=element_text(size=14)) +
  scale_x_continuous(limits=c(-20, 10)) + scale_y_continuous(limits=c(-20, 20))


### mid region / region2 (where three clusters overlap)
midregion_SNPs <- SNPslist %>% filter(av>=32.785 & av<=48.712) %>% pull(., contig_ID)
midregion_genos <- ANG_female %>% select(all_of(midregion_SNPs))

midregion_genind <- df2genind(midregion_genos, ploidy=2, sep="/")

midregion_scaled <- scaleGen(midregion_genind, NA.method="mean")

midregion_PCA <- dudi.pca(midregion_scaled, cent=FALSE, scale=FALSE, scannf=TRUE)

midregion_PCs <- midregion_PCA$li

midregion_var <- round((midregion_PCA$eig/(sum(midregion_PCA$eig)))*100, 2)
midregion_var[1:2]

midregion_PCs <- rownames_to_column(midregion_PCs, var="snail_ID")
midregion_PCs <- midregion_PCs %>% left_join(., ANG_female_info, by="snail_ID")

ggplot(midregion_PCs, aes(x=Axis1, y=Axis2, colour=ecotype)) + 
  geom_point(alpha=0.4, size=5) + 
  theme_bw() +
  xlab("PC1") + ylab("PC2") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(text=element_text(size=14))




#####
#### carry out PCA on the same regions but in males
#####


ANGmale <- read.csv("ANGmale.csv")

# transpose dataframe as adegenet need SNPs in columns
ANGmale <- ANGmale %>% unite(contig_ID, c("CHROM", "POS"), sep="_")
testm <- ANGmale %>% select(-contig_ID)
testmtranspose <- as_tibble(t(testm))
namesm <- colnames(testm)
testmtranspose <- data.frame(namesm, testmtranspose)
names(testmtranspose)

# get list of SNPs to be colnames
SNPslist_names <- SNPslist[,1]
snailsm <- testmtranspose %>% select(namesm)
testmtranspose <- testmtranspose %>% select(-namesm)
colnames(testmtranspose) <- SNPslist_names
ANG_male <- cbind(snailsm, testmtranspose)
ANG_male <- ANG_male %>% rename("snail_ID"="namesm")
# make snail names the rownames - so genind can take them
rownames(ANG_male) <- ANG_male$snail_ID
ANG_male <- ANG_male %>% select(-snail_ID)

# change genotype coding to adegenet format
ANG_male[ANG_male==0] <- "0/0"
ANG_male[ANG_male==1] <- "0/1"
ANG_male[ANG_male==2] <- "1/1"

# make a dataframe of ecotype info for the snails (in case I want to colour the PCAs)
ANG_male_info <- read.csv("ANG_snail_info.csv") %>% filter(sex_dissection=="male") %>% 
  select(-sex_dissection, -DistAlongPath) %>% arrange(snail_ID)

# filter the genotypes to only the ones in the focal cluster
cluster2317_0.89_SNPs <- SNPslist %>% filter(av>=32.98692 & av<=43.801) %>% pull(., contig_ID)
malecluster1_genos <- ANG_male %>% select(all_of(cluster2317_0.89_SNPs))

# make genind object
malecluster1_genind <- df2genind(malecluster1_genos, ploidy=2, sep="/")
summary(malecluster1_genind)

# scale
malecluster1_scaled <- scaleGen(malecluster1_genind, NA.method="mean")

# PCA
malecluster1_PCA <- dudi.pca(malecluster1_scaled, cent=FALSE, scale=FALSE, scannf=TRUE)

# save PCs 
malecluster1_PCs <- malecluster1_PCA$li

# add grouping info to add colour to plot - first make rownames a col for the PCs
malecluster1_PCs <- rownames_to_column(malecluster1_PCs, var="snail_ID")
malecluster1_PCs <- malecluster1_PCs %>% left_join(., ANG_male_info, by="snail_ID")

# plot with colour 
ggplot(malecluster1_PCs, aes(x=Axis1, y=Axis2, colour=ecotype)) + 
  geom_point(alpha=0.4, size=5) + 
  theme_bw() +
  xlab("PC1") + ylab("PC2") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(text=element_text(size=14))

# also a zoom in of the clustered section
ggplot(malecluster1_PCs, aes(x=Axis1, y=Axis2, colour=ecotype)) + 
  geom_point(alpha=0.4, size=5) + 
  theme_bw() +
  xlab("PC1") + ylab("PC2") +
  #guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(text=element_text(size=14)) +
  scale_x_continuous(limits=c(0, 10)) + scale_y_continuous(limits=c(-5, 5))



################
### next cluster
### 2318_0.89

cluster2318_0.89_SNPs <- SNPslist %>% filter(av>=33.624 & av<=48.712) %>% pull(., contig_ID)
malecluster2_genos <- ANG_male %>% select(all_of(cluster2318_0.89_SNPs))

malecluster2_genind <- df2genind(malecluster2_genos, ploidy=2, sep="/")

malecluster2_scaled <- scaleGen(malecluster2_genind, NA.method="mean")

malecluster2_PCA <- dudi.pca(malecluster2_scaled, cent=FALSE, scale=FALSE, scannf=TRUE)

malecluster2_PCs <- malecluster2_PCA$li

malecluster2_PCs <- rownames_to_column(malecluster2_PCs, var="snail_ID")
malecluster2_PCs <- malecluster2_PCs %>% left_join(., ANG_male_info, by="snail_ID")

ggplot(malecluster2_PCs, aes(x=Axis1, y=Axis2, colour=ecotype)) + 
  geom_point(alpha=0.4, size=5) + 
  theme_bw() +
  xlab("PC1") + ylab("PC2") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(text=element_text(size=14))

ggplot(malecluster2_PCs, aes(x=Axis1, y=Axis2, colour=ecotype)) + 
  geom_point(alpha=0.4, size=5) + 
  theme_bw() +
  xlab("PC1") + ylab("PC2") +
  #guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(text=element_text(size=14)) +
  scale_x_continuous(limits=c(-10, 10)) + scale_y_continuous(limits=c(-10, 10))



################
### next cluster
### 2398_0.88

cluster2398_0.88_SNPs <- SNPslist %>% filter(av>=48.712 & av<=60.235) %>% pull(., contig_ID)
malecluster3_genos <- ANG_male %>% select(all_of(cluster2398_0.88_SNPs))

malecluster3_genind <- df2genind(malecluster3_genos, ploidy=2, sep="/")

malecluster3_scaled <- scaleGen(malecluster3_genind, NA.method="mean")

malecluster3_PCA <- dudi.pca(malecluster3_scaled, cent=FALSE, scale=FALSE, scannf=TRUE)

malecluster3_PCs <- malecluster3_PCA$li


malecluster3_PCs <- rownames_to_column(malecluster3_PCs, var="snail_ID")
malecluster3_PCs <- malecluster3_PCs %>% left_join(., ANG_male_info, by="snail_ID") %>%
  rename(Ecotype=ecotype)
malecluster3_PCs[is.na(malecluster3_PCs)] <- "Hybrid"
malecluster3_PCs$Ecotype <- factor(malecluster3_PCs$Ecotype, levels=c("Crab", "Wave", "Hybrid"))
malecluster3_var <- round((malecluster3_PCA$eig/(sum(malecluster3_PCA$eig)))*100, 2)
malecluster3_var[1:2]

ggplot(malecluster3_PCs, aes(x=Axis1, y=Axis2, colour=Ecotype)) + 
  geom_point(alpha=0.4, size=5) + 
  theme_bw() +
  xlab("PC1 (19.1%)") + ylab("PC2 (3.5%)") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(text=element_text(size=14)) +
  scale_colour_manual(labels=c("Crab", "Wave", "Hybrid"), values=c("#F8766D", "#00BFC4", "grey70"))




################
### next cluster
### 3414_0.69

cluster3414_0.69_SNPs <- SNPslist %>% filter(av>=0 & av<=32.785) %>% pull(., contig_ID)
malecluster4_genos <- ANG_male %>% select(all_of(cluster3414_0.69_SNPs))

malecluster4_genind <- df2genind(malecluster4_genos, ploidy=2, sep="/")

malecluster4_scaled <- scaleGen(malecluster4_genind, NA.method="mean")

malecluster4_PCA <- dudi.pca(malecluster4_scaled, cent=FALSE, scale=FALSE, scannf=TRUE)

malecluster4_PCs <- malecluster4_PCA$li

malecluster4_PCs <- rownames_to_column(malecluster4_PCs, var="snail_ID")
malecluster4_PCs <- malecluster4_PCs %>% left_join(., ANG_male_info, by="snail_ID") %>%
  rename(Ecotype=ecotype)
malecluster4_PCs[is.na(malecluster4_PCs)] <- "Hybrid"
malecluster4_PCs$Ecotype <- factor(malecluster4_PCs$Ecotype, levels=c("Crab", "Wave", "Hybrid"))
malecluster4_var <- round((malecluster4_PCA$eig/(sum(malecluster4_PCA$eig)))*100, 2)
malecluster4_var[1:2]

ggplot(malecluster4_PCs, aes(x=Axis1, y=Axis2, colour=Ecotype)) + 
  geom_point(alpha=0.4, size=5) + 
  theme_bw() +
  xlab("PC1 (17.5%)") + ylab("PC2 (3.0%)") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(text=element_text(size=14)) +
  scale_colour_manual(labels=c("Crab", "Wave", "Hybrid"), values=c("#F8766D", "#00BFC4", "grey70"))

ggplot(malecluster4_PCs, aes(x=Axis1, y=Axis2, colour=ecotype)) + 
  geom_point(alpha=0.4, size=5) + 
  theme_bw() +
  xlab("PC1") + ylab("PC2") +
  #guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(text=element_text(size=14)) +
  scale_x_continuous(limits=c(-12, -2)) + scale_y_continuous(limits=c(-3, 3))


################
### next cluster
### 3508_0.67

cluster3508_0.67_SNPs <- SNPslist %>% filter(av>=33.904 & av<=43.801) %>% pull(., contig_ID)
malecluster5_genos <- ANG_male %>% select(all_of(cluster3508_0.67_SNPs))

malecluster5_genind <- df2genind(malecluster5_genos, ploidy=2, sep="/")

malecluster5_scaled <- scaleGen(malecluster5_genind, NA.method="mean")

malecluster5_PCA <- dudi.pca(malecluster5_scaled, cent=FALSE, scale=FALSE, scannf=TRUE)

malecluster5_PCs <- malecluster5_PCA$li


malecluster5_PCs <- rownames_to_column(malecluster5_PCs, var="snail_ID")
malecluster5_PCs <- malecluster5_PCs %>% left_join(., ANG_male_info, by="snail_ID")

ggplot(malecluster5_PCs, aes(x=Axis1, y=Axis2, colour=ecotype)) + 
  geom_point(alpha=0.4, size=5) + 
  theme_bw() +
  xlab("PC1") + ylab("PC2") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(text=element_text(size=14))

ggplot(malecluster5_PCs, aes(x=Axis1, y=Axis2, colour=ecotype)) + 
  geom_point(alpha=0.4, size=5) + 
  theme_bw() +
  xlab("PC1") + ylab("PC2") +
  #guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(text=element_text(size=14)) +
  scale_x_continuous(limits=c(0, 4)) + scale_y_continuous(limits=c(-3, 5))


#### mid reion / region2
midregion_SNPs <- SNPslist %>% filter(av>=32.785 & av<=48.712) %>% pull(., contig_ID)
midregion_male_genos <- ANG_male %>% select(all_of(midregion_SNPs))

midregion_male_genind <- df2genind(midregion_male_genos, ploidy=2, sep="/")

midregion_male_scaled <- scaleGen(midregion_male_genind, NA.method="mean")

midregion_male_PCA <- dudi.pca(midregion_male_scaled, cent=FALSE, scale=FALSE, scannf=TRUE)

midregion_male_PCs <- midregion_male_PCA$li
midregion_male_var <- round((midregion_male_PCA$eig/(sum(midregion_male_PCA$eig)))*100, 2)
midregion_male_var[1:2]

midregion_male_PCs <- rownames_to_column(midregion_male_PCs, var="snail_ID")
midregion_male_PCs <- midregion_male_PCs %>% left_join(., ANG_male_info, by="snail_ID")

ggplot(midregion_male_PCs, aes(x=Axis1, y=Axis2, colour=ecotype)) + 
  geom_point(alpha=0.4, size=5) + 
  theme_bw() +
  xlab("PC1") + ylab("PC2") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(text=element_text(size=14))



#
#
#
##### PCAs of the same regions but for both sexes together
#
#
#

ANG_snails <- read.csv("ANG_genotype_table.csv")

# transpose dataframe as adegenet need SNPs in columns
ANG_snails <- ANG_snails %>% rename("contig_ID"="snail_ID")
test <- ANG_snails %>% select(-contig_ID)
testtranspose <- as_tibble(t(test))
names <- colnames(test)
testtranspose <- data.frame(names, testtranspose)
names(testtranspose)

# get list of SNPs to be colnames
SNPslist_names <- SNPslist[,1]
snails <- testtranspose %>% select(names)
testtranspose <- testtranspose %>% select(-names)
colnames(testtranspose) <- SNPslist_names
ANG_snails <- cbind(snails, testtranspose)
ANG_snails <- ANG_snails %>% rename("snail_ID"="names")
# make snail names the rownames - so genind can take them
#rownames(ANG_snails) <- ANG_snails$snail_ID
#ANG_snails <- ANG_snails %>% select(-snail_ID)

# change genotype coding to adegenet format
ANG_snails[ANG_snails==0] <- "0/0"
ANG_snails[ANG_snails==1] <- "0/1"
ANG_snails[ANG_snails==2] <- "1/1"


# make a dataframe of ecotype/sex info for the snails to use to colour the PCAs 
dissect2 <- read.csv("ANG_dissections_20150723.csv")
dissect2 <- dissect2 %>% filter(sex_dissection=="female" | sex_dissection=="male") %>%
  select(snail_ID, sex_dissection)
path <- read.csv("ANG_SnailDataLitPath_20160702.csv")
path <- path %>% select(snail_ID, DistAlongPath)

ANG_snail_info_new <- full_join(dissect2, path, by="snail_ID")

ANG_snail_info_new <- ANG_snail_info_new %>% 
  mutate(ecotype=ifelse(DistAlongPath>88, "Wave",
                 ifelse(DistAlongPath<68, "Crab", "Hybrid"))) %>% arrange(snail_ID)
ANG_snail_info_new <- ANG_snail_info_new %>% 
  mutate(group=ifelse(sex_dissection=="male" & ecotype=="Wave", "W_M",
                      ifelse(sex_dissection=="male" & ecotype=="Crab", "C_M",
                             ifelse(sex_dissection=="female" & ecotype=="Wave", "W_F",
                                    ifelse(sex_dissection=="female" & ecotype=="Crab", "C_F", NA)))))
ANG_snail_info <- ANG_snail_info %>% 
  mutate(group=ifelse(sex_dissection=="male" & ecotype=="Wave", "W_M",
                      ifelse(sex_dissection=="male" & ecotype=="Crab", "C_M",
                             ifelse(sex_dissection=="female" & ecotype=="Wave", "W_F",
                                    ifelse(sex_dissection=="female" & ecotype=="Crab", "C_F", NA)))))


### cluster 2318_0.89 - cluster 2
### will do all and part region, as before
# filter the genotypes to only the ones in the focal cluster region
cluster2318_0.89_SNPs <- SNPslist %>% filter(av>=33.624 & av<=48.712) %>% pull(., contig_ID)
cluster2318_0.89_allsnails_genos <- ANG_snails %>% select(snail_ID, all_of(cluster2318_0.89_SNPs))
cluster2318_0.89_allsnails_genos <- cluster2318_0.89_allsnails_genos %>% slice(-c(381:382))
rownames(cluster2318_0.89_allsnails_genos) <- cluster2318_0.89_allsnails_genos$snail_ID
cluster2318_0.89_allsnails_genos <- cluster2318_0.89_allsnails_genos %>% select(-snail_ID)

# make genind object
cluster2318_0.89_allsnails_genind <- df2genind(cluster2318_0.89_allsnails_genos, ploidy=2, sep="/")
#summary(cluster2317_0.89_allsnails_genind)

# scale
cluster2318_0.89_allsnails_scaled <- scaleGen(cluster2318_0.89_allsnails_genind, NA.method="mean")

# PCA
cluster2318_0.89_allsnails_PCA <- dudi.pca(cluster2318_0.89_allsnails_scaled, cent=FALSE, scale=FALSE, scannf=TRUE)

# save PCs 
cluster2318_0.89_allsnails_PCs <- cluster2318_0.89_allsnails_PCA$li


# add grouping info to add colour to plot - first make rownames a col for the PCs
cluster2318_0.89_allsnails_PCs <- rownames_to_column(cluster2318_0.89_allsnails_PCs, var="snail_ID")
cluster2318_0.89_allsnails_PCs <- cluster2318_0.89_allsnails_PCs %>% 
  left_join(., ANG_snail_info_new, by="snail_ID")

# plot with colour 
ggplot(cluster2318_0.89_allsnails_PCs, aes(x=Axis1, y=Axis2, colour=group)) + 
  geom_point(alpha=0.4, size=5) + 
  theme_bw() +
  xlab("PC1") + ylab("PC2") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(text=element_text(size=14))

# also a zoom in of the clustered section
#ggplot(cluster2318_0.89_allsnails_PCs, aes(x=Axis1, y=Axis2, colour=ecotype)) + 
#  geom_point(alpha=0.4, size=5) + 
#  theme_bw() +
#  xlab("PC1") + ylab("PC2") +
#  #guides(colour = guide_legend(override.aes = list(alpha = 1))) +
#  theme(text=element_text(size=14)) +
#  scale_x_continuous(limits=c(-1, 1)) + scale_y_continuous(limits=c(-1, 1))


#### cluster 2317_0.89 (cluster1)
# filter the genotypes to only the ones in the focal cluster region
cluster2317_0.89_SNPs <- SNPslist %>% filter(av>=32.98692 & av<=43.801) %>% pull(., contig_ID)
cluster2317_0.89_allsnails_genos <- ANG_snails %>% select(snail_ID, all_of(cluster2317_0.89_SNPs))
cluster2317_0.89_allsnails_genos <- cluster2317_0.89_allsnails_genos %>% slice(-c(381:382))
rownames(cluster2317_0.89_allsnails_genos) <- cluster2317_0.89_allsnails_genos$snail_ID
cluster2317_0.89_allsnails_genos <- cluster2317_0.89_allsnails_genos %>% select(-snail_ID)

# make genind object
cluster2317_0.89_allsnails_genind <- df2genind(cluster2317_0.89_allsnails_genos, ploidy=2, sep="/")
summary(cluster2317_0.89_allsnails_genind)

# scale
cluster2317_0.89_allsnails_scaled <- scaleGen(cluster2317_0.89_allsnails_genind, NA.method="mean")

# PCA
cluster2317_0.89_allsnails_PCA <- dudi.pca(cluster2317_0.89_allsnails_scaled, cent=FALSE, scale=FALSE, scannf=TRUE)

# save PCs 
cluster2317_0.89_allsnails_PCs <- cluster2317_0.89_allsnails_PCA$li


# add grouping info to add colour to plot - first make rownames a col for the PCs
cluster2317_0.89_allsnails_PCs <- rownames_to_column(cluster2317_0.89_allsnails_PCs, var="snail_ID")
cluster2317_0.89_allsnails_PCs <- cluster2317_0.89_allsnails_PCs %>% 
  left_join(., ANG_snail_info_new, by="snail_ID")

# plot with colour 
ggplot(cluster2317_0.89_allsnails_PCs, aes(x=Axis1, y=Axis2, colour=group)) + 
  geom_point(alpha=0.4, size=5) + 
  theme_bw() +
  xlab("PC1") + ylab("PC2") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(text=element_text(size=14))

# also a zoom in of the clustered section
#ggplot(cluster2317_0.89_allsnails_PCs, aes(x=Axis1, y=Axis2, colour=ecotype)) + 
#  geom_point(alpha=0.4, size=5) + 
#  theme_bw() +
#  xlab("PC1") + ylab("PC2") +
#  #guides(colour = guide_legend(override.aes = list(alpha = 1))) +
#  theme(text=element_text(size=14)) +
#  scale_x_continuous(limits=c(-1, 1)) + scale_y_continuous(limits=c(-1, 1))


### cluster 2398_0.88
### will do all and part region, as before
# filter the genotypes to only the ones in the focal cluster region
cluster2398_0.88_SNPs <- SNPslist %>% filter(av>=48.712 & av<=60.235) %>% pull(., contig_ID)
cluster2398_0.88_allsnails_genos <- ANG_snails %>% select(snail_ID, all_of(cluster2398_0.88_SNPs))
cluster2398_0.88_allsnails_genos <- cluster2398_0.88_allsnails_genos %>% slice(-c(381:382))
rownames(cluster2398_0.88_allsnails_genos) <- cluster2398_0.88_allsnails_genos$snail_ID
cluster2398_0.88_allsnails_genos <- cluster2398_0.88_allsnails_genos %>% select(-snail_ID)

# make genind object
cluster2398_0.88_allsnails_genind <- df2genind(cluster2398_0.88_allsnails_genos, ploidy=2, sep="/")
#summary(cluster2317_0.89_allsnails_genind)

# scale
cluster2398_0.88_allsnails_scaled <- scaleGen(cluster2398_0.88_allsnails_genind, NA.method="mean")

# PCA
cluster2398_0.88_allsnails_PCA <- dudi.pca(cluster2398_0.88_allsnails_scaled, cent=FALSE, scale=FALSE, scannf=TRUE)

# save PCs 
cluster2398_0.88_allsnails_PCs <- cluster2398_0.88_allsnails_PCA$li


# add grouping info to add colour to plot - first make rownames a col for the PCs
cluster2398_0.88_allsnails_PCs <- rownames_to_column(cluster2398_0.88_allsnails_PCs, var="snail_ID")
cluster2398_0.88_allsnails_PCs <- cluster2398_0.88_allsnails_PCs %>% 
  left_join(., ANG_snail_info_new, by="snail_ID")
cluster2398_0.88_allsnails_var <- round((cluster2398_0.88_allsnails_PCA$eig/(sum(cluster2398_0.88_allsnails_PCA$eig)))*100, 2)
cluster2398_0.88_allsnails_var[1:2]

# plot with colour 
ggplot(cluster2398_0.88_allsnails_PCs, aes(x=Axis1, y=Axis2, colour=group, fill=group, shape=sex_dissection)) + 
  geom_point(alpha=0.4, size=4) + 
  scale_shape_manual(values=c(23, 19)) +
  theme_bw() +
  xlab("PC1 (20.3%)") + ylab("PC2 (3.0%)") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(text=element_text(size=14)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# also a zoom in of the clustered section
#ggplot(cluster2398_0.88_allsnails_PCs, aes(x=Axis1, y=Axis2, colour=group)) + 
#  geom_point(alpha=0.4, size=5) + 
#  theme_bw() +
#  xlab("PC1") + ylab("PC2") +
#  #guides(colour = guide_legend(override.aes = list(alpha = 1))) +
#  theme(text=element_text(size=14)) +
#  scale_x_continuous(limits=c(-1, 1)) + scale_y_continuous(limits=c(-1, 1))



### cluster 3414_0.69
### will do all and part region, as before
# filter the genotypes to only the ones in the focal cluster region
cluster3414_0.69_SNPs <- SNPslist %>% filter(av>=0 & av<=32.785) %>% pull(., contig_ID)
cluster3414_0.69_allsnails_genos <- ANG_snails %>% select(snail_ID, all_of(cluster3414_0.69_SNPs))
cluster3414_0.69_allsnails_genos <- cluster3414_0.69_allsnails_genos %>% slice(-c(381:382))
rownames(cluster3414_0.69_allsnails_genos) <- cluster3414_0.69_allsnails_genos$snail_ID
cluster3414_0.69_allsnails_genos <- cluster3414_0.69_allsnails_genos %>% select(-snail_ID)

# make genind object
cluster3414_0.69_allsnails_genind <- df2genind(cluster3414_0.69_allsnails_genos, ploidy=2, sep="/")
#summary(cluster3414_0.69_allsnails_genind)

# scale
cluster3414_0.69_allsnails_scaled <- scaleGen(cluster3414_0.69_allsnails_genind, NA.method="mean")

# PCA
cluster3414_0.69_allsnails_PCA <- dudi.pca(cluster3414_0.69_allsnails_scaled, cent=FALSE, scale=FALSE, scannf=TRUE)

# save PCs 
cluster3414_0.69_allsnails_PCs <- cluster3414_0.69_allsnails_PCA$li


# add grouping info to add colour to plot - first make rownames a col for the PCs
cluster3414_0.69_allsnails_PCs <- rownames_to_column(cluster3414_0.69_allsnails_PCs, var="snail_ID")
cluster3414_0.69_allsnails_PCs <- cluster3414_0.69_allsnails_PCs %>% 
  left_join(., ANG_snail_info_new, by="snail_ID")
cluster3414_0.69_allsnails_var <- round((cluster3414_0.69_allsnails_PCA$eig/(sum(cluster3414_0.69_allsnails_PCA$eig)))*100, 2)
cluster3414_0.69_allsnails_var[1:2]

# plot with colour 
ggplot(cluster3414_0.69_allsnails_PCs, aes(x=Axis1, y=Axis2, colour=group, fill=group, shape=sex_dissection)) + 
  geom_point(alpha=0.4, size=4) + 
  scale_shape_manual(values=c(23, 19)) +
  theme_bw() +
  xlab("PC1 (17.6%)") + ylab("PC2 (3.1%)") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(text=element_text(size=14)) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())


# also a zoom in of the clustered section
#ggplot(cluster3414_0.69_allsnails_PCs, aes(x=Axis1, y=Axis2, colour=group)) + 
#  geom_point(alpha=0.4, size=5) + 
#  theme_bw() +
#  xlab("PC1") + ylab("PC2") +
#  #guides(colour = guide_legend(override.aes = list(alpha = 1))) +
#  theme(text=element_text(size=14)) +
#  scale_x_continuous(limits=c(-1, 1)) + scale_y_continuous(limits=c(-1, 1))



### cluster 3508_0.67
### will do all and part region, as before
# filter the genotypes to only the ones in the focal cluster region
cluster3508_0.67_SNPs <- SNPslist %>% filter(av>=33.904 & av<=43.801) %>% pull(., contig_ID)
cluster3508_0.67_allsnails_genos <- ANG_snails %>% select(snail_ID, all_of(cluster3508_0.67_SNPs))
cluster3508_0.67_allsnails_genos <- cluster3508_0.67_allsnails_genos %>% slice(-c(381:382))
rownames(cluster3508_0.67_allsnails_genos) <- cluster3508_0.67_allsnails_genos$snail_ID
cluster3508_0.67_allsnails_genos <- cluster3508_0.67_allsnails_genos %>% select(-snail_ID)

# make genind object
cluster3508_0.67_allsnails_genind <- df2genind(cluster3508_0.67_allsnails_genos, ploidy=2, sep="/")
#summary(cluster3508_0.67_allsnails_genind)

# scale
cluster3508_0.67_allsnails_scaled <- scaleGen(cluster3508_0.67_allsnails_genind, NA.method="mean")

# PCA
cluster3508_0.67_allsnails_PCA <- dudi.pca(cluster3508_0.67_allsnails_scaled, cent=FALSE, scale=FALSE, scannf=TRUE)

# save PCs 
cluster3508_0.67_allsnails_PCs <- cluster3508_0.67_allsnails_PCA$li

# add grouping info to add colour to plot - first make rownames a col for the PCs
cluster3508_0.67_allsnails_PCs <- rownames_to_column(cluster3508_0.67_allsnails_PCs, var="snail_ID")
cluster3508_0.67_allsnails_PCs <- cluster3508_0.67_allsnails_PCs %>% 
  left_join(., ANG_snail_info_new, by="snail_ID")

# plot with colour 
ggplot(cluster3508_0.67_allsnails_PCs, aes(x=Axis1, y=Axis2, colour=group)) + 
  geom_point(alpha=0.4, size=5) + 
  theme_bw() +
  xlab("PC1") + ylab("PC2") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(text=element_text(size=14))

# also a zoom in of the clustered section
#ggplot(cluster3508_0.67_allsnails_PCs, aes(x=Axis1, y=Axis2, colour=group)) + 
#  geom_point(alpha=0.4, size=5) + 
#  theme_bw() +
#  xlab("PC1") + ylab("PC2") +
#  #guides(colour = guide_legend(override.aes = list(alpha = 1))) +
#  theme(text=element_text(size=14)) +
#  scale_x_continuous(limits=c(-1, 1)) + scale_y_continuous(limits=c(-1, 1))


###
## mid region (where three clusters overlap)
midregion_SNPs <- SNPslist %>% filter(av>=32.785 & av<=48.801) %>% pull(., contig_ID)
midregion_allsnails_genos <- ANG_snails %>% select(snail_ID, all_of(midregion_SNPs))
midregion_allsnails_genos <- midregion_allsnails_genos %>% slice(-c(381:382))
rownames(midregion_allsnails_genos) <- midregion_allsnails_genos$snail_ID
midregion_allsnails_genos <- midregion_allsnails_genos %>% select(-snail_ID)

# make genind object
midregion_allsnails_genind <- df2genind(midregion_allsnails_genos, ploidy=2, sep="/")
#summary(cluster3508_0.67_allsnails_genind)

# scale
midregion_allsnails_scaled <- scaleGen(midregion_allsnails_genind, NA.method="mean")

# PCA
midregion_allsnails_PCA <- dudi.pca(midregion_allsnails_scaled, cent=FALSE, scale=FALSE, scannf=TRUE)

# save PCs 
midregion_allsnails_PCs <- midregion_allsnails_PCA$li

# add grouping info to add colour to plot - first make rownames a col for the PCs
midregion_allsnails_PCs <- rownames_to_column(midregion_allsnails_PCs, var="snail_ID")
midregion_allsnails_PCs <- midregion_allsnails_PCs %>% 
  left_join(., ANG_snail_info_new, by="snail_ID")

# plot with colour 
ggplot(midregion_allsnails_PCs, aes(x=Axis1, y=Axis2, colour=group)) + 
  geom_point(alpha=0.4, size=5) + 
  theme_bw() +
  xlab("PC1") + ylab("PC2") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(text=element_text(size=14))




############# check regions 2a and 2b - split middle region at 43.801cM

## region 2a females
region2a_SNPs <- SNPslist %>% filter(av>=32.98 & av<=43.801) %>% pull(., contig_ID)
region2a_genos <- ANG_female %>% select(all_of(region2a_SNPs))

region2a_genind <- df2genind(region2a_genos, ploidy=2, sep="/")

region2a_scaled <- scaleGen(region2a_genind, NA.method="mean")

region2a_PCA <- dudi.pca(region2a_scaled, cent=FALSE, scale=FALSE, scannf=TRUE)

region2a_PCs <- region2a_PCA$li

region2a_PCs <- rownames_to_column(region2a_PCs, var="snail_ID")
region2a_PCs <- region2a_PCs %>% left_join(., ANG_female_info, by="snail_ID") %>%
  rename(Ecotype=ecotype)
region2a_PCs[is.na(region2a_PCs)] <- "Hybrid"
region2a_PCs$Ecotype <- factor(region2a_PCs$Ecotype, levels=c("Crab", "Wave", "Hybrid"))
region2a_var <- round((region2a_PCA$eig/(sum(region2a_PCA$eig)))*100, 2)
region2a_var[1:2]

ggplot(region2a_PCs, aes(x=Axis1, y=Axis2, colour=Ecotype)) + 
  geom_point(alpha=0.4, size=5) + 
  theme_bw() +
  xlab("PC1 (13.7%)") + ylab("PC2 (9.2%)") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(text=element_text(size=14)) +
  scale_colour_manual(labels=c("Crab", "Wave", "Hybrid"), values=c("#F8766D", "#00BFC4", "grey70"))


## region 2a males
region2a_SNPs <- SNPslist %>% filter(av>=32.98 & av<=43.801) %>% pull(., contig_ID)
region2a_m_genos <- ANG_male %>% select(all_of(region2a_SNPs))

region2a_m_genind <- df2genind(region2a_m_genos, ploidy=2, sep="/")

region2a_m_scaled <- scaleGen(region2a_m_genind, NA.method="mean")

region2a_m_PCA <- dudi.pca(region2a_m_scaled, cent=FALSE, scale=FALSE, scannf=TRUE)

region2a_m_PCs <- region2a_m_PCA$li

region2a_m_PCs <- rownames_to_column(region2a_m_PCs, var="snail_ID")
region2a_m_PCs <- region2a_m_PCs %>% left_join(., ANG_male_info, by="snail_ID") %>%
  rename(Ecotype=ecotype)
region2a_m_PCs[is.na(region2a_m_PCs)] <- "Hybrid"
region2a_m_PCs$Ecotype <- factor(region2a_m_PCs$Ecotype, levels=c("Crab", "Wave", "Hybrid"))
region2a_m_var <- round((region2a_m_PCA$eig/(sum(region2a_m_PCA$eig)))*100, 2)
region2a_m_var[1:2]

ggplot(region2a_m_PCs, aes(x=Axis1, y=Axis2, colour=Ecotype)) + 
  geom_point(alpha=0.4, size=5) + 
  theme_bw() +
  xlab("PC1 (17.7%)") + ylab("PC2 (4.5%)") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(text=element_text(size=14)) +
  scale_colour_manual(labels=c("Crab", "Wave", "Hybrid"), values=c("#F8766D", "#00BFC4", "grey70"))

## region 2a all snails
region2a_SNPs <- SNPslist %>% filter(av>=32.98 & av<=43.801) %>% pull(., contig_ID)
region2a_allsnails_genos <- ANG_snails %>% select(snail_ID, all_of(region2a_SNPs))
region2a_allsnails_genos <- region2a_allsnails_genos %>% slice(-c(381:382))
rownames(region2a_allsnails_genos) <- region2a_allsnails_genos$snail_ID
region2a_allsnails_genos <- region2a_allsnails_genos %>% select(-snail_ID)

region2a_allsnails_genind <- df2genind(region2a_allsnails_genos, ploidy=2, sep="/")

region2a_allsnails_scaled <- scaleGen(region2a_allsnails_genind, NA.method="mean")

region2a_allsnails_PCA <- dudi.pca(region2a_allsnails_scaled, cent=FALSE, scale=FALSE, scannf=TRUE)

region2a_allsnails_PCs <- region2a_allsnails_PCA$li

region2a_allsnails_PCs <- rownames_to_column(region2a_allsnails_PCs, var="snail_ID")
region2a_allsnails_PCs <- region2a_allsnails_PCs %>% left_join(., ANG_snail_info_new, by="snail_ID")
region2a_allsnails_var <- round((region2a_allsnails_PCA$eig/(sum(region2a_allsnails_PCA$eig)))*100, 2)
region2a_allsnails_var[1:2]

ggplot(region2a_allsnails_PCs, aes(x=Axis1, y=Axis2, colour=group, fill=group, shape=sex_dissection)) + 
  geom_point(alpha=0.4, size=4) + 
  scale_shape_manual(values=c(23, 19)) +
  theme_bw() +
  xlab("PC1 (18.0%)") + ylab("PC2 (6.7%)") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(text=element_text(size=14)) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

## region 2b females
region2b_SNPs <- SNPslist %>% filter(av>=43.801 & av<=48.712) %>% pull(., contig_ID)
region2b_genos <- ANG_female %>% select(all_of(region2b_SNPs))

region2b_genind <- df2genind(region2b_genos, ploidy=2, sep="/")

region2b_scaled <- scaleGen(region2b_genind, NA.method="mean")

region2b_PCA <- dudi.pca(region2b_scaled, cent=FALSE, scale=FALSE, scannf=TRUE)

region2b_PCs <- region2b_PCA$li

region2b_PCs <- rownames_to_column(region2b_PCs, var="snail_ID")
region2b_PCs <- region2b_PCs %>% left_join(., ANG_female_info, by="snail_ID") %>%
  rename(Ecotype=ecotype)
region2b_PCs[is.na(region2b_PCs)] <- "Hybrid"
region2b_PCs$Ecotype <- factor(region2b_PCs$Ecotype, levels=c("Crab", "Wave", "Hybrid"))
region2b_var <- round((region2b_PCA$eig/(sum(region2b_PCA$eig)))*100, 2)
region2b_var[1:2]

ggplot(region2b_PCs, aes(x=Axis1, y=Axis2, colour=Ecotype)) + 
  geom_point(alpha=0.4, size=5) + 
  theme_bw() +
  xlab("PC1 (12.3%)") + ylab("PC2 (3.2%)") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(text=element_text(size=14)) +
  scale_colour_manual(labels=c("Crab", "Wave", "Hybrid"), values=c("#F8766D", "#00BFC4", "grey70"))


## region 2b males
region2b_SNPs <- SNPslist %>% filter(av>=43.801 & av<=48.712) %>% pull(., contig_ID)
region2b_m_genos <- ANG_male %>% select(all_of(region2b_SNPs))

region2b_m_genind <- df2genind(region2b_m_genos, ploidy=2, sep="/")

region2b_m_scaled <- scaleGen(region2b_m_genind, NA.method="mean")

region2b_m_PCA <- dudi.pca(region2b_m_scaled, cent=FALSE, scale=FALSE, scannf=TRUE)

region2b_m_PCs <- region2b_m_PCA$li

region2b_m_PCs <- rownames_to_column(region2b_m_PCs, var="snail_ID")
region2b_m_PCs <- region2b_m_PCs %>% left_join(., ANG_male_info, by="snail_ID") %>%
  rename(Ecotype=ecotype)
region2b_m_PCs[is.na(region2b_m_PCs)] <- "Hybrid"
region2b_m_PCs$Ecotype <- factor(region2b_m_PCs$Ecotype, levels=c("Crab", "Wave", "Hybrid"))
region2b_m_var <- round((region2b_m_PCA$eig/(sum(region2b_m_PCA$eig)))*100, 2)
region2b_m_var[1:2]

ggplot(region2b_m_PCs, aes(x=Axis1, y=Axis2, colour=Ecotype)) + 
  geom_point(alpha=0.4, size=5) + 
  theme_bw() +
  xlab("PC1 (7.8%)") + ylab("PC2 (5.4%)") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(text=element_text(size=14)) +
  scale_colour_manual(labels=c("Crab", "Wave", "Hybrid"), values=c("#F8766D", "#00BFC4", "grey70"))

## region 2b all snails
region2b_SNPs <- SNPslist %>% filter(av>=43.801 & av<=48.712) %>% pull(., contig_ID)
region2b_allsnails_genos <- ANG_snails %>% select(snail_ID, all_of(region2b_SNPs))
region2b_allsnails_genos <- region2b_allsnails_genos %>% slice(-c(381:382))
rownames(region2b_allsnails_genos) <- region2b_allsnails_genos$snail_ID
region2b_allsnails_genos <- region2b_allsnails_genos %>% select(-snail_ID)

region2b_allsnails_genind <- df2genind(region2b_allsnails_genos, ploidy=2, sep="/")

region2b_allsnails_scaled <- scaleGen(region2b_allsnails_genind, NA.method="mean")

region2b_allsnails_PCA <- dudi.pca(region2b_allsnails_scaled, cent=FALSE, scale=FALSE, scannf=TRUE)

region2b_allsnails_PCs <- region2b_allsnails_PCA$li

region2b_allsnails_PCs <- rownames_to_column(region2b_allsnails_PCs, var="snail_ID")
region2b_allsnails_PCs <- region2b_allsnails_PCs %>% left_join(., ANG_snail_info_new, by="snail_ID")
region2b_allsnails_var <- round((region2b_allsnails_PCA$eig/(sum(region2b_allsnails_PCA$eig)))*100, 2)
region2b_allsnails_var[1:2]

ggplot(region2b_allsnails_PCs, aes(x=Axis1, y=Axis2, colour=group, fill=group, shape=sex_dissection)) + 
  geom_point(alpha=0.4, size=4) + 
  scale_shape_manual(values=c(23, 19)) +
  theme_bw() +
  xlab("PC1 (12.6%)") + ylab("PC2 (4.6%)") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(text=element_text(size=14)) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())


