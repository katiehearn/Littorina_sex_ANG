rm(list=ls())
# set working directory

library(genetics)
library(dplyr)
library(LDna)


# reads ANG genotypes, gets a random sample of SNPs (1 per contig),
# then calculates LD matrix using the <genetics> package


# read genotype file
geno <- read.csv("ANGfemale.csv", header=T, stringsAsFactors = F)
# make SNP ID the rowname
geno <- geno %>% unite(SNP, c("CHROM", "POS"), sep="_")
rownames(geno) <- geno$SNP

# read snail data
snail <- read.csv("ANG_SnailDataLitPath_20160702.csv")
snail <- snail[is.na(snail$DistAlongPath)==F,]
snail <- filter(snail, snail_ID %in% colnames(geno)) # keep only snails with genotype information
snails.nm <- c(as.character(snail$snail_ID))
geno <- select(geno, one_of(snails.nm)) # keep only genotypes with snail information
#make SNP names 
geno$SNP <- rownames(geno)

# get genotypes in columns in <genetics> package format
# NB does not use fit information so alleles not identified
# as crab and wave, and LD is not directional
snail.gen <- data.frame(colnames(geno)[1:(ncol(geno)-1)])
colnames(snail.gen) <- "snail_ID"

for (snp in 1:length(geno$SNP)){
  temp <- as.character(geno[snp,1:(ncol(geno)-1)])
  temp[temp=="0"]<-"AA"
  temp[temp=="1"]<-"AB"
  temp[temp=="2"]<-"BB"
  temp[temp=="NA"]<-NA
  
  snail.gen[,snp+1] <- temp
  colnames(snail.gen)[snp+1]<-geno$SNP[snp]
}


# remove genotype columns where there are fewer than 500 genotypes
# a check; shouldn't remove anything as should have previously been filtered
###################################################

gn <- 2
for (g in 2:length(snail.gen)){
  if (sum(is.na(snail.gen[,g]))<500){
    snail.gen[,gn] <- snail.gen[,g]
    gn <- gn+1
  }
}
gn <- gn-1  # Marina's correction!
snail.gen <- snail.gen[,1:gn]



###########################################################
# create a genotype matrix from the columns of snail.gen,
# checking that each SNP is polymorphic

snail_gen <- data.frame(snail.gen$snail_ID)
colnames(snail_gen)[1]<- "snail_ID"

i <- 0
cn <- rep(NA,gn-1)

# first get polymorphic loci
for (g in 2:gn){
  
  print (g)
  snail_gen$l1 <- genotype(snail.gen[,g],sep=1)
  af <- summary(snail_gen$l1)[[2]][1,2] # allele frequency
  #only keep the locus if it is polymorphic
  #save contig and position for loci retained
  if (af>0 & af<1){
    names(snail_gen)[names(snail_gen)=="l1"] <- colnames(snail.gen)[g]
    i <- i+1 # count up if this SNP is in a new contig
    cn[i] <- colnames(snail.gen)[g]
  }
  else{
    snail_gen <- select(snail_gen,-l1)
  }
  
}



# get LD matrix (lddc is a list with several LD statistics: D, D', r)

lddc <- LD(snail_gen)
matr<-lddc$`R^2`

ldna <- LDnaRaw(matr)
# 'ldna' object is the LD matrix- save the 'matr' object as generating this is the slowest step
#write.csv(matr, "ldna_ANG_female.csv")

#### here, investigate clusters by varying phi and number of edges and saving the outputs for comparison

par(mfcol=c(2,3))
clusters1 <- extractClusters(ldna, min.edges = 100, phi = 10, rm.COCs = TRUE, plot.tree=TRUE)

clusters1

summary <- summaryLDna(ldna, clusters1, matr)  
summary
# save/copy summary to somewhere after each iteration of edges and phi so can compare later
#write.csv(summary, file="summary_edges100_phi10.csv")


 

