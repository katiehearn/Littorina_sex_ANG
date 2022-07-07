rm(list=ls())
library(tidyverse)


####
### Plot genetic maps, by family and sex, coloured by inversion region
### Including conversion of inversion breakpoints from CxC map positions to CxW map positions
### then replot- colour points by inversion region and use shape for sex
####


###### import data - CxC data; CxW data, then reformat data
# set working directory
CC_female_data <- read.table("chr12_female.3_markerName", sep="\t")
CC_male_data <- read.table("chr12_male.3_markerName", sep="\t")
mapLG12 <- read.table("mapposhernan.txt", header=T) %>% filter(LG==12) %>% select("contig", "pos", "Mpos_majority") %>%
  rename("av"="Mpos_majority") %>% unite(contig_pos, c("contig", "pos"), sep="_")
family_81_f <- read.table("marker.Female.informative_fam81_1_LG12.txt", header=F, sep="\t") %>% select(-V5)
colnames(family_81_f) <- c("contig", "pos", "male_pos", "female_pos")
family_81_m <- read.table("marker.Male.informative_fam81_1_LG12.txt", header=F, sep="\t") %>% select(-V5)
colnames(family_81_m) <- c("contig", "pos", "male_pos", "female_pos")
family_82_f <- read.table("marker.Female.informative_fam82_1_LG12.txt", header=F, sep="\t") %>% select(-V5)
colnames(family_82_f) <- c("contig", "pos", "male_pos", "female_pos")
family_82_m <- read.table("marker.Male.informative_fam82_1_LG12.txt", header=F, sep="\t") %>% select(-V5)
colnames(family_82_m) <- c("contig", "pos", "male_pos", "female_pos")
family_83_f <- read.table("marker.Female.informative_fam83_1_LG12.txt", header=F, sep="\t") %>% select(-V5)
colnames(family_83_f) <- c("contig", "pos", "male_pos", "female_pos")
family_83_m <- read.table("marker.Male.informative_fam83_1_LG12.txt", header=F, sep="\t") %>% select(-V5)
colnames(family_83_m) <- c("contig", "pos", "male_pos", "female_pos")
family_91_f <- read.table("marker.Female.informative_fam91_2_LG12.txt", header=F, sep="\t") %>% select(-V5)
colnames(family_91_f) <- c("contig", "pos", "male_pos", "female_pos")
family_91_m <- read.table("marker.Male.informative_fam91_2_LG12.txt", header=F, sep="\t") %>% select(-V5)
colnames(family_91_m) <- c("contig", "pos", "male_pos", "female_pos")
family_91_f <- read.table("marker.Female.informative_fam91_2_LG12.txt", header=F, sep="\t") %>% select(-V5)
colnames(family_91_f) <- c("contig", "pos", "male_pos", "female_pos")
family_91_m <- read.table("marker.Male.informative_fam91_2_LG12.txt", header=F, sep="\t") %>% select(-V5)
colnames(family_91_m) <- c("contig", "pos", "male_pos", "female_pos")
family_92_f <- read.table("marker.Female.informative_fam92_2_LG12.txt", header=F, sep="\t") %>% select(-V5)
colnames(family_92_f) <- c("contig", "pos", "male_pos", "female_pos")
family_92_m <- read.table("marker.Male.informative_fam92_2_LG12.txt", header=F, sep="\t") %>% select(-V5)
colnames(family_92_m) <- c("contig", "pos", "male_pos", "female_pos")
family_93_f <- read.table("marker.Female.informative_fam93_LG12.txt", header=F, sep="\t") %>% select(-V5)
colnames(family_93_f) <- c("contig", "pos", "male_pos", "female_pos")
family_93_m <- read.table("marker.Male.informative_fam93_LG12.txt", header=F, sep="\t") %>% select(-V5)
colnames(family_93_m) <- c("contig", "pos", "male_pos", "female_pos")

## unite contig and pos columns, add family column, full join data frames to make single one per sex to work with
family_81_f <- family_81_f %>% unite(contig_pos, c("contig", "pos"), sep="_") %>% mutate(family="f81") %>% select(-male_pos)
family_81_m <- family_81_m %>% unite(contig_pos, c("contig", "pos"), sep="_") %>% mutate(family="f81") %>% select(-female_pos)
family_82_f <- family_82_f %>% unite(contig_pos, c("contig", "pos"), sep="_") %>% mutate(family="f82") %>% select(-male_pos)
family_82_m <- family_82_m %>% unite(contig_pos, c("contig", "pos"), sep="_") %>% mutate(family="f82") %>% select(-female_pos)
family_83_f <- family_83_f %>% unite(contig_pos, c("contig", "pos"), sep="_") %>% mutate(family="f83") %>% select(-male_pos)
family_83_m <- family_83_m %>% unite(contig_pos, c("contig", "pos"), sep="_") %>% mutate(family="f83") %>% select(-female_pos)
family_91_f <- family_91_f %>% unite(contig_pos, c("contig", "pos"), sep="_") %>% mutate(family="f91") %>% select(-male_pos)
family_91_m <- family_91_m %>% unite(contig_pos, c("contig", "pos"), sep="_") %>% mutate(family="f91") %>% select(-female_pos)
family_92_f <- family_92_f %>% unite(contig_pos, c("contig", "pos"), sep="_") %>% mutate(family="f92") %>% select(-male_pos)
family_92_m <- family_92_m %>% unite(contig_pos, c("contig", "pos"), sep="_") %>% mutate(family="f92") %>% select(-female_pos)
family_93_f <- family_93_f %>% unite(contig_pos, c("contig", "pos"), sep="_") %>% mutate(family="f93") %>% select(-male_pos)
family_93_m <- family_93_m %>% unite(contig_pos, c("contig", "pos"), sep="_") %>% mutate(family="f93") %>% select(-female_pos)

families_f <- bind_rows(family_81_f, family_82_f) %>% bind_rows(., family_83_f) %>% bind_rows(., family_91_f) %>%
  bind_rows(., family_92_f) %>% bind_rows(family_93_f)
families_m <- bind_rows(family_81_m, family_82_m) %>% bind_rows(., family_83_m) %>% bind_rows(., family_91_m) %>%
  bind_rows(., family_92_m) %>% bind_rows(family_93_m)

# add headers to CxC dataframes and remove unneeded columns
CC_female_data <- CC_female_data %>% select(-V5, -V6)
colnames(CC_female_data) <- c("contig", "pos", "male_pos", "female_pos")
CC_male_data <- CC_male_data %>% select(-V5, -V6)
colnames(CC_male_data) <- c("contig", "pos", "male_pos", "female_pos")

# make dataset of markers that are in both the male and female CxC dataset
CC_female_common <- CC_female_data %>% unite(contig_pos, c("contig", "pos"), sep="_")
CC_male_common <- CC_male_data %>% unite(contig_pos, c("contig", "pos"), sep="_")
CC_female_common <- CC_female_common %>% filter(contig_pos %in% CC_male_common$contig_pos)
CC_male_common <- CC_male_common %>% filter(contig_pos %in% CC_female_common$contig_pos)
# 505 common SNPs (out of 537 in female map and 940 in male)- so a good amount


###### data frames to use: CC_female_common, CC_male_common, families_f, families_m

## add inversion breakpoint information to CxC data

# SNPs in CxC map are contig_pos column in female or male data frame (same SNPs in either)

# make data frame of CC SNPs and what inversion they're in - need to import combined map to do this
# then join this to CxW datasets
# then if CxW data is in order, should be able to label everything according to inversion
# need to do separately for each sex

CCf <- CC_female_common %>% rename("f_female_pos"="female_pos",
                                   "f_male_pos"="male_pos")
CCm <- CC_male_common %>% rename("m_female_pos"="female_pos",
                                 "m_male_pos"="male_pos")
CC_common <- full_join(CCf, CCm, by="contig_pos")

CC_SNPs_inv <- CC_common %>% left_join(., mapLG12, by="contig_pos") %>% 
  mutate(inv=ifelse(av<=32.785, "LGC121",
                    ifelse(av>=32.98692 & av<=43.801, "LGC122",
                           ifelse(av>=43.801 & av<=48.712, "LGC123", 
                                  ifelse(av>=48.712 & av<=60.235, "LGC124", NA)))))
CC_SNPs_inv <- CC_SNPs_inv %>% 
  mutate(new_inv=ifelse(m_male_pos<=10.269, "LGC121",
                    ifelse(m_male_pos>=11.108 & m_male_pos<=32.301, "LGC122",
                           ifelse(m_male_pos>=34.613 & m_male_pos<=42.123, "LGC123", 
                                  ifelse(m_male_pos>=42.794 , "LGC124", NA)))))  

## CxC inv pos sorted; just need to plot


## add inversion breakpoint information to CxW data

# add inversion labels from CxC dataframe using common contigs rather than common SNPs
  
map_contigs <- mapLG12 %>% separate(contig_pos, c("contig", "pos"), sep="_") %>% select(-pos)
CW_f_SNPs_inv_contig <- families_f %>% separate(contig_pos, c("contig", "pos"), sep="_") %>% 
  left_join(., map_contigs, by="contig") %>%
  mutate(inv=ifelse(av<=32.785, "LGC121",
                    ifelse(av>=32.98692 & av<=43.801, "LGC122",
                           ifelse(av>=43.801 & av<=48.712, "LGC123", 
                                  ifelse(av>=48.712 & av<=60.235, "LGC124", NA)))))
CW_f_SNPs_inv_contig %>% arrange(family, inv, female_pos) %>% View()

CW_m_SNPs_inv_contig <- families_m %>% separate(contig_pos, c("contig", "pos"), sep="_") %>% 
  left_join(., map_contigs, by="contig") %>%
  mutate(inv=ifelse(av<=32.785, "LGC121",
                    ifelse(av>=32.98692 & av<=43.801, "LGC122",
                           ifelse(av>=43.801 & av<=48.712, "LGC123", 
                                  ifelse(av>=48.712 & av<=60.235, "LGC124", NA)))))
CW_m_SNPs_inv_contig %>% arrange(family, inv, male_pos) %>% View()


### add SNP index and plot

## CxC
CC_SNPs_inv <- CC_SNPs_inv %>% arrange(f_male_pos) %>% mutate(index=row_number())
# sex is colour; inv is shape
ggplot(CC_SNPs_inv, aes(x=index, y=f_female_pos, shape=new_inv)) + 
  geom_point(colour="goldenrod2", alpha=0.5, size=4) +
  geom_point(aes(y=f_male_pos), colour="green4", alpha=0.5, size=4) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
# sex is shape; inv is colour
ggplot(CC_SNPs_inv, aes(x=index, y=f_female_pos, colour=new_inv)) + 
  geom_point(shape=5, alpha=0.8, size=4, show.legend = FALSE) +
  geom_path(colour="black", alpha=0.6) +
  geom_point(aes(y=f_male_pos), shape=19, alpha=0.6, size=4, show.legend = FALSE) +
  geom_path(aes(y=f_male_pos), colour="black", alpha=0.6) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


## CxW
CW_f_SNPs_inv_contig <- CW_f_SNPs_inv_contig %>% group_by(family) %>% arrange(female_pos) %>% mutate(index=row_number())
CW_m_SNPs_inv_contig <- CW_m_SNPs_inv_contig %>% group_by(family) %>% arrange(male_pos)%>% mutate(index=row_number())

ggplot(CW_f_SNPs_inv_contig, aes(x=index, y=female_pos, shape=inv)) + 
  geom_point(colour="goldenrod2", alpha=0.5, size=3) +
  geom_point(data=CW_m_SNPs_inv_contig, aes(x=index, y=male_pos), colour="green4", alpha=0.5, size=3) +
  facet_wrap(~family) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggplot(CW_f_SNPs_inv_contig, aes(x=index, y=female_pos, colour=inv)) + 
  geom_point(shape=25, alpha=0.5, size=3) +
  geom_point(data=CW_m_SNPs_inv_contig, aes(x=index, y=male_pos), shape=21, alpha=0.5, size=3) +
  facet_wrap(~family) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
# SNPs within a map pos are not in inversion order yet

# need to add in inversion labels for rest of CW SNPs, based on outermost SNPs, per family
# plot each family separately so they can have different scale axes

CW_f_fam81 <- CW_f_SNPs_inv_contig %>% filter(family=="f81") %>% 
  mutate(new_inv=ifelse(inv=="LGC121", "LGC121",
                        ifelse(inv=="LGC122", "LGC122",
                               ifelse(inv=="LGC123", "LGC123",
                                      ifelse(inv=="LGC124", "LGC124",
                                             ifelse(female_pos<14.799, "LGC121",
                                            ifelse(female_pos>14.799 & female_pos<52.923, "LGC122", 
                                             ifelse(female_pos>52.923, "LGC124", NA))))))))
CW_f_fam82 <- CW_f_SNPs_inv_contig %>% filter(family=="f82") %>% 
  mutate(new_inv=ifelse(inv=="LGC121", "LGC121",
                        ifelse(inv=="LGC122", "LGC122",
                               ifelse(inv=="LGC123", "LGC123",
                                      ifelse(inv=="LGC124", "LGC124",
                                             ifelse(female_pos<18.323, "LGC121",
                                             ifelse(female_pos>=20.823 & female_pos<43.460, "LGC123", 
                                             ifelse(female_pos>43.460, "LGC124", NA))))))))
CW_f_fam83 <- CW_f_SNPs_inv_contig %>% filter(family=="f83") %>% 
  mutate(new_inv=ifelse(inv=="LGC121", "LGC121",
                        ifelse(inv=="LGC122", "LGC122",
                               ifelse(inv=="LGC123", "LGC123",
                                      ifelse(inv=="LGC124", "LGC124",
                                             ifelse(female_pos<14.799, "LGC121",
                                             ifelse(female_pos>14.799 & female_pos<52.923, "LGC122",
                                             ifelse(female_pos>61.847, "LGC124", NA))))))))
CW_f_fam91 <- CW_f_SNPs_inv_contig %>% filter(family=="f91") %>% 
  mutate(new_inv=ifelse(inv=="LGC121", "LGC121",
                        ifelse(inv=="LGC122", "LGC122",
                               ifelse(inv=="LGC123", "LGC123",
                                      ifelse(inv=="LGC124", "LGC124",
                                             ifelse(female_pos<=66.576, "LGC121",
                                             ifelse(female_pos>=71.125 & female_pos<92.497, "LGC122",
                                             ifelse(female_pos>92.947, "LGC124", NA))))))))
CW_f_fam92 <- CW_f_SNPs_inv_contig %>% filter(family=="f92") %>% 
  mutate(new_inv=inv)

CW_f_fam93 <- CW_f_SNPs_inv_contig %>% filter(family=="f93") %>% 
  mutate(new_inv=ifelse(inv=="LGC121", "LGC121",
                        ifelse(inv=="LGC122", "LGC122",
                               ifelse(inv=="LGC123", "LGC123",
                                      ifelse(inv=="LGC124", "LGC124",
                                             ifelse(female_pos>=39.610 & female_pos<=69.662, "LGC121",
                                             ifelse(female_pos>0.781 & female_pos<=38.828, "LGC122", NA)))))))

CW_m_fam81 <- CW_m_SNPs_inv_contig %>% filter(family=="f81") %>% 
  mutate(new_inv=ifelse(inv=="LGC121", "LGC121",
                        ifelse(inv=="LGC122", "LGC122",
                               ifelse(inv=="LGC123", "LGC123",
                                      ifelse(inv=="LGC124", "LGC124",
                                             ifelse(male_pos<4.766, "LGC121",
                                                    ifelse(male_pos>4.766 & male_pos<28.593, "LGC122", 
                                                           ifelse(male_pos>28.593 & male_pos<42.890, "LGC123", 
                                                                  ifelse(male_pos>42.890, "LGC124", NA)))))))))

CW_m_fam82 <- CW_m_SNPs_inv_contig %>% filter(family=="f82") %>% 
  mutate(new_inv=ifelse(inv=="LGC121", "LGC121",
                        ifelse(inv=="LGC122", "LGC122",
                               ifelse(inv=="LGC123", "LGC123",
                                      ifelse(inv=="LGC124", "LGC124",
                                             ifelse(male_pos>91.923 & male_pos<112.059, "LGC121",
                                                    ifelse(male_pos>56.652 & male_pos<86.790, "LGC122", NA)))))))

CW_m_fam83 <- CW_m_SNPs_inv_contig %>% filter(family=="f83") %>% 
  mutate(new_inv=ifelse(inv=="LGC121", "LGC121",
                        ifelse(inv=="LGC122", "LGC122",
                               ifelse(inv=="LGC123", "LGC123",
                                      ifelse(inv=="LGC124", "LGC124", NA)))))

CW_m_fam91 <- CW_m_SNPs_inv_contig %>% filter(family=="f91") %>% 
  mutate(new_inv=ifelse(inv=="LGC121", "LGC121",
                        ifelse(inv=="LGC122", "LGC122",
                               ifelse(inv=="LGC123", "LGC123",
                                      ifelse(inv=="LGC124", "LGC124", NA)))))

CW_m_fam92 <- CW_m_SNPs_inv_contig %>% filter(family=="f92") %>% 
  mutate(new_inv=ifelse(inv=="LGC121", "LGC121",
                        ifelse(inv=="LGC122", "LGC122",
                               ifelse(inv=="LGC123", "LGC123",
                                      ifelse(inv=="LGC124", "LGC124", 
                                             ifelse(male_pos>40.681 & male_pos<55.948, "LGC121", NA))))))

CW_m_fam93 <- CW_m_SNPs_inv_contig %>% filter(family=="f93") %>% 
  mutate(new_inv=ifelse(inv=="LGC121", "LGC121",
                        ifelse(inv=="LGC122", "LGC122",
                               ifelse(inv=="LGC123", "LGC123",
                                      ifelse(inv=="LGC124", "LGC124", 
                                             ifelse(male_pos>40.681 & male_pos<55.948, "LGC121", NA))))))


## include only contigs that are in the CxC maps (remove NA coloured SNPs)
## index markers according to contig order in CxC map (so SNPs with no recombination are put in order of inversion)
## flip those ones where the order is reversed (LGC12.1 on the right)
##    ones that need reversing are: f82_m, f83_m, f92_f, f92_m, f93_f

## first reverse the order on those that need it
CW_m_fam82 <- CW_m_fam82 %>% select(-index) %>% arrange(desc(male_pos)) %>% mutate(index=row_number()) %>%
  mutate(new_male_pos=-(male_pos-112.059)) %>% select(-male_pos) %>% rename("male_pos"="new_male_pos")
CW_m_fam83 <- CW_m_fam83 %>% select(-index) %>% arrange(desc(male_pos)) %>% mutate(index=row_number()) %>%
  mutate(new_male_pos=-(male_pos-62.882)) %>% select(-male_pos) %>% rename("male_pos"="new_male_pos")
CW_m_fam92 <- CW_m_fam92 %>% select(-index) %>% arrange(desc(male_pos)) %>% mutate(index=row_number()) %>%
  mutate(new_male_pos=-(male_pos-55.948)) %>% select(-male_pos) %>% rename("male_pos"="new_male_pos")
CW_f_fam92 <- CW_f_fam92 %>% select(-index) %>% arrange(desc(female_pos)) %>% mutate(index=row_number()) %>%
  mutate(new_female_pos=-(female_pos-85.571)) %>% select(-female_pos) %>% rename("female_pos"="new_female_pos")
CW_f_fam93 <- CW_f_fam93 %>% select(-index) %>% arrange(desc(female_pos)) %>% mutate(index=row_number()) %>%
  mutate(new_female_pos=-(female_pos-69.662)) %>% select(-female_pos) %>% rename("female_pos"="new_female_pos")

## index markers according to position on CxC map so non-recomb areas aren't mixed up
CW_f_fam81 <- CW_f_fam81 %>% arrange(female_pos, av) %>% mutate(new_index=row_number())
CW_m_fam81 <- CW_m_fam81 %>% arrange(male_pos, av) %>% mutate(new_index=row_number())
CW_f_fam82 <- CW_f_fam82 %>% arrange(female_pos, av) %>% mutate(new_index=row_number())
CW_m_fam82 <- CW_m_fam82 %>% arrange(male_pos, av) %>% mutate(new_index=row_number())
CW_f_fam83 <- CW_f_fam83 %>% arrange(female_pos, av) %>% mutate(new_index=row_number())
CW_m_fam83 <- CW_m_fam83 %>% arrange(male_pos, av) %>% mutate(new_index=row_number())
CW_f_fam91 <- CW_f_fam91 %>% arrange(female_pos, av) %>% mutate(new_index=row_number())
CW_m_fam91 <- CW_m_fam91 %>% arrange(male_pos, av) %>% mutate(new_index=row_number())
CW_f_fam92 <- CW_f_fam92 %>% arrange(female_pos, av) %>% mutate(new_index=row_number())
CW_m_fam92 <- CW_m_fam92 %>% arrange(male_pos, av) %>% mutate(new_index=row_number())
CW_f_fam93 <- CW_f_fam93 %>% arrange(female_pos, av) %>% mutate(new_index=row_number())
CW_m_fam93 <- CW_m_fam93 %>% arrange(male_pos, av) %>% mutate(new_index=row_number())


## remove NAs and plot
CW_f_fam81_noNAs <- CW_f_fam81 %>% filter(!is.na(new_inv)) %>% arrange(new_index) %>% 
  mutate(new_index2=row_number()) %>% select(-new_index) %>% rename("new_index"="new_index2")
CW_m_fam81_noNAs <- CW_m_fam81 %>% filter(!is.na(new_inv)) %>% arrange(new_index) %>% 
  mutate(new_index2=row_number()) %>% select(-new_index) %>% rename("new_index"="new_index2")
CW_f_fam82_noNAs <- CW_f_fam82 %>% filter(!is.na(new_inv)) %>% arrange(new_index) %>% 
  mutate(new_index2=row_number()) %>% select(-new_index) %>% rename("new_index"="new_index2")
CW_m_fam82_noNAs <- CW_m_fam82 %>% filter(!is.na(new_inv)) %>% arrange(new_index) %>% 
  mutate(new_index2=row_number()) %>% select(-new_index) %>% rename("new_index"="new_index2")
CW_f_fam83_noNAs <- CW_f_fam83 %>% filter(!is.na(new_inv)) %>% arrange(new_index) %>% 
  mutate(new_index2=row_number()) %>% select(-new_index) %>% rename("new_index"="new_index2")
CW_m_fam83_noNAs <- CW_m_fam83 %>% filter(!is.na(new_inv)) %>% arrange(new_index) %>% 
  mutate(new_index2=row_number()) %>% select(-new_index) %>% rename("new_index"="new_index2")
CW_f_fam91_noNAs <- CW_f_fam91 %>% filter(!is.na(new_inv)) %>% arrange(new_index) %>% 
  mutate(new_index2=row_number()) %>% select(-new_index) %>% rename("new_index"="new_index2")
CW_m_fam91_noNAs <- CW_m_fam91 %>% filter(!is.na(new_inv)) %>% arrange(new_index) %>% 
  mutate(new_index2=row_number()) %>% select(-new_index) %>% rename("new_index"="new_index2")
CW_f_fam92_noNAs <- CW_f_fam92 %>% filter(!is.na(new_inv)) %>% arrange(new_index) %>% 
  mutate(new_index2=row_number()) %>% select(-new_index) %>% rename("new_index"="new_index2")
CW_m_fam92_noNAs <- CW_m_fam92 %>% filter(!is.na(new_inv)) %>% arrange(new_index) %>% 
  mutate(new_index2=row_number()) %>% select(-new_index) %>% rename("new_index"="new_index2")
CW_f_fam93_noNAs <- CW_f_fam93 %>% filter(!is.na(new_inv)) %>% arrange(new_index) %>% 
  mutate(new_index2=row_number()) %>% select(-new_index) %>% rename("new_index"="new_index2")
CW_m_fam93_noNAs <- CW_m_fam93 %>% filter(!is.na(new_inv)) %>% arrange(new_index) %>% 
  mutate(new_index2=row_number()) %>% select(-new_index) %>% rename("new_index"="new_index2")


## plots
 
ggplot(CW_f_fam81_noNAs, aes(x=new_index, y=female_pos, colour=new_inv)) + 
  geom_point(shape=25, alpha=0.5, size=3) +
  geom_point(data=CW_m_fam81_noNAs, aes(x=new_index, y=male_pos, colour=new_inv), shape=21, alpha=0.5, size=3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggplot(CW_f_fam82_noNAs, aes(x=new_index, y=female_pos, colour=new_inv)) + 
  geom_point(shape=25, alpha=0.5, size=3) +
  geom_point(data=CW_m_fam82_noNAs, aes(x=new_index, y=male_pos, colour=new_inv), shape=21, alpha=0.5, size=3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggplot(CW_f_fam83_noNAs, aes(x=new_index, y=female_pos, colour=new_inv)) + 
  geom_point(shape=25, alpha=0.5, size=3) +
  geom_point(data=CW_m_fam83_noNAs, aes(x=new_index, y=male_pos, colour=new_inv), shape=21, alpha=0.5, size=3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggplot(CW_f_fam91_noNAs, aes(x=new_index, y=female_pos, colour=new_inv)) + 
  geom_point(shape=25, alpha=0.5, size=3) +
  geom_point(data=CW_m_fam91_noNAs, aes(x=new_index, y=male_pos, colour=new_inv), shape=21, alpha=0.5, size=3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggplot(CW_f_fam92_noNAs, aes(x=new_index, y=female_pos, colour=new_inv)) + 
  geom_point(shape=25, alpha=0.5, size=3) +
  geom_point(data=CW_m_fam92_noNAs, aes(x=new_index, y=male_pos, colour=new_inv), shape=21, alpha=0.5, size=3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggplot(CW_f_fam93_noNAs, aes(x=new_index, y=female_pos, colour=new_inv)) + 
  geom_point(shape=25, alpha=0.5, size=3) +
  geom_point(data=CW_m_fam93_noNAs, aes(x=new_index, y=male_pos, colour=new_inv), shape=21, alpha=0.5, size=3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


### remove the few misplaced markers at start and ends and redo index
# just manually check dataframe; very few to remove

CW_f_fam81_noNAs_cleaned <- CW_f_fam81_noNAs %>% slice(-1, -2, -449) %>% mutate(new_index2=row_number()) %>%
  select(-new_index) %>% rename("new_index"="new_index2")
CW_m_fam81_noNAs_cleaned <- CW_m_fam81_noNAs %>% slice(-438, -439)
CW_f_fam82_noNAs_cleaned <- CW_f_fam82_noNAs %>% slice(-1) %>% mutate(new_index2=row_number()) %>%
  select(-new_index) %>% rename("new_index"="new_index2")
CW_m_fam82_noNAs_cleaned <- CW_m_fam82_noNAs %>% slice(-391)
CW_f_fam83_noNAs_cleaned <- CW_f_fam83_noNAs %>% slice(-1) %>% mutate(new_index2=row_number()) %>%
  select(-new_index) %>% rename("new_index"="new_index2")
CW_m_fam83_noNAs_cleaned <- CW_m_fam83_noNAs %>% slice(-503)
CW_f_fam91_noNAs_cleaned <- CW_f_fam91_noNAs %>% slice(-1) %>% mutate(new_index2=row_number()) %>%
  select(-new_index) %>% rename("new_index"="new_index2")
CW_m_fam91_noNAs_cleaned <- CW_m_fam91_noNAs %>% slice(-c(1:20)) %>% mutate(new_index2=row_number()) %>%
  select(-new_index) %>% rename("new_index"="new_index2")
CW_f_fam93_noNAs_cleaned <- CW_f_fam93_noNAs %>% slice(-1, -4) %>% mutate(new_index2=row_number()) %>%
  select(-new_index) %>% rename("new_index"="new_index2")
CW_m_fam93_noNAs_cleaned <- CW_m_fam93_noNAs %>% slice(-1) %>% mutate(new_index2=row_number()) %>%
  select(-new_index) %>% rename("new_index"="new_index2")

### replot

ggplot(CW_f_fam81_noNAs_cleaned, aes(x=new_index, y=female_pos, colour=new_inv)) + 
  geom_point(shape=5, alpha=0.8, size=3, show.legend = FALSE) +
  geom_path(colour="black", alpha=0.6) +
  geom_point(data=CW_m_fam81_noNAs_cleaned, aes(x=new_index, y=male_pos, colour=new_inv), shape=19, alpha=0.6, size=3, show.legend = FALSE) +
  geom_path(data=CW_m_fam81_noNAs_cleaned, aes(x=new_index, y=male_pos), colour="black", alpha=0.6) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

ggplot(CW_f_fam82_noNAs_cleaned, aes(x=new_index, y=female_pos, colour=new_inv)) + 
  geom_point(shape=5, alpha=0.8, size=3, show.legend = FALSE) +
  geom_path(colour="black", alpha=0.6) +
  geom_point(data=CW_m_fam82_noNAs_cleaned, aes(x=new_index, y=male_pos, colour=new_inv), shape=19, alpha=0.6, size=3, show.legend = FALSE) +
  geom_path(data=CW_m_fam82_noNAs_cleaned, aes(x=new_index, y=male_pos), colour="black", alpha=0.6) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

ggplot(CW_f_fam83_noNAs_cleaned, aes(x=new_index, y=female_pos, colour=new_inv)) + 
  geom_point(shape=5, alpha=0.8, size=3, show.legend = FALSE) +
  geom_path(colour="black", alpha=0.6) +
  geom_point(data=CW_m_fam83_noNAs_cleaned, aes(x=new_index, y=male_pos, colour=new_inv), shape=19, alpha=0.6, size=3, show.legend = FALSE) +
  geom_path(data=CW_m_fam83_noNAs_cleaned, aes(x=new_index, y=male_pos), colour="black", alpha=0.6) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

ggplot(CW_f_fam91_noNAs_cleaned, aes(x=new_index, y=female_pos, colour=new_inv)) + 
  geom_point(shape=5, alpha=0.8, size=3, show.legend = FALSE) +
  geom_path(colour="black", alpha=0.6) +
  geom_point(data=CW_m_fam91_noNAs_cleaned, aes(x=new_index, y=male_pos, colour=new_inv), shape=19, alpha=0.6, size=3, show.legend = FALSE) +
  geom_path(data=CW_m_fam91_noNAs_cleaned, aes(x=new_index, y=male_pos), colour="black", alpha=0.6) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

ggplot(CW_f_fam92_noNAs, aes(x=new_index, y=female_pos, colour=new_inv)) + 
  geom_point(shape=5, alpha=0.8, size=3, show.legend = FALSE) +
  geom_path(colour="black", alpha=0.6) +
  geom_point(data=CW_m_fam92_noNAs, aes(x=new_index, y=male_pos, colour=new_inv), shape=19, alpha=0.6, size=3, show.legend = FALSE) +
  geom_path(data=CW_m_fam92_noNAs, aes(x=new_index, y=male_pos), colour="black", alpha=0.6) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggplot(CW_f_fam93_noNAs_cleaned, aes(x=new_index, y=female_pos, colour=new_inv)) + 
  geom_point(shape=5, alpha=0.8, size=3, show.legend = FALSE) +
  geom_path(colour="black", alpha=0.6) +
  geom_point(data=CW_m_fam93_noNAs_cleaned, aes(x=new_index, y=male_pos, colour=new_inv), shape=19, alpha=0.6, size=3, show.legend = FALSE) +
  geom_path(data=CW_m_fam93_noNAs_cleaned, aes(x=new_index, y=male_pos), colour="black", alpha=0.6) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


