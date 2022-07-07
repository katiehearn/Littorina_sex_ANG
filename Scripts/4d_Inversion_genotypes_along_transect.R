####
## again, to be run after LD cluster investigation script

###############
#### Plot PCA cluster identity along transect
#### Then calculate inversion (PC cluster) frequency in windows along transect


## Generate dataframes of PCA values and snail info for each cluster and sex


clust3_females_PCandinfo <- cluster2398_0.88_PCs %>% left_join(., ANG_snail_info_new, by="snail_ID")
clust4_females_PCandinfo <- cluster3414_0.69_PCs %>% left_join(., ANG_snail_info_new, by="snail_ID")
region2a_females_PCandinfo <- region2a_PCs %>% left_join(., ANG_snail_info_new, by="snail_ID")
region2b_females_PCandinfo <- region2b_PCs %>% left_join(., ANG_snail_info_new, by="snail_ID")

clust3_males_PCandinfo <- malecluster3_PCs %>% left_join(., ANG_snail_info_new, by="snail_ID")
clust4_males_PCandinfo <- malecluster4_PCs %>% left_join(., ANG_snail_info_new, by="snail_ID")
region2a_males_PCandinfo <- region2a_m_PCs %>% left_join(., ANG_snail_info_new, by="snail_ID")
region2b_males_PCandinfo <- region2b_m_PCs %>% left_join(., ANG_snail_info_new, by="snail_ID")


## assign cluster identities and plot along transect (Dist vs cluster) separately for males and females

### Region1 females
# check PCA 
ggplot(clust4_females_PCandinfo, aes(x=Axis1, y=Axis2, colour=ecotype.x)) + geom_point(alpha=0.5)
region1_females <- clust4_females_PCandinfo #%>% select(-groups)
# define groups
region1_females <- region1_females %>% mutate(cluster=ifelse(Axis1<1, 1,
                                                             ifelse(Axis1>75, 3, 2)))
region1_females$cluster <- as.factor(region1_females$cluster)
ggplot(region1_females, aes(x=Axis1, y=Axis2, colour=cluster)) + geom_point(alpha=0.4, size=3) + theme_bw() 
# plot along transect
ggplot(region1_females, aes(x=DistAlongPath, y=cluster, colour=cluster)) + geom_point(shape=21, size=3, alpha=0.8, stroke=1.5) + theme_bw()

### Region2a females
ggplot(region2a_females_PCandinfo, aes(x=Axis1, y=Axis2, colour=ecotype.x)) + geom_point(alpha=0.5)
region2a_females <- region2a_females_PCandinfo #%>% select(-groups)
region2a_females <- region2a_females %>% 
  mutate(cluster=ifelse(Axis1<(-40), 1,
                        ifelse(Axis1>40, 3, 
                               ifelse(between(Axis1,-40,40), 2, NA))))
region2a_females$cluster <- as.factor(region2a_females$cluster)
ggplot(region2a_females, aes(x=Axis1, y=Axis2, colour=cluster)) + geom_point(alpha=0.4, size=3) + theme_bw()
ggplot(region2a_females, aes(x=DistAlongPath, y=cluster, colour=cluster)) + 
  geom_point(shape=21, alpha=0.8, stroke=1.5, size=3) + theme_bw()

### Region2b females
ggplot(region2b_females_PCandinfo, aes(x=Axis1, y=Axis2, colour=ecotype.x)) + geom_point(alpha=0.5)
region2b_females <- region2b_females_PCandinfo #%>% select(-groups)
region2b_females <- region2b_females %>% 
  mutate(cluster=ifelse(Axis1<(-40), 1,
                        ifelse(Axis1>20, 3, 
                               ifelse(between(Axis1,-30,15), 2, NA))))
region2b_females$cluster <- as.factor(region2b_females$cluster)
ggplot(region2b_females, aes(x=Axis1, y=Axis2, colour=cluster)) + geom_point(alpha=0.4, size=3) + theme_bw()
ggplot(region2b_females, aes(x=DistAlongPath, y=cluster, colour=cluster)) + 
  geom_point(shape=21, alpha=0.8, stroke=1.5, size=3) + theme_bw()

### Region3 females
ggplot(clust3_females_PCandinfo, aes(x=Axis1, y=Axis2, colour=ecotype.x)) + geom_point(alpha=0.5)
region3_females <- clust3_females_PCandinfo #%>% select(-groups)
region3_females <- region3_females %>% 
  mutate(cluster=ifelse(Axis1<(-30), 1,
                        ifelse(Axis1>40, 3, 
                               ifelse(between(Axis1,-15,15), 2, NA))))
region3_females$cluster <- as.factor(region3_females$cluster)
ggplot(region3_females, aes(x=Axis1, y=Axis2, colour=cluster)) + geom_point(alpha=0.4, size=3) + theme_bw()
ggplot(region3_females, aes(x=DistAlongPath, y=cluster, colour=cluster)) + geom_point(shape=21, alpha=0.8, stroke=1.5, size=3) + theme_bw()

##### now males
## remember to check and try to match up cluster labelling with female clusters

### Region1 males
# check PCA 
ggplot(clust4_males_PCandinfo, aes(x=Axis1, y=Axis2, colour=ecotype.x)) + geom_point(alpha=0.5)
region1_males <- clust4_males_PCandinfo #%>% select(-groups)
# define groups
region1_males <- region1_males %>% mutate(cluster=ifelse(Axis1<0, 1,
                                                             ifelse(Axis1>75, 3, 
                                                                    ifelse(between(Axis1, 25, 65), 2, NA))))
region1_males$cluster <- as.factor(region1_males$cluster)
ggplot(region1_males, aes(x=Axis1, y=Axis2, colour=cluster)) + geom_point(alpha=0.4, size=3) + theme_bw()
# plot along transect
ggplot(region1_males, aes(x=DistAlongPath, y=cluster, colour=cluster)) + geom_point(shape=21, size=3, alpha=0.8, stroke=1.5) + theme_bw()

### Region2a males
ggplot(region2a_males_PCandinfo, aes(x=Axis1, y=Axis2, colour=ecotype.x)) + geom_point(alpha=0.5)
region2a_males <- region2a_males_PCandinfo #%>% select(-groups)
region2a_males <- region2a_males %>% 
  mutate(cluster=ifelse(Axis1<(-75), 3,
                        ifelse(Axis1>0, 1, 
                               ifelse(between(Axis1,-60,-10), 2, NA))))
region2a_males$cluster <- as.factor(region2a_males$cluster)
ggplot(region2a_males, aes(x=Axis1, y=Axis2, colour=cluster)) + geom_point(alpha=0.4, size=3) + theme_bw()
ggplot(region2a_males, aes(x=DistAlongPath, y=cluster, colour=cluster)) + 
  geom_point(shape=21, alpha=0.8, stroke=1.5, size=3) + theme_bw()


### Region2b males
ggplot(region2b_males_PCandinfo, aes(x=Axis1, y=Axis2, colour=ecotype.x)) + geom_point(alpha=0.5)
region2b_males <- region2b_males_PCandinfo #%>% select(-groups)
region2b_males <- region2b_males %>% 
  mutate(cluster=ifelse(Axis1<20, 3,
                        ifelse(Axis1>80, 1, 
                               ifelse(between(Axis1,20,80), 2, NA))))
region2b_males$cluster <- as.factor(region2b_males$cluster)
ggplot(region2b_males, aes(x=Axis1, y=Axis2, colour=cluster)) + geom_point(alpha=0.4, size=3) + theme_bw()
ggplot(region2b_males, aes(x=DistAlongPath, y=cluster, colour=cluster)) + 
  geom_point(shape=21, alpha=0.8, stroke=1.5, size=3) + theme_bw()


### Region3 males
# cluster1
ggplot(clust3_males_PCandinfo, aes(x=Axis1, y=Axis2, colour=ecotype.x)) + geom_point(alpha=0.5)
region3_males <- clust3_males_PCandinfo #%>% select(-groups)
region3_males <- region3_males %>% 
  mutate(cluster=ifelse(Axis1<(-60), 3,
                        ifelse(Axis1>15, 1, 
                               ifelse(between(Axis1,-45,0), 2, NA))))
region3_males$cluster <- as.factor(region3_males$cluster)
ggplot(region3_males, aes(x=Axis1, y=Axis2, colour=cluster)) + geom_point(alpha=0.4, size=3) + theme_bw()
ggplot(region3_males, aes(x=DistAlongPath, y=cluster, colour=cluster)) + geom_point(shape=21, size=3, alpha=0.8, stroke=1.5) + theme_bw()



##############
##############

### now calculate frequencies of each inversion haplotype in windows along transect
### like doing allele frequencies


### region1 females

# windows of 25 snails, shifting along 1 snails at a time

# loop that slices the dataframe from x to x+24 , with x being increased by 1 (or more) each repetition
#   then calculates the inv frequencies
#   then with sections being saved into a list

# set up empty list, seq of snails to start at, and x counter
region1_females_pathwins <- list()
seq=seq(1, 185, by=5)
n <- 1
# make loop to cut into list of windows
for(i in seq) {
  nam <- paste(n)
  window <- region1_females %>% arrange(DistAlongPath) %>% dplyr::slice(n:(n+24))
  n <- n+5
  region1_females_pathwins[[nam]] <- window
}

# next do loop of calculating inv frequencies (set up empty list to save into first) and mean distalongpath
windownames <- as.character(seq)
region1_females_invfreqs <- list()
no <- 0
windows_dists <- list()

for(i in region1_females_pathwins) {
  no <- no+1
  windownam <- windownames[no]
  window_freqs <- i %>% group_by(cluster) %>% summarise(freq=n())
  window_dist <- i %>% summarise(mean_dist=mean(DistAlongPath))
  region1_females_invfreqs[[windownam]] <- window_freqs
  windows_dists[[windownam]] <- window_dist
}


# next need to collapse into a data frame (collapse both freqs and the mean dist)
region1_females_invfreqscol <- do.call(rbind.data.frame, region1_females_invfreqs)
region1_females_invfreqscol <- tibble::rownames_to_column(region1_females_invfreqscol, "window") %>%
  separate(window, c("window", "rep"), sep="\\.") %>%
  select(-rep)
region1_females_invfreqscol$window <- as.numeric(region1_females_invfreqscol$window)
# now spread and calculate inv freqs/props
region1_females_invfreqscol <- region1_females_invfreqscol %>% 
  spread(cluster, freq) %>% arrange(window)
region1_females_invfreqscol[is.na(region1_females_invfreqscol)] <- 0
region1_females_invfreqscol <- region1_females_invfreqscol %>%
  mutate(prop_x=(`1`+0.5*`2`)/25,
         prop_y=(`3`+0.5*`2`)/25)
region1_females_meandist <- do.call(rbind.data.frame, windows_dists)
region1_females_meandist <- tibble::rownames_to_column(region1_females_meandist, "window")
region1_females_meandist$window <- as.numeric(region1_females_meandist$window)

# add distance info to  windows and plot
region1_females_invfreqscol <- region1_females_invfreqscol %>% 
  left_join(., region1_females_meandist, by="window")
ggplot(region1_females_invfreqscol, aes(x=mean_dist, y=prop_x)) + 
  geom_point(shape=21, size=3, stroke=2, alpha=0.8, colour="#F8766D") + 
  geom_point(aes(y=prop_y), shape=21, size=3, stroke=2, alpha=0.8, colour="#619CFF") + 
  scale_y_continuous(limits=c(0,1)) + theme_bw()


### region2a females
region2a_females_pathwins <- list()
seq=seq(1, 185, by=5)
n <- 1
for(i in seq) {
  nam <- paste(n)
  window <- region2a_females %>% arrange(DistAlongPath) %>% slice(n:(n+24))
  n <- n+5
  region2a_females_pathwins[[nam]] <- window
}
windownames <- as.character(seq)
invfreqs <- list()
no <- 0
windows_dists <- list()
for(i in region2a_females_pathwins) {
  no <- no+1
  windownam <- windownames[no]
  window_freqs <- i %>% group_by(cluster) %>% summarise(freq=n())
  window_dist <- i %>% summarise(mean_dist=mean(DistAlongPath))
  invfreqs[[windownam]] <- window_freqs
  windows_dists[[windownam]] <- window_dist
}
invfreqscol <- do.call(rbind.data.frame, invfreqs)
invfreqscol <- tibble::rownames_to_column(invfreqscol, "window") %>%
  separate(window, c("window", "rep"), sep="\\.") %>%
  select(-rep)
invfreqscol$window <- as.numeric(invfreqscol$window)
invfreqscol <- invfreqscol %>% 
  spread(cluster, freq) %>% arrange(window)
invfreqscol[is.na(invfreqscol)] <- 0

region2a_females_invfreqscol <- invfreqscol %>%
  mutate(prop_x=(`1`+0.5*`2`)/(`1`+`2`+`3`),
         prop_y=(`3`+0.5*`2`)/(`1`+`2`+`3`))
meandist <- do.call(rbind.data.frame, windows_dists)
meandist <- tibble::rownames_to_column(meandist, "window")
meandist$window <- as.numeric(meandist$window)

region2a_females_invfreqscol <- region2a_females_invfreqscol %>% 
  left_join(., meandist, by="window")
ggplot(region2a_females_invfreqscol, aes(x=mean_dist, y=prop_x)) + 
  geom_point(shape=21, alpha=0.8, size=3, stroke=2, colour="#F8766D") + 
  geom_point(aes(y=prop_y), shape=21, alpha=0.8,size=3, stroke=2,  colour="#619CFF") +
  scale_y_continuous(limits=c(0,1)) + theme_bw()

### region2b females
region2b_females_pathwins <- list()
seq=seq(1, 185, by=5)
n <- 1
for(i in seq) {
  nam <- paste(n)
  window <- region2b_females %>% arrange(DistAlongPath) %>% slice(n:(n+24))
  n <- n+5
  region2b_females_pathwins[[nam]] <- window
}
windownames <- as.character(seq)
invfreqs <- list()
no <- 0
windows_dists <- list()
for(i in region2b_females_pathwins) {
  no <- no+1
  windownam <- windownames[no]
  window_freqs <- i %>% group_by(cluster) %>% summarise(freq=n())
  window_dist <- i %>% summarise(mean_dist=mean(DistAlongPath))
  invfreqs[[windownam]] <- window_freqs
  windows_dists[[windownam]] <- window_dist
}
invfreqscol <- do.call(rbind.data.frame, invfreqs)
invfreqscol <- tibble::rownames_to_column(invfreqscol, "window") %>%
  separate(window, c("window", "rep"), sep="\\.") %>%
  select(-rep)
invfreqscol$window <- as.numeric(invfreqscol$window)
invfreqscol <- invfreqscol %>% 
  spread(cluster, freq) %>% arrange(window)
invfreqscol[is.na(invfreqscol)] <- 0

region2b_females_invfreqscol <- invfreqscol %>%
  mutate(prop_x=(`1`+0.5*`2`)/(`1`+`2`+`3`),
         prop_y=(`3`+0.5*`2`)/(`1`+`2`+`3`))
meandist <- do.call(rbind.data.frame, windows_dists)
meandist <- tibble::rownames_to_column(meandist, "window")
meandist$window <- as.numeric(meandist$window)

region2b_females_invfreqscol <- region2b_females_invfreqscol %>% 
  left_join(., meandist, by="window")
ggplot(region2b_females_invfreqscol, aes(x=mean_dist, y=prop_x)) + 
  geom_point(shape=21, alpha=0.8, size=3, stroke=2, colour="#F8766D") + 
  geom_point(aes(y=prop_y), shape=21, alpha=0.8,size=3, stroke=2,  colour="#619CFF") +
  scale_y_continuous(limits=c(0,1)) + theme_bw()


### region3 females
region3_females_pathwins <- list()
seq=seq(1, 185, by=5)
n <- 1
for(i in seq) {
  nam <- paste(n)
  window <- region3_females %>% arrange(DistAlongPath) %>% slice(n:(n+24))
  n <- n+5
  region3_females_pathwins[[nam]] <- window
}
windownames <- as.character(seq)
invfreqs <- list()
no <- 0
windows_dists <- list()
for(i in region3_females_pathwins) {
  no <- no+1
  windownam <- windownames[no]
  window_freqs <- i %>% group_by(cluster) %>% summarise(freq=n())
  window_dist <- i %>% summarise(mean_dist=mean(DistAlongPath))
  invfreqs[[windownam]] <- window_freqs
  windows_dists[[windownam]] <- window_dist
}
invfreqscol <- do.call(rbind.data.frame, invfreqs)
invfreqscol <- tibble::rownames_to_column(invfreqscol, "window") %>%
  separate(window, c("window", "rep"), sep="\\.") %>%
  select(-rep)
invfreqscol$window <- as.numeric(invfreqscol$window)
invfreqscol <- invfreqscol %>% 
  spread(cluster, freq) %>% arrange(window)
invfreqscol[is.na(invfreqscol)] <- 0

region3_females_invfreqscol <- invfreqscol %>%
  mutate(prop_x=(`1`+0.5*`2`)/(`1`+`2`+`3`),
         prop_y=(`3`+0.5*`2`)/(`1`+`2`+`3`))
meandist <- do.call(rbind.data.frame, windows_dists)
meandist <- tibble::rownames_to_column(meandist, "window")
meandist$window <- as.numeric(meandist$window)

region3_females_invfreqscol <- region3_females_invfreqscol %>% 
  left_join(., meandist, by="window")
ggplot(region3_females_invfreqscol, aes(x=mean_dist, y=prop_x)) + 
  geom_point(shape=21, alpha=0.8, size=3, stroke=2, colour="#F8766D") + 
  geom_point(aes(y=prop_y), shape=21, alpha=0.8,size=3, stroke=2,  colour="#619CFF") +
  scale_y_continuous(limits=c(0,1)) + theme_bw()





####### now males 
## (remember to change the seq due to different number of males)
region1_males_pathwins <- list()
seq=seq(1, 146, by=5)
n <- 1
for(i in seq) {
  nam <- paste(n)
  window <- region1_males %>% arrange(DistAlongPath) %>% slice(n:(n+24))
  n <- n+5
  region1_males_pathwins[[nam]] <- window
}
windownames <- as.character(seq)
invfreqs <- list()
no <- 0
windows_dists <- list()
for(i in region1_males_pathwins) {
  no <- no+1
  windownam <- windownames[no]
  window_freqs <- i %>% group_by(cluster) %>% summarise(freq=n())
  window_dist <- i %>% summarise(mean_dist=mean(DistAlongPath))
  invfreqs[[windownam]] <- window_freqs
  windows_dists[[windownam]] <- window_dist
}
invfreqscol <- do.call(rbind.data.frame, invfreqs)
invfreqscol <- tibble::rownames_to_column(invfreqscol, "window") %>%
  separate(window, c("window", "rep"), sep="\\.") %>%
  select(-rep)
invfreqscol$window <- as.numeric(invfreqscol$window)
invfreqscol <- invfreqscol %>% 
  spread(cluster, freq) %>% arrange(window)
invfreqscol[is.na(invfreqscol)] <- 0

region1_males_invfreqscol <- invfreqscol %>%
  mutate(prop_x=(`1`+0.5*`2`)/(`1`+`2`+`3`),
         prop_y=(`3`+0.5*`2`)/(`1`+`2`+`3`))
meandist <- do.call(rbind.data.frame, windows_dists)
meandist <- tibble::rownames_to_column(meandist, "window")
meandist$window <- as.numeric(meandist$window)

region1_males_invfreqscol <- region1_males_invfreqscol %>% 
  left_join(., meandist, by="window")
ggplot(region1_males_invfreqscol, aes(x=mean_dist, y=prop_x)) + 
  geom_point(shape=21, alpha=0.8, size=3, stroke=2, colour="#F8766D") + 
  geom_point(aes(y=prop_y), shape=21, alpha=0.8, size=3, stroke=2, colour="#619CFF") +
  scale_y_continuous(limits=c(0,1)) + theme_bw()


### region2a males
region2a_males_pathwins <- list()
seq=seq(1, 146, by=5)
n <- 1
for(i in seq) {
  nam <- paste(n)
  window <- region2a_males %>% arrange(DistAlongPath) %>% slice(n:(n+24))
  n <- n+5
  region2a_males_pathwins[[nam]] <- window
}
windownames <- as.character(seq)
invfreqs <- list()
no <- 0
windows_dists <- list()
for(i in region2a_males_pathwins) {
  no <- no+1
  windownam <- windownames[no]
  window_freqs <- i %>% group_by(cluster) %>% summarise(freq=n())
  window_dist <- i %>% summarise(mean_dist=mean(DistAlongPath))
  invfreqs[[windownam]] <- window_freqs
  windows_dists[[windownam]] <- window_dist
}
invfreqscol <- do.call(rbind.data.frame, invfreqs)
invfreqscol <- tibble::rownames_to_column(invfreqscol, "window") %>%
  separate(window, c("window", "rep"), sep="\\.") %>%
  select(-rep)
invfreqscol$window <- as.numeric(invfreqscol$window)
invfreqscol <- invfreqscol %>% 
  spread(cluster, freq) %>% arrange(window)
invfreqscol[is.na(invfreqscol)] <- 0

region2a_males_invfreqscol <- invfreqscol %>%
  mutate(prop_x=(`1`+0.5*`2`)/(`1`+`2`+`3`),
         prop_y=(`3`+0.5*`2`)/(`1`+`2`+`3`))
meandist <- do.call(rbind.data.frame, windows_dists)
meandist <- tibble::rownames_to_column(meandist, "window")
meandist$window <- as.numeric(meandist$window)

region2a_males_invfreqscol <- region2a_males_invfreqscol %>% 
  left_join(., meandist, by="window")
ggplot(region2a_males_invfreqscol, aes(x=mean_dist, y=prop_x)) + 
  geom_point(shape=21, size=3, stroke=2, alpha=0.8, colour="#F8766D") + 
  geom_point(aes(y=prop_y),  size=3, stroke=2, shape=21, alpha=0.8, size=3, stroke=2, colour="#619CFF") +
  scale_y_continuous(limits=c(0,1)) + theme_bw()


### region2b males
region2b_males_pathwins <- list()
seq=seq(1, 146, by=5)
n <- 1
for(i in seq) {
  nam <- paste(n)
  window <- region2b_males %>% arrange(DistAlongPath) %>% slice(n:(n+24))
  n <- n+5
  region2b_males_pathwins[[nam]] <- window
}
windownames <- as.character(seq)
invfreqs <- list()
no <- 0
windows_dists <- list()
for(i in region2b_males_pathwins) {
  no <- no+1
  windownam <- windownames[no]
  window_freqs <- i %>% group_by(cluster) %>% summarise(freq=n())
  window_dist <- i %>% summarise(mean_dist=mean(DistAlongPath))
  invfreqs[[windownam]] <- window_freqs
  windows_dists[[windownam]] <- window_dist
}
invfreqscol <- do.call(rbind.data.frame, invfreqs)
invfreqscol <- tibble::rownames_to_column(invfreqscol, "window") %>%
  separate(window, c("window", "rep"), sep="\\.") %>%
  select(-rep)
invfreqscol$window <- as.numeric(invfreqscol$window)
invfreqscol <- invfreqscol %>% 
  spread(cluster, freq) %>% arrange(window)
invfreqscol[is.na(invfreqscol)] <- 0

region2b_males_invfreqscol <- invfreqscol %>%
  mutate(prop_x=(`1`+0.5*`2`)/(`1`+`2`+`3`),
         prop_y=(`3`+0.5*`2`)/(`1`+`2`+`3`))
meandist <- do.call(rbind.data.frame, windows_dists)
meandist <- tibble::rownames_to_column(meandist, "window")
meandist$window <- as.numeric(meandist$window)

region2b_males_invfreqscol <- region2b_males_invfreqscol %>% 
  left_join(., meandist, by="window")
ggplot(region2b_males_invfreqscol, aes(x=mean_dist, y=prop_x)) + 
  geom_point(shape=21, size=3, stroke=2, alpha=0.8, colour="#F8766D") + 
  geom_point(aes(y=prop_y),  size=3, stroke=2, shape=21, alpha=0.8, size=3, stroke=2, colour="#619CFF") +
  scale_y_continuous(limits=c(0,1)) + theme_bw()


### region3 males
region3_males_pathwins <- list()
seq=seq(1, 146, by=5)
n <- 1
for(i in seq) {
  nam <- paste(n)
  window <- region3_males %>% arrange(DistAlongPath) %>% slice(n:(n+24))
  n <- n+5
  region3_males_pathwins[[nam]] <- window
}
windownames <- as.character(seq)
invfreqs <- list()
no <- 0
windows_dists <- list()
for(i in region3_males_pathwins) {
  no <- no+1
  windownam <- windownames[no]
  window_freqs <- i %>% group_by(cluster) %>% summarise(freq=n())
  window_dist <- i %>% summarise(mean_dist=mean(DistAlongPath))
  invfreqs[[windownam]] <- window_freqs
  windows_dists[[windownam]] <- window_dist
}
invfreqscol <- do.call(rbind.data.frame, invfreqs)
invfreqscol <- tibble::rownames_to_column(invfreqscol, "window") %>%
  separate(window, c("window", "rep"), sep="\\.") %>%
  select(-rep)
invfreqscol$window <- as.numeric(invfreqscol$window)
invfreqscol <- invfreqscol %>% 
  spread(cluster, freq) %>% arrange(window)
invfreqscol[is.na(invfreqscol)] <- 0

region3_males_invfreqscol <- invfreqscol %>%
  mutate(prop_x=(`1`+0.5*`2`)/(`1`+`2`+`3`),
         prop_y=(`3`+0.5*`2`)/(`1`+`2`+`3`))
meandist <- do.call(rbind.data.frame, windows_dists)
meandist <- tibble::rownames_to_column(meandist, "window")
meandist$window <- as.numeric(meandist$window)

region3_males_invfreqscol <- region3_males_invfreqscol %>% 
  left_join(., meandist, by="window")
ggplot(region3_males_invfreqscol, aes(x=mean_dist, y=prop_x)) + 
  geom_point(shape=21, size=3, stroke=2, alpha=0.8, colour="#F8766D") + 
  geom_point(aes(y=prop_y),  size=3, stroke=2, shape=21, alpha=0.8, size=3, stroke=2, colour="#619CFF") +
  scale_y_continuous(limits=c(0,1)) + theme_bw()




#### save inversion frequency files and PCA data incl. individuals' inversion genotypes 
## females, males, each of 4 regions. One by one
# set working directory
#write.csv(region3_males, file="region3_males_PCandcluster.csv", row.names=F)
#write.csv(region3_males_invfreqscol, file="region3_males_invfreqs.csv", row.names=F)



