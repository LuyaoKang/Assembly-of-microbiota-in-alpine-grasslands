# Biogeographic pattern of soil microbiome
## Variation and Distance-Decay relationships
deg2rad <- function(deg) return(deg*pi/180)
# Calculates the geodesic distance between two points specified by 
# radian latitude/longitude using the Haversine formula
# Ouputs distance between sites 1 and 2 as meters
gcd.hf <- function(long1, lat1, long2, lat2) {
  R <- 6371 # Earth mean radius [km]
  delta.long <- (long2 - long1)
  delta.lat <- (lat2 - lat1)
  a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
  c <- 2 * asin(min(1,sqrt(a)))
  d = (R * c)*1000
  return(d) # Distance in meters
}


# Fxn to calculate matrix of distances between each two sites
# INPUT: a data frame in which longs are in first column and lats in second column
# OUTPUT: a distance matrix (class dist) between all pairwise sites
# Output distances are in meters
CalcDists <- function(longlats) {
  name <- list(rownames(longlats), rownames(longlats))
  n <- nrow(longlats)
  z <- matrix(0, n, n, dimnames = name)
  for (i in 1:n) {
    for (j in 1:n) z[i, j] <- gcd.hf(long1 = deg2rad(longlats[i, 1]), 
                                     lat1 = deg2rad(longlats[i, 2]), long2 = deg2rad(longlats[j, 1]), 
                                     lat2 = deg2rad(longlats[j, 2]))
  }
  z <- as.dist(z)
  return(z)
}

library(vegan)
library(ggplot2)
######bacteria#########
# calculate the disimilarity matrix of environment factors based on euclid distance
env.table<-sample_data(bac.phylo)
geo_dist<-data.frame(env.table[,1:2])
geo_dist<-CalcDists(geo_dist)/100000
Geography.dist <- as.numeric(geo_dist)
#determine the euclidean distance of environmental variables
env_vars<-env.table[,-(1:4)]
#standardize the env data
env_vars<-decostand(as.data.frame(env_vars), 'standardize')
env_dist<-vegdist(env_vars, 'euclidean',upper=F)
env.dist<-as.numeric(env_dist)
#DDR
ddr_fun <- function(phylo){
  #otu_tab <- decostand(otu_tab, 'hellinger', MARGIN = 1)
  dist_mat <- vegdist(t(otu_table(phylo)), 'bray', upper=F)
  com_dissimilarty <- as.numeric(dist_mat)*100
  com_similarty <- (1-as.numeric(dist_mat))*100
  dist_tab<-as.data.frame(cbind(com_similarty, com_dissimilarty,
                                Geography.dist, env.dist))
  colnames(dist_tab) <-c ('com_similarty', 'com_dissimilarty',
                          'Geography', 'environment')
  list <- list(dist_tab, dist_mat)
  return(list)
}
bacteria.dis.tab <- ddr_fun(bac.phylo)
fungi.dis.tab <- ddr_fun(fungi.phylo)
protist.dis.tab <- ddr_fun(protist.phylo)
animal.dis.tab <- ddr_fun(animal.phylo)

#ordinary least-squares regressions
bacteria.fit<-lm(com_similarty~Geography,data=bacteria.dis.tab[[1]])
summary(bacteria.fit)
bacteria.ln.fit<-lm(com_similarty/100~log(Geography*100+1),data=bacteria.dis.tab[[1]])
summary(bacteria.ln.fit)

fungi.fit<-lm(com_similarty~Geography,data=fungi.dis.tab[[1]])
summary(fungi.fit)
fungi.ln.fit<-lm(com_similarty/100~log(Geography*100+1),data=fungi.dis.tab[[1]])
summary(fungi.ln.fit)

protist.fit<-lm(com_similarty~Geography,data=protist.dis.tab[[1]])
summary(protist.fit)
protist.ln.fit<-lm(com_similarty/100~log(Geography*100+1),data=protist.dis.tab[[1]])
summary(protist.ln.fit)

animal.fit<-lm(com_similarty~Geography,data=animal.dis.tab[[1]])
summary(animal.fit)
animal.ln.fit<-lm(com_similarty/100~log(Geography*100+1),data=animal.dis.tab[[1]])
summary(animal.ln.fit)

#mantel test
mantel(geo_dist,bacteria.dis.tab[[2]],method="spearman")

mantel(geo_dist,fungi.dis.tab[[2]],method="spearman")

mantel(geo_dist,protist.dis.tab[[2]],method="spearman")

mantel(geo_dist,animal.dis.tab[[2]],method="spearman")

#compare the difference between multi regression slope
library(emmeans)
dat<-rbind(cbind(taxa=c(rep('Bacteria',nrow(bacteria.dis.tab[[1]]))),bacteria.dis.tab[[1]]),
           cbind(taxa=c(rep('Fungi',nrow(fungi.dis.tab[[1]]))),fungi.dis.tab[[1]]),
           cbind(taxa=c(rep('Protists',nrow(protist.dis.tab[[1]]))),protist.dis.tab[[1]]),
           cbind(taxa=c(rep('Animal',nrow(animal.dis.tab[[1]]))),animal.dis.tab[[1]]))

##comparison of slope of DDRs
m.interaction <- lm(com_similarty ~ Geography*taxa, data = dat)
anova(m.interaction)

# Obtain slopes
m.interaction$coefficients
geo.lst <- emtrends(m.interaction, "taxa", var="Geography")

# Compare slopes
pairs(geo.lst)

#compare the variations of composition
library(plyr)
ddply(dat, c("taxa"), summarise,
      mean = mean(com_dissimilarty/100), sd = sd(com_dissimilarty/100),
      se = sd(com_dissimilarty/100)/sqrt(length(com_dissimilarty)))

mode_ln<-aov(com_dissimilarty/100~taxa,dat)
summary(mode_ln)
TukeyHSD(mode_ln)$taxa

##plot
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 4
cols = gg_color_hue(n)
dat$taxa<-factor(dat$taxa,ordered=T,levels = c('Bacteria','Fungi', 'Protists','Animal'))

p_dissimilarty<-ggplot(dat,aes(taxa,com_dissimilarty/100))+
  geom_violin(trim=T,aes(fill=taxa))+
  scale_fill_manual(values= cols)+
  geom_boxplot(width=0.1,fill="white")+
  xlab('Taxa')+ylab('Community dissimilarity')+
  theme(axis.title = element_text(size=9),axis.text = element_text(size=7),
        panel.grid = element_blank(),legend.position='none')

p_distance<-ggplot(dat,aes(Geography,com_similarty)) + 
  geom_point(shape=19,alpha=0.1,aes(colour = taxa))+
  scale_color_manual(values=cols)+
  geom_smooth(method="lm", size=1, se=T,colour='black') +
  facet_wrap( .~ taxa , scales="free",ncol = 2) +
  ylab('Community similarity')+xlab('Geographic distance(100km)')+
  theme(axis.title = element_text(size=9),axis.text = element_text(size=7),
        strip.text = element_text(size = 9),legend.position='none',panel.grid = element_blank())

p_ln.distance<-ggplot(dat,aes(log(Geography*100+1),com_similarty/100)) + 
  geom_point(shape=19,alpha=0.1,aes(colour = taxa))+
  scale_color_manual(values=cols)+
  geom_smooth(method="lm", size=1, se=T,colour='black') +
  facet_wrap( .~ taxa , scales="free",ncol = 2) +
  ylab('Community similarity')+xlab('Ln(Geographic distance(km)+1)')+
  theme(axis.title = element_text(size=9),axis.text = element_text(size=7),
        strip.text = element_text(size = 9),legend.position='none',panel.grid = element_blank())

library(cowplot)
combine1<-plot_grid(p_dissimilarty,p_distance,
                    labels = c('(a)','(b)'), ncol = 2, label_x = .03,label_y = 1,
                    rel_widths = c(2, 3),hjust = 0, label_size=8)
combine1

## Test sampling effort
library(vegan)
library(ggplot2)
## calculate the disimilarity matrix of environment factors based on euclid distance
env.table<-sample_data(bac.phylo)
geo_dist<-data.frame(env.table[,1:2])
geo_dist<-CalcDists(geo_dist)/100000
Geography.dist <- as.numeric(geo_dist)

# This function was modified according Meyer et al. 2018 ISME-J, DOI: 10.1038/s41396-018-0103-3, 
# which performs a single rarefaction of a community matrix, calculates the slope 
# of the distance-decay relationship of that community, then repeats iteratively. It was
# designed to explore the relationship between sampling effort and biogeographic patterns.

# The user provides: 
# 'otu' - a community matrix with samples as rows and columns as taxa, 
# 'depth'- the level to which the community should be rarefied,  
# 'geo'- a matrix of geographic distance between samples (with sample names that match 'otu'),
# 'method'- the desired community dissimilarity metric (see Vegan package, function 'vegdist' for metric options),
# 'niter' -  the number of iterations that a community matrix will be rarefied then regressed against distance.

sampling_effort_dd <- function(otu, depth, geo, method, niter){
  require(vegan)
  require(labdsv)
  slope_df <- data.frame('Slope_est'=as.numeric(), 'Std_er'=as.numeric(), 't_val'=as.numeric(), 'p'=as.numeric())
  for(i in 1:niter){
    otu_rare <-suppressWarnings(rrarefy(otu, sample=depth))
    otu_rare <- otu_rare[which(rowSums(otu_rare) >= (0.95 *depth)),] 
    geo <- geo[which(row.names(geo) %in% row.names(otu_rare)),
               which(row.names(geo) %in% row.names(otu_rare))]
    
    bray <- as.matrix(vegdist(otu_rare, method=method))
    bray[upper.tri(bray)] <- 0 
    bray_long <- dematrify(as.data.frame(bray))
    bray_long$Bray <- bray_long$abundance
    bray_long$abundance <- NULL
    bray_long$Geodist <- NA
    
    geo[upper.tri(geo)] <- 0 
    geo_long <- dematrify(as.data.frame(geo))
    dist_combined <- bray_long
    dist_combined$Geodist <- geo_long$abundance
    dist_combined$Bray <- (1- dist_combined$Bray)*100
    dist_combined$Geodist <- dist_combined$Geodist
    lmcoef <- data.frame(summary(lm(dist_combined$Bray ~ dist_combined$Geodist))$coefficient)
    slope_df[i,] <- lmcoef[2,]
  }
  return(slope_df)
}

## calculate the geographic distance based on euclid distance
env.table<-sample_data(bac.phylo)
geo_dist<-data.frame(env.table[,1:2])
geo_dist<-as.matrix(CalcDists(geo_dist)/100000)

#calculates the slope of the distance-decay relationship with 100 iterations at defferent rarefication levels.
#bacteria
otu.table <- read.delim('./data/bacteria/otutab.txt',sep="\t", row.names=1, header = T)
bac.otu.table <- as.matrix(t(otu.table))
slop_bac_rare10000 <- sampling_effort_dd(bac.otu.table, 10000, geo_dist, "bray", 100)
slop_bac_rare7000 <- sampling_effort_dd(bac.otu.table, 7000, geo_dist, "bray", 100)
slop_bac_rare3500 <- sampling_effort_dd(bac.otu.table, 3500, geo_dist, "bray", 100)

#fungi
otu.table <- read.delim('./data/fungi/otutab.txt',sep="\t", row.names=1, header = T)
fungi.otu.table <- as.matrix(t(otu.table))
slop_fungi_rare2000 <- sampling_effort_dd(fungi.otu.table, 2000, geo_dist, "bray", 100)
slop_fungi_rare1500 <- sampling_effort_dd(fungi.otu.table, 1500, geo_dist, "bray", 100)
slop_fungi_rare800 <- sampling_effort_dd(fungi.otu.table, 800, geo_dist, "bray", 100)

#protist
otu.table <- read.delim('./data/protist/otutab.txt',sep="\t", row.names=1, header = T)
protist.otu.table <- as.matrix(t(otu.table))
slop_protist_rare550 <- sampling_effort_dd(protist.otu.table, 550, geo_dist, "bray", 100)
slop_protist_rare400 <- sampling_effort_dd(protist.otu.table, 400, geo_dist, "bray", 100)
slop_protist_rare250 <- sampling_effort_dd(protist.otu.table, 250, geo_dist, "bray", 100)

#animal
otu.table <- read.delim('./data/animal/otutab.txt',sep="\t", row.names=1, header = T)
animal.otu.table <- as.matrix(t(otu.table))
slop_animal_rare280 <- sampling_effort_dd(animal.otu.table, 280, geo_dist, "bray", 100)
slop_animal_rare150 <- sampling_effort_dd(animal.otu.table, 150, geo_dist, "bray", 100)
slop_animal_rare100 <- sampling_effort_dd(animal.otu.table, 100, geo_dist, "bray", 100)

#melt the data
bac.slop <- data.frame(tax = rep('Bacteria', 300), 
                       depth = factor(c(rep(10000, 100), rep(7000, 100), rep(3500, 100))),
                       Slope_est = c(slop_bac_rare10000$Slope_est, slop_bac_rare7000$Slope_est,
                                     slop_bac_rare3500$Slope_est))
mode <- aov(Slope_est~depth, data = bac.slop)
TukeyHSD(mode)
fungi.slop <- data.frame(tax = rep('Fungi', 300), 
                         depth = factor(c(rep(2000, 100), rep(1500, 100), rep(800, 100))),
                         Slope_est = c(slop_fungi_rare2000$Slope_est, slop_fungi_rare1500$Slope_est,
                                       slop_fungi_rare800$Slope_est))
mode <- aov(Slope_est~depth, data = fungi.slop)
TukeyHSD(mode)
protist.slop <- data.frame(tax = rep('Protist', 300), 
                           depth = factor(c(rep(550, 100), rep(400, 100), rep(250, 100))),
                           Slope_est = c(slop_protist_rare550$Slope_est, slop_protist_rare400$Slope_est,
                                         slop_protist_rare250$Slope_est))
mode <- aov(Slope_est~depth, data = protist.slop)
TukeyHSD(mode)
animal.slop <- data.frame(tax = rep('Animal', 300), 
                          depth = factor(c(rep(280, 100), rep(150, 100), rep(100, 100))),
                          Slope_est = c(slop_animal_rare280$Slope_est, slop_animal_rare150$Slope_est,
                                        slop_animal_rare100$Slope_est))
mode <- aov(Slope_est~depth, data = animal.slop)
TukeyHSD(mode)
slop_df <- data.frame(rbind(bac.slop, fungi.slop, protist.slop, animal.slop))
slop_df$depth <- as.character(slop_df$depth)
slop_df$depth <- factor(slop_df$depth, ordered = T,
                        levels = c('10000', '7000', '3500',
                                   '2000', '1500', '800',
                                   '550', '400', '250',
                                   '280', '150', '100'))
slop_df$tax <- factor(slop_df$tax, ordered = T, 
                      levels = c('Bacteria', 'Fungi', 'Protist', 'Animal'))
#boxplot
library(ggplot2)
g_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 4
cols = gg_color_hue(n)
mytheme<-theme_bw()+
  theme(axis.title = element_text(size = 9),
        axis.text = element_text(size=7, colour = 'black'),
        panel.grid = element_blank(),
        #strip.background = element_rect(fill = 'skyblue'),
        strip.text = element_text(size = 9))

p_rare_slope<-ggplot(slop_df,aes(x=depth,y=Slope_est))+
  geom_boxplot(width=0.3,aes(fill=tax))+
  scale_fill_manual(values=cols)+
  ylab('DD Slope')+ xlab('Rarefication level') +
  #scale_x_discrete(limits=c('Bacteria','Fungi','Protist','Animal'))+
  #facet_wrap(variable~., scales="free",ncol = 3)+
  mytheme
p_rare_slope