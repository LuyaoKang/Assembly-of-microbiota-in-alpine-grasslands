# NULL model analysis
#function for extract the most 1000 abundant ASVs
filter.abun <- function(phylo, N) {
  require(phyloseq)
  comun <- otu_table(phylo)
  rowmean <-sapply(1:nrow(comun), function(x) mean(comun[x, ]))
  comun <- comun[order(rowmean, decreasing = TRUE), ]
  taxa_list <- rownames(comun)[1:N]
  otu_tab <- comun[rownames(comun) %in% taxa_list, ]
  otu_tab <- data.frame(t(otu_tab))
  return(otu_tab)
}

filter.tree <- function(phylo, otu_tab) {
  require(phyloseq)
  require(ape)
  tree <- phy_tree(phylo)
  tips<-colnames(otu_tab)
  tree<-keep.tip(tree, tips)
}

#extract the most 1000 abundant ASVs
#bacteria
filter.bac.tab <- filter.abun(bac.phylo, 1000)
ncol(filter.bac.tab)
sum(filter.bac.tab)/sum(otu_table(bac.phylo))
bacteria.tree <- filter.tree(bac.phylo, filter.bac.tab)
#fungi
filter.fungi.tab <- filter.abun(fungi.phylo, 1000)
ncol(filter.fungi.tab)
sum(filter.fungi.tab)/sum(otu_table(fungi.phylo))
fungi.tree <- filter.tree(fungi.phylo, filter.fungi.tab)
#protist
filter.protist.tab <- filter.abun(protist.phylo, 1000)
ncol(filter.protist.tab)
sum(filter.protist.tab)/sum(otu_table(protist.phylo))
protist.tree <- filter.tree(protist.phylo, filter.protist.tab)
#animal
filter.animal.tab <- filter.abun(animal.phylo, nrow(otu_table(animal.phylo)))
ncol(filter.animal.tab)
sum(filter.animal.tab)/sum(otu_table(animal.phylo))
animal.tree <- filter.tree(animal.phylo, filter.animal.tab)


#R codes for null model analysis modified according to Stegen et al. (2013) ####
#phylo: Phylogenetic tree of each OTU
#comun: A community table with samples as rows and OTUs as columns. 
#Beta_NTI
Beta_NTI<-function(tree, comun){
  require(picante)
  ## make sure the names on the phylogeny are ordered the same as the names in otu table
  match.phylo.comun = match.phylo.data(tree, t(comun));
  ## calculate empirical betaMNTD
  beta.mntd.weighted = as.matrix(comdistnt(t(match.phylo.comun$data), 
                                           cophenetic(match.phylo.comun$phy), abundance.weighted = T));
  identical(colnames(match.phylo.comun$data), colnames(beta.mntd.weighted)); # just a check, should be TRUE
  identical(colnames(match.phylo.comun$data), rownames(beta.mntd.weighted)); # just a check, should be TRUE
  # calculate randomized betaMNTD
  beta.reps = 9; # number of randomizations
  rand.weighted.bMNTD.comp = array(c(-9), dim = c(ncol(match.phylo.comun$data), ncol(match.phylo.comun$data), beta.reps));
  for (rep in 1:beta.reps) {
    rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(match.phylo.comun$data),
                                                          taxaShuffle(cophenetic(match.phylo.comun$phy)),
                                                          abundance.weighted = T, exclude.conspecifics = F));
    print(c(date(), rep));
  }
  weighted.bNTI = matrix(c(NA), nrow = ncol(match.phylo.comun$data), ncol = ncol(match.phylo.comun$data));
  dim(weighted.bNTI);
  for (columns in 1:(ncol(match.phylo.comun$data)-1)) {
    for (rows in (columns+1):ncol(match.phylo.comun$data)) {
      
      rand.vals = rand.weighted.bMNTD.comp[rows, columns,];
      weighted.bNTI[rows, columns] = (beta.mntd.weighted[rows, columns] - mean(rand.vals)) / sd(rand.vals);
      rm("rand.vals");
    };
  };
  rownames(weighted.bNTI) = colnames(match.phylo.comun$data);
  colnames(weighted.bNTI) = colnames(match.phylo.comun$data);
  results<-as.dist(weighted.bNTI);
  return(results)
}

#RC_bray
raup_crick= function(comun, reps=9){
  require(ecodist) 
  ## count number of sites and total species richness across all plots (gamma)
  n_sites<-nrow(comun)
  gamma<-ncol(comun)
  ##build a site by site matrix for the results, with the names of the sites in the row and col names:
  results<-matrix(data = NA, nrow = n_sites, ncol = n_sites, dimnames = list(row.names(comun), row.names(comun)))
  ##make the comun matrix into a new, pres/abs. matrix:
  ceiling(comun/max(comun))->comun.inc
  ##create an occurrence vector- used to give more weight to widely distributed species in the null model:
  occur <- apply(comun.inc, MARGIN = 2, FUN = sum)
  ##create an abundance vector- used to give more weight to abundant species in the second step of the null model:
  abundance <- apply(comun, MARGIN = 2, FUN = sum)
  ##make_null:
  ##looping over each pairwise community combination:
  for(null.one in 1:(nrow(comun)-1)){
    for(null.two in (null.one+1):nrow(comun)){
      null_bray_curtis<-NULL
      for(i in 1:reps){
        ##two empty null communities of size gamma:
        com1 <- rep(0, gamma)
        com2 <- rep(0, gamma)
        ##add observed number of species to com1, weighting by species occurrence frequencies:
        com1[sample(1:gamma, sum(comun.inc[null.one,]), replace = FALSE, prob = occur)]<-1
        com1.samp.sp = sample(which(com1>0), (sum(comun[null.one,])-sum(com1)), replace = TRUE, prob = abundance[which(com1>0)]);
        com1.samp.sp = cbind(com1.samp.sp,1); # head(com1.samp.sp);
        com1.sp.counts = as.data.frame(tapply(com1.samp.sp[,2], com1.samp.sp[,1], FUN = sum)); colnames(com1.sp.counts) = 'counts'; # head(com1.sp.counts);
        com1.sp.counts$sp = as.numeric(rownames(com1.sp.counts)); # head(com1.sp.counts);
        com1[com1.sp.counts$sp] = com1[com1.sp.counts$sp] + com1.sp.counts$counts; # com1;
        #sum(com1) - sum(spXsite[null.one,]); ## this should be zero if everything work properly
        rm('com1.samp.sp','com1.sp.counts');			
        ##same for com2:
        com2[sample(1:gamma, sum(comun.inc[null.two,]), replace = FALSE, prob = occur)]<-1
        com2.samp.sp = sample(which(com2>0),(sum(comun[null.two,])-sum(com2)),replace = TRUE, prob = abundance[which(com2>0)]);
        com2.samp.sp = cbind(com2.samp.sp,1); # head(com2.samp.sp);
        com2.sp.counts = as.data.frame(tapply(com2.samp.sp[,2], com2.samp.sp[,1], FUN = sum)); colnames(com2.sp.counts) = 'counts'; # head(com2.sp.counts);
        com2.sp.counts$sp = as.numeric(rownames(com2.sp.counts)); # head(com2.sp.counts);
        com2[com2.sp.counts$sp] = com2[com2.sp.counts$sp] + com2.sp.counts$counts; # com2;
        # sum(com2) - sum(spXsite[null.two,]); ## this should be zero if everything work properly
        rm('com2.samp.sp','com2.sp.counts');
        null.comun = rbind(com1, com2); # null.comun;
        ##calculate null bray curtis
        null_bray_curtis[i] = distance(null.comun, method = 'bray-curtis');
      }; # end reps loop
      ## empirically observed bray curtis
      obs.bray = distance(comun[c(null.one, null.two),], method = 'bray-curtis');
      ##how many null observations is the observed value tied with?
      num_exact_matching_in_null = sum(null_bray_curtis == obs.bray);
      ##how many null values are smaller than the observed *dissimilarity*?
      num_less_than_in_null = sum(null_bray_curtis < obs.bray);
      rc = ((num_less_than_in_null +(num_exact_matching_in_null)/2)/reps)
      ##modification of raup crick standardizes the metric to range from -1 to 1 instead of 0 to 1
      rc = (rc-.5)*2
      results[null.two,null.one] = round(rc, digits=2); ##store the metric in the results matrix
      print(c(null.one, null.two,date()));
    }; ## end null.two loop
  }; ## end null.one loop
  
  results <- as.dist(results)
  return(results)
}


#calculate the βNTI and RC_bray, this procedure is excessively time consuming. 
#In our study, it spent about half a month to finish this procedure.

#Beta_NTI_bacteria <- Beta_NTI(bacteria.tree, filter.bac.tab)
#Beta_NTI_bacteria_1000 <- as.matrix(Beta_NTI_bacteria)
#write.csv(Beta_NTI_bacteria_1000,"./result/null_model/bacteria_bNTI_1000.csv", quote = F);

#Beta_NTI_fungi <- Beta_NTI(fungi.tree, filter.fun.tab)
#Beta_NTI_fungi_1000 <- as.matrix(Beta_NTI_fungi)
#write.csv(Beta_NTI_fungi_1000,"./result/null_model/fungi_bNTI_1000.csv", quote = F);

#protist
#Beta_NTI_protists <- Beta_NTI(animal.tree, filter.anim.tab)
#Beta_NTI_protists_1000 <- as.matrix(Beta_NTI_protists)
#write.csv(Beta_NTI_protists_1000,"./result/null_model/protists_bNTI_1000.csv", quote = F);

#Beta_NTI_animal <- Beta_NTI(animal.tree, filter.anim.tab)
#Beta_NTI_animal_1000 <- as.matrix(Beta_NTI_animal)
#write.csv(b_NTI_animal_1000,"./result/null_mode/animal_bNTI_1000.csv", quote = F);

#calculate the RC_bray
#rcbray_bacteria <- raup_crick(filter.bac.tab)
#rcbray_bacteria_1000<-as.matrix(rcbray_bacteria)
#write.csv(rcbray_bacteria_1000,"./result/null_model/raup_crick_bacteria_1000.csv", quote = F);

#rcbray_fungi <- raup_crick(filter.fun.tab)
#rcbray_fungi_1000<-as.matrix(rcbray_fungi)
#write.csv(rcbray_fungi,"./result/null_model/raup_crick_fungi_1000.csv",quote = F);

#rcbray_protists<-raup_crick(filter.proti.tab)
#rcbray_protists_1000<-as.matrix(rcbray_protists)
#write.csv(rcbray_protists_1000,"./result/null_model/raup_crick_protists_1000.csv", quote = F);

#rcbray_animal <- raup_crick(filter.anim.tab)
#rcbray_animal_1000<-as.matrix(rcbray_animal)
#write.csv(rcbray_animal_1000,"./result/null_model/animal/raup_crick_animal_1000.csv", quote = F);


#read βNTI and RC_bray matrixes
Beta_NTI_bacteria <- read.csv("./result/null_model/bacteria_bNTI_1000.csv", header = T, row.names = 1, stringsAsFactors = F);
rcbray_bacteria <- read.csv("./result/null_model/raup_crick_bacteria_1000.csv", header = T, row.names = 1, stringsAsFactors = F);

Beta_NTI_fungi <- read.csv("./result/null_model/fungi_bNTI_1000.csv", header = T, row.names = 1, stringsAsFactors = F);
rcbray_fungi <- read.csv("./result/null_model/raup_crick_fungi_1000.csv", header = T, row.names = 1, stringsAsFactors = F);

Beta_NTI_protists <- read.csv("./result/null_model/protist_bNTI_1000.csv", header = T, row.names = 1, stringsAsFactors = F);
rcbray_protists <- read.csv("./result/null_model/raup_crick_protists_1000.csv", header = T, row.names = 1, stringsAsFactors = F);

Beta_NTI_animal <- read.csv("./result/null_model/animal_bNTI_1000.csv", header = T, row.names = 1, stringsAsFactors = F);
rcbray_animal <- read.csv("./result/null_model/raup_crick_animal_1000.csv", header = T, row.names = 1, stringsAsFactors = F);
#determine the relative contribution of each process in shaping community
null.index <- function(bNTI.df, rcbray.df){
  bNTI.dist <- as.dist(bNTI.df)
  RC.dist <- as.dist(rcbray.df)
  library(NST)
  null.value <- cbind(dist.3col(bNTI.dist), dist.3col(RC.dist)[3])
  colnames(null.value) <- c('name1', 'name2', 'bNTI', 'RC_bray')
  heterogeneous_selection <- nrow(subset(null.value, bNTI > 2))
  homogeneous_selection <- nrow(subset(null.value, bNTI < (-2)))
  dispersal_limitation <- nrow(subset(null.value,bNTI < 2 & bNTI > (-2) & RC_bray > 0.95))
  homogenizing_dispersal <- nrow(subset(null.value,bNTI < 2 & bNTI > (-2) & RC_bray < (-0.95)))
  undominated <- nrow(subset(null.value, bNTI < 2 & bNTI > (-2) & RC_bray < 0.95 & RC_bray > (-0.95)))
  null.pro <- data.frame(homogeneous_selection, heterogeneous_selection, dispersal_limitation,
                       homogenizing_dispersal, undominated)
  process <- colnames(null.pro)
  Proportion <- as.numeric(null.pro[1, ]/sum(null.pro[1, ]))
  df <- data.frame(process,Proportion)
  return(df)
}

#determine the relative contribution of each process in shaping community
bacteria.null <- null.index(Beta_NTI_bacteria, rcbray_bacteria)
fungi.null <- null.index(Beta_NTI_fungi, rcbray_fungi)
protists.null <- null.index(Beta_NTI_protists, rcbray_protists)
animal.null <- null.index(Beta_NTI_animal, rcbray_animal)

#transform the wide dataframe to long format
null.3col <- function(bNTI.df){
  bnti_mat <- as.matrix(bNTI.df)
  bNTI.mat <- matrix(0,nrow = 147, ncol = 147)
  bNTI.mat[lower.tri(bNTI.mat)] <- bnti_mat[lower.tri(bnti_mat, diag=TRUE)]
  bNTI.dist <- as.dist(bNTI.df)
  library(NST)
  null.value<-dist.3col(bNTI.dist)
  colnames(null.value)<-c('name1','name2','bNTI')
  return(null.value)
}

bNTI_bacteria <- null.3col(Beta_NTI_bacteria)
bNTI_fungi <- null.3col(Beta_NTI_fungi)
bNTI_protists <- null.3col(Beta_NTI_protists)
bNTI_animal <- null.3col(Beta_NTI_animal)

#density plot of βNTI
#set colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 4
cols = gg_color_hue(n)

bNTI_all<-data.frame(taxa=c(rep('Bacteria', nrow(bNTI_bacteria)),
                            rep('Fungi', nrow(bNTI_fungi)),
                            rep('Protists', nrow(bNTI_protists)),
                            rep('Animal', nrow(bNTI_animal))),
                     rbind(bNTI_bacteria, bNTI_fungi,bNTI_protists, bNTI_animal))
bNTI_all$taxa <- factor(bNTI_all$taxa, ordered=T, levels = c('Bacteria', 'Fungi', 'Protists', 'Animal'))
#plot for the frequence distribution of bNTI
bNTI_freq_plot <- ggplot(bNTI_all, aes(bNTI, fill = taxa)) + 
  geom_histogram(color = "black", position= "identity", bins = 20) +
  facet_wrap( ~ taxa, scales = 'free', ncol = 2) +
  scale_fill_manual(values=c(cols[1], cols[2], cols[3], cols[4])) +
  ylab('Frequency') + xlab('βNTI') +
  theme_bw() +
  theme(axis.title = element_text(size=9), axis.text = element_text(size=7, colour = 'black'),
        strip.text = element_text(size = 9), legend.position='none', panel.grid = element_blank())
bNTI_freq_plot

#pie plot for the proportion of each process
doughnut.plot<-function(data){
  # load library
  library(ggplot2)
  library(dplyr)
  # Compute percentages
  data$process <- factor(data$process, ordered = T, 
                         levels = c('heterogeneous_selection', 'homogeneous_selection',
                                    'dispersal_limitation', 'homogenizing_dispersal', 'undominated'))
  data$Proportion <- round(data$Proportion*100, 2)
  
  count.data <- data %>%
    arrange(desc(process)) %>%
    mutate(lab.ypos = cumsum(Proportion) - 0.5*Proportion)
  mycols <- c("#CD534CFF", "#EFC000FF", "#0073C2FF", "#5d8aa8", '#009E73')
  # Make the plot
  p<-ggplot(count.data, aes(x = 2, y = Proportion, fill = process)) +
    geom_bar(stat = "identity", color = "white") +
    coord_polar(theta = "y", start = 0)+
    geom_text(aes(y = lab.ypos, label = Proportion), color = "black")+
    scale_fill_manual(values = mycols) +
    theme(legend.title = element_text(size = 8),
          legend.text = element_text(size = 7))+
    theme_void()+
    xlim(0.5, 2.5)
  return(p)
}

B.doug <- doughnut.plot(bacteria.null)
F.doug <- doughnut.plot(fungi.null)
P.doug <- doughnut.plot(protists.null)
A.doug <- doughnut.plot(animal.null)

library(cowplot)
plot_grid(B.doug, F.doug, P.doug, A.doug,
          labels = "auto", ncol = 2, nrow = 2, label_x = .01, 
          label_y = 1.005, hjust = 0, label_size = 17, align = "v")
