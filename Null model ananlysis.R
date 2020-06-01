#Beta_NTI
Beta_NTI<-function(tree,otu_table){
  require(picante)
  ## make sure the names on the phylogeny are ordered the same as the names in otu table
  match.phylo.otu = match.phylo.data(tree, t(otu_table));
  ## calculate empirical betaMNTD
  beta.mntd.weighted = as.matrix(comdistnt(t(match.phylo.otu$data),cophenetic(match.phylo.otu$phy),abundance.weighted=T));
  identical(colnames(match.phylo.otu$data),colnames(beta.mntd.weighted)); # just a check, should be TRUE
  identical(colnames(match.phylo.otu$data),rownames(beta.mntd.weighted)); # just a check, should be TRUE
  # calculate randomized betaMNTD
  beta.reps = 999; # number of randomizations
  rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(match.phylo.otu$data),ncol(match.phylo.otu$data),beta.reps));
  for (rep in 1:beta.reps) {
    rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(match.phylo.otu$data),
                                                          taxaShuffle(cophenetic(match.phylo.otu$phy)),
                                                          abundance.weighted=T,exclude.conspecifics = F));
    print(c(date(),rep));
  }
  weighted.bNTI = matrix(c(NA),nrow=ncol(match.phylo.otu$data),ncol=ncol(match.phylo.otu$data));
  dim(weighted.bNTI);
  for (columns in 1:(ncol(match.phylo.otu$data)-1)) {
    for (rows in (columns+1):ncol(match.phylo.otu$data)) {
      
      rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
      weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
      rm("rand.vals");
    };
  };
  rownames(weighted.bNTI) = colnames(match.phylo.otu$data);
  colnames(weighted.bNTI) = colnames(match.phylo.otu$data);
  results<-as.dist(weighted.bNTI);
  return(results)
}


#RC_bray
raup_crick= function(spXsite, reps=999){
  require(ecodist)
  ## count number of sites and total species richness across all plots (gamma)
  n_sites<-nrow(spXsite)
  gamma<-ncol(spXsite)
  ##build a site by site matrix for the results, with the names of the sites in the row and col names:
  results<-matrix(data=NA, nrow=n_sites, ncol=n_sites, dimnames=list(row.names(spXsite), row.names(spXsite)))
  ##make the spXsite matrix into a new, pres/abs. matrix:
  ceiling(spXsite/max(spXsite))->spXsite.inc
  ##create an occurrence vector- used to give more weight to widely distributed species in the null model:
  occur<-apply(spXsite.inc, MARGIN=2, FUN=sum)
  ##create an abundance vector- used to give more weight to abundant species in the second step of the null model:
  abundance<-apply(spXsite, MARGIN=2, FUN=sum)
  ##make_null:
  ##looping over each pairwise community combination:
  for(null.one in 1:(nrow(spXsite)-1)){
    for(null.two in (null.one+1):nrow(spXsite)){
      null_bray_curtis<-NULL
      for(i in 1:reps){
        ##two empty null communities of size gamma:
        com1<-rep(0,gamma)
        com2<-rep(0,gamma)
        ##add observed number of species to com1, weighting by species occurrence frequencies:
        com1[sample(1:gamma, sum(spXsite.inc[null.one,]), replace=FALSE, prob=occur)]<-1
        com1.samp.sp = sample(which(com1>0),(sum(spXsite[null.one,])-sum(com1)),replace=TRUE,prob=abundance[which(com1>0)]);
        com1.samp.sp = cbind(com1.samp.sp,1); # head(com1.samp.sp);
        com1.sp.counts = as.data.frame(tapply(com1.samp.sp[,2],com1.samp.sp[,1],FUN=sum)); colnames(com1.sp.counts) = 'counts'; # head(com1.sp.counts);
        com1.sp.counts$sp = as.numeric(rownames(com1.sp.counts)); # head(com1.sp.counts);
        com1[com1.sp.counts$sp] = com1[com1.sp.counts$sp] + com1.sp.counts$counts; # com1;
        #sum(com1) - sum(spXsite[null.one,]); ## this should be zero if everything work properly
        rm('com1.samp.sp','com1.sp.counts');			
        ##same for com2:
        com2[sample(1:gamma, sum(spXsite.inc[null.two,]), replace=FALSE, prob=occur)]<-1
        com2.samp.sp = sample(which(com2>0),(sum(spXsite[null.two,])-sum(com2)),replace=TRUE,prob=abundance[which(com2>0)]);
        com2.samp.sp = cbind(com2.samp.sp,1); # head(com2.samp.sp);
        com2.sp.counts = as.data.frame(tapply(com2.samp.sp[,2],com2.samp.sp[,1],FUN=sum)); colnames(com2.sp.counts) = 'counts'; # head(com2.sp.counts);
        com2.sp.counts$sp = as.numeric(rownames(com2.sp.counts)); # head(com2.sp.counts);
        com2[com2.sp.counts$sp] = com2[com2.sp.counts$sp] + com2.sp.counts$counts; # com2;
        # sum(com2) - sum(spXsite[null.two,]); ## this should be zero if everything work properly
        rm('com2.samp.sp','com2.sp.counts');
        null.spXsite = rbind(com1,com2); # null.spXsite;
        ##calculate null bray curtis
        null_bray_curtis[i] = distance(null.spXsite,method='bray-curtis');
      }; # end reps loop
      ## empirically observed bray curtis
      obs.bray = distance(spXsite[c(null.one,null.two),],method='bray-curtis');
      ##how many null observations is the observed value tied with?
      num_exact_matching_in_null = sum(null_bray_curtis==obs.bray);
      ##how many null values are smaller than the observed *dissimilarity*?
      num_less_than_in_null = sum(null_bray_curtis<obs.bray);
      rc = ((num_less_than_in_null +(num_exact_matching_in_null)/2)/reps)
      ##modification of raup crick standardizes the metric to range from -1 to 1 instead of 0 to 1
      rc = (rc-.5)*2
      results[null.two,null.one] = round(rc,digits=2); ##store the metric in the results matrix
      print(c(null.one,null.two,date()));
    }; ## end null.two loop
  }; ## end null.one loop
  
  results<-as.dist(results)
  return(results)
}

#determine the relative importance of each ecological process based on the βNTI and RC_bray indexes. bNTI.df and rcbray.df 
#indicate the bNTI and raup_crick (based on bray-curties index) dataframe, respectively.

null.index<-function(bNTI.df,rcbray.df){
  #coerce bNTI.df to a distance matrix
  bnti_mat <- as.matrix(bNTI.df)
  bNTI.mat<-matrix(0,nrow = nrow(bNTI.df),ncol = ncol(bNTI.df))
  bNTI.mat[lower.tri(bNTI.mat)] <- bnti_mat[lower.tri(bnti_mat, diag=TRUE)]
  bNTI.dist<-as.dist(bNTI.mat)
  #coerce rcbray.df to a distance matrix
  rcbray.mat <- as.matrix(rcbray.df)
  RC.mat<-matrix(0,nrow = nrow(rcbray.df),ncol = ncol(rcbray.df))
  RC.mat[lower.tri(RC.mat)] <- rcbray.mat [lower.tri(rcbray.mat, diag=TRUE)]
  RC.dist<-as.dist(RC.mat)
  #Transform distance matrix to a 3-column matrix in which the first 2 columns indicate the pairwised sites number.
  library(NST)
  null.value<-cbind(dist.3col(bNTI.dist),dist.3col(RC.dist)[3])
  colnames(null.value)<-c('name1','name2','bNTI','RC_bray')
  #partition the pairwise comparisons of βNTI and RC_bray that were assigned to different ecological processes.
  homogeneous_selection<-nrow(subset(null.value,bNTI>2))
  heterogeneous_selection<-nrow(subset(null.value,bNTI<(-2)))
  dispersal_limitation<-nrow(subset(null.value,bNTI<2 & bNTI>(-2) & RC_bray>0.95))
  homogenizing_dispersal<-nrow(subset(null.value,bNTI<2 & bNTI>(-2) & RC_bray<(-0.95)))
  undominated<-nrow(subset(null.value,bNTI<2 & bNTI>(-2) & RC_bray<0.95 & RC_bray>(-0.95)))
  #determine the relative importance of each ecological process.
  null.pro<-data.frame(homogeneous_selection,heterogeneous_selection,dispersal_limitation,
                       homogenizing_dispersal,undominated)
  process<-colnames(null.pro)
  Proportion<-as.numeric(null.pro[1,]/sum(null.pro[1,]))
  df<-data.frame(process,count)
  return(df)
}
