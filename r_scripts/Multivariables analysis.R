# Multivariables analysis
## variation partitioning with two-way PERMANOVA
#Regression approach with distance as response and raw predictors
#reference:https://media.nature.com/original/nature-assets/ismej/journal/v12/n2/extref/ismej2017183x2.txt
DistResponseReg <-function (DistMatrix,X) {
  DistMatrix <- as.matrix(DistMatrix)
  X <- as.matrix(X)
  n <- nrow(DistMatrix)
  p <- ncol(X)
  row.wt = rep(1, nrow(DistMatrix))
  col.wt = rep(1, ncol(DistMatrix))
  st <- sum(col.wt)
  sr <- sum(row.wt)
  row.wt <- row.wt/sr
  col.wt <- col.wt/st
  DistMatrix <- -0.5*(DistMatrix*DistMatrix)
  row.mean <- apply(row.wt * DistMatrix, 2, sum)
  col.mean <- apply(col.wt *t(DistMatrix), 2, sum)
  col.mean <- col.mean - sum(row.mean * col.wt)
  DistMatrix <- sweep(DistMatrix, 2, row.mean)
  G <- t(sweep(t(DistMatrix), 2, col.mean))
  H<-X%*%solve(t(X)%*%X)%*%t(X)
  I<-diag(n)
  predicted <- H%*%G%*%H
  residuals <- (I-H)%*%G%*%(I-H)
  MS_regression<-sum(diag(predicted))/p
  MS_residual<-sum(diag(residuals))/(n-p)
  F<-MS_regression/MS_residual
  MS_Total=sum(diag(G))/n;
  RsqAdj=1-MS_residual/MS_Total;
  
  result <- list(Rsq_Adj=RsqAdj,F_value=F,res_matrix=residuals,pred_matrix=predicted)
  return(result)
  
}

# variation partitioning with distance matrix as response
# and X (e.g., environment) and W (e.g., space or connectivity)
VarPartDistResponse<-function (Dist_Matrix,X,W,Number_Permutations) {
  X <- as.matrix(X)
  X <- apply(X, 2, scale)
  W <- as.matrix(W)
  W <- apply(W, 2, scale)
  Dist_Matrix <- as.matrix(Dist_Matrix)
  Number_Predictors_X <- ncol(X)
  Number_Predictors_W <- ncol(W)
  n<-nrow(X)
  
  XW <- as.matrix(cbind(X,W));
  result <- DistResponseReg(Dist_Matrix,XW)
  abc <- result$Rsq_Adj
  Fabc <- result$F_value
  
  result <- DistResponseReg(Dist_Matrix,X)
  ab <- result$Rsq_Adj
  Fab <- result$F_value
  residuals_X <- result$res_matrix 
  predicted_X <- result$pred_matrix
  
  result <- DistResponseReg(Dist_Matrix,W)
  bc <- result$Rsq_Adj
  Fbc <- result$F_value
  residuals_W <- result$res_matrix 
  predicted_W <- result$pred_matrix
  
  # unique fraction of contribution related to X
  a <- abc - bc
  # unique fraction of contribution related to W
  c <- abc- ab
  # common fraction of contribution between X and W
  b <- abc - a - c
  # residual fraction
  d <- 1-abc
  
  Fa=(a/Number_Predictors_X)/(d/(n-Number_Predictors_X-Number_Predictors_W));
  Fc=(c/Number_Predictors_W)/(d/(n-Number_Predictors_X-Number_Predictors_W));
  
  Prob_abc=1/Number_Permutations; Prob_ab=1/Number_Permutations; Prob_bc=1/Number_Permutations; Prob_a=1/Number_Permutations; Prob_c=1/Number_Permutations;
  
  # permutations test
  for (i in 1:(Number_Permutations-1)) {
    # testing fraction a; notice that we permute the residual values in W and not in X
    permuted_rows=sample(n,replace=FALSE)
    # permuting the residual matrix, which is from the distance, and hence the need to permute
    # rows and columns in the same way, hence the use of permuted_rows for columns and rows below
    # testing fraction a
    # Yperm=predicted_W+residuals_W[permuted_rows,permuted_rows] # implement permutation of residuals in the future
    
    # testing fraction a
    result <- DistResponseReg(Dist_Matrix,XW[permuted_rows,])
    abcRnd <- result$Rsq_Adj
    FabcRnd <- result$F_value
    result <- DistResponseReg(Dist_Matrix,W[permuted_rows,])
    bcRnd <- result$Rsq_Adj
    FbcRnd <- result$F_value
    aRnd=abcRnd-bcRnd;
    dRnd=1-abcRnd;
    FaRnd=(aRnd/Number_Predictors_X)/(dRnd/(n-Number_Predictors_X-Number_Predictors_W));
    if (FaRnd >= Fa) {Prob_a<-Prob_a+1/Number_Permutations}
    
    # testing fraction c
    result <- DistResponseReg(Dist_Matrix,X[permuted_rows,])
    abRnd <- result$Rsq_Adj
    FabRnd <- result$F_value
    cRnd=abcRnd-abRnd;
    FcRnd=(cRnd/Number_Predictors_W)/(dRnd/(n-Number_Predictors_X-Number_Predictors_W));
    if (FcRnd >= Fc) {Prob_c<-Prob_c+1/Number_Permutations}
    
    # testing abc
    if (FabcRnd >= Fabc) {Prob_abc<-Prob_abc+1/Number_Permutations}
    # testing ab
    if (FabRnd >= Fab) {Prob_ab<-Prob_ab+1/Number_Permutations}
    # testing bc
    if (FbcRnd >= Fbc) {Prob_bc<-Prob_bc+1/Number_Permutations}
    
  }
  
  result <- mat.or.vec(7,2)
  
  result[1,1] <- abc
  result[2,1] <- ab
  result[3,1] <- bc
  result[4,1] <- a
  result[5,1] <- c
  result[6,1] <- b
  result[7,1] <- d
  
  result[1,2] <- Prob_abc
  result[2,2] <- Prob_ab
  result[3,2] <- Prob_bc
  result[4,2] <- Prob_a
  result[5,2] <- Prob_c
  result[6,2] <- NA
  result[7,2] <- NA
  
  colnames(result) <- c("Estimate","p-value")
  rownames(result) <- c("abc","ab","bc","a","c","b","d")
  result<-as.data.frame(result)
  return(result)
}

# model selection for linear model with distance as response
#reference:
SelectionDistResponseReg<-function (Dist_Matrix,X,Number_Permutations,alpha) {
  X <- as.matrix(X)
  Number_Predictors<-ncol(X)
  n<-nrow(X)
  Dist_Matrix <- as.matrix(Dist_Matrix)
  X <- apply(X, 2, scale)
  result<-DistResponseReg(Dist_Matrix,X);
  F_Observed<-result$F_value
  # global test with all predictors
  Prob_Global<-1/Number_Permutations
  for (i in 1:(Number_Permutations-1)) {
    X_permuted<-X[sample(n,replace=FALSE),]
    result<-DistResponseReg(Dist_Matrix,X_permuted)
    F_Random<-result$F_value
    if (F_Random >= F_Observed) {Prob_Global<-Prob_Global+1/Number_Permutations}
  }
  
  # set of these as NA in the global model is not significant
  Variables_In_Model <- NA
  Rsq_Final_Model <- NA
  if (Prob_Global < alpha) {
    # calculate contribution of each predictor separately as they're orthgonal
    F_Ind_X<-mat.or.vec(Number_Predictors,1)
    for (i in 1:Number_Predictors) {
      F_Ind_X[i]<-DistResponseReg(Dist_Matrix,X[,i])[2]
    }
    # start selection
    Variables_In_Model<-as.matrix(which.max(F_Ind_X))
    F_Ind_X[Variables_In_Model[1]]<-NA
    #Original_Columns=as.matrix(1:n)
    #Original_Columns[Variable_In_Model]<-NA
    found <- FALSE
    while (found == FALSE) { 
      # this could be made faster by not calculating again the F for the first variable entered, but this version is more general
      # contrast current model with the largest variable contribution not entered in the model
      candidate_model <- cbind(X[,Variables_In_Model],X[,which.max(F_Ind_X)])
      result_candidate <- DistResponseReg(Dist_Matrix,candidate_model)
      F_candidate_Obs <- result_candidate[2]
      # test wheter the entered variable improves fit
      Prob_F <- 1/Number_Permutations
      for (i in 1:(Number_Permutations-1)) {
        candidate_predictor <- which.max(F_Ind_X)
        candidate_model <- cbind(X[,Variables_In_Model],X[sample(n,replace=FALSE),candidate_predictor])
        result_candidate <- DistResponseReg(Dist_Matrix,candidate_model);
        F_candidate_Rnd <- result_candidate$F_value
        if (F_candidate_Rnd >= F_candidate_Obs) {Prob_F<-Prob_F+1/Number_Permutations}
      }
      if (Prob_F > alpha) {found <- TRUE} else {
        F_Ind_X[candidate_predictor] <- NA
        Variables_In_Model <- append(Variables_In_Model,candidate_predictor)
      }
      
      # nrow(na.omit(OriginalColumns))
    }
    result<-DistResponseReg(Dist_Matrix,X[,Variables_In_Model])
    Rsq_Final_Model<-result$Rsq_Adj
  }
  
  result <- list(Global_P=Prob_Global,Rsq_Selected_Model=Rsq_Final_Model,Selected_Variables=Variables_In_Model)
  return(result)
}

#Construct a variation partition function 
vp <- function(phylo){
  set.seed(999)
  require(vegan)
  siteXspe <- t(otu_table(phylo))
  siteXspe <- as.matrix(siteXspe)
  env <- sample_data(phylo)
  env_vars<-env[,-c(1:4)]
  Geography_vars<-env[,c(1:2)]
  otu_tab <- decostand(siteXspe, 'hellinger')
  commu.dist<-vegdist(otu_tab, 'bray',upper=F)
  pcnm_vars <- (pcnm(dist(Geography_vars)))$vectors
  # check collinearities among all variables
  env_df <- as.matrix(env_vars)
  env_df <- decostand(env_df,"standardize")
  pcnm_df <- as.matrix(pcnm_vars)
  all.vars <- data.frame(env_df, pcnm_df)
  ord <- dbrda(otu_tab~., all.vars)
  ord.vif <- vif.cca(ord)
  #remove the varibles with vif >20
  all.vars <- all.vars[,names(ord.vif[ord.vif < 20])]
  env_mat <- env_df[,c(intersect(colnames(env_df),colnames(all.vars)))]
  pcnm_mat <- pcnm_df[,c(intersect(colnames(pcnm_df),colnames(all.vars)))]
  #variables selection
  env_sel_vars<-env_mat[,SelectionDistResponseReg(commu.dist, env_mat,999,0.05)$Selected_Variables]
  pcnm_sel_vars<-pcnm_mat[,SelectionDistResponseReg(commu.dist, pcnm_mat,999,0.05)$Selected_Variables]
  #conduct VPA
  vp<-VarPartDistResponse(commu.dist,env_sel_vars,pcnm_sel_vars,999)
  results<-list(vp,env_sel_vars,pcnm_sel_vars)
  return(results)
}

vp.bacteria<-vp(bac.phylo)
vp.fungi<-vp(fungi.phylo)
vp.protists<-vp(protist.phylo)
vp.animal<-vp(animal.phylo)

bacteria_vp<-vp.bacteria[[1]]
bacteria_env_df<-vp.bacteria[[2]]
bacteria_geo_df<-vp.bacteria[[3]]
colnames(bacteria_env_df)
colnames(bacteria_geo_df)

fungi_vp<-vp.fungi[[1]]
fungi_env_df<-vp.fungi[[2]]
fungi_geo_df<-vp.fungi[[3]]
colnames(fungi_env_df)
colnames(fungi_geo_df)

protists_vp<-vp.protists[[1]]
protists_env_df<-vp.protists[[2]]
protists_geo_df<-vp.protists[[3]]
colnames(protists_env_df)
colnames(protists_geo_df)

animal_vp<-vp.animal[[1]]
animal_env_df<-vp.animal[[2]]
animal_geo_df<-vp.animal[[3]]
colnames(animal_env_df)
colnames(animal_geo_df)

#plot
#first, variable partition plot
theme_set(theme_grey())
library(reshape)
library(ggplot2)
Sorting_effect=c(bacteria_vp[4,1],fungi_vp[4,1],protists_vp[4,1],animal_vp[4,1])*100
dispersal_effect=c(bacteria_vp[5,1],fungi_vp[5,1],protists_vp[5,1],animal_vp[5,1])*100
sorting.dispersal=c(bacteria_vp[6,1],fungi_vp[6,1],protists_vp[6,1],animal_vp[6,1])*100
sort.dispersal.ration=c(Sorting_effect/dispersal_effect)
vp.dat<-data.frame(taxa=c('Bacteria','Fungi','Protists','Animal'),Sorting_effect,
                   dispersal_effect,sorting.dispersal,sort.dispersal.ration)
vp.dat$taxa<-factor(vp.dat$taxa,ordered = T,
                    levels=c('Bacteria','Fungi','Protists','Animal'))

melted <- melt(vp.dat, id.vars=c("taxa"))
vp.dat <- melted[1:12, ]
sdr.dat <- melted[13:16, ]

# Stacked plot for vpa
#set color
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 4
cols = gg_color_hue(4)
vp.plot<-ggplot(vp.dat, aes(fill= variable, y=value, x=taxa)) +
  geom_bar(position="stack", stat="identity", width = 0.75)+
  #scale_fill_manual(values=cols)+
  xlab('Taxa')+ylab('Explained variation (%)')+
  scale_x_discrete(limits=c('Bacteria','Fungi','Protists','Animal'))+
  theme(axis.title = element_text(size=9),
        axis.text = element_text(size=7, colour = 'black'),
        strip.text = element_text(size = 9),legend.position=c(0.70, 0.85),
        panel.grid = element_blank())

#barplot for Sorting/dispersal effect ratio
sdr.plot<-ggplot(sdr.dat,aes(x=taxa, y=value)) +
  geom_bar(stat = 'identity',aes(fill='a1caf1'), width = 0.75)+
  scale_fill_manual(values=cols)+
  xlab('Taxa')+ylab('Sorting/dispersal effect ratio')+
  scale_x_discrete(limits=c('Bacteria','Fungi','Protists','Animal'))+
  theme(axis.title = element_text(size=9),
        axis.text = element_text(size=7, colour = 'black'),
        strip.text = element_text(size = 9),legend.position ='none',
        panel.grid = element_blank())
library(cowplot)
plot_grid(vp.plot, sdr.plot,
          labels = c("(a)", "(b)"), ncol = 2, nrow = 1, label_x = .01,
          label_y = 1.005,hjust = 0, label_size=8,align = "v")
