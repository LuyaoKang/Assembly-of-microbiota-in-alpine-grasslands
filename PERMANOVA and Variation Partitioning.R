##############################################
############### 2-way PERMANOVA ##############
##############################################

# Regression approach with distance as response and raw predictors
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
  return(result)
}
