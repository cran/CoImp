#############################################################
# ALGORITMO DI IMPUTAZIONE DATI MAR VIA FUNZIONE COPULA     #
#                                                           #
# INPUT:data matrix, vector of probabilities, K records more #
#similar to the missing one, distance measure,    #
#number of uniform numbers to numerically construct the empirical conditional #
#copula, set.seed for reproducibility  ###


### NPCoImp
### A NON PARAMETRIC COPULA BASED IMPUTATION METHOD
##
##  The authors of this software are
##  F. Marta L. Di Lascio, and
##  Aurora Gatto, Copyright (c) 2024

##  Permission to use, copy, modify, and distribute this software for any
##  purpose without fee is hereby granted, provided that this entire notice
##  is included in all copies of any software which is or includes a copy
##  or modification of this software and in all copies of the supporting
##  documentation for such software.
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 2 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  A copy of the GNU General Public License is available at
##  http://www.r-project.org/Licenses/

## ***************************************************************************************************

setClass("NPCoImp",
         representation(Imputed.matrix  = "matrix"
                        ,Selected.alpha  = "vector"
                        ,numFlat      = "numeric"
                        ),
         prototype = list(Imputed.matrix = matrix(0,0,0)
                        ,Selected.alpha = vector()
                        ,numFlat    = numeric()
                        )
         )

## ***************************************************************************************************

setMethod(
    f="show",
    signature="NPCoImp",
    definition=function(object){
        out <- object
        cat (" Main output of the function NPCoImp \n")
        cat (" -------------------------------------------------------------------------- \n")
        cat (" Imputed data matrix : \n")
        print(out@"Imputed.matrix")
        cat (" -------------------------------------------------------------------------- \n")
        cat (" Best quantile for the imputation : \n")
        print(out@"Selected.alpha")
        cat (" -------------------------------------------------------------------------- \n")
        cat ("Number of flat conditional empirical copulas : \n")
        print(out@"numFlat")
        cat (" -------------------------------------------------------------------------- \n")
    }
)

## ***************************************************************************************************

NPCoImp <- function(X, Psi=seq(0.05,0.45,by=0.05), K=7, method="gower", N=1000, seed=set.seed(7)){
  if(is.matrix(X)==FALSE){stop("Please, provide a matrix data object in input")}
  if(length(Psi)<1){stop("Please, provide at least one value for Psi")}
  if(max(Psi)>=0.5){stop("The (a)symmetry is investigated around 0.5, thus Psi should be less than 0.5")} 
  if(min(Psi)<0){stop("The minimum value for Psi should be greater than 0")} 
  if(length(Psi)==1 && Psi==0){cat("Median of the conditional empirical copula is used to impute data (Psi=0)","\n")}
  if(K%%1!=0){
    K <- round(K, digits=0)
    cat("Note that K has been rounded since the provided value is not integer.","\n")
  }
  if(K==0){stop("Please, provide a value different from zero for K")}
  if(N<500){cat("Note that the conditional empirical copula is poorly estimated. To improve the estimation use a larger N.","\n")}
  #
  lPsi <- length(Psi)
  x.miss <- X
  u.miss <- pobs(x.miss)
  ind.complete <- complete.cases(x.miss)
  r.miss <- which(ind.complete==FALSE)
  x.complete <- x.miss[ind.complete,]
  X.imp <- x.miss
  #     
  cat("It is thinking.. - Do you remember Commodor 64?","\n")
  #
  h <- 0
  alpha <- vector(length=length(r.miss)) 
  for(m in 1:length(r.miss)){
    miss  <- r.miss[m] 
    x.na    <- x.miss[miss,]
    j.cond <- which(!is.na(x.na))
    j.na     <- which(is.na(x.na))
    u.na    <- u.miss[miss,]
    set.seed(seed)
    u.        <- replicate(N, runif(length(j.na)))
    if(length(j.na)==1){u.<- t(u.)}
    cec      <- vector()
    for(i in 1:ncol(u.)){
      nan <- u.na
      nan[which(is.na(u.na))] <- u.[,i]
      cec[i] <- condEmpCop(nan,x.complete,j.cond)
    }
    if(all(cec==0)){
      print("The empirical conditional copula is flat.")
      h <- h+1
      X.imp[miss,]  <- x.miss[miss,] 
      alpha_fin <- 0
    }else{
      fin <- cbind(t(u.),cec)
      col.cec <- length(j.na)+1
      fin.ord <- fin[order(fin[,col.cec]),]
      cec.ord <- fin.ord[,col.cec]
      if(sum(cec>0.5)==0 | sum(cec<0.5)==0){         
        print("No Fu>0.5 or Fu<0.5")
        Falpha <- ncol(u.)  
        u.imp <- fin.ord[Falpha,-col.cec] 
        alpha_fin <- fin.ord[Falpha,col.cec] # max values of the computed values of cec
      }else{
        Falpha1 <- vector(length=lPsi) 
        Falpha2 <- vector(length=lPsi) 
        for(i in 1: lPsi){
          Falpha1[i] <- cec.ord[which(cec.ord>=(0.5-Psi[i]))[1]]
          Falpha2[i] <- cec.ord[which(cec.ord>=(0.5+Psi[i]))[1]]
        }
        cec_diff <- round((Falpha1-(1-Falpha2)),2)
        cec_diff <- cec_diff[which(!is.na(cec_diff))]
        lcecdiff <- length(cec_diff)
        if(sum(cec_diff)==0){
          alpha_fin <- 0.5
          u.imp <- fin.ord[which(cec.ord>=alpha_fin)[1],-col.cec] 
        }else{
          if(sum(cec_diff)>0){ 
            alpha_fin <- (0.5-Psi)[which((cec_diff==max(cec_diff)))[1]]
            u.imp <- mean(fin.ord[which(fin.ord[,col.cec]==fin.ord[which((cec.ord>=alpha_fin))[1],col.cec]),-col.cec]) 
          }else{
            if(sum(cec_diff)<0){ 
              alpha_fin <- (0.5+Psi)[which((abs(cec_diff)==max(abs(cec_diff))))[1]]
              u.imp <-  mean(fin.ord[which(fin.ord[,col.cec]==fin.ord[which((cec.ord>=alpha_fin))[1],col.cec]),-col.cec])  # quantile alpha_fin
            }
          }
        }
      }
      u.na[j.na] <- u.imp      
      if(method=="gower"){
        distances <- as.matrix(daisy(rbind(u.na, u.miss[-r.miss,]), metric=method))[-1, 1] 
      }else{
        if(method=="kendall"){
          distances <- as.matrix(as.dist(sqrt(1-cor(t(rbind(u.na, u.miss[-r.miss,])),method=method)^2)))[-1, 1] 
        }else{
          if(method=="kendall2"){
            distances <- as.matrix(as.dist(1-abs(cor(t(rbind(u.na, u.miss[-r.miss,])),method="kendall"))))[-1, 1] 
          }else{
            distances <- as.matrix(dist(rbind(u.na, u.miss[-r.miss,]), method=method))[-1, 1] 
          }
        }   
      }
      smallest_indices <- order(distances)[1:K] 
      w_miss <-   x.complete[smallest_indices, ] 
      x.na[j.na] <- colMeans(w_miss)[j.na]  
      x.imp <- x.na
      X.imp[miss,]  <- x.imp
    }
    alpha[m] <- alpha_fin
  }
    out <- new("NPCoImp")
    out@Imputed.matrix   <- X.imp;
    out@Selected.alpha  <- alpha;
    out@numFlat   <- h;
    return(out);    
}
