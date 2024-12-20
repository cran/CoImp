#####################################################################
# ALGORITHM OF THE NON PARAMETRIC COPULA BASED IMPUTATION METHOD    
#                                                           
# INPUT: data matrix, vector of probabilities, K records to be selected for the
#        imputation, distance measure, kind of empirical copula function 
#
#
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
                        ,Selected.quantile.alpha  = "vector"
                        ,numFlat      = "numeric"
                        ),
         prototype = list(Imputed.matrix = matrix(0,0,0)
                        ,Selected.quantile.alpha = vector()
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
        cat (" Best quantile for the imputation and the corresponding alpha: \n")
        print(out@"Selected.quantile.alpha")
        cat (" -------------------------------------------------------------------------- \n")
        cat ("Number of flat conditional empirical copulas : \n")
        print(out@"numFlat")
        cat (" -------------------------------------------------------------------------- \n")
    }
)

## ***************************************************************************************************

NPCoImp <- function(X, Psi=seq(0.05,0.45,by=0.05), smoothing="beta", K=7, method="gower"){
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
  if(length(Psi)<9){cat("Note that the imputation through the conditional empirical beta copula could be poorly accurate. To improve it use a larger Psi.","\n")}
  #
  lPsi <- length(Psi)
  x.miss <- X
  u.miss <- pobs(x.miss)
  ind.complete <- complete.cases(x.miss)
  r.miss <- which(ind.complete==FALSE)
  x.complete <- x.miss[ind.complete,]
  X.imp <- x.miss
  #     
  cec.alpha <- vector()
  qfin      <- vector()
  #
  h <- 0
  alpha <- vector(length=length(r.miss)) 
  for(m in 1:length(r.miss)){
      miss  <- r.miss[m] 
      x.na    <- x.miss[miss,]
      j.cond <- which(!is.na(x.na))
      j.na     <- which(is.na(x.na))
      u.na    <- u.miss[miss,]    
      cecL   <- vector()
      cecU   <- vector()
      for(j in 1:lPsi){       
          nanL <- u.na
          nanL[j.na] <- rep(0.5-Psi[j], length(j.na))
          nanU <- u.na
          nanU[j.na] <- rep(0.5+Psi[j], length(j.na))
          cecL[j] <- condEmpCop(nanL,x.complete,j.cond,smoothing=smoothing)
          cecU[j] <- condEmpCop(nanU,x.complete,j.cond,smoothing=smoothing) 
      }
      cecASY <- c(cecL,cecU) 
      if(all(cecASY==0)){ 
          print("It is not possible to assess the (a)symmetry since the empirical conditional copula is equal to zero in the considered values. Hence, it is not possible to impute.")
          h <- h+1
          X.imp[miss,]  <- x.miss[miss,] 
          q_fin <- cec_alpha <- u.imp <- NA
      }else{
          if(sum(cecASY>0.5)==0){      
              print("No C(u)>0.5") 
              q_fin <- (0.5-Psi)[which((cecASY==max(cecASY)))[1]]
              cec_alpha <- max(cecASY)
              u.imp <- q_fin
          }else{
              if(sum(cecASY<0.5)==0){  
                  print("No C(u)<0.5") 
                  q_fin <- (0.5+Psi)[which((cecASY==min(cecASY)))[1]]
                  cec_alpha <- min(cecASY)
                  u.imp <- q_fin
              }else{
                  cec_diff <- round(cecL-(1-cecU),2) 
                  cec_diff <- cec_diff[which(!is.na(cec_diff))]
                  lcecdiff <- length(cec_diff)
                  if(sum(cec_diff)==0){
                      q_fin <- 0.5 
                      cec_alpha <- 0.5
                      u.imp <- q_fin
                  }else{
                      if(sum(cec_diff)>0){ 
                          q_fin <- (0.5-Psi)[which((cec_diff==max(cec_diff)))[1]] 
                          nan <- u.na
                          nan[which(is.na(u.na))] <- q_fin 
                          cec_alpha <- condEmpCop(nan,x.complete,j.cond,smoothing=smoothing) 
                          u.imp <- q_fin
                      }else{
                          if(sum(cec_diff)<0){ 
                              q_fin <- (0.5+Psi)[which((abs(cec_diff)==max(abs(cec_diff))))[1]] 
                              nan <- u.na
                              nan[which(is.na(u.na))] <- q_fin 
                              cec_alpha <- condEmpCop(nan,x.complete,j.cond,smoothing=smoothing)
                              u.imp <- q_fin            
                          }
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
          cec.alpha[m] <- cec_alpha
          qfin[m] <- q_fin            
      }
  }
  out <- new("NPCoImp")
  out@Imputed.matrix <- X.imp;
  out@Selected.quantile.alpha <- cbind(qfin, cec.alpha);
  out@numFlat <- h;
  return(out);
}
