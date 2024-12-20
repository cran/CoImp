#############################################################
# CONDITIONAL EMPIRICAL COPULA FUNCTION     
#                                                           
# INPUT: evaluation points of the empirical copula, complete sample data matrix, 
#        vector of indices of conditioning variables, kind of empirical copula
#
### NPCoImp
### A NON PARAMETRIC COPULA BASED IMPUTATION METHOD
##
##  The authors of this software are
##  F. Marta L. Di Lascio, and
##  Aurora Gatto, Copyright (c) 2024
##
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


condEmpCop <- function(u, X, j.cond, smoothing){
  # u evaluation points of the empirical copula
  # X complete sample data (if pobs(X) is given in input, nothing changes since C.n
  #     calls F.n that applies pobs and apply twice pobs do not determing anythings)
  # j.cond vector of indices of conditioning variables
  if(any(u<0, 1<u)) stop("'u' must be in [0,1]")
  if(!is.matrix(X)) stop("'X' should be a data matrix")
  if(any(is.na(X))){
    print("Complete sample data are used")
    X <- X[complete.cases(X),]
  }
  n.marg <- ncol(X)
  if(length(u)!=n.marg) stop("'X' and u have incompatible dimensions")
  if(length(j.cond)>=n.marg) stop("'X' and 'j.cond' have incompatible dimensions")
  if(length(j.cond)==1){
    res <- C.n(u,X=X, smoothing=smoothing)/u[j.cond]
  }else{
    num <- C.n(u,X=X, smoothing=smoothing)
    den <- C.n(u[j.cond],X=X[,j.cond], smoothing=smoothing)
    if(num==0 & den==0){
      rum <- runif(1, 0, 0.000001)
      res <- (num+rum)/(den+rum)
    }else{
      if(den==0){
        rum <- runif(1, 0, 0.000001)
        res <- num/(den+rum)
      }else{
        res <- num/den
      }
    }
  }
  return(res)
}
