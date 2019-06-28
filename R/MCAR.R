
### CoImp
### A COPULA BASED IMPUTATION METHOD
##
##  The authors of this software are
##  Francesca Marta Lilja Di Lascio, and
##  Simone Giannerini, Copyright (c) 2013

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

setClass("MCAR",
         representation(db.missing = "matrix"
                        ),
         prototype = list(db.missing = matrix()
                        )
         )


## ***************************************************************************************************

MCAR <- function(db.complete, perc.miss = 0.3, setseed = 13, mcols = NULL, ...){
                # introduce MCAR in a dataset
                # a subset of cols where introduce MCAR values can be specified in "mcols"
                if(!is.matrix(db.complete))
                    stop("the data matrix in entry should be a matrix")
                if(perc.miss<=0)
                    stop("the missing percentage should be positive")
                if(sum(is.na(db.complete))!=0)
                    stop("the data matrix in entry should be complete")
                #
                if(is.null(mcols)==TRUE){
                    db.missing <- db.complete
                }else{
                    db.missing <- db.complete[,mcols]
                }
                n.marg <- ncol(db.missing)
                n      <- nrow(db.missing)
                if((n*perc.miss)<1)
                    stop("there are no missing to be introduced; perc.miss should be increased")
                set.seed(setseed)
                idMiss <- sample(1:n, n*perc.miss)                                    # sample missing cases
                nMiss  <- length(idMiss)
                mMax   <- n.marg-1                                                    # maximum num of missing variables for each record
                set.seed(setseed)
                howmanyMiss <- sapply(idMiss, function(x) sample(1:mMax, 1))          # num of missing to be introduced in each selected id (idMiss)
                misscols    <- lapply(howmanyMiss, function(x) sample(1:n.marg, x))   # variables (in num=howmanyMiss) to be missed for each id
                for(i in 1:nMiss){
                    for (j in misscols[[i]]){
                        db.missing[idMiss[i],j] <- NA
                    }
                }
                db.missing.fin <- db.complete
                if(is.null(mcols)==TRUE){
                    db.missing.fin <- db.missing
                }else{
                    db.missing.fin[,mcols] <- db.missing
                }
                #
                out       <- new("MCAR")
                out@db.missing  <- db.missing.fin;
                return(out);
}
