#############################################################
# ALGORITMO DI IMPUTAZIONE DATI MAR VIA FUNZIONE COPULA     #
#                                                           #
# INPUT: data matrix, number of replications HITorMISS, par #
# for barplot and copula models to be used                  #
#                                                           #
### CoImp
### A COPULA BASED IMPUTATION METHOD
##
##  The authors of this software are
##  F. Marta L. Di Lascio, and
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

setClass("CoImp",
         representation(Missing.data.matrix    = "matrix"
                        ,Perc.miss             = "matrix"
                        ,Estimated.Model       = "list"
                        ,Estimation.Method     = "character"
                        ,Index.matrix.NA       = "matrix"
                        ,Smooth.param          = "vector"
                        ,Imputed.data.matrix   = "matrix"
                        ,Estimated.Model.Imp   = "list"
                        ,Estimation.Method.Imp = "character"
                        ),
         prototype = list(Missing.data.matrix  = matrix(0,0,0)
                        ,Perc.miss             = matrix(0,0,0)
                        ,Estimated.Model       = list()
                        ,Estimation.Method     = character()
                        ,Index.matrix.NA       = matrix(0,0,0)
                        ,Smooth.param          = vector()
                        ,Imputed.data.matrix   = matrix(0,0,0)
                        ,Estimated.Model.Imp   = list()
                        ,Estimation.Method.Imp = character()
                        )
         )

 ## ***************************************************************************************************

setMethod(
    f="plot",
    signature(x = "CoImp",y = "missing"),
    definition=function(x, y, plot.legend = TRUE, args.legend = list(y = 110, cex = 0.8), plot.name = NULL, ...){
        par(mfrow=c(1,1));
        bar.plot(tab = x@Perc.miss, legend.plot = plot.legend, args.legend = args.legend, ...)
        X <- x@Imputed.data.matrix
        smoothing <- x@Smooth.param
        n.marg <- ncol(X)
        fit0_nn <- list()
        for(i in 1:n.marg){
            fit0_nn[[i]] <- locfit::locfit( ~ locfit::lp(X[,i], nn=smoothing[i], deg=1))      # alpha=0.9 (deprecated) or nn=0.5
        }
        fit0 <- fit0_nn
        alpha <- smoothing
        f.x  <- list()
        for(j in 1:n.marg){
            f.x[[j]] <- function(x){predict(fit0[[j]], newdata=x)}
            knorm <- try(integrate(f.x[[j]], lower=0, upper=Inf),silent=TRUE)
            if(knorm[[1]]==0 | inherits(knorm, "try-error")==TRUE){knorm[[1]] <- 0.001}
            f.x[[j]] <- function(x){predict(fit0[[j]], newdata=x)/as.numeric(knorm[[1]])}
        }
        #
        if(n.marg<=3){
            dev.new();
            par(mfrow=c(1,ceiling(n.marg)))
        }else{
            if(n.marg >3 && n.marg<=6){
                dev.new();
                par(mfrow=c(2,ceiling(n.marg/2)))
            }else{
                dev.new();
                par(mfrow=c(3,ceiling(n.marg/3)))
            }
        }
        par(mai=c(0.5,0.5,0.3,0.5))
        for(i in 1:n.marg){
            j       <- i
            minimo  <- min(X[complete.cases(X[,i]),i])
            massimo <- max(X[complete.cases(X[,i]),i])
            opt1     <- f.x[[i]](optimize(f.x[[i]], c(minimo,massimo), maximum=TRUE)$maximum)
            his  <- hist(X[,i], plot=FALSE)
            opt2 <- max(his$density)
            opt  <- max(opt1,opt2)
            plot(his, freq=FALSE, ylim=c(0,opt), main="", xlab="", ylab="")
            plot(f.x[[i]], lwd=2, xlim=c(minimo-1, massimo+1),ylim=c(0,opt), col="blue", add=TRUE)
            if(is.null(plot.name)){
                title(main = paste("Variable X", i,sep=""))
            }else{
                title(main = plot.name[i])
            }
        }
    }
)

## ***************************************************************************************************

setMethod(
    f="show",
    signature="CoImp",
    definition=function(object){
        out <- object
        cat (" Main output of the function CoImp \n")
        cat (" -------------------------------------------------------------------------- \n")
        cat (" Percentage of missing and available data : \n")
        print(out@"Perc.miss")
        cat (" -------------------------------------------------------------------------- \n")
        cat (" Imputed data matrix : \n")
        print(out@"Imputed.data.matrix")
        cat (" -------------------------------------------------------------------------- \n")
    }
)

## ***************************************************************************************************

CoImp <- function(X, n.marg = ncol(X), x.up = NULL, x.lo = NULL, q.up = rep(0.85, n.marg),
            q.lo = rep(0.15, n.marg), type.data = "continuous", smoothing =
            rep(0.5, n.marg), plot = TRUE, model = list(normalCopula(0.5,
            dim=n.marg),  claytonCopula(10, dim=n.marg), gumbelCopula(10,
            dim=n.marg), frankCopula(10, dim=n.marg), tCopula(0.5,
            dim=n.marg,...), rotCopula(claytonCopula(10,dim=n.marg),
            flip=rep(TRUE,n.marg)),...), start. = NULL, ...){
    #
    if(!is.matrix(X))
        stop("Tha data should be a matrix")
    if(n.marg<=1)
        stop("The data matrix should contain at least two variables")
    if(nrow(X)<1)
        stop("The data matrix should contain at least one observation")
    if(any(is.na(X))==FALSE)
        stop("The data matrix should contain at least one missing value")
    if(sum(complete.cases(X))==0)
        stop("The data matrix should contain at least one complete record")
    if(type.data != "discrete" & type.data != "continuous")
        stop("The data must be either continuous or discrete")
    if(type.data == "discrete")
        warning("The variables are treated as continuous and rounded off.")
    #
    if(!is.null(x.up) & !is.null(q.up))
        stop("Specify either x.up xor q.up")
    if(!is.null(x.lo) & !is.null(q.lo))
        stop("Specify either x.lo xor q.lo")
    if(!is.null(q.up) & (any(q.up<=0) || any(q.up>=1)))
        stop("q.up must lie in (0, 1)")
    if(!is.null(q.lo) & (any(q.lo<=0) || any(q.lo>=1)))
        stop("q.lo must lie in (0, 1)")
    if(is.null(x.up) & is.null(q.up))
        stop("Specify either x.up xor q.up")
    if(is.null(x.lo) & is.null(q.lo))
        stop("Specify either x.lo xor q.lo")
    #
    fit0_nn <- list()
    for(i in 1:n.marg){
        fit0_nn[[i]] <- locfit::locfit( ~ locfit::lp(X[,i], nn=smoothing[i], deg=1))      # alpha=0.9 (deprecated) or nn=0.5
    }
    fit0 <- fit0_nn
    alpha <- smoothing
    names(alpha) <- c(paste("x",1:n.marg,sep=""))
    f.x  <- list()
    F.x  <- list()
    for(j in 1:n.marg){
        f.x[[j]] <- function(x){predict(fit0[[j]], newdata=x)}
        knorm <- try(integrate(f.x[[j]], lower=0, upper=Inf),silent=TRUE)
        if(knorm[[1]]==0 | inherits(knorm, "try-error")==TRUE){knorm[[1]] <- 0.001}
        f.x[[j]] <- function(x){predict(fit0[[j]], newdata=x)/as.numeric(knorm[[1]])}
        F.x[[j]] <- function(x){ifelse(!is.na(x),integrate(f.x[[j]], lower=0, upper=x)$value,return(NA))}
    }
    ind.miss     <- which(is.na(X), arr.ind = TRUE)
    if(nrow(ind.miss)>1){ind.miss <-  ind.miss[order(ind.miss[,1]),]}
    ind.cols.na  <- split(ind.miss[,2],ind.miss[,1])
    num.rows.na  <- length(ind.cols.na)
    #
    n.mod        <- length(model)
    loglik       <- double(length=n.mod)
    metodo.fin   <- double(length=n.mod)
    metodo.optim.fin <- double(length=n.mod)
    for(i in 1:n.mod){
        metodo   <- "ml"
        metodo.c <- "BFGS"
        udat.na  <- fit.margin(dataset=t(X), param=list(dimc=n.marg))
        udat     <- udat.na[complete.cases(udat.na),]
        fitc     <- try(fitCopula(data=udat, copula=model[[i]], start=start., method=metodo, optim.method=metodo.c), silent=TRUE)
        if(inherits(fitc, "try-error")==TRUE){
            metodo <- c("mpl", "itau", "irho")
            repeat{
                if(length(metodo)==0 || inherits(fitc, "try-error")==FALSE)
                    break
                fitc <- try(fitCopula(data = udat, copula = model[[i]],
                        start = start., method = metodo[[1]]), silent = TRUE)
                metodo <- setdiff(metodo, metodo[[1]])
            }
        }
        if(inherits(fitc, "try-error")==TRUE){
            metodo.c <- c("Nelder-Mead", "CG", "L-BFGS-B", "SANN")
            repeat{
                if(length(metodo.c)==0 || inherits(fitc, "try-error")==FALSE)
                    break
                metodo <- c("ml","mpl","itau","irho")
                repeat{
                    if(length(metodo)==0 || inherits(fitc, "try-error")==FALSE)
                        break
                    fitc <- try(fitCopula(data = udat, copula = model[[i]],
                            start = start., method = metodo[[1]], optim.method=metodo.c[[1]]), silent = TRUE)
                    metodo <- setdiff(metodo, metodo[[1]])
                }
                metodo.c <- setdiff(metodo.c, metodo.c[[1]])
            }
        }
        if (inherits(fitc, "try-error") || is.nan(suppressWarnings(loglikCopula(param=fitc@estimate, u=udat, copula=model[[i]])))) {
            loglik[i] <- -10000
        }else{
            loglik[i] <- suppressWarnings(loglikCopula(param=fitc@estimate, u=udat, copula=model[[i]]))
        }
        if(length(metodo)==0){
            metodo.fin[i] <- 0
        }else{
            metodo.fin[i] <- metodo[[1]]
        }
        if(length(metodo.c)==0){
            metodo.optim.fin[i] <- 0
        }else{
            metodo.optim.fin[i] <- metodo.c[[1]]
        }
    }
    best <- which(loglik==max(loglik[which(!is.na(loglik))]))[[1]]
    mod.fin.base <- model[[best]]
    metodo.fin.base <- metodo.fin[[best]]
    metodo.optim.fin.base <- metodo.optim.fin[[best]]
    if(metodo.fin.base=="ml"){
        udat.na <- fit.margin(dataset=t(X), param=list(dimc=n.marg))
        udat <- udat.na[complete.cases(udat.na),]
    }else{
        udat.na <- fit.margin2(dataset=t(X), param=list(dimc=n.marg))
        udat <- udat.na[complete.cases(udat.na),]
    }
    mod.fin <- try(fitCopula(data=udat, copula=mod.fin.base, start=start., method=metodo.fin.base, optim.method=metodo.optim.fin.base),silent=TRUE)
    if (inherits(mod.fin, "try-error")) {
        stop("Imputation failed")
    }else{
        if(inherits(mod.fin.base, what="rotExplicitCopula")  | inherits(mod.fin.base, what="rotCopula")){
            mod.fin.base@copula@parameters <- mod.fin@estimate
        }else{
            mod.fin.base@parameters <- mod.fin@estimate
        }
    }
    #
    dati.fin <- X
    x.min  <- apply(X, 2, min, na.rm=TRUE)
    x.max  <- apply(X, 2, max, na.rm=TRUE)
    if(is.null(x.up)){
        x.up <- diag(apply(X, 2, quantile, p=q.up, na.rm=TRUE))
    }
    if(is.null(x.lo)){
        x.lo <- diag(apply(X, 2, quantile, p=q.lo, na.rm=TRUE))
    }
    for(i in 1:num.rows.na){
        cat("\r Number of imputed rows: ", i, "\n");
        cols.na <- ind.cols.na[[i]]
        cols.no.na <- seq(1:n.marg)[-cols.na]
        rows.na <- as.numeric(names(ind.cols.na))[i]
        y <- X[rows.na,]
        lcn <- length(cols.na)
        #
        meanX   <- colMeans(X,na.rm=TRUE)
        medianX <- apply(X=X, FUN=median, MARGIN=2, na.rm=TRUE)
        #
        zz <- y
        fcond <- function(x){
            zz[cols.na] <- x
            fcond.mod(x, y=zz, ind=cols.na, model=mod.fin.base, distr=F.x, dens=f.x)
        }
        fcond2 <- Vectorize(fcond)
        #
        if(lcn==n.marg)
            stop("It cannot impute a full NA record")
        #
        a  <- x.min[cols.na]
        b  <- x.max[cols.na]
        aa <- x.lo[cols.na]
        bb <- x.up[cols.na]

        if(lcn==1){
            maxf <- suppressWarnings(optimize(fcond2, interval=c(a, b), maximum=TRUE)$max)
            if(maxf>=aa & maxf<=bb){
                umissM <- hitormiss(FUN=fcond, p=lcn, h=fcond(maxf), a=aa, b=bb)
            }else{
                umissM <- hitormiss(FUN=fcond, p=lcn, h=fcond(maxf), a=a, b=b)
            }
        }else{
            maxf <- suppressWarnings(try(optim(par=meanX[cols.na], fn=fcond, lower=a, upper=b, control = list(fnscale=-1))$par,silent=TRUE))
            if(inherits(maxf, "try-error")==TRUE){
                maxf <- suppressWarnings(optim(par=medianX[cols.na], fn=fcond, lower=a, upper=b, control = list(fnscale=-1))$par)
            }
            if(all(maxf>=aa) & all(maxf<=bb)){
                umissM <- suppressWarnings(hitormiss(FUN=fcond,p=lcn,h=fcond(maxf),a=aa,b=bb))
            }else{
                umissM <- suppressWarnings(hitormiss(FUN=fcond,p=lcn,h=fcond(maxf),a=a,b=b))
            }
        }
        dati.fin[rows.na,cols.na] <- umissM
    }
    #
    loglik.imp <- double(length=n.mod)
    metodo.fin.imp <- double(length=n.mod)
    metodo.optim.fin.imp <- double(length=n.mod)
    if(type.data=="discrete"){
        dati.fin <- round(dati.fin,0)
    }
    for(i in 1:n.mod){
        metodo <- "ml"
        metodo.c <- "BFGS"
        udat.na.imp <- fit.margin(dataset=t(dati.fin), param=list(dimc=n.marg))
        udat.imp <- udat.na.imp[complete.cases(udat.na.imp),]
        fitc.imp <- try(fitCopula(data=udat.imp, copula=model[[i]], start=start., method=metodo, optim.method=metodo.c), silent=TRUE)
        if(inherits(fitc.imp, "try-error")==TRUE){
            metodo <- c("mpl","itau","irho")
            repeat{
                if(length(metodo)==0 || inherits(fitc.imp, "try-error")==FALSE)
                    break
                fitc.imp <- try(fitCopula(data = udat.imp, copula = model[[i]],
                        start = start., method = metodo[[1]]), silent = TRUE)
                metodo <- setdiff(metodo, metodo[[1]])
            }
        }
        if(inherits(fitc.imp, "try-error")==TRUE){
            metodo.c <- c("Nelder-Mead", "CG", "L-BFGS-B", "SANN")
            repeat{
                if(length(metodo.c)==0 || inherits(fitc.imp, "try-error")==FALSE)
                    break
                metodo <- c("ml","mpl","itau","irho")
                repeat{
                    if(length(metodo)==0 || inherits(fitc.imp, "try-error")==FALSE)
                        break
                    fitc.imp <- try(fitCopula(data = udat.imp, copula = model[[i]],
                            start = start., method = metodo[[1]], optim.method=metodo.c[[1]]), silent = TRUE)
                    metodo <- setdiff(metodo, metodo[[1]])
                }
                metodo.c <- setdiff(metodo.c, metodo.c[[1]])
            }
        }
        if (inherits(fitc.imp, "try-error") || is.nan(suppressWarnings(loglikCopula(param=fitc.imp@estimate, u=udat.imp, copula=model[[i]])))) {
            loglik.imp[i] <- -10000
        }else{
            loglik.imp[i] <- suppressWarnings(loglikCopula(param=fitc.imp@estimate, u=udat.imp, copula=model[[i]]))
        }
        if(length(metodo)==0){
            metodo.fin.imp[i] <- 0
        }else{
            metodo.fin.imp[i] <- metodo[[1]]
        }
        if(length(metodo.c)==0){
            metodo.optim.fin.imp[i] <- 0
        }else{
            metodo.optim.fin.imp[i] <- metodo.c[[1]]
        }
    }
    best.imp <- which(loglik.imp==max(loglik.imp[which(!is.na(loglik.imp))]))[[1]]
    mod.fin.base.imp <- model[[best.imp]]
    metodo.fin.base.imp <- metodo.fin.imp[[best.imp]]
    metodo.optim.fin.base.imp <- metodo.optim.fin.imp[[best.imp]]
    if(metodo.fin.base.imp=="ml"){
        udat.na.imp <- fit.margin(dataset=t(dati.fin), param=list(dimc=n.marg))
        udat.imp <- udat.na.imp[complete.cases(udat.na.imp),]
    }else{
        udat.na.imp <- fit.margin2(dataset=t(dati.fin), param=list(dimc=n.marg))
        udat.imp <- udat.na.imp[complete.cases(udat.na.imp),]
    }
    mod.fin.imp <- try(fitCopula(data=udat.imp, copula=mod.fin.base.imp, start=start., method=metodo.fin.base.imp, optim.method=metodo.optim.fin.base.imp),silent=TRUE)
    if (inherits(mod.fin.imp, "try-error")) {
        stop("Imputation failed")
    }else{
        if(inherits(mod.fin.base.imp,what="rotExplicitCopula") | inherits(mod.fin.base.imp,what="rotCopula")){
            mod.fin.base.imp@copula@parameters <- mod.fin.imp@estimate
        }else{
            mod.fin.base.imp@parameters <- mod.fin.imp@estimate
        }
    }
    #
    perc.miss <- round(colMeans(is.na(X)*100),2)
    perc.data <- round(rbind(100-perc.miss, perc.miss),3)
    rownames(perc.data) <- c("Data","Missing")
    ifelse(is.null(colnames(X)), colnames(perc.data) <- paste("X",c(1:n.marg),sep=""), colnames(perc.data) <- colnames(X))
    #
    if(is.null(colnames(X))) colnames(X) <- paste("X",c(1:n.marg),sep="")
    ifelse(is.null(colnames(X)), colnames(dati.fin) <- paste("X",c(1:n.marg),sep=""), colnames(dati.fin) <- colnames(X))
    #
    if(inherits(mod.fin.base, what="rotExplicitCopula")  | inherits(mod.fin.base, what="rotCopula")){
         mod.fin.base <- mod.fin.base@copula
    }
    if(inherits(mod.fin.base.imp, what="rotExplicitCopula")  | inherits(mod.fin.base.imp, what="rotCopula")){
        mod.fin.base.imp <- mod.fin.base.imp@copula
    }
    mod.pre  <- list(model = mod.fin.base@class, dimension = mod.fin.base@dimension,  parameter = mod.fin.base@parameters, number = best)
    mod.post <- list(model = mod.fin.base.imp@class, dimension = mod.fin.base.imp@dimension, parameter = mod.fin.base.imp@parameters, number = best.imp)
    #
    out <- new("CoImp")
    out@Missing.data.matrix   <- X;
    out@Perc.miss             <- perc.data;
    out@Estimated.Model       <- mod.pre;
    out@Estimation.Method     <- metodo.fin.base;
    out@Index.matrix.NA       <- ind.miss;
    out@Smooth.param          <- alpha;
    out@Imputed.data.matrix   <- dati.fin;
    out@Estimated.Model.Imp   <- mod.post;
    out@Estimation.Method.Imp <- metodo.fin.base.imp;
    if(plot==TRUE) plot(out)
    return(out);
}
