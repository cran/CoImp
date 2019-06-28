\name{PerfMeasure}
\alias{PerfMeasure}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Performance measures for evaluating the goodness of an imputed database}
\description{
    Set of measures useful to evaluate the goodness of the used imputation method.
}
\usage{
PerfMeasure(db.complete, db.imputed, db.missing, n.marg = 2, model = 
                        list(normalCopula(0.5, dim=n.marg), claytonCopula(10, 
                        dim=n.marg), gumbelCopula(10, dim=n.marg), frankCopula(10, 
                        dim=n.marg), tCopula(0.5, dim=n.marg,...), 
                        rotCopula(claytonCopula(10,dim=n.marg),flip=rep(TRUE,n.marg)),
                        ...), ...)
}                    
                    
 %- maybe also 'usage' for other objects documented here.
\arguments{
  \item{db.complete}{the complete data matrix.}
  \item{db.imputed}{the imputed data matrix.}
  \item{db.missing}{the data matrix with \code{NA} data.}
  \item{n.marg}{the number of variables in db.complete.}
  \item{model}{a list of copula models to be used for the imputation. See the Details section.
        This should be one of \code{normal} and \code{t} (with \code{dispstr} as in the \code{copula} package), \code{frank}, \code{clayton}, \code{gumbel}, and \code{rotated copulas}.
        As in \code{fitCopula}, itau fitting coerced tCopula to 'df.fixed=TRUE'.}
  \item{\dots}{further parameters for \code{\link{fitCopula}}.}
}
\details{
%%    \describe{
%%    \item{Usage for continuous variables:}{.\cr
%%    \code{PerfMeasure(db.complete, db.imputed, db.missing, n.marg = 2, 
%%      model = list(normalCopula(0.5, dim=n.marg),  
 %%    claytonCopula(10, dim=n.marg), gumbelCopula(10, dim=n.marg), 
%%     frankCopula(10, dim=n.marg), tCopula(0.5, dim=n.marg,...),
%%    rotCopula(claytonCopula(10,dim=n.marg),flip=rep(TRUE,n.marg)),...))}}
%%    }

    PerfMeasure computes some measures useful for evaluating the goodness of the used imputation method.
    PerfMeasure requires in input the imputed, the complete and the missing data matrix and gives in output
    five different measures of performance. See below for details
}
\value{An object of S4 class "PerfMeasure", which is a list with the following elements:
    \item{\code{MARE}}{Object of class \code{"numeric"}. The mean (on the replications performed) of the absolute relative error between the imputed and the corresponding original value.}
    \item{RB}{Object of class \code{"numeric"}. The relative bias of the estimator for the dependence parameter.}
    \item{RRMSE}{Object of class \code{"numeric"}. The relative root mean squared error of the estimator for the dependence parameter.}
    %\item{MS}{Object of class \code{"numeric"}. The average difference (on the replications performed) between the dependence of
    %the imputed variables and that of the complete ones computed in terms of \code{\link{Srho.ts}}.}
    \item{TID}{Object of class \code{"vector"}. Upper and lower tail dependence indexes for bivariate copulas. Original function is in \code{\link{tailIndex}}.}
}

\references{
Di Lascio, F.M.L. Giannerini, S. and Reale A. (201x) "A multivariate technique based on conditional copula specification for the imputation of complex dependent data". Working paper.

Di Lascio, F.M.L., Giannerini, S. and Reale, A. (2015) "Exploring Copulas for the Imputation of Complex Dependent Data". Statistical Methods & Applications, 24(1), p. 159-175. DOI 10.1007/s10260-014-0287-2.

Di Lascio, F.M.L., Giannerini, S. and Reale, A. (2014) "Imputation of complex dependent data by conditional copulas: analytic versus semiparametric approach", Book of proceedings of the 21st International Conference on Computational Statistics (COMPSTAT 2014), p. 491-497. ISBN 9782839913478.

Bianchi, G. Di Lascio, F.M.L. Giannerini, S. Manzari, A. Reale, A. and Ruocco, G. (2009) "Exploring copulas for the imputation of missing nonlinearly dependent data". Proceedings of the VII Meeting Classification and Data Analysis Group of the Italian Statistical Society (Cladag), Editors: Salvatore Ingrassia and Roberto Rocci, Cleup, p. 429-432. ISBN: 978-88-6129-406-6.
}

\author{
Francesca Marta Lilja Di Lascio <marta.dilascio@unibz.it>,

Simone Giannerini <simone.giannerini@unibo.it>}

\examples{
\dontrun{
# generate data from a 4-variate Gumbel copula with different margins

set.seed(11)
n.marg <- 4
theta  <- 5
copula <- frankCopula(theta, dim = n.marg)
mymvdc <- mvdc(copula, c("norm", "gamma", "beta","gamma"), list(list(mean=7, sd=2),
 list(shape=3, rate=2), list(shape1=4, shape2=1), list(shape=4, rate=3)))
n      <- 20
x.samp <- rMvdc(n, mymvdc)

# randomly introduce univariate and multivariate missing

perc.mis    <- 0.3
set.seed(11)
miss.row    <- sample(1:n, perc.mis*n, replace=TRUE)
miss.col    <- sample(1:n.marg, perc.mis*n, replace=TRUE)
miss        <- cbind(miss.row,miss.col)
x.samp.miss <- replace(x.samp,miss,NA)

# impute missing values

imp <- CoImp(x.samp.miss, n.marg=n.marg, smoothing=rep(0.6,n.marg), plot=TRUE,
            type.data="continuous");
imp

# apply PerfMeasure to the imputed data set

pm <- PerfMeasure(db.complete=x.samp, db.missing=x.samp.miss,
            db.imputed=imp@"Imputed.data.matrix", n.marg=4)
pm
str(pm)
}

}

% Add one or more standard keywords, see file 'KEYWORDS' in the % R documentation directory.
\keyword{imputation}
\keyword{copula}
\keyword{multivariate}
