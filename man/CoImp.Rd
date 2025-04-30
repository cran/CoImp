\name{CoImp}
\alias{CoImp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Copula-Based Imputation Method}
\description{
    Imputation method based on conditional copula functions.
}
\usage{
CoImp(X, n.marg = ncol(X), x.up = NULL, x.lo = NULL, q.up = rep(0.85, n.marg), 
            q.lo = rep(0.15, n.marg), type.data = "continuous", smoothing = 
            rep(0.5, n.marg), plot = TRUE, model = list(normalCopula(0.5, 
            dim=n.marg),  claytonCopula(10, dim=n.marg), gumbelCopula(10, 
            dim=n.marg), frankCopula(10, dim=n.marg), tCopula(0.5, 
            dim=n.marg,...), rotCopula(claytonCopula(10,dim=n.marg), 
            flip=rep(TRUE,n.marg)),...), start. = NULL, ...)
} 
 
 %- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a data matrix with missing values. Missing values should be denoted
        with \code{NA}.}
  \item{n.marg}{the number of variables in X.}
  \item{x.up}{a numeric vector of length n.marg with the upper value of each margin used in the Hit or Miss method. Specify either x.up xor q.up.}
  \item{x.lo}{a numeric vector of length n.marg with the lower value of each margin used in the Hit or Miss method. Specify either x.lo xor q.lo.}
  \item{q.up}{a numeric vector of length n.marg with the probability of the quantile used to define x.up for each margin. Specify either x.up xor q.up.}
  \item{q.lo}{a numeric vector of length n.marg with the probability of the quantile used to define x.lo for each margin. Specify either x.lo xor q.lo.}   
  \item{type.data}{the nature of the variables in X: \code{discrete} or \code{continuous}.}
  \item{smoothing}{values for the nearest neighbour component of the smoothing parameter of the \code{\link[locfit]{lp}} function.}
  \item{plot}{logical: if \code{TRUE} plots the estimated marginal densities and a bar plot of the percentages of missing and available data for each margin.}  
  \item{model}{a list of copula models to be used for the imputation, see the Details section.
        This should be one of \code{normal} and \code{t} (with \code{dispstr} as in the \code{copula} package), \code{frank}, \code{clayton}, \code{gumbel}, and \code{rotated copulas}.
        As in \code{\link[copula]{fitCopula}}, itau fitting coerced tCopula to 'df.fixed=TRUE'.}
  \item{start.}{a numeric vector of starting values for the parameter optimization via \code{\link{optim}}.}
  \item{\dots}{further parameters for \code{\link[copula]{fitCopula}}, \code{\link[locfit]{lp}} and further graphical arguments.}
}
\details{
    CoImp is an imputation method based on conditional copula functions that allows to impute missing observations according to the multivariate dependence structure of the generating process without any assumptions on the margins. This method can be used independently from the dimension and the kind (monotone or non monotone) of the missing patterns.

    Brief description of the approach:

   \enumerate{
    \item estimate both the margins and the copula model on available data by means of the semi-parametric sequential two-step inference for margins;
    \item derive conditional density functions of the missing variables given non-missing ones through the corresponding conditional copulas obtained by using the Bayes' rule;
    \item impute missing values by drawing observations from the conditional density functions derived at the previous step. The Monte Carlo method used is the Hit or Miss.
}
    The estimation approach for the copula fit is semiparametric: a range of nonparametric margins and parametric copula models can be selected by the user.
}
\value{An object of S4 class "CoImp", which is a list with the following elements:
     \item{Missing.data.matrix}{the original missing data matrix to be imputed.}
     \item{Perc.miss}{the matrix of the percentage of missing and available data.}
     \item{Estimated.Model}{the estimated copula model on the available data.}
     \item{Estimation.Method}{the estimation method used for the copula \code{Estimated.Model}.}
     \item{Index.matrix.NA}{matrix indices of the missing data.}
     \item{Smooth.param}{the smoothing parameter alpha selected on the basis of the AIC.}
     \item{Imputed.data.matrix}{the imputed data matrix.}
     \item{Estimated.Model.Imp}{the estimated copula model on the imputed data matrix.}
     \item{Estimation.Method.Imp}{the estimation method used for the copula \code{Estimated.Model.Imp}.}
}

\references{
Di Lascio, F.M.L., Giannerini, S. and Reale, A. (2015) "Exploring Copulas for the Imputation of Complex Dependent Data". Statistical Methods & Applications, 24(1), p. 159-175. DOI 10.1007/s10260-014-0287-2.

Di Lascio, F.M.L., Giannerini, S. and Reale, A. (2014) "Imputation of complex dependent data by conditional copulas: analytic versus semiparametric approach", Book of proceedings of the 21st International Conference on Computational Statistics (COMPSTAT 2014), p. 491-497. ISBN 9782839913478.

Bianchi, G. Di Lascio, F.M.L. Giannerini, S. Manzari, A. Reale, A. and Ruocco, G. (2009) "Exploring copulas for the imputation of missing nonlinearly dependent data". Proceedings of the VII Meeting Classification and Data Analysis Group of the Italian Statistical Society (Cladag), Editors: Salvatore Ingrassia and Roberto Rocci, Cleup, p. 429-432. ISBN: 978-88-6129-406-6.
}

\author{
F. Marta L. Di Lascio <marta.dilascio@unibz.it>,
Simone Giannerini <simone.giannerini@unibo.it>
}
\seealso{
    \code{\link[CoImp]{NPCoImp}},
    \code{\link[CoImp]{MCAR}},
    \code{\link[CoImp]{MAR}},
    \code{\link[copula]{fitCopula}},
    \code{\link[locfit]{lp}}.
}
\examples{

## generate data from a 4-variate Frank copula with different margins

set.seed(21)
n.marg <- 4
theta  <- 5
copula <- frankCopula(theta, dim = n.marg)
mymvdc <- mvdc(copula, c("norm", "gamma", "beta","gamma"), list(list(mean=7, sd=2),
list(shape=3, rate=2), list(shape1=4, shape2=1), list(shape=4, rate=3)))
n      <- 20
x.samp <- copula::rMvdc(n, mymvdc)

# randomly introduce univariate and multivariate missing

perc.mis    <- 0.3
set.seed(11)
miss.row    <- sample(1:n, perc.mis*n, replace=TRUE)
miss.col    <- sample(1:n.marg, perc.mis*n, replace=TRUE)
miss        <- cbind(miss.row,miss.col)
x.samp.miss <- replace(x.samp,miss,NA)

# impute missing values

imp <- CoImp(x.samp.miss, n.marg=n.marg, smoothing = rep(0.6,n.marg), plot=TRUE,
       type.data="continuous", model=list(normalCopula(0.5, dim=n.marg),
       frankCopula(10, dim=n.marg), gumbelCopula(10, dim=n.marg)));

# methods show and plot

show(imp)
plot(imp)

\dontrun{
## generate data from a 3-variate Clayton copula and introduce missing by
## using the MCAR function and try to impute through a rotated copula

set.seed(11)
n.marg <- 3
theta  <- 5
copula <- claytonCopula(theta, dim = n.marg)
mymvdc <- mvdc(copula, c("beta", "beta", "beta"), list(list(shape1=4, shape2=1),
            list(shape1=.5, shape2=.5), list(shape1=2, shape2=3)))
n      <- 50
x.samp <- copula::rMvdc(n, mymvdc)

# randomly introduce MCAR univariate and multivariate missing

perc.miss <- 0.15
setseed   <- set.seed(13)
x.samp.miss <- MCAR(x.samp, perc.miss, setseed)
x.samp.miss <- x.samp.miss@"db.missing"

# impute missing values

imp <- CoImp(x.samp.miss, n.marg=n.marg, smoothing = c(0.45,0.2,0.5), plot=TRUE,
        q.lo=rep(0.1,n.marg), q.up=rep(0.9,n.marg), model=list(claytonCopula(0.5,
        dim=n.marg),  rotCopula(claytonCopula(0.5,dim=n.marg))));

# methods show and plot

show(imp)
plot(imp)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the % R documentation directory.
\keyword{imputation}
\keyword{copula}
\keyword{multivariate}
