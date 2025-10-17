\name{NPCoImp}
\alias{NPCoImp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Nonparametric Copula-Based Imputation Method}
\description{
    Imputation method based on empirical copula-based conditional cumulative probability functions.
}
\usage{
NPCoImp(X, Psi=seq(0.05,0.45,by=0.05), smoothing="beta", K=7, method="gower", ...)
} 
\arguments{
  \item{X}{a data matrix with missing values. Missing values should be denoted with \code{NA}.}
  \item{Psi}{vector of probabilities to evaluate the radial symmetry/asymmetry of the empirical copula-based conditional cumulative probability function and find the best lower-orthant quantile for the imputation (see below for details).}
  \item{smoothing}{the character string specifying the type of smoothing of the empirical copula. Default is "beta" (empirical beta copula) but also "none" (the original empirical copula) can be used.}
  \item{K}{the number of data matrix rows more similar to the missing one that are used for the imputation.}
  \item{method}{the distance measure used for the imputation, among Euclidean, Manhattan, Canberra, Gower, and two based on the Kendall-correlation coefficient (see below for details).}  
  \item{...}{further parameters for daisy, e.g. weights.}
}
\details{
    NPCoImp is a nonparametric imputation method based on the empirical copula-based conditional cumulative probability functions. To choose the best lower-orthant quantile for the imputation it evaluates the radial (a)symmetry of the empirical copula-based conditional cumulative probability functions and it uses the K pseudo-observations more similar to the missing one. The NPCoImp allows the imputation of missing observations according to the multivariate dependence structure of the data generating process without any assumptions on the margins. This method can be used independently from the dimension and the kind (monotone or non monotone) of the missing patterns. 
    Brief description of the approach:
   \enumerate{
    \item estimate the empirical copula-based conditional cumulative probability function of the missing observation(s) given the available ones; 
    \item evaluate the radial (a)symmetry of the empirical copula-based cumulative probability distribution function around 0.5 (see the paper in the references for details); 
    \item select the lower-orthant quantile of the empirical copula-based conditional cumulative probability function on the basis of its radial (a)symmetry (see the paper in the references for details);
    \item select the K pseudo-observations closest to the imputed one and the corresponding original observations;
    \item impute missing values by replacing them from the average of the original observations derived at the previous step;
    \item calculate the empirical copula-based conditional cumulative probability function of the lower-orthant quantile used for imputing.
    }
}
\value{An object of S4 class "NPCoImp", which is a list with the following elements:
     \item{Imputed.matrix}{the imputed data matrix.}
     \item{Selected.alpha}{the (conditional) probability of the lower-orthant quantile selected for the imputation.}
     \item{numFlat}{the number of possible flat empirical copula-based conditional cumulative probability functions, i.e. when the copula is always zero.}
}

\references{
Di Lascio, F.M.L, Gatto A. (202x) "A nonparametric copula-based imputation method". Under review.
}

\author{
F. Marta L. Di Lascio <marta.dilascio@unibz.it>,
Aurora Gatto <aurora.gatto@unibz.it>
}
\seealso{
    \code{\link[CoImp]{CoImp}},
    \code{\link[CoImp]{MCAR}},
    \code{\link[CoImp]{MAR}}.
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

perc.mis    <- 0.25
set.seed(14)
miss.row    <- sample(1:n, perc.mis*n, replace=TRUE)
miss.col    <- sample(1:n.marg, perc.mis*n, replace=TRUE)
miss        <- cbind(miss.row,miss.col)
x.samp.miss <- replace(x.samp,miss,NA)
x.samp.miss
probs <- seq(0.05,0.45,by=0.1)
ndist <- 7
dist.meth <- "gower"  

# impute missing values
NPimp <- NPCoImp(X=x.samp.miss, Psi=probs, smoothing="beta", K=ndist, 
                    method=dist.meth)

# methods show

show(NPimp)

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
probs <- seq(0.05,0.45,by=0.05)
ndist <- 7
dist.meth <- "gower" 

# impute missing values

NPimp2 <- NPCoImp(X=x.samp.miss, Psi=probs, smoothing="beta", K=ndist, 
                    method=dist.meth)

# methods show and plot

show(NPimp2)
}
}

