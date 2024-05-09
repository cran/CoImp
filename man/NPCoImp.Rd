\name{NPCoImp}
\alias{NPCoImp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Non-Parametric Copula-Based Imputation Method}
\description{
    Imputation method based on empirical conditional copula functions.
}
\usage{
NPCoImp(X, Psi=seq(0.05,0.45,by=0.05), K=7, method="gower", N=1000, seed=set.seed(7))
} 
\arguments{
  \item{X}{a data matrix with missing values. Missing values should be denoted with \code{NA}.}
  \item{Psi}{vector of probabilities to assess the symmetry/asymmetry of the empirical conditional copula (ecc) function and find the best quantile for the imputation (see below for details).}
  \item{K}{the number of data matrix rows more similar to the missing one that are used for the imputation.}
  \item{method}{the distance measure used for the imputation, among Euclidean, Manhattan, Canberra, Gower, and two based on the Kendall-correlation coefficient (see below for details).}  
  \item{N}{the number of uniform numbers to numerically construct the empirical conditional copula.}
  \item{seed}{the set.seed for the numerical inversion Monte Carlo method used for constructing the empirical conditional copula.}
}
\details{
    NPCoImp is a non-parametric imputation method based on the empirical conditional copula function and the numerical inversion Monte Carlo method. To choose the best quantile for the imputation it assesses the (a)symmetry of the empirical conditional copula and it uses the K pseudo-observations more similar to the missing one. The NPCoImp allows the imputation of missing observations according to the multivariate dependence structure of the data generating process without any assumptions on the margins. This method can be used independently from the dimension and the kind (monotone or non monotone) of the missing patterns. 
    Brief description of the approach:
   \enumerate{
    \item numerically reconstruct the empirical conditional copula of the missing observation(s); 
    \item evaluate the (a)symmetry of the empirical conditional copula around 0.5. Inspired by Patil et al. (2012), we compute C(0.5-Psi|available obs)-( 1 - C(0.5+Psi|available obs)) where Psi in ]0,0.5[ and say that C(.|available obs) is symmetric (negative asymmetric or positive asymmetric) if the above equation is null (positive or negative, respectively);
    \item select the quantile of the empirical conditional copula on the basis of its (a)symmetry. Therefore:
    \itemize{
    \item symmetry: we impute through the median of the empirical conditional copula; 
    \item negative asymmetry: we impute with a quantile on the left tail of the ecc (see the paper in the references for details);
    \item positive asymmetry: we impute with a quantile on the right tail of the ecc (see the paper in the references for details);
    }
    \item select the K pseudo-observations closest to the imputed one and the corresponding original observations;
    \item impute missing values by replacing them from the average of the original observations derived at the previous step.
    }
}
\value{An object of S4 class "NPCoImp", which is a list with the following elements:
     \item{Imputed.matrix}{the imputed data matrix.}
     \item{Selected.alpha}{the selected alpha values. The best quantile for the imputation.}
     \item{numFlat}{the number of possible flat empirical conditional copulas, i.e. when ecc is always zero.}
}

\references{
Di Lascio, F.M.L, Gatto A. (202x) "A non-parametric conditional copula-based imputation method". Working paper.
}

\author{
F. Marta L. Di Lascio <marta.dilascio@unibz.it>,
Aurora Gatto <aurora.gatto@unibz.it>
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
probs <- seq(0.05,0.45,by=0.05)
ndist <- 7
dist.meth <- "gower"  

# impute missing values
NPimp <- NPCoImp(X=x.samp.miss, Psi=probs, K=ndist, method=dist.meth, 
                 N=1000, seed=set.seed(7))

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

NPimp2 <- NPCoImp(X=x.samp.miss, Psi=probs, K=ndist, method=dist.meth, N=1000)
# methods show and plot

show(NPimp2)
}
}

