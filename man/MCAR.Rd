\name{MCAR}
\alias{MCAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generation of multivariate MCAR data}
\description{
    Introduction of artificial missing completely at random (MCAR) data in a given data set. Missing values are multivariate and have generic pattern.
}
\usage{
MCAR(db.complete, perc.miss = 0.3, setseed = 13, mcols = NULL, ...)

}
 %- maybe also 'usage' for other objects documented here.
\arguments{
  \item{db.complete}{the complete data matrix.}
  \item{perc.miss}{the percentage of missing value to be generated.}
  \item{setseed}{the seed for the generation of the missing values.}
  \item{mcols}{the index of the columns in which to introduce MCAR values.}
  \item{\dots}{further parameters for \code{\link[copula]{fitCopula}}.}
}
\details{
%    \describe{
%    \item{Usage for continuous variables:}{.\cr
%    \code{MCAR <- function(db.complete, perc.miss = 0.3, setseed = 13, mcols = NULL ...)}}
%    }

    MCAR introduce artificial missing completely at random values in a given complete data set. Missing values are multivariate and have generic pattern.
}
\value{An object of S4 class "MCAR", which is a list with the following element:
    \item{db.missing}{Object of class \code{"matrix"}. A data set with artificial multivariate MCAR.}
}

\references{
Di Lascio, F.M.L., Giannerini, S. and Reale, A. (2015) "Exploring Copulas for the Imputation of Complex Dependent Data". Statistical Methods & Applications, 24(1), p. 159-175. DOI 10.1007/s10260-014-0287-2.

Di Lascio, F.M.L., Giannerini, S. and Reale, A. (2014) "Imputation of complex dependent data by conditional copulas: analytic versus semiparametric approach", Book of proceedings of the 21st International Conference on Computational Statistics (COMPSTAT 2014), p. 491-497. ISBN 9782839913478.

Bianchi, G. Di Lascio, F.M.L. Giannerini, S. Manzari, A. Reale, A. and Ruocco, G. (2009) "Exploring copulas for the imputation of missing nonlinearly dependent data". Proceedings of the VII Meeting Classification and Data Analysis Group of the Italian Statistical Society (Cladag), Editors: Salvatore Ingrassia and Roberto Rocci, Cleup, p. 429-432. ISBN: 978-88-6129-406-6.
}

\author{
F. Marta L. Di Lascio <marta.dilascio@unibz.it>,

Simone Giannerini <simone.giannerini@unibo.it>}
\seealso{
    \code{\link{CoImp}},
    \code{\link{NPCoImp}},
    \code{\link[copula]{fitCopula}}.
}
\examples{

# generate data from a 4-variate Gumbel copula with different margins

set.seed(11)
n.marg <- 4
theta  <- 5
copula <- frankCopula(theta, dim = n.marg)
mymvdc <- mvdc(copula, c("norm", "gamma", "beta","gamma"), list(list(mean=7, sd=2),
list(shape=3, rate=2), list(shape1=4, shape2=1), list(shape=4, rate=3)))
n      <- 30
x.samp <- rMvdc(n, mymvdc)

# apply MCAR by introducing 30\% of missing data

mcar   <- MCAR(db.complete = x.samp, perc.miss = 0.3, setseed = 11)

mcar

# same example as above but introducing missing only in the first and third column

mcar2   <- MCAR(db.complete = x.samp, perc.miss = 0.3, setseed = 11, mcols=c(1,3))

mcar2

}

% Add one or more standard keywords, see file 'KEYWORDS' in the % R documentation directory.
\keyword{imputation}
\keyword{copula}
\keyword{multivariate}
