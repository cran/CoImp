\name{MAR}
\alias{MAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generation of multivariate missing at random (MAR) data}
\description{
    Introduction of artificial missing at random (MAR) data in a given data set. Missing values are multivariate and have generic pattern.
}
\usage{
MAR(db.complete, perc.miss = 0.3, setseed = 13, mcols = NULL, ...)

}
 %- maybe also 'usage' for other objects documented here.
\arguments{
  \item{db.complete}{the complete data matrix.}
  \item{perc.miss}{the percentage of missing values to be generated.}
  \item{setseed}{the seed for the generation of the missing values.}
  \item{mcols}{the index of the columns in which to introduce MAR values.}
  \item{\dots}{further parameters for \code{\link[copula]{fitCopula}}.}
}
\details{
%%    \describe{
%%    \item{Usage for continuous variables:}{.\cr
%%    \code{MCAR <- function(db.complete, perc.miss = 0.3, setseed = 13, ...)}}
%%    }

    MAR introduce artificial missing at random values in a given complete data set. Missing values are univariate and multivariate and have generic pattern.
}
\value{An object of S4 class "MAR", which is a list with the following element:
    \item{perc.record.missing}{Object of class \code{"numeric"}. A percentage value.}
    \item{db.missing}{Object of class \code{"matrix"}. A data set with artificial multivariate MAR.}
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
    \code{\link[CoImp]{CoImp}},
    \code{\link[CoImp]{NPCoImp}},
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

# apply MAR by introducing 30\% of missing data

mar    <- MAR(db.complete = x.samp, perc.miss = 0.3, setseed = 11)

mar
}

% Add one or more standard keywords, see file 'KEYWORDS' in the % R documentation directory.
\keyword{imputation}
\keyword{copula}
\keyword{multivariate}
