\name{CoImp}
\alias{CoImp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Copula-Based Imputation Method}
\description{
    Imputation method based on conditional copula functions.
}
\usage{
CoImp(X, n.marg = 2, type.data = "continuous", smoothing = rep(0.5,n.marg),
      plot.marg = TRUE, plot.bar = TRUE, plot.legend = TRUE, args.legend =
      list(y = 110, cex = 0.8), model = list(normalCopula(0.5, dim=n.marg,
      dispstr="ex"), claytonCopula(10, dim=n.marg),gumbelCopula(10, dim=n.marg),
      frankCopula(10, dim=n.marg)), ...)
}
 %- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a data matrix with missing values. Missing values should be denoted
        with \code{NA}.}
  \item{n.marg}{the number of variables in X.}
  \item{type.data}{the nature of the variables in X: \code{discrete} or \code{continuous}.}
  \item{smoothing}{values for the nearest neighbour component of the smoothing parameter of the \code{lp} function.}
  \item{plot.marg}{logical: if \code{TRUE} plots the estimated marginal densities.}
  \item{plot.bar}{logical: if \code{TRUE} shows a bar plot of the percentages of missing and available data for each margin.}
  \item{plot.legend}{logical: see \code{\link{barplot}}.}
  \item{args.legend}{list of additional arguments to pass to \code{\link{legend}}.}
  \item{model}{a list of copula models to be used for the imputation, see the Details section.
        This should be one of \code{normal}, \code{frank}, \code{clayton} and \code{gumbel}.}
  \item{\dots}{further parameters for \code{\link{fitCopula}}, \code{\link{lp}} and further graphical arguments.}
}
\details{
%    \describe{
%    \item{Usage for continuous variables:}{.\cr
%    \code{CoImp(X, n.marg = 2, type.data = "continuous", smoothing = rep(0.5,n.marg),
%        plot.marg = TRUE, plot.bar = TRUE, plot.legend = TRUE, args.legend = list(y = 110,
%        cex = 0.8), model = list(normalCopula(0.5, dim=n.marg, dispstr="ex"),
%        claytonCopula(10, dim=n.marg),gumbelCopula(10, dim=n.marg),
%        frankCopula(10, dim=n.marg)), ...)}}
%    }

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
Di Lascio, F.M.L. Giannerini, S. and Reale A. (201x) "A multivariate technique based on conditional copula specification for the imputation of complex dependent data". Working paper.

Di Lascio, F.M.L. Giannerini, S. and Reale A. (201x) "Exploring Copulas for the Imputation of Complex Dependent Data". Under review.

Bianchi, G. Di Lascio, F.M.L. Giannerini, S. Manzari, A. Reale, A. and Ruocco, G. (2009) "Exploring copulas for the imputation of missing nonlinearly dependent data". Proceedings of the VII Meeting Classification and Data Analysis Group of the Italian Statistical Society (Cladag), Editors: Salvatore Ingrassia and Roberto Rocci, Cleup, p. 429-432. ISBN: 978-88-6129-406-6.
}

\author{
Francesca Marta Lilja Di Lascio <marta.dilascio@unibz.it>,

Simone Giannerini <simone.giannerini@unibo.it>}

\examples{

# generate data from a 4-variate Gumbel copula with different margins

set.seed(11)
n.marg <- 4
theta  <- 5
copula <- frankCopula(theta, dim = n.marg)
mymvdc <- mvdc(copula, c("norm", "gamma", "beta","gamma"), list(list(mean=7, sd=2),
list(shape=3, rate=2), list(shape1=4, shape2=1), list(shape=4, rate=3)))
n      <- 20
x.samp <- copula::rMvdc(n, mymvdc)

# randomly introduce univariate and multivariate missing

perc.mis    <- 0.5
set.seed(11)
miss.row    <- sample(1:n, perc.mis*n, replace=TRUE)
miss.col    <- sample(1:n.marg, perc.mis*n, replace=TRUE)
miss        <- cbind(miss.row,miss.col)
x.samp.miss <- replace(x.samp,miss,NA)

# impute missing values

imp <- CoImp(x.samp.miss, n.marg=n.marg, smoothing = rep(0.6,n.marg), TRUE, TRUE, TRUE,
            type.data="continuous");

# methods show and plot

show(imp)
plot(imp)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the % R documentation directory.
\keyword{imputation}
\keyword{copula}
\keyword{multivariate}
