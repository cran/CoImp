\name{CoImp-class}
\docType{class}
\alias{CoImp-class}
\alias{show,CoImp-method}
\alias{plot,CoImp,missing-method}

\title{Class "CoImp"}
\description{A class for \code{CoImp} and its extensions}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CoImp", ...)}.
}
\section{Slots}{
     \describe{
    \item{\code{Missing.data.matrix}:}{Object of class \code{"matrix"}. Original missing data matrix to be imputed.}
    \item{\code{Perc.miss}:}{Object of class \code{"matrix"}. Missing and available data percentage for each variable.}
    \item{\code{Estimated.Model}:}{Object of class \code{"list"}. The list contains:
        \tabular{ll}{
            \code{model}    \tab the copula model selected and estimated on the complete cases. \cr
            \code{dimension}\tab the dimension of the \code{model}.\cr
            \code{parameter}\tab the estimated dependence parameter of the \code{model}.\cr
            \code{number}   \tab the index of the estimated model in the list of models given in input.\cr
            }
        }
    \item{\code{Estimation.Method}:}{Object of class \code{"character"}. The estimation method used for the copula model in \code{Estimated.Model}. Allowed methods are in \code{\link{fitCopula}}.}
    \item{\code{Index.matrix.NA}:}{Object of class \code{"matrix"}. Matrix of row and column indexes of missing data.}
    \item{\code{Smooth.param}:}{Object of class \code{"numeric"}. The values of the nearest neighbor component of the smoothing parameter of the \code{lp} function.}
    \item{Imputed.data.matrix}{Object of class \code{"matrix"}. The imputed data matrix.}
    \item{Estimated.Model.Imp}{Object of class \code{"list"}. The list contains:
        \tabular{ll}{
            \code{model}    \tab the copula model selected and estimated on the imputed cases. \cr
            \code{dimension}\tab the dimension of the \code{model}.\cr
            \code{parameter}\tab the estimated dependence parameter of the \code{model}.\cr
            \code{number}   \tab the index of the estimated model in the list of models given in input.\cr
            }
        }
    \item{Estimation.Method.Imp}{Object of class \code{"character"}.The estimation method used for the copula model in \code{Estimated.Model.Imp}. Allowed methods are in \code{\link{fitCopula}}.}
  }
}

\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "CoImp", y = "missing")}: ...}
    \item{show}{\code{signature(object = "CoImp")}: ...}
	 }
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
    See Also \code{\link{CoImp}}, \code{\link{MCAR}}, \code{\link{MAR}}, \code{\link{lp}} and \code{\linkS4class{copula}}.
}

\examples{
showClass("CoImp")
}

\keyword{classes}
