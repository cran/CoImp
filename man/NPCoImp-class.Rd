\name{NPCoImp-class}
\docType{class}
\alias{NPCoImp-class}
\alias{show,NPCoImp-method}

\title{Class "NPCoImp"}
\description{A class for \code{NPCoImp} and its extensions}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("NPCoImp", ...)}.
}
\section{Slots}{
     \describe{
    \item{\code{Imputed.matrix}}{Object of class \code{"matrix"}. The imputed data matrix.}
    \item{\code{Selected.alpha}}{Object of class \code{"vector"}. The (conditional) probability of the lower-orthant quantile selected for the imputation.}
    \item{\code{numFlat}}{Object of class \code{"numeric"}. The number of possible flat empirical copula-based conditional probability functions, i.e. when the function cannot be empirically estimated.}
}
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "NPCoImp")}: ...} 
	 }
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
showClass("NPCoImp")
}
\keyword{classes}
