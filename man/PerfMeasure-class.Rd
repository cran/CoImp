\name{PerfMeasure-class}
\alias{PerfMeasure-class}
\alias{show,PerfMeasure-method}

\title{Class "PerfMeasure"}
\description{A class for \code{PerfMeasure} and its extensions}

\section{Objects from the Class}{
    Objects can be created by calls of the form \code{new("PerfMeasure", ...)}.
}

\section{Slots}{
    \describe{
        \item{\code{MARE}:}{Object of class \code{"numeric"}. The mean (on the replications performed) of the absolute relative error between the imputed and the corresponding original value.}
        \item{\code{RB}:}{Object of class \code{"numeric"}. The relative bias of the estimator for the dependence parameter.}
        \item{\code{RRMSE}:}{Object of class \code{"numeric"}. The relative root mean squared error of the estimator for the dependence parameter.}
        %\item{\code{MS}:}{Object of class \code{"numeric"}. The average difference (on the replications performed) between the dependence of
        %    the imputed variables and that of the complete ones computed in terms of \code{\link{Srho.ts}}.}
        \item{\code{TID}:}{Object of class \code{"vector"}. Upper and lower tail dependence indexes for bivariate copulas. Original function is in \code{\link[copula]{tailIndex}}.}
    }
}

\section{Methods}{
    \describe{
        \item{show}{\code{signature(object = "PerfMeasure")}: ... }
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
    \code{\link{CoImp}},
    \code{\link{NPCoImp}},
    \code{\link{MCAR}},
    \code{\link{MAR}}.
}
\examples{
showClass("PerfMeasure")
}
\keyword{classes}
