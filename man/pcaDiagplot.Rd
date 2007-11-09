\name{pcaDiagplot}
\alias{pcaDiagplot}
\title{ Diagnostics plot for PCA }
\description{
Score distances and of the orthogonal distances are computed and plotted.
}
\usage{
pcaDiagplot(X, X.pca, a = 2, quantile = 0.975, plot = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{ numeric data frame or matrix }
  \item{X.pca}{ PCA object resulting e.g. from \code{\link{princomp}} }
  \item{a}{ number of principal components }
  \item{quantile}{ quantile for the critical cut-off values }
  \item{plot}{ if TRUE a plot is generated }
  \item{\dots}{ additional graphics parameters, see \code{\link{par}} }
}
\details{
The score distance measures the outlyingness of the onjects within the PCA
space using Mahalanobis distances. The orthogonal distance measures the 
distance of the objects orthogonal to the PCA space. Cut-off values for both
distance measures help to distinguish between outliers and regular observations.
}
\value{
  \item{SDist }{Score distances}
  \item{ODist }{Orthogonal distances}
  \item{critSD }{critical cut-off value for the score distances}
  \item{critOD }{critical cut-off value for the orthogonal distances}
}
\references{ K. Varmuza and P. Filzmoser: Introduction to Multivariate Statistical Analysis
in Chemometrics. CRC Press. To appear. }
\author{ P. Filzmoser <P.Filzmoser@tuwien.ac.at}
\seealso{\code{\link{princomp}}}
\examples{
data(glass)
require(robustbase)
glass.mcd <- covMcd(glass)
rpca <- princomp(glass,covmat=glass.mcd)
res <- pcaDiagplot(glass,rpca,a=2)
}
\keyword{ multivariate }
\keyword{ robust }
