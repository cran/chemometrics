\name{pcaVarexpl}
\alias{pcaVarexpl}
\title{ PCA diagnostics for variables}
\description{
Diagnostics of PCA to see the explained variance for each variable.
}
\usage{
pcaVarexpl(X, a, center = TRUE, scale = TRUE, plot = TRUE, ...)
}
\arguments{
  \item{X}{ numeric data frame or matrix }
  \item{a}{ number of principal components }
  \item{center}{ centring of X (FALSE or TRUE) }
  \item{scale}{ scaling of X (FALSE or TRUE) }
  \item{plot}{ if TRUE make plot with explained variance }
  \item{\dots}{ additional graphics parameters, see \code{\link{par}} }
}
\details{
For a desired number of principal components the percentage of explained
variance is computed for each variable and plotted.
}
\value{
  \item{ExplVar }{explained variance for each variable}
}
\references{ K. Varmuza and P. Filzmoser: Introduction to Multivariate Statistical Analysis
in Chemometrics. CRC Press. To appear. }
\author{ P. Filzmoser <P.Filzmoser@tuwien.ac.at} }
\seealso{\code{\link{princomp}}}
\examples{
data(glass)
res <- pcaVarexpl(glass,a=2)
}
\keyword{ multivariate }
