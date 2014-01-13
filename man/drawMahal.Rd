\name{drawMahal}
\alias{drawMahal}
\title{ Draws ellipses according to Mahalanobis distances }
\description{
For 2-dimensional data a scatterplot is made. Additionally, ellipses
corresponding to certain Mahalanobis distances and quantiles of
the data are drawn.
}
\usage{
drawMahal(x, center, covariance, quantile = c(0.975, 0.75, 0.5, 0.25), m = 1000, 
lwdcrit = 1, ...)
}
\arguments{
  \item{x}{ numeric data frame or matrix with 2 columns }
  \item{center}{ vector of length 2 with multivariate center of x }
  \item{covariance}{ 2 by 2 covariance matrix of x }
  \item{quantile}{ vector of quantiles for the Mahalanobis distance }
  \item{m}{ number of points where the ellipses should pass through }
  \item{lwdcrit}{ line width of the ellipses }
  \item{\dots}{ additional graphics parameters, see \code{\link{par}} }
}
\details{
For multivariate normally distributed data, a fraction of 1-quantile of data
should be outside the ellipses. For center and covariance also robust
estimators, e.g. from the MCD estimator, can be supplied.
}
\value{
A scatterplot with the ellipses is generated. 
}
\references{ K. Varmuza and P. Filzmoser: Introduction to Multivariate Statistical 
Analysis in Chemometrics. CRC Press, Boca Raton, FL, 2009. }
\author{ Peter Filzmoser <P.Filzmoser@tuwien.ac.at> }
\seealso{\code{\link[robustbase]{covMcd}}}
\examples{
data(glass)
data(glass.grp)
x=glass[,c(2,7)]
require(robustbase)
x.mcd=covMcd(x)
drawMahal(x,center=x.mcd$center,covariance=x.mcd$cov,quantile=0.975,pch=glass.grp)
}
\keyword{ multivariate }
\keyword{ robust }
