\name{Moutlier}
\alias{Moutlier}
\title{ Plots classical and robust Mahalanobis distances }
\description{
For multivariate outlier detection the Mahalanobis distance can be used.
Here a plot of the classical and the robust (based on the MCD) 
Mahalanobis distance is drawn.
}
\usage{
Moutlier(X, quantile = 0.975, plot = TRUE, ...)
}
\arguments{
  \item{X}{ numeric data frame or matrix }
  \item{quantile}{ cut-off value (quantile) for the Mahalanobis distance }
  \item{plot}{ if TRUE a plot is generated }
  \item{\dots}{ additional graphics parameters, see \code{\link{par}} }
}
\details{
For multivariate normally distributed data, a fraction of 1-quantile of data
can be declared as potential multivariate outliers. These would be identified
with the Mahalanobis distance based on classical mean and covariance.
For deviations from multivariate normality center and covariance have to 
be estimated in a robust way, e.g. by the MCD estimator. The resulting
robust Mahalanobis distance is suitable for outlier detection. Two plots
are generated, showing classical and robust Mahalanobis distance versus
the observation numbers. 
}
\value{
  \item{md }{Values of the classical Mahalanobis distance }
  \item{rd }{Values of the robust Mahalanobis distance }
  \item{cutoff }{Value with the outlier cut-off}
  ...
}
\references{ K. Varmuza and P. Filzmoser: Introduction to Multivariate Statistical Analysis
in Chemometrics. CRC Press. To appear. }
\author{ Peter Filzmoser <P.Filzmoser@tuwien.ac.at> }
\seealso{\code{\link[robustbase]{covMcd}}}
\examples{
data(glass)
data(glass.grp)
x=glass[,c(2,7)]
require(robustbase)
res <- Moutlier(glass,quantile=0.975,pch=glass.grp)
}
\keyword{ multivariate }
\keyword{ robust }
