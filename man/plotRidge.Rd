\name{plotRidge}
\alias{plotRidge}
\title{Plot results of Ridge regression}
\description{
Two plots from Ridge regression are generated: The MSE resulting from
Generalized Cross Validation (GCV) versus the Ridge parameter lambda,
and the regression coefficients versus lambda. The optimal choice for
lambda is indicated.
}
\usage{
plotRidge(formula, data, lambda = seq(0.5, 50, by = 0.05), ...)
}
\arguments{
  \item{formula}{formula, like y~X, i.e., dependent~response variables}
  \item{data}{data frame to be analyzed}
  \item{lambda}{possible values for the Ridge parameter to evaluate}
  \item{\dots}{additional plot arguments}
}
\details{
For all values provided in lambda the results for Ridge regression are computed.
The function \code{\link[MASS]{lm.ridge}} is used for cross-validation and 
Ridge regression.
}
\value{
  \item{predicted}{predicted values for the optimal lambda}
  \item{lambdaopt}{optimal Ridge parameter lambda from GCV}
}
\references{ K. Varmuza and P. Filzmoser: Introduction to Multivariate Statistical Analysis
in Chemometrics. CRC Press. To appear. }
\author{ Peter Filzmoser <P.Filzmoser@tuwien.ac.at> }
\seealso{\code{\link[MASS]{lm.ridge}}, \code{\link{plotRidge}}}
\examples{
data(PAC)
res=plotRidge(y~X,data=PAC,lambda=seq(1,20,by=0.5))
}
\keyword{ multivariate }
