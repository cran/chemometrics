\name{lassocoef}
\alias{lassocoef}
\title{Plot Lasso coefficients}
\description{
Plots the coefficients of Lasso regression
}
\usage{
lassocoef(formula, data, sopt, plot.opt = TRUE, ...)
}
\arguments{
  \item{formula}{formula, like y~X, i.e., dependent~response variables}
  \item{data}{data frame to be analyzed}
  \item{sopt}{optimal fraction from Lasso regression, see details}
  \item{plot.opt}{if TRUE a plot will be generated}
  \item{\dots}{additional plot arguments}
}
\details{
Using the function \code{\link{lassoCV}} for cross-validation, the optimal
fraction sopt can be determined. Besides a plot for the Lasso coefficients
for all values of fraction, the optimal fraction is taken to compute the
number of coefficients that are exactly zero.
}
\value{
  \item{coefficients}{regression coefficients for the optimal Lasso parameter}
  \item{sopt}{optimal value for fraction}
  \item{numb.zero}{number of zero coefficients for optimal fraction}
  \item{numb.nonzero}{number of nonzero coefficients for optimal fraction}
  \item{ind}{index of fraction with optimal choice for fraction}
}
\references{ K. Varmuza and P. Filzmoser: Introduction to Multivariate Statistical Analysis
in Chemometrics. CRC Press, Boca Raton, FL, 2009. }
\author{ Peter Filzmoser <P.Filzmoser@tuwien.ac.at> }
\seealso{\code{\link[lars]{cv.lars}}, \code{\link{lassoCV}}}
\examples{
data(PAC)
res=lassocoef(y~X,data=PAC,sopt=0.3)
}
\keyword{ multivariate }
