\name{lassoCV}
\alias{lassoCV}
\title{CV for Lasso regression}
\description{
Performs cross-validation (CV) for Lasso regression and plots the results
in order to select the optimal Lasso parameter.
}
\usage{
lassoCV(formula, data, K = 10, fraction = seq(0, 1, by = 0.05), trace = FALSE, plot.opt = TRUE, sdfact = 2, legpos = "topright", ...)
}
\arguments{
  \item{formula}{formula, like y~X, i.e., dependent~response variables}
  \item{data}{data frame to be analyzed}
  \item{K}{the number of segments to use for CV}
  \item{fraction}{fraction for Lasso parameters to be used for evaluation, see details}
  \item{trace}{if 'TRUE', intermediate results are printed}
  \item{plot.opt}{if TRUE a plot will be generated that shows optimal choice for "fraction"}
  \item{sdfact}{factor for the standard error for selection of the optimal parameter, see details}
  \item{legpos}{position of the legend in the plot}
  \item{\dots}{additional plot arguments}
}
\details{
The parameter "fraction" is the sum of absolute values of the regression coefficients for a particular
Lasso parameter on the sum of absolute values of the regression coefficients for the maximal
possible value of the Lasso parameter (unconstrained case), see also \code{\link[lars]{lars}}.
The optimal fraction is chosen according to the following criterion:
Within the CV scheme, the mean of the SEPs is computed, as well as their standard errors. Then one 
searches for the minimum of the mean SEPs and adds sdfact*standarderror. The optimal fraction
is the smallest fraction that is below this bound.
}
\value{
  \item{cv}{CV curve at each value of fraction}
  \item{cv.error}{standard errors for each value of fraction}
  \item{SEP}{SEP value for each value of fraction}
  \item{ind}{index of fraction with optimal choice for fraction}
  \item{sopt}{optimal value for fraction}
  \item{fraction}{all values considered for fraction}
}
\references{ K. Varmuza and P. Filzmoser: Introduction to Multivariate Statistical Analysis
in Chemometrics. CRC Press, Boca Raton, FL, 2009. }
\author{ Peter Filzmoser <P.Filzmoser@tuwien.ac.at> }
\seealso{\code{\link[lars]{cv.lars}}, \code{\link{lassocoef}}}
\examples{
data(PAC)
res=lassoCV(y~X,data=PAC,K=5,fraction=seq(0.1,0.5,by=0.1))
}
\keyword{ multivariate }
