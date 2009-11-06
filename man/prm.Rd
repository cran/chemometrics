\name{prm}
\alias{prm}
\title{Robust PLS}
\description{
Robust PLS by partial robust M-regression.
}
\usage{
prm(X, y, a, fairct = 4, opt = "l1m",usesvd=FALSE)
}
\arguments{
  \item{X}{predictor matrix}
  \item{y}{response variable}
  \item{a}{number of PLS components}
  \item{fairct}{tuning constant, by default fairct=4}
  \item{opt}{if "l1m" the mean centering is done by the l1-median,
       otherwise if "median" the coordinate-wise median is taken}
  \item{usesvd}{if TRUE, SVD will be used if X has more columns than rows}
}
\details{
M-regression is used to robustify PLS, with initial weights based
on the FAIR weight function.
}
\value{
  \item{coef}{vector with regression coefficients}
  \item{intercept}{coefficient for intercept}
  \item{wy}{vector of length(y) with residual weights}
  \item{wt}{vector of length(y) with weights for leverage}
  \item{w}{overall weights}
  \item{scores}{matrix with PLS X-scores}
  \item{loadings}{matrix with PLS X-loadings}
  \item{fitted.values}{vector with fitted y-values}
  \item{mx}{column means of X}
  \item{my}{mean of y}
}
\references{S. Serneels, C. Croux, P. Filzmoser, and P.J. Van Espen. 
  Partial robust M-regression. Chemometrics and Intelligent Laboratory Systems, 
  Vol. 79(1-2), pp. 55-64, 2005.}
\author{ Peter Filzmoser <P.Filzmoser@tuwien.ac.at> }
\seealso{\code{\link[pls]{mvr}}}
\examples{
data(PAC)
res <- prm(PAC$X,PAC$y,a=5)
}
\keyword{ multivariate }
