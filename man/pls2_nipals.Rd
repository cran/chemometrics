\name{pls2_nipals}
\alias{pls2_nipals}
\title{PLS2 by NIPALS}
\description{
NIPALS algorithm for PLS2 regression (y is multivariate)
}
\usage{
pls2_nipals(X, Y, a, it = 50, tol = 1e-08, scale = FALSE)
}
\arguments{
  \item{X}{original X data matrix}
  \item{Y}{original Y-data matrix}
  \item{a}{number of PLS components}
  \item{it}{number of iterations}
  \item{tol}{tolerance for convergence}
  \item{scale}{if TRUE the X and y data will be scaled in addition to centering,
   if FALSE only mean centering is performed}
}
\details{
The NIPALS algorithm is the originally proposed algorithm for PLS. Here, the
Y-data matrix is multivariate.
}
\value{
  \item{P}{matrix with loadings for X}
  \item{T}{matrix with scores for X}
  \item{Q}{matrix with loadings for Y}
  \item{U}{matrix with scores for Y}
  \item{D}{D-matrix within the algorithm}
  \item{W}{weights for X}
  \item{C}{weights for Y}
  \item{B}{final regression coefficients}
}
\references{ K. Varmuza and P. Filzmoser: Introduction to Multivariate Statistical Analysis
in Chemometrics. CRC Press, Boca Raton, FL, 2009. }
\author{ Peter Filzmoser <P.Filzmoser@tuwien.ac.at> }
\seealso{\code{\link[pls]{mvr}}, \code{\link{pls1_nipals}}}
\examples{
data(cereal)
res <- pls2_nipals(cereal$X,cereal$Y,a=5)
}
\keyword{ multivariate }
