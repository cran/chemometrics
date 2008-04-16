\name{pls1_nipals}
\alias{pls1_nipals}
\title{PLS1 by NIPALS}
\description{
NIPALS algorithm for PLS1 regression (y is univariate)
}
\usage{
pls1_nipals(X, y, a, it = 50, tol = 1e-08, scale = FALSE)
}
\arguments{
  \item{X}{original X data matrix}
  \item{y}{original y-data}
  \item{a}{number of PLS components}
  \item{it}{number of iterations}
  \item{tol}{tolerance for convergence}
  \item{scale}{if TRUE the X and y data will be scaled in addition to centering,
   if FALSE only mean centering is performed}
}
\details{
The NIPALS algorithm is the originally proposed algorithm for PLS. Here, the
y-data are only allowed to be univariate. This simplifies the algorithm.
}
\value{
  \item{P}{matrix with loadings for X}
  \item{T}{matrix with scores for X}
  \item{W}{weights for X}
  \item{C}{weights for Y}
  \item{b}{final regression coefficients}
}
\references{ K. Varmuza and P. Filzmoser: Introduction to Multivariate Statistical Analysis
in Chemometrics. CRC Press. To appear. }
\author{ Peter Filzmoser <P.Filzmoser@tuwien.ac.at> }
\seealso{\code{\link[pls]{mvr}}, \code{\link{pls2_nipals}}}
\examples{
data(PAC)
res <- pls1_nipals(PAC$X,PAC$y,a=5)
}
\keyword{ multivariate }

