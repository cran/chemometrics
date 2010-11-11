\name{nipals}
\alias{nipals}
\title{ PCA calculation with the NIPALS algorithm }
\description{
NIPALS is an algorithm for computing PCA scores and loadings. 
}
\usage{
nipals(X, a, it = 10, tol = 1e-04)
}
\arguments{
  \item{X}{ numeric data frame or matrix }
  \item{a}{ maximum number of principal components to be computed }
  \item{it}{ maximum number of iterations }
  \item{tol}{ tolerance limit for convergence of the algorithm }
}
\details{
The NIPALS algorithm is well-known in chemometrics. It is an
algorithm for computing PCA scores and loadings. The advantage
is that the components are computed one after the other, and
one could stop at a desired number of components.
}
\value{
  \item{T }{matrix with the PCA scores}
  \item{P }{matrix with the PCA loadings}
}
\references{ K. Varmuza and P. Filzmoser: Introduction to Multivariate Statistical Analysis
in Chemometrics. CRC Press, Boca Raton, FL, 2009. }
\author{ Peter Filzmoser <P.Filzmoser@tuwien.ac.at> }
\seealso{\code{\link{princomp}}}
\examples{
data(glass)
res <- nipals(glass,a=2)
}
\keyword{ multivariate }
