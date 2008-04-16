\name{pls_eigen}
\alias{pls_eigen}
\title{Eigenvector algorithm for PLS}
\description{
Computes the PLS solution by eigenvector decompositions.
}
\usage{
pls_eigen(X, Y, a)
}
\arguments{
  \item{X}{X input data, centered (and scaled)}
  \item{Y}{Y input data, centered (and scaled)}
  \item{a}{number of PLS components}
}
\details{
The X loadings (P) and scores (T) are found by the eigendecomposition
of X'YY'X. The Y loadings (Q) and scores (U) come from the
eigendecomposition of Y'XX'Y. The resulting P and Q are orthogonal.
The first score vectors are the same as for standard PLS, subsequent 
score vectors different.
}
\value{
  \item{P}{matrix with loadings for X}
  \item{T}{matrix with scores for X}
  \item{Q}{matrix with loadings for Y}
  \item{U}{matrix with scores for Y}
}
\references{ K. Varmuza and P. Filzmoser: Introduction to Multivariate Statistical Analysis
in Chemometrics. CRC Press. To appear. }
\author{ Peter Filzmoser <P.Filzmoser@tuwien.ac.at> }
\seealso{\code{\link[pls]{mvr}}}
\examples{
data(cereal)
res <- pls_eigen(cereal$X,cereal$Y,a=6)
}
\keyword{ multivariate }

