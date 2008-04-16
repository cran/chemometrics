\name{clr}
\alias{clr}
\title{ centered logratio transformation }
\description{
A data transformation according to the centered logratio transformation is done.
}
\usage{
clr(X)
}
\arguments{
  \item{X}{ numeric data frame or matrix }
}
\details{
The clr transformation is one possibility to transform compositional data
to a real space. Afterwards, the transformed data can be analyzed in the
usual way.
}
\value{
   Returns the transformed data matrix with the same dimension as X.
}
\references{ K. Varmuza and P. Filzmoser: Introduction to Multivariate Statistical Analysis
in Chemometrics. CRC Press. To appear. }
\author{ Peter Filzmoser <P.Filzmoser@tuwien.ac.at> }
\seealso{ \code{\link{alr}},\code{\link{ilr}} }
\examples{
data(glass)
glass_clr <- clr(glass)
}
\keyword{ manip }
