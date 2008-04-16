\name{ilr}
\alias{ilr}
\title{ isometric logratio transformation }
\description{
A data transformation according to the isometric logratio transformation is done.
}
\usage{
ilr(X)
}
\arguments{
  \item{X}{ numeric data frame or matrix }
}
\details{
The ilr transformation is one possibility to transform compositional data
to a real space. Afterwards, the transformed data can be analyzed in the
usual way.
}
\value{
   Returns the transformed data matrix with one dimension less than X.
}
\references{ K. Varmuza and P. Filzmoser: Introduction to Multivariate Statistical Analysis
in Chemometrics. CRC Press. To appear. }
\author{ Peter Filzmoser <P.Filzmoser@tuwien.ac.at> }
\seealso{ \code{\link{alr}},\code{\link{clr}} }
\examples{
data(glass)
glass_ilr <- ilr(glass)
}
\keyword{ manip }
