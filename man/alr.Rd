\name{alr}
\alias{alr}
\title{additive logratio transformation}
\description{
A data transformation according to the additive logratio transformation is done.
}
\usage{
alr(X, divisorvar)
}
\arguments{
  \item{X}{ numeric data frame or matrix }
  \item{divisorvar}{ number of the column of X for the variable to divide with }
}
\details{
The alr transformation is one possibility to transform compositional data
to a real space. Afterwards, the transformed data can be analyzed in the
usual way.
}
\value{
   Returns the transformed data matrix with one variable (divisor variable) less.
}
\references{ K. Varmuza and P. Filzmoser: Introduction to Multivariate Statistical Analysis
in Chemometrics. CRC Press. To appear. }
\author{ Peter Filzmoser <P.Filzmoser@tuwien.ac.at> }
\seealso{ \code{\link{clr}},\code{\link{ilr}} }
\examples{
data(glass)
glass_alr <- alr(glass,1)
}
\keyword{ manip }
