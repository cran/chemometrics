\name{plotcompprm}
\alias{plotcompprm}
\title{Component plot for repeated DCV of PRM}
\description{
Generate plot showing optimal number of components for
Repeated Double Cross-Validation of Partial Robust M-regression}
\usage{
plotcompprm(prmdcvobj, ...)
}
\arguments{
  \item{prmdcvobj}{object from repeated double-CV of PRM, see \code{\link{prm_dcv}} }
  \item{\dots}{additional plot arguments}
}
\details{
After running repeated double-CV for PRM, this plot helps to decide on the final
number of components.
}
\value{
  \item{optcomp}{optimal number of components}
  \item{compdistrib}{frequencies for the optimal number of components}
}
\references{ K. Varmuza and P. Filzmoser: Introduction to Multivariate Statistical Analysis
in Chemometrics. CRC Press, Boca Raton, FL, 2009. }
\author{ Peter Filzmoser <P.Filzmoser@tuwien.ac.at> }
\seealso{\code{\link{prm}}}
\examples{
require(pls)
data(yarn)
res <- prm_dcv(yarn$NIR,yarn$density,a=4,repl=2)
plot2 <- plotcompprm(res)
}
\keyword{ multivariate }
