\name{plotcompmvr}
\alias{plotcompmvr}
\title{Component plot for repeated DCV}
\description{
Generate plot showing optimal number of components for
Repeated Double Cross-Validation}
\usage{
plotcompmvr(mvrdcvobj, ...)
}
\arguments{
  \item{mvrdcvobj}{object from repeated double-CV, see \code{\link{mvr_dcv}} }
  \item{\dots}{additional plot arguments}
}
\details{
After running repeated double-CV, this plot helps to decide on the final
number of components.
}
\value{
  \item{optcomp}{optimal number of components}
  \item{compdistrib}{frequencies for the optimal number of components}
}
\references{ K. Varmuza and P. Filzmoser: Introduction to Multivariate Statistical Analysis
in Chemometrics. CRC Press. To appear. }
\author{ Peter Filzmoser <P.Filzmoser@tuwien.ac.at> }
\seealso{\code{\link[pls]{mvr}}}
\examples{
require(pls)
data(yarn)
res <- mvr_dcv(density~NIR,ncomp=10,data=yarn,method="simpls",repl=10)
plot2 <- plotcompmvr(res)
}
\keyword{ multivariate }
