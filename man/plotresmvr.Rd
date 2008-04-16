\name{plotresmvr}
\alias{plotresmvr}
\title{Plot residuals from repeated DCV}
\description{
Generate plot showing residuals for Repeated Double Cross Validation
}
\usage{
plotresmvr(mvrdcvobj, optcomp, y, X, method = "simpls", ...)
}
\arguments{
  \item{mvrdcvobj}{object from repeated double-CV, see \code{\link{mvr_dcv}}}
  \item{optcomp}{optimal number of components}
  \item{y}{data from response variable}
  \item{X}{data with explanatory variables}
  \item{method}{the multivariate regression method to be used, see \code{\link[pls]{mvr}}}
  \item{\dots}{additional plot arguments}
}
\details{
After running repeated double-CV, this plot visualizes the residuals.
}
\value{
  A plot is generated.
}
\references{ K. Varmuza and P. Filzmoser: Introduction to Multivariate Statistical Analysis
in Chemometrics. CRC Press. To appear. }
\author{ Peter Filzmoser <P.Filzmoser@tuwien.ac.at> }
\seealso{\code{\link[pls]{mvr}}}
\examples{
require(pls)
data(yarn)
res <- mvr_dcv(density~NIR,ncomp=10,data=yarn,method="simpls",repl=10)
plot4 <- plotresmvr(res,opt=7,yarn$density,yarn$NIR,method="simpls")
}
\keyword{ multivariate }

