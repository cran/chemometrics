\name{plotpredprm}
\alias{plotpredprm}
\title{Plot predictions from repeated DCV of PRM}
\description{
Generate plot showing predicted values for Repeated Double Cross Validation
of Partial Robust M-regression
}
\usage{
plotpredprm(prmdcvobj, optcomp, y, X, ...)
}
\arguments{
  \item{prmdcvobj}{object from repeated double-CV of PRM, see \code{\link{prm_dcv}}}
  \item{optcomp}{optimal number of components}
  \item{y}{data from response variable}
  \item{X}{data with explanatory variables}
  \item{\dots}{additional plot arguments}
}
\details{
After running repeated double-CV for PRM, this plot visualizes the predicted values.
The result is compared with predicted values obtained via usual CV of PRM.
}
\value{
  A plot is generated.
}
\references{ K. Varmuza and P. Filzmoser: Introduction to Multivariate Statistical Analysis
in Chemometrics. CRC Press, Boca Raton, FL, 2009. }
\author{ Peter Filzmoser <P.Filzmoser@tuwien.ac.at> }
\seealso{\code{\link{prm}}}
\examples{
require(pls)
data(yarn)
res <- prm_dcv(yarn$NIR,yarn$density,a=4,repl=2)
plot3 <- plotpredprm(res,opt=res$afinal,yarn$density,yarn$NIR)
}
\keyword{ multivariate }
