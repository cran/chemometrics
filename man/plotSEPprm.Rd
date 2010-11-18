\name{plotSEPprm}
\alias{plotSEPprm}
\title{Plot trimmed SEP from repeated DCV of PRM}
\description{
Generate plot showing trimmed SEP values for Repeated Double Cross Validation for
Partial RObust M-Regression (PRM)
}
\usage{
plotSEPprm(prmdcvobj, optcomp, y, X, complete = TRUE, ...)
}
\arguments{
  \item{prmdcvobj}{object from repeated double-CV of PRM, see \code{\link{prm_dcv}}}
  \item{optcomp}{optimal number of components}
  \item{y}{data from response variable}
  \item{X}{data with explanatory variables}
  \item{complete}{if TRUE the trimmed SEPcv values are drawn and computed from
      \code{\link{prm_cv}} for the same range of components as included in the prmdcvobj 
      object; if FALSE only optcomp components are computed and their results are displayed}
  \item{\dots}{additional arguments ofr \code{\link{prm_cv}}}
}
\details{
After running repeated double-CV for PRM, this plot visualizes the distribution of the SEP values.
While the gray lines represent the resulting trimmed SEP values from repreated double CV,
the black line is the result for standard CV with PRM, and it is usually too optimistic.
}
\value{
  \item{SEPdcv}{all trimmed SEP values from repeated double-CV}
  \item{SEPcv}{trimmed SEP values from usual CV}
}
\references{ K. Varmuza and P. Filzmoser: Introduction to Multivariate Statistical Analysis
in Chemometrics. CRC Press, Boca Raton, FL, 2009. }
\author{ Peter Filzmoser <P.Filzmoser@tuwien.ac.at> }
\seealso{\code{\link{prm}}}
\examples{
require(pls)
data(yarn)
res <- prm_dcv(yarn$NIR,yarn$density,a=4,repl=2)
plot1 <- plotSEPprm(res,opt=res$afinal,yarn$density,yarn$NIR)
}
\keyword{ multivariate }

