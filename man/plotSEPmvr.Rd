\name{plotSEPmvr}
\alias{plotSEPmvr}
\title{Plot SEP from repeated DCV}
\description{
Generate plot showing SEP values for Repeated Double Cross Validation 
}
\usage{
plotSEPmvr(mvrdcvobj, optcomp, y, X, method = "simpls", complete = TRUE, ...)
}
\arguments{
  \item{mvrdcvobj}{object from repeated double-CV, see \code{\link{mvr_dcv}}}
  \item{optcomp}{optimal number of components}
  \item{y}{data from response variable}
  \item{X}{data with explanatory variables}
  \item{method}{the multivariate regression method to be used, see \code{\link[pls]{mvr}}}
  \item{complete}{if TRUE the SEPcv values are drawn and computed for the same range of
      components as included in the mvrdcvobj object; if FALSE only optcomp components
      are computed and their results are displayed}
  \item{\dots}{additional plot arguments}
}
\details{
After running repeated double-CV, this plot visualizes the distribution of the SEP values.
}
\value{
  \item{SEPdcv}{all SEP values from repeated double-CV}
  \item{SEPcv}{SEP values from classical CV}
}
\references{ K. Varmuza and P. Filzmoser: Introduction to Multivariate Statistical Analysis
in Chemometrics. CRC Press, Boca Raton, FL, 2009. }
\author{ Peter Filzmoser <P.Filzmoser@tuwien.ac.at> }
\seealso{\code{\link[pls]{mvr}}}
\examples{
require(pls)
data(yarn)
res <- mvr_dcv(density~NIR,ncomp=10,data=yarn,method="simpls",repl=10)
plot1 <- plotSEPmvr(res,opt=7,yarn$density,yarn$NIR,method="simpls")
}
\keyword{ multivariate }

