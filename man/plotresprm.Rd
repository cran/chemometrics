\name{plotresprm}
\alias{plotresprm}
\title{Plot residuals from repeated DCV of PRM}
\description{
Generate plot showing residuals for Repeated Double Cross Validation for 
Partial Robust M-regression
}
\usage{
plotresprm(prmdcvobj, optcomp, y, X, ...)
}
\arguments{
  \item{prmdcvobj}{object from repeated double-CV of PRM, see \code{\link{prm_dcv}}}
  \item{optcomp}{optimal number of components}
  \item{y}{data from response variable}
  \item{X}{data with explanatory variables}
  \item{\dots}{additional plot arguments}
}
\details{
After running repeated double-CV for PRM, this plot visualizes the residuals.
The result is compared with predicted values obtained via usual CV of PRM.
}
\value{
  A plot is generated.
}
\references{ K. Varmuza and P. Filzmoser: Introduction to Multivariate Statistical 
Analysis in Chemometrics. CRC Press, Boca Raton, FL, 2009. }
\author{ Peter Filzmoser <P.Filzmoser@tuwien.ac.at> }
\seealso{\code{\link{prm}}}
\examples{
data(NIR)
X <- NIR$xNIR[1:30,]      # first 30 observations - for illustration
y <- NIR$yGlcEtOH[1:30,1] # only variable Glucose
NIR.Glc <- data.frame(X=X, y=y)
res <- prm_dcv(X,y,a=4,repl=2)
plot4 <- plotresprm(res,opt=res$afinal,y,X)
}
\keyword{ multivariate }

