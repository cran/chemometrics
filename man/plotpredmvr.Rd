\name{plotpredmvr}
\alias{plotpredmvr}
\title{Plot predictions from repeated DCV}
\description{
Generate plot showing predicted values for Repeated Double Cross Validation
}
\usage{
plotpredmvr(mvrdcvobj, optcomp, y, X, method = "simpls", ...)
}
\arguments{
  \item{mvrdcvobj}{object from repeated double-CV, see \code{\link{mvr_dcv}}}
  \item{optcomp}{optimal number of components}
  \item{y}{data from response variable}
  \item{X}{data with explanatory variables}
  \item{method}{the multivariate regression method to be used, see 
  \code{\link[pls]{mvr}}}
  \item{\dots}{additional plot arguments}
}
\details{
After running repeated double-CV, this plot visualizes the predicted values.
}
\value{
  A plot is generated.
}
\references{ K. Varmuza and P. Filzmoser: Introduction to Multivariate Statistical 
Analysis in Chemometrics. CRC Press, Boca Raton, FL, 2009. }
\author{ Peter Filzmoser <P.Filzmoser@tuwien.ac.at> }
\seealso{\code{\link[pls]{mvr}}}
\examples{
data(NIR)
X <- NIR$xNIR[1:30,]      # first 30 observations - for illustration
y <- NIR$yGlcEtOH[1:30,1] # only variable Glucose
NIR.Glc <- data.frame(X=X, y=y)
res <- mvr_dcv(y~.,data=NIR.Glc,ncomp=10,method="simpls",repl=10)
plot3 <- plotpredmvr(res,opt=7,y,X,method="simpls")
}
\keyword{ multivariate }
