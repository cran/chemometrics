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
plot2 <- plotcompmvr(res)
}
\keyword{ multivariate }
