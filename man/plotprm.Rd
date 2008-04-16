\name{plotprm}
\alias{plotprm}
\title{Plot results from robust PLS}
\description{
The predicted values and the residuals are shown for robust PLS using the optimal number
of components.
}
\usage{
plotprm(prmobj, y, ...)
}
\arguments{
  \item{prmobj}{resulting object from CV of robust PLS, see \code{\link{prm_cv}}}
  \item{y}{vector with values of response variable}
  \item{\dots}{additional plot arguments}
}
\details{
Robust PLS based on partial robust M-regression is available at \code{\link{prm}}.
Here the function \code{\link{prm_cv}} has to be used first, applying cross-validation
with robust PLS. Then the result is taken by this routine and two plots are generated
for the optimal number of PLS components: The measured versus the predicted y, and
the predicted y versus the residuals.
}
\value{
  A plot is generated.
}
\references{ K. Varmuza and P. Filzmoser: Introduction to Multivariate Statistical Analysis
in Chemometrics. CRC Press. To appear. }
\author{ Peter Filzmoser <P.Filzmoser@tuwien.ac.at> }
\seealso{\code{\link{prm}}}
\examples{
data(cereal)
set.seed(123)
res <- prm_cv(cereal$X,cereal$Y[,1],a=10,segments=4,plot.opt=FALSE)
plotprm(res,cereal$Y[,1])
}
\keyword{ multivariate }
