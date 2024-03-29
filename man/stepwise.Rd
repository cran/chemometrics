\name{stepwise}
\alias{stepwise}
\title{Stepwise regression}
\description{
Stepwise regression, starting from the empty model, with scope to the full model
}
\usage{
stepwise(formula, data, k, startM, maxTime = 1800, direction = "both", 
writeFile = FALSE, maxsteps = 500, ...)
}
\arguments{
  \item{formula}{formula, like y~X, i.e., dependent~response variables}
  \item{data}{data frame to be analyzed}
  \item{k}{sensible values are log(nrow(x)) for BIC or 2 for AIC; if not 
  provided -> BIC}
  \item{startM}{optional, the starting model; provide a binary vector}
  \item{maxTime}{maximal time to be used for algorithm}
  \item{direction}{either "forward" or "backward" or "both"}
  \item{writeFile}{if TRUE results are shown on the screen}
  \item{maxsteps}{maximum number of steps}
  \item{\dots}{additional plot arguments}
}
\details{
This function is similar to the function \code{\link[stats]{step}} for stepwise 
regression. It is especially designed for cases where the number of regressor 
variables is much higher than the number of objects. The formula for the full model 
(scope) is automatically generated.
}
\value{
  \item{usedTime}{time that has been used for algorithm}
  \item{bic}{BIC values for different models}
  \item{models}{matrix with no. of models rows and no. of variables columns, and 
  0/1 entries defining the models}
}
\references{ K. Varmuza and P. Filzmoser: Introduction to Multivariate Statistical 
Analysis in Chemometrics. CRC Press, Boca Raton, FL, 2009. }
\author{ Leonhard Seyfang and (marginally) Peter Filzmoser <P.Filzmoser@tuwien.ac.at> }
\seealso{\code{\link[stats]{step}}}
\examples{
data(NIR)
X <- NIR$xNIR[1:30,]      # first 30 observations - for illustration
y <- NIR$yGlcEtOH[1:30,1] # only variable Glucose
NIR.Glc <- data.frame(X=X, y=y)
res <- stepwise(y~.,data=NIR.Glc,maxsteps=2)   
}
\keyword{ multivariate }
