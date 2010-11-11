\name{stepwise}
\alias{stepwise}
\title{Stepwise regression}
\description{
Stepwise regression, starting from the empty model, with scope to the full model
}
\usage{
stepwise(formula, data, k, startM, maxTime = 1800, direction = "both", writeFile = FALSE,
   resname = "stepres00", maxsteps = 500, ...)
}
\arguments{
  \item{formula}{formula, like y~X, i.e., dependent~response variables}
  \item{data}{data frame to be analyzed}
  \item{k}{sensible values are log(nrow(x)) for BIC or 2 for AIC; if not provided -> BIC}
  \item{startM}{optional, the starting model; provide a binary vector}
  \item{maxTime}{maximal time to be used for algorithm}
  \item{direction}{either "forward" or "backward" or "both"}
  \item{writeFile}{if TRUE results are stored in the file "resname"}
  \item{resname}{filename where results are stored, only if writeFile is TRUE}
  \item{maxsteps}{maximum number of steps}
  \item{\dots}{additional plot arguments}
}
\details{
This function is similar to the function \code{\link[stats]{step}} for stepwise regression. 
It is especially designed for cases where the number of regressor variables is much 
higher than the number of objects. The formula for the full model (scope) is automatically 
generated.
}
\value{
  \item{usedTime}{time that has been used for algorithm}
  \item{bic}{BIC values for different models}
  \item{models}{matrix with no. of models rows and no. of variables columns, and 0/1 entries defining the models}
}
\references{ K. Varmuza and P. Filzmoser: Introduction to Multivariate Statistical Analysis
in Chemometrics. CRC Press, Boca Raton, FL, 2009. }
\author{ Leonhard Seyfang and (marginally) Peter Filzmoser <P.Filzmoser@tuwien.ac.at> }
\seealso{\code{\link[stats]{step}}}
\examples{
data(yarn)
res=stepwise(density~NIR,data=yarn,maxsteps=2) # increase maxsteps to get a good solution
}
\keyword{ multivariate }
