\name{lmCV}
\alias{lmCV}
\title{Repeated Cross Validation for lm}
\description{
Repeated Cross Validation for multiple linear regression: a cross-validation
is performed repeadedly, and standard evaluation measured are returned.
}
\usage{
lmCV(formula, data, repl = 100, segments = 4, segment.type = c("random", "consecutive", "interleaved"), length.seg, trace = FALSE, ...)
}
\arguments{
  \item{formula}{formula, like y~X, i.e., dependent~response variables}
  \item{data}{data set including y and X}
  \item{repl}{number of replication for Cross Validation}
  \item{segments}{number of segments used for splitting into training and test data}
  \item{segment.type}{"random", "consecutive", "interleaved" splitting into training and test data}
  \item{length.seg}{number of parts for training and test data, overwrites segments}
  \item{trace}{if TRUE intermediate results are reported}
  \item{\dots}{additional plotting arguments}
}
\details{
Repeating the cross-validation with allow for a more careful evaluation.
}
\value{
  \item{residuals}{matrix of size length(y) x repl with residuals}
  \item{predicted}{matrix of size length(y) x repl with predicted values}
  \item{SEP}{Standard Error of Prediction computed for each column of "residuals"}
  \item{SEPm}{mean SEP value}
  \item{RMSEP}{Root MSEP value computed for each column of "residuals"}
  \item{RMSEPm}{mean RMSEP value}
}
\references{ K. Varmuza and P. Filzmoser: Introduction to Multivariate Statistical Analysis
in Chemometrics. CRC Press, Boca Raton, FL, 2009. }
\author{ Peter Filzmoser <P.Filzmoser@tuwien.ac.at> }
\seealso{\code{\link[pls]{mvr}}}
\examples{
data(ash)
set.seed(100)
res=lmCV(SOT~.,data=ash,repl=10)
hist(res$SEP)
}
\keyword{ multivariate }
