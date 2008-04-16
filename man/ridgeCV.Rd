\name{ridgeCV}
\alias{ridgeCV}
\title{Repeated CV for Ridge regression}
\description{
Performs repeated cross-validation (CV) to evaluate the result of Ridge regression 
where the optimal Ridge parameter lambda was chosen on a fast evaluation scheme.
}
\usage{
ridgeCV(formula, data, lambdaopt, repl = 5, segments = 10, 
   segment.type = c("random", "consecutive", "interleaved"), length.seg, trace = FALSE,
   plot.opt = TRUE, ...)
}
\arguments{
  \item{formula}{formula, like y~X, i.e., dependent~response variables}
  \item{data}{data frame to be analyzed}
  \item{lambdaopt}{optimal Ridge parameter lambda}
  \item{repl}{number of replications for the CV}
  \item{segments}{the number of segments to use for CV,
           or a list with segments (see \code{\link[pls]{mvrCv}})}
  \item{segment.type}{the type of segments to use.  Ignored if 'segments' is a list}
  \item{length.seg}{Positive integer.  The length of the segments to use.  If
          specified, it overrides 'segments' unless 'segments' is a list}
  \item{trace}{logical; if 'TRUE', the segment number is printed for each segment}
  \item{plot.opt}{if TRUE a plot will be generated that shows the predicted
          versus the observed y-values}
  \item{\dots}{additional plot arguments}
}
\details{
Generalized Cross Validation (GCV) is used by the function \code{\link[MASS]{lm.ridge}} 
to get a quick answer for the optimal Ridge parameter. This function should make a
careful evaluation once the optimal parameter lambda has been selected. Measures
for the prediction quality are computed and optionally plots are shown.
}
\value{
  \item{residuals}{matrix of size length(y) x repl with residuals}
  \item{predicted}{matrix of size length(y) x repl with predicted values}
  \item{SEP}{Standard Error of Prediction computed for each column of "residuals"}
  \item{SEPm}{mean SEP value}
  \item{sMAD}{MAD of Prediction computed for each column of "residuals"}
  \item{sMADm}{mean of MAD values}
  \item{RMSEP}{Root MSEP value computed for each column of "residuals"}
  \item{RMSEPm}{mean RMSEP value}
}
\references{ K. Varmuza and P. Filzmoser: Introduction to Multivariate Statistical Analysis
in Chemometrics. CRC Press. To appear. }
\author{ Peter Filzmoser <P.Filzmoser@tuwien.ac.at> }
\seealso{\code{\link[MASS]{lm.ridge}}, \code{\link{plotRidge}}}
\examples{
data(PAC)
res=ridgeCV(y~X,data=PAC,lambdaopt=4.3,repl=5,segments=5)
}
\keyword{ multivariate }
