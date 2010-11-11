\name{mvr_dcv}
\alias{mvr_dcv}
\title{Repeated double-cross-validation for PLS and PCR}
\description{
Performs a careful evaluation by repeated double-CV for multivariate regression methods, like PLS and PCR.
}
\usage{
mvr_dcv(formula, ncomp, data, subset, na.action, 
  method = c("kernelpls", "widekernelpls", "simpls", "oscorespls", "svdpc"), 
  scale = FALSE, repl = 100, sdfact = 2, 
  segments0 = 4, segment0.type = c("random", "consecutive", "interleaved"), length.seg0, 
  segments = 10, segment.type = c("random", "consecutive", "interleaved"), length.seg, 
  trace = FALSE, plot.opt = FALSE, selstrat = "hastie", ...)
}
\arguments{
  \item{formula}{formula, like y~X, i.e., dependent~response variables}
  \item{ncomp}{number of PLS components}
  \item{data}{data frame to be analyzed}
  \item{subset}{optional vector to define a subset}
  \item{na.action}{a function which indicates what should happen when the data
          contain missing values}
  \item{method}{the multivariate regression method to be used, see \code{\link[pls]{mvr}}}
  \item{scale}{numeric vector, or logical.  If numeric vector, X is scaled
          by dividing each variable with the corresponding element of
          'scale'.  If 'scale' is 'TRUE', X is scaled by dividing each
          variable by its sample standard deviation.  If
          cross-validation is selected, scaling by the standard
          deviation is done for every segment.}
  \item{repl}{Number of replicattion for the double-CV}
  \item{sdfact}{factor for the multiplication of the standard deviation for
          the determination of the optimal number of components}
  \item{segments0}{the number of segments to use for splitting into training and test data, 
           or a list with segments (see \code{\link[pls]{mvrCv}})}
  \item{segment0.type}{the type of segments to use.  Ignored if 'segments0' is a list}
  \item{length.seg0}{Positive integer.  The length of the segments to use.  If
          specified, it overrides 'segments' unless 'segments0' is a list}
  \item{segments}{the number of segments to use for selecting the optimal number if components,
           or a list with segments (see \code{\link[pls]{mvrCv}})}
  \item{segment.type}{the type of segments to use.  Ignored if 'segments' is a list}
  \item{length.seg}{Positive integer.  The length of the segments to use.  If
          specified, it overrides 'segments' unless 'segments' is a list}
  \item{trace}{logical; if 'TRUE', the segment number is printed for each segment}
  \item{plot.opt}{if TRUE a plot will be generated that shows the selection of the
          optimal number of components for each step of the CV}
  \item{selstrat}{method that defines how the optimal number of components is selected,
          should be one of "diffnext", "hastie", "relchange"; see details}
  \item{\dots}{additional parameters}
}
\details{
In this cross-validation (CV) scheme, the optimal number of components is determined by an additional
CV in the training set, and applied to the test set. The procedure is repeated repl times.
There are different strategies for determining the optimal number of components (parameter selstrat):
"diffnext" compares MSE+sdfact*sd(MSE) among the neighbors, and if the MSE falls outside this bound, 
this is the optimal number. "hastie" searches for the number of components with the minimum of the
mean MSE's. The optimal number of components is the model with the smallest number of components which 
is still in the range of the MSE+sdfact*sd(MSE), where MSE and sd are taken from the minimum.
"relchange" is a strategy where the relative change is combined with "hastie": First the minimum of the 
mean MSE's is searched, and MSE's of larger components are omitted. For this selection, the relative 
change in MSE compared to the min, and relative to the max, is computed. If this change
is very small (e.g. smaller than 0.005), these components are omitted.
Then the "hastie" strategy is applied for the remaining MSE's.
}
\value{
  \item{resopt}{array [nrow(Y) x ncol(Y) x repl] with residuals using optimum number of components}
  \item{predopt}{array [nrow(Y) x ncol(Y) x repl] with predicted Y using optimum number of components}
  \item{optcomp}{matrix [segments0 x repl] optimum number of components for each training set}
  \item{pred}{array [nrow(Y) x ncol(Y) x ncomp x repl] with predicted Y for all numbers of components}
  \item{SEPopt}{SEP over all residuals using optimal number of components}
  \item{sIQRopt}{spread of inner half of residuals as alternative robust spread measure to the SEPopt}
  \item{sMADopt}{MAD of residuals as alternative robust spread measure to the SEPopt}
  \item{MSEPopt}{MSEP over all residuals using optimal number of components}
  \item{afinal}{final optimal number of components}
  \item{SEPfinal}{vector of length ncomp with final SEP values; use the element afinal for the 
		optimal SEP}
}
\references{ K. Varmuza and P. Filzmoser: Introduction to Multivariate Statistical Analysis
in Chemometrics. CRC Press, Boca Raton, FL, 2009. }
\author{ Peter Filzmoser <P.Filzmoser@tuwien.ac.at> }
\seealso{\code{\link[pls]{mvr}}}
\examples{
require(pls)
data(yarn)
res <- mvr_dcv(density~NIR,ncomp=10,data=yarn,method="simpls",repl=10)
}
\keyword{ multivariate }
