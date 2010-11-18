\name{prm_dcv}
\alias{prm_dcv}
\title{Repeated double-cross-validation for robust PLS}
\description{
Performs a careful evaluation by repeated double-CV for robust PLS, called
PRM (partial robust M-estimation).
}
\usage{
prm_dcv(X,Y,a=10,repl=10,segments0=4,segments=7,segment0.type="random",
  segment.type="random",sdfact=2,fairct=4,trim=0.2,opt="median",plot.opt=FALSE, ...) 
}
\arguments{
  \item{X}{predictor matrix}
  \item{Y}{response variable}
  \item{a}{number of PLS components}
  \item{repl}{Number of replicattion for the double-CV}
  \item{segments0}{the number of segments to use for splitting into training and 
      test data, or a list with segments (see \code{\link[pls]{mvrCv}})}
  \item{segments}{the number of segments to use for selecting the optimal number if components,
           or a list with segments (see \code{\link[pls]{mvrCv}})}
  \item{segment0.type}{the type of segments to use.  Ignored if 'segments0' is a list}
  \item{segment.type}{the type of segments to use.  Ignored if 'segments' is a list}
  \item{sdfact}{factor for the multiplication of the standard deviation for
          the determination of the optimal number of components, see \code{\link{mvr_dcv}}}
  \item{fairct}{tuning constant, by default fairct=4}
  \item{trim}{trimming percentage for the computation of the SEP}
  \item{opt}{if "l1m" the mean centering is done by the l1-median,
       otherwise if "median", by the coordinate-wise median}
  \item{plot.opt}{if TRUE a plot will be generated that shows the selection of the
          optimal number of components for each step of the CV}
  \item{\dots}{additional parameters}
}
\details{
In this cross-validation (CV) scheme, the optimal number of components is determined by an additional
CV in the training set, and applied to the test set. The procedure is repeated repl times.
The optimal number of components is the model with the smallest number of components which 
is still in the range of the MSE+sdfact*sd(MSE), where MSE and sd are taken from the minimum.
}
\value{
  \item{b}{estimated regression coefficients}
  \item{intercept}{estimated regression intercept}
  \item{resopt}{array [nrow(Y) x ncol(Y) x repl] with residuals using optimum number of components}
  \item{predopt}{array [nrow(Y) x ncol(Y) x repl] with predicted Y using optimum number of components}
  \item{optcomp}{matrix [segments0 x repl] optimum number of components for each training set}
  \item{residcomp}{array [nrow(Y) x ncomp x repl] with residuals using optimum number of components}
  \item{pred}{array [nrow(Y) x ncol(Y) x ncomp x repl] with predicted Y for all numbers of components}
  \item{SEPall}{matrix [ncomp x repl] with SEP values}
  \item{SEPtrim}{matrix [ncomp x repl] with trimmed SEP values}
  \item{SEPcomp}{vector of length ncomp with trimmed SEP values; use the element afinal for the 
		optimal trimmed SEP}
  \item{afinal}{final optimal number of components}
  \item{SEPopt}{trimmed SEP over all residuals using optimal number of components}
}
\references{ K. Varmuza and P. Filzmoser: Introduction to Multivariate Statistical Analysis
in Chemometrics. CRC Press, Boca Raton, FL, 2009. }
\author{ Peter Filzmoser <P.Filzmoser@tuwien.ac.at> }
\seealso{\code{\link[pls]{mvr}}}
\examples{
require(pls)
data(yarn)
res <- prm_dcv(yarn$NIR,yarn$density,a=3,repl=2)
}
\keyword{ multivariate }
