\name{prm_cv}
\alias{prm_cv}
\title{Cross-validation for robust PLS}
\description{
Cross-validation (CV) is carried out with robust PLS based on partial robust M-regression.
A plot with the choice for the optimal number of components is generated.
This only works for univariate y-data.
}
\usage{
prm_cv(X, y, a, fairct = 4, opt = "median", subset = NULL, segments = 10, segment.type = "random", 
       trim = 0.2, sdfact = 2, plot.opt = TRUE)
}
\arguments{
  \item{X}{predictor matrix}
  \item{y}{response variable}
  \item{a}{number of PLS components}
  \item{fairct}{tuning constant, by default fairct=4}
  \item{opt}{if "l1m" the mean centering is done by the l1-median,
       otherwise by the coordinate-wise median}
  \item{subset}{optional vector defining a subset of objects}
  \item{segments}{the number of segments to use or a list with segments (see \code{\link[pls]{mvrCv}})}
  \item{segment.type}{the type of segments to use.  Ignored if 'segments' is a list}
  \item{trim}{trimming percentage for the computation of the SEP}
  \item{sdfact}{factor for the multiplication of the standard deviation for
          the determination of the optimal number of components, see \code{\link{mvr_dcv}}}
  \item{plot.opt}{if TRUE a plot will be generated that shows the selection of the
          optimal number of components for each step of the CV, see \code{\link{mvr_dcv}}}
}
\details{
A function for robust PLS based on partial robust M-regression is available at \code{\link{prm}}.
The optimal number of robust PLS components is chosen according to the following criterion:
Within the CV scheme, the mean of the trimmed SEPs SEPtrimave is computed for each number of components,
as well as their standard errors SEPtrimse. Then one searches for the minimum of the SEPtrimave values
and adds sdfact*SEPtrimse. The optimal number of components is the most parsimonious model that
is below this bound.
}
\value{
  \item{predicted}{matrix with length(y) rows and a columns with predicted values}
  \item{SEPall}{vector of length a with SEP values for each number of components}
  \item{SEPtrim}{vector of length a with trimmed SEP values for each number of components}
  \item{SEPj}{matrix with segments rows and a columns with SEP values within the CV for each number of components}
  \item{SEPtrimj}{matrix with segments rows and a columns with trimmed SEP values within the CV for each number of components}
  \item{optcomp}{final optimal number of PLS components}
  \item{SEPop}{trimmed SEP value for final optimal number of PLS components}
}
\references{ K. Varmuza and P. Filzmoser: Introduction to Multivariate Statistical Analysis
in Chemometrics. CRC Press. To appear. }
\author{ Peter Filzmoser <P.Filzmoser@tuwien.ac.at> }
\seealso{\code{\link{prm}}}
\examples{
data(cereal)
set.seed(123)
res <- prm_cv(cereal$X,cereal$Y[,1],a=10,segments=4,plot.opt=TRUE)
}
\keyword{ multivariate }
