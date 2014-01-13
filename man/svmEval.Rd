\name{svmEval}
\alias{svmEval}
\title{Support Vector Machine evaluation by CV}
\description{
Evaluation for Support Vector Machines (SVM) by cross-validation
}
\usage{
svmEval(X, grp, train, kfold = 10, gamvec = seq(0, 10, by = 1), kernel = "radial", 
degree = 3, plotit = TRUE, legend = TRUE, legpos = "bottomright", ...)
}
\arguments{
  \item{X}{standardized complete X data matrix (training and test data)}
  \item{grp}{factor with groups for complete data (training and test data)}
  \item{train}{row indices of X indicating training data objects}
  \item{kfold}{number of folds for cross-validation}
  \item{gamvec}{range for gamma-values, see \code{\link[e1071]{svm}}}
  \item{kernel}{kernel to be used for SVM, should be one of "radial", "linear", 
  "polynomial", "sigmoid", default to "radial", see \code{\link[e1071]{svm}}}
  \item{degree}{degree of polynome if kernel is "polynomial", default to 3, see 
  \code{\link[e1071]{svm}}}
  \item{plotit}{if TRUE a plot will be generated}
  \item{legend}{if TRUE a legend will be added to the plot}
  \item{legpos}{positioning of the legend in the plot}
  \item{\dots}{additional plot arguments}
}
\details{
The data are split into a calibration and a test data set (provided by "train").
Within the calibration set "kfold"-fold CV is performed by applying the
classification method to "kfold"-1 parts and evaluation for the last part.
The misclassification error is then computed for the training data, for the
CV test data (CV error) and for the test data.
}
\value{
  \item{trainerr}{training error rate}
  \item{testerr}{test error rate}
  \item{cvMean}{mean of CV errors}
  \item{cvSe}{standard error of CV errors}
  \item{cverr}{all errors from CV}
  \item{gamvec}{range for gamma-values, taken from input}
}
\references{ K. Varmuza and P. Filzmoser: Introduction to Multivariate Statistical 
Analysis in Chemometrics. CRC Press, Boca Raton, FL, 2009. }
\author{ Peter Filzmoser <P.Filzmoser@tuwien.ac.at> }
\seealso{\code{\link[e1071]{svm}}}
\examples{
data(fgl,package="MASS")
grp=fgl$type
X=scale(fgl[,1:9])
k=length(unique(grp))
dat=data.frame(grp,X)
n=nrow(X)
ntrain=round(n*2/3)
require(e1071)
set.seed(143)
train=sample(1:n,ntrain)
ressvm=svmEval(X,grp,train,gamvec=c(0,0.05,0.1,0.2,0.3,0.5,1,2,5),
  legpos="topright")
title("Support vector machines")
}
\keyword{ multivariate }

