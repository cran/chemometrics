\name{knnEval}
\alias{knnEval}
\title{kNN evaluation by CV}
\description{
Evaluation for k-Nearest-Neighbors (kNN) classification by cross-validation
}
\usage{
knnEval(X, grp, train, kfold = 10, knnvec = seq(2, 20, by = 2), plotit = TRUE, 
    legend = TRUE, legpos = "bottomright", ...)
}
\arguments{
  \item{X}{standardized complete X data matrix (training and test data)}
  \item{grp}{factor with groups for complete data (training and test data)}
  \item{train}{row indices of X indicating training data objects}
  \item{kfold}{number of folds for cross-validation}
  \item{knnvec}{range for k for the evaluation of kNN}
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
  \item{knnvec}{range for k for the evaluation of kNN, taken from input}
}
\references{ K. Varmuza and P. Filzmoser: Introduction to Multivariate Statistical Analysis
in Chemometrics. CRC Press, Boca Raton, FL, 2009. }
\author{ Peter Filzmoser <P.Filzmoser@tuwien.ac.at> }
\seealso{\code{\link[class]{knn}}}
\examples{
data(fgl,package="MASS")
grp=fgl$type
X=scale(fgl[,1:9])
k=length(unique(grp))
dat=data.frame(grp,X)
n=nrow(X)
ntrain=round(n*2/3)
require(class)
set.seed(123)
train=sample(1:n,ntrain)
resknn=knnEval(X,grp,train,knnvec=seq(1,30,by=1),legpos="bottomright")
title("kNN classification")
}
\keyword{ multivariate }

