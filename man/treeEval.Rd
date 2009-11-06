\name{treeEval}
\alias{treeEval}
\title{Classification tree evaluation by CV}
\description{
Evaluation for classification trees by cross-validation
}
\usage{
treeEval(X, grp, train, kfold = 10, cp = seq(0.01, 0.1, by = 0.01), plotit = TRUE, 
   legend = TRUE, legpos = "bottomright", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{standardized complete X data matrix (training and test data)}
  \item{grp}{factor with groups for complete data (training and test data)}
  \item{train}{row indices of X indicating training data objects}
  \item{kfold}{number of folds for cross-validation}
  \item{cp}{range for tree complexity parameter, see \code{\link[rpart]{rpart}}}
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
  \item{cp}{range for tree complexity parameter, taken from input}
}
\references{ K. Varmuza and P. Filzmoser: Introduction to Multivariate Statistical Analysis
in Chemometrics. CRC Press, Boca Raton, FL, 2009. }
\author{ Peter Filzmoser <P.Filzmoser@tuwien.ac.at> }
\seealso{\code{\link[rpart]{rpart}}}
\examples{
data(fgl,package="MASS")
grp=fgl$type
X=scale(fgl[,1:9])
k=length(unique(grp))
dat=data.frame(grp,X)
n=nrow(X)
ntrain=round(n*2/3)
require(rpart)
set.seed(123)
train=sample(1:n,ntrain)
par(mar=c(4,4,3,1))
restree=treeEval(X,grp,train,cp=c(0.01,0.02:0.05,0.1,0.15,0.2:0.5,1))
title("Classification trees")

}
\keyword{ multivariate }
