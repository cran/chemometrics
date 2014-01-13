\name{nnetEval}
\alias{nnetEval}
\title{Neural network evaluation by CV}
\description{
Evaluation for Artificial Neural Network (ANN) classification by cross-validation
}
\usage{
nnetEval(X, grp, train, kfold = 10, decay = seq(0, 10, by = 1), size = 30, 
maxit = 100, plotit = TRUE, legend = TRUE, legpos = "bottomright", ...)
}
\arguments{
  \item{X}{standardized complete X data matrix (training and test data)}
  \item{grp}{factor with groups for complete data (training and test data)}
  \item{train}{row indices of X indicating training data objects}
  \item{kfold}{number of folds for cross-validation}
  \item{decay}{weight decay, see \code{\link[nnet]{nnet}}, 
       can be a vector with several values - but then "size" can be only one value}
  \item{size}{number of hidden units, see \code{\link[nnet]{nnet}}, 
       can be a vector with several values - but then "decay" can be only one value}
  \item{maxit}{maximal number of iterations for ANN, see \code{\link[nnet]{nnet}}}
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
  \item{decay}{value(s) for weight decay, taken from input}
  \item{size}{value(s) for number of hidden units, taken from input}
}
\references{ K. Varmuza and P. Filzmoser: Introduction to Multivariate Statistical 
Analysis in Chemometrics. CRC Press, Boca Raton, FL, 2009. }
\author{ Peter Filzmoser <P.Filzmoser@tuwien.ac.at> }
\seealso{\code{\link[nnet]{nnet}}}
\examples{
data(fgl,package="MASS")
grp=fgl$type
X=scale(fgl[,1:9])
k=length(unique(grp))
dat=data.frame(grp,X)
n=nrow(X)
ntrain=round(n*2/3)
require(nnet)
set.seed(123)
train=sample(1:n,ntrain)
resnnet=nnetEval(X,grp,train,decay=c(0,0.01,0.1,0.15,0.2,0.3,0.5,1),
   size=20,maxit=20)
}
\keyword{ multivariate }

