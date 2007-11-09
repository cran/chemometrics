\name{pcaCV}
\alias{pcaCV}
\title{ Determine the number of PCA components with repeated cross validation }
\description{
By splitting data into training and test data repeatedly the number of
principal components can be determined by inspecting the distribution
of the explained variances.
}
\usage{
pcaCV(X, amax, center = TRUE, scale = TRUE, repl = 50, segments = 4, segment.type = c("random", "consecutive", "interleaved"), length.seg, trace = FALSE, plot.opt = TRUE, ...)
}
\arguments{
  \item{X}{ numeric data frame or matrix }
  \item{amax}{ maximum number of components for evaluation }
  \item{center}{ should the data be centered? TRUE or FALSE  }
  \item{scale}{ should the data be scaled? TRUE or FALSE  }
  \item{repl}{ number of replications of the CV procedure }
  \item{segments}{ number of segments for CV }
  \item{segment.type}{ "random", "consecutive", "interleaved" splitting into training and test data }
  \item{length.seg}{ number of parts for training and test data, overwrites segments }
  \item{trace}{ if TRUE intermediate results are reported }
  \item{plot.opt}{ if TRUE the results are shown by boxplots }
  \item{\dots}{ additional graphics parameters, see \code{\link{par}} }
}
\details{
For cross validation the data are split into a number of segments, PCA is
computed (using 1 to amax components) for all but one segment, and the 
scores of the segment left out are calculated. This is done in turn, by 
omitting each segment one time. Thus, a complete score matrix results for
each desired number of components, and the error martrices of fit can be
computed. A measure of fit is the explained variance, which is computed
for each number of components. Then the whole procedure is repeated
(repl times), which results in repl numbers of explained variance for
1 to amax components, i.e. a matrix. The matrix is presented by boxplots,
where each boxplot summarized the explained variance for a certain number
of principal components.
}
\value{
  \item{ExplVar }{matrix with explained variances, repl rows, and amax columns}
  \item{MSEP }{matrix with MSEP values, repl rows, and amax columns}
}
\references{ K. Varmuza and P. Filzmoser: Introduction to Multivariate Statistical Analysis
in Chemometrics. CRC Press. To appear. }
\author{ P. Filzmoser <P.Filzmoser@tuwien.ac.at} }
\seealso{\code{\link{princomp}}}
\examples{
data(glass)
require(robustbase)
res <- pcaCV(glass,segments=4,repl=100,cex.lab=1.2,ylim=c(0,1),las=1)
}
\keyword{ multivariate }
