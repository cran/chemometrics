\name{sd_trim}
\alias{sd_trim}
\title{ Trimmed standard deviation }
\description{
The trimmed standard deviation as a robust estimator of scale is computed.
}
\usage{
sd_trim(x,trim=0.2,const=TRUE)
}
\arguments{
  \item{x}{ numeric vector, data frame or matrix }
  \item{trim}{ trimming proportion; should be between 0 and 0.5 }
  \item{const}{ if TRUE, the appropriate consistency correction is done }
}
\details{
The trimmed standard deviation is defined as the average trimmed sum of
squared deviations around the trimmed mean. A consistency factor for
normal distribution is included. However, this factor is only available
now for trim equal to 0.1 or 0.2. For different trimming percentages the
appropriate constant needs to be used.
If the input is a data matrix, the trimmed standard deviation of the
columns is computed.
}
\value{
   Returns the trimmed standard deviations of the vector x, or in case of
a matrix, of the columns of x.
}
\references{ K. Varmuza and P. Filzmoser: Introduction to Multivariate Statistical 
Analysis in Chemometrics. CRC Press, Boca Raton, FL, 2009. }
\author{ Peter Filzmoser <P.Filzmoser@tuwien.ac.at> }
\seealso{ \code{\link{sd}},\code{\link{mean}} }
\examples{
x <- c(rnorm(100),100) # outlier 100 is included
sd(x) # classical standard deviation
sd_trim(x) # trimmed standard deviation
}
\keyword{ manip }
