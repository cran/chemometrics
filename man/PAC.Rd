\name{PAC}
\alias{PAC}
\docType{data}
\title{GC retention indices}
\description{
For 209 objects an X-data set (467 variables) and a y-data set (1 variable)
is available. The data describe GC-retention indices of polycyclic aromatic 
compounds (y) which have been modeled by molecular descriptors (X).
}
\usage{data(PAC)}
\format{
  A data frame with 209 objects and 2 list elements:
  \describe{
    \item{\code{y}}{numeric vector with length 209}
    \item{\code{X}}{matrix with 209 rows and 467 columns}
  }
}
\details{
The data can be used for linear and non-linear models.
}
\source{
K. Varmuza and P. Filzmoser: Introduction to Multivariate Statistical Analysis
in Chemeometrics. CRC Press. To appear. }
\references{
K. Varmuza and P. Filzmoser: Introduction to Multivariate Statistical Analysis
in Chemeometrics. CRC Press. To appear. }
\examples{
data(PAC)
names(PAC)
}
\keyword{datasets}
