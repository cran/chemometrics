\name{cereal}
\alias{cereal}
\docType{data}
\title{Data from cereals}
\description{
For 15 cereals an X and Y data set, measured on the same objects, is available.
The X data are 145 infrared spectra, and the Y data are 6 chemical/technical
properties (Heating value, C, H, N, Starch, Ash). Also the scaled Y data are 
included (mean 0, variance 1 for each column).
The cereals come from 5 groups B=Barley, M=Maize, R=Rye, T=Triticale, W=Wheat.
}
\usage{data(cereal)}
\format{
  A data frame with 15 objects and 3 list elements:
  \describe{
    \item{\code{X}}{matrix with 15 rows and 145 columns}
    \item{\code{Y}}{matrix with 15 rows and 6 columns}
    \item{\code{Ysc}}{matrix with 15 rows and 6 columns}
  }
}
\details{
The data set can be used for PLS2.
}
\source{
K. Varmuza and P. Filzmoser: Introduction to Multivariate Statistical Analysis
in Chemeometrics. CRC Press. To appear. }
\references{
K. Varmuza and P. Filzmoser: Introduction to Multivariate Statistical Analysis
in Chemeometrics. CRC Press. To appear. }
\examples{
data(cereal)
names(cereal)
}
\keyword{datasets}
