\name{hyptis}
\alias{hyptis}
\docType{data}
\title{Hyptis data set}
\description{
30 objects (Wild growing, flowering Hyptis suaveolens) and 7 variables
(chemotypes), and 2 variables that explain the grouping (4 groups).
}
\usage{data(hyptis)}
\format{
  A data frame with 30 observations on the following 9 variables.
  \describe{
    \item{\code{Sabinene}}{a numeric vector}
    \item{\code{Pinene}}{a numeric vector}
    \item{\code{Cineole}}{a numeric vector}
    \item{\code{Terpinene}}{a numeric vector}
    \item{\code{Fenchone}}{a numeric vector}
    \item{\code{Terpinolene}}{a numeric vector}
    \item{\code{Fenchol}}{a numeric vector}
    \item{\code{Location}}{a factor with levels \code{East-high} \code{East-low} \code{North} \code{South}}
    \item{\code{Group}}{a numeric vector with the group information}
  }
}
\details{
This data set can be used for cluster analysis.
}
\source{
  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
P. Grassi, M.J. Nunez, K. Varmuza, and C. Franz:
Chemical polymorphism of essential oils of Hyptis suaveolens from El Salvador.
Flavour \& Fragrance, 20, 131-135, 2005.
K. Varmuza and P. Filzmoser: Introduction to Multivariate Statistical Analysis
in Chemeometrics. CRC Press. To appear. }
}
\examples{
data(hyptis)
str(hyptis)
}
\keyword{datasets}
