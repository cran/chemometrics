\name{glass}
\alias{glass}
\docType{data}
\title{ glass vessels data }
\description{
13 different measurements for 180 archaeological glass vessels from
different groups are included.
}
\usage{data(glass)}
\format{
  The format is:
 num [1:180, 1:13] 13.9 14.2 14.7 14.8 14.1 ...
 - attr(*, "dimnames")=List of 2
  ..$ : NULL
  ..$ : chr [1:13] "Na2O" "MgO" "Al2O3" "SiO2" ...
}
\details{
This is a matrix with 180 objects and 13 columns.
}
\source{
Janssen, K.H.A., De Raedt, I., Schalm, O., Veeckman, J.: Microchim. Acta 15 (suppl.) 
(1998) 253-267. Compositions of 15th - 17th century archaeological glass vessels 
excavated in Antwerp.
}
\references{
K. Varmuza and P. Filzmoser: Introduction to Multivariate Statistical Analysis
in Chemeometrics. CRC Press. To appear. }
}
\examples{
data(glass)
str(glass)
}
\keyword{datasets}
