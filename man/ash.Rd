\name{ash}
\alias{ash}
\docType{data}
\title{ ash data }
\description{
Data from 99 ash samples originating from different biomass, measured
on 9 variables; 8 log-transformed variables are added.
}
\usage{data(ash)}
\format{
  A data frame with 99 observations on the following 17 variables.
  \describe{
    \item{\code{SOT}}{a numeric vector}
    \item{\code{P2O5}}{a numeric vector}
    \item{\code{SiO2}}{a numeric vector}
    \item{\code{Fe2O3}}{a numeric vector}
    \item{\code{Al2O3}}{a numeric vector}
    \item{\code{CaO}}{a numeric vector}
    \item{\code{MgO}}{a numeric vector}
    \item{\code{Na2O}}{a numeric vector}
    \item{\code{K2O}}{a numeric vector}
    \item{\code{log(P2O5)}}{a numeric vector}
    \item{\code{log(SiO2)}}{a numeric vector}
    \item{\code{log(Fe2O3)}}{a numeric vector}
    \item{\code{log(Al2O3)}}{a numeric vector}
    \item{\code{log(CaO)}}{a numeric vector}
    \item{\code{log(MgO)}}{a numeric vector}
    \item{\code{log(Na2O)}}{a numeric vector}
    \item{\code{log(K2O)}}{a numeric vector}
  }
}
\details{
The dependent variable Softening Temperature (SOT) of ash should be 
modeled by the elemental composition of the ash data. Data from 99
ash samples - originating from different biomass - comprise the
experimental SOT (630-1410 centigrades), and the experimentally
determined eight mass concentrations the listed elements. Since the
distribution of the elements is skweed, the log-transformed variables
have been added. 
}
\source{
K. Varmuza and P. Filzmoser: Introduction to Multivariate Statistical Analysis
in Chemometrics. CRC Press, Boca Raton, FL, 2009. }
\references{
K. Varmuza and P. Filzmoser: Introduction to Multivariate Statistical Analysis
in Chemometrics. CRC Press, Boca Raton, FL, 2009. }
\examples{
data(ash)
str(ash)
}
\keyword{datasets}
