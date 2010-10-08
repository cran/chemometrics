\name{NIR}
\alias{NIR}
\docType{data}
\title{
NIR data
}
\description{
For 166 alcoholic fermentation mashes of different feedstock (rye, wheat and corn)
we have 235 variables (X) containing the first derivatives of near infrared spectroscopy
(NIR) absorbance values at 1115-2285 nm, and two variables (Y) containing the 
concentration of glucose and ethanol (in g/L).
}
\usage{data(NIR)}
\format{
A data frame with 166 objects and 2 list elements:
  \describe{
    \item{\code{xNIR}}{data frame with 166 rows and 235 columns}
    \item{\code{yGlcEtOH}}{data frame with 166 rows and 2 columns}
  }
}
\details{
The data can be used for linear and non-linear models.
}
\source{
B. Liebmann, A. Friedl, and K. Varmuza. Determination of glucose and ethanol in
bioethanol production by near infrared spectroscopy and chemometrics.
Anal. Chim. Acta, 642:171-178, 2009.
}
\references{
B. Liebmann, A. Friedl, and K. Varmuza. Determination of glucose and ethanol in
bioethanol production by near infrared spectroscopy and chemometrics.
Anal. Chim. Acta, 642:171-178, 2009.
}
\examples{
data(NIR)
str(NIR)
}
\keyword{datasets}
