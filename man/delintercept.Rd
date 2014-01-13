\name{delintercept}
\alias{delintercept}
\title{Delete intercept from model matrix}
\description{
A utility function to delete any intercept column from a model
matrix, and adjust the assign attribute correspondingly. 
}
\usage{
delintercept(mm)
}
\arguments{
  \item{mm}{Model matrix}
}
\value{
  A model matrix without intercept column.
}
\author{ B.-H. Mevik and Ron Wehrens }
\seealso{\code{\link[pls]{delete.intercept}}}
\keyword{ multivariate }
