\name{plotsom}
\alias{plotsom}
\title{ Plot SOM results }
\description{
Plot results of Self Organizing Maps (SOM).
}
\usage{
plotsom(obj, grp, type = c("num", "bar"), margins = c(3,2,2,2), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{ result object from \code{\link[som]{som}} }
  \item{grp}{ numeric vector or factor with group information }
  \item{type}{ type of presentation for output, see details }
  \item{margins}{ plot margins for output, see \code{\link{par}} }
  \item{\dots}{ additional graphics parameters, see \code{\link{par}} }
}
\details{
The results of Self Organizing Maps (SOM) are plotted either in a table
with numbers (type="num") or with barplots (type="bar"). 
There is a limitation to at most 9 groups. A summary table is returned.
}
\value{
  \item{sumtab }{Summary table}
}
\references{ K. Varmuza and P. Filzmoser: Introduction to Multivariate Statistical Analysis
in Chemometrics. CRC Press. To appear. }
\author{ Peter Filzmoser <P.Filzmoser@tuwien.ac.at> }
\seealso{\code{\link[som]{som}}}
\examples{
data(glass)
require(som)
Xs <- scale(glass)
Xn <- Xs/sqrt(apply(Xs^2,1,sum))
X_SOM <- som(Xn,xdim=4,ydim=4) # 4x4 fields
data(glass.grp)
res <- plotsom(X_SOM,glass.grp,type="bar")
}
\keyword{ multivariate }
\keyword{ robust }
