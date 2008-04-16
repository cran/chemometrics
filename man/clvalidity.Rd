\name{clvalidity}
\alias{clvalidity}
\title{compute and plot cluster validity}
\description{
A cluster validity measure based on within- and between-sum-of-squares
is computed and plotted for the methods k-means, fuzzy c-means, and
model-based clustering.
}
\usage{
clvalidity(x, clnumb = c(2:10))
}
\arguments{
  \item{x}{input data matrix}
  \item{clnumb}{range for the desired number of clusters}
}
\details{
The validity measure for a number $k$ of clusters is 
$\sum_j W_j$ divided by $\sum_{j<l} B_{jl}$ with
$W_j$ is the sum of squared distances of the objects in each cluster 
cluster to its center, and $B_{jl}$ is the squared distance between
the cluster centers of cluster $j$ and $l$.
}
\value{
  \item{validity}{vector with validity measure for the desired numbers of clusters}
}
\references{ K. Varmuza and P. Filzmoser: Introduction to Multivariate Statistical Analysis
in Chemometrics. CRC Press. To appear. }
\author{ Peter Filzmoser <P.Filzmoser@tuwien.ac.at> }
\seealso{\code{\link{princomp}}}
\examples{
data(glass)
require(robustbase)
res <- pcaCV(glass,segments=4,repl=100,cex.lab=1.2,ylim=c(0,1),las=1)
}
\keyword{ multivariate }

