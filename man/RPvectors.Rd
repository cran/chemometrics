\name{RPvectors}
\alias{RPvectors}
\title{ Generating random projection directions}
\description{
A matrix with pandom projection (RP) directions (columns) is generated
according to a chosen distributions; optionally the random vectors are
orthogonalized.
}
\usage{
RPvectors(a, m, ortho = "none", distr = "uniform", par_unif = c(-1, 1), 
par_norm = c(0, 1), par_eq = c(-1, 0, 1), par_uneq = c(-sqrt(3), 0, sqrt(3)), 
par_uneqprob = c(1/6, 2/3, 1/6))
}
\arguments{
  \item{a}{
number of generated vectors (>=1)
}
  \item{m}{
dimension of generated vectors (>=2)
}
  \item{ortho}{
orthogonalization of vectors:
"none" ... no orthogonalization (default);
"onfly" ... orthogonalization on the fly after each generated vector;
"end" ... orthogonalization at the end, after the whole random matrix was generated
}
  \item{distr}{
distribution of generated random vector components:
"uniform" ... uniformly distributed in range par_unif (see below);
                   default U[-1, +1];
"normal" ... normally distributed with parameters par_norm (see below);
                   typical N(0, 1);
"randeq" ... random selection of values par_eq (see below) with equal probabilities;
                   typically -1, 0, +1;
"randuneq" ... random selection of values par_uneq (see below) with
probabilties par_uneqprob (see below);
                   typical -(3)^0.5 with probability 1/6;
                             0      with probability 2/3;
                           +(3)^0.5 with probability 1/6
}
  \item{par_unif}{
parameters for range for distr=="uniform"; default to c(-1,1)
}
  \item{par_norm}{
parameters for mean and sdev for distr=="normal"; default to c(0,1)
}
  \item{par_eq}{
values for distr=="randeq" which are replicated; default to c(-1,0,1)
}
  \item{par_uneq}{
values for distr=="randuneq" which are replicated with probabilties par_uneqprob;
default to c(-sqrt(3),0,sqrt(3))
}
  \item{par_uneqprob}{
probabilities for distr=="randuneq" to replicate values par_uneq;
default to c(1/6,2/3,1/6)
}
}
\details{
The generated random projections can be used for dimension reduction of multivariate
data. Suppose we have a data matrix X with n rows and m columns. Then the call
B <- RPvectors(a,m) will produce a matrix B with the random directions in its
columns. The matrix product X times t(B) results in a matrix of lower dimension a.
There are several options to generate the projection directions, like orthogonal
directions, and different distributions with different parameters to generate 
the random numbers.
Random Projection (RP) can have comparable performance for dimension reduction
like PCA, but gives a big advantage in terms of computation time. 
}
\value{
The value returned is the matrix B with a columns of length m, representing
the random vectors
}
\references{
 K. Varmuza, P. Filzmoser, and B. Liebmann.
 Random projection experiments with chemometric data.
 Journal of Chemometrics. To appear.
}
\author{
 Peter Filzmoser <P.Filzmoser@tuwien.ac.at>
}

\examples{
B <- RPvectors(a=5,m=10)
res <- t(B)%*%B # the result is an orthonormal random projection matrix
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{  multivariate }
