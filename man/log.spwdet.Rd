% Copyright 2002 by Roger S. Bivand
\name{log.spwdet}
\alias{log.spwdet}
\alias{spwdet}
\title{Sparse computation of Jacobian for spatial SAR estimation}
\description{
  This function calls compiled code for computing the Jacobian of the
likelihood of spatial SAR models:

\deqn{\log(\det[I - \rho W])}

where $W$ is the n by n spatial weights matrix. For smaller n, \code{eigenw()}
provides a practical solution, but for larger n, finding the eigenvalues of
$W$ can be timeconsuming. Since W is in general sparse (but not for weights
based on inverse distances), and become sparser with increasing n, sparse
matrix techniques may be used to find the determinant (Pace and Barry 1997).
The underlying code used here is the same as that used by the SPLUS Spatial 
Statistics module, and was written by Kundert and Sangiovanni-Vincentelli 
(1988), and is available from netlib.
}
\usage{
log.spwdet(sparseweights, rho)
spwdet(sparseweights, rho)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{sparseweights}{a \code{spatial neighbour} object representation 
typically from \code{listw2sn}}
  \item{rho}{the spatial parameter whose value is being estimated}
}

\value{
  \code{log.spwdet} return the log of the determinant of I - rho*W
}
\references{Pace, R. K. and Barry, R. 1997 Quick computation of spatial 
autoregressive estimators. \emph{Geographical Analysis}, 29, 232--247; 
Kaluzny, S. P., Vega, S. C., Cardoso, T. P. and Shelly, A. A. 1996 
\emph{S+SPATIALSTATS users manual version 1.0}; Kundert, K. and 
Sangiovanni-Vincentelli, A. 1988 Sparse user's guide, version 1.3a. 
Department of Electrical Engineering and Computer Science, University 
of California, Berkeley.
}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{eigenw}}}

\examples{
library(spweights)
library(sptests)
small.nb <- cell2nb(7, 7)
rho <- seq(0.1, 0.8, 0.02)
res <- numeric(length(rho))
system.time(sn <- listw2sn(nb2listw(small.nb)))
system.time(for(i in 1:length(rho)) res[i] <- log.spwdet(sn, rho[i]))
cbind(rho, res)
system.time(eig <- eigenw(nb2listw(small.nb)))
system.time(for(i in 1:length(rho)) res[i] <- log(prod(1 - rho[i]*eig)))
cbind(rho, res)
system.time({I <- diag(length(small.nb)); W <- nb2mat(small.nb)})
system.time(for(i in 1:length(rho)) res[i] <- log(det(I - rho[i]*W)))
cbind(rho, res)
larger.nb <- cell2nb(15, 15)
system.time(sn <- listw2sn(nb2listw(larger.nb)))
system.time(for(i in 1:length(rho)) res[i] <- log.spwdet(sn, rho[i]))
cbind(rho, res)
system.time(eig <- eigenw(nb2listw(larger.nb)))
system.time(for(i in 1:length(rho)) res[i] <- log(prod(1 - rho[i]*eig)))
cbind(rho, res)
system.time({I <- diag(length(larger.nb)); W <- nb2mat(larger.nb)})
system.time(for(i in 1:length(rho)) res[i] <- log(det(I - rho[i]*W)))
cbind(rho, res)
}
\keyword{spatial}
