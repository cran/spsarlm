% Copyright 2002 by Roger S. Bivand
\name{lagsarlm}
\alias{lagsarlm}
\alias{sar.lag.mixed.f}
%\alias{sar.lag.mixed.f.s}
%\alias{dosparse}
\title{Spatial simultaneous autoregressive lag model estimation}
\description{
  Maximum likelihood estimation of spatial simultaneous autoregressive
lag and mixed models of the form:

\deqn{y = \rho W y + X \beta + \varepsilon}

where $\rho$ is found by \code{optimize()} first and $\beta$ and other
parameters by generalized least squares subsequently. In the mixed model,
the spatially lagged independent variables are added to X.
}
\usage{
lagsarlm(formula, data=list(), listw, type="lag", method="eigen", quiet=T,
  zero.policy=FALSE)
%sar.lag.mixed.f.s(rho, sn, e.a, e.b, e.c, n, quiet)
sar.lag.mixed.f(rho, eig, e.a, e.b, e.c, n, quiet)
%dosparse(listw, y, x, wy, K, quiet)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{formula}{a symbolic description of the model to be fit. The details 
of model specification are given for \code{lm()}}
  \item{data}{an optional data frame containing the variables in the model. 
By default the variables are taken from the environment which the function 
is called.}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{type}{default "lag", may be set to "mixed"}
  \item{method}{"eigen" (default) - the Jacobian is computed as the product 
of (1 - rho*eigenvalue) using \code{eigenw} %, and "sparse" - computes the 
%determinant of the sparse matrix (I - rho*W) directly using \code{log.spwdet}.
}
  \item{quiet}{default=T; if F, reports function values during optimization.}
  \item{zero.policy}{if TRUE assign zero to the lagged value of zones without 
neighbours, if FALSE (default) assign NA - causing \code{lagsarlm()} to terminate with an error}
  \item{rho}{value of the spatial parameter}
  \item{eig}{eigenvalues of the full spatial weights matrix from \code{eigenw}}
  \item{y}{dependent variable}
  \item{wy}{spatially lagged dependent variable}
  \item{x}{independent variables}
  \item{n}{length of y (and eig)}
  \item{e.a}{term used in computing likelihood}
  \item{e.b}{term used in computing likelihood}
  \item{e.c}{term used in computing likelihood}
  \item{K}{1 if no intercept, 2 if intercept present in x}
%  \item{sn}{sparse spatial neighbour object from \code{listw2sn}}
}
\details{
  
}
\value{
  A list object of class \code{sarlm}
  \item{type}{"lag" or "mixed"}
  \item{rho}{simultaneous autoregressive lag coefficient}
  \item{coefficients}{GLS coefficient estimates}
  \item{rest.se}{asymptotic standard errors if ase=TRUE}
  \item{LL}{log likelihood value at computed optimum}
  \item{s2}{GLS residual variance}
  \item{SSE}{sum of squared GLS errors}
  \item{parameters}{number of parameters estimated}
  \item{lm.model}{the \code{lm} object returned when estimating for $\rho=0$}
  \item{method}{the method used to calculate the Jacobian}
  \item{call}{the call used to create this object}
  \item{residuals}{GLS residuals}
  \item{lm.target}{the \code{lm} object returned for the GLS fit}
  \item{fitted.values}{GLS fitted values}
  \item{se.fit}{The GLS standard errors of the fitted values (not taking into
account the standard error of $\lambda$)}
  \item{ase}{TRUE if method=eigen}
  \item{LLs}{if ase=FALSE (for method="sparse"), the log likelihood values of
models estimated dropping each of the independent variables in turn, used
in the summary function as a substitute for variable coefficient
significance tests}
  \item{rho.se}{if ase=TRUE, the asymptotic standard error of $\rho$}
  \item{LMtest}{if ase=TRUE, the Lagrange Multiplier test for the absence
of spatial autocorrelation in the lag model residuals}
  \item{zero.policy}{zero.policy for this model}

The sar.lag.mixed.* functions return the value of the log likelihood function
at $\rho$.

}
\references{Cliff, A. D., Ord, J. K. 1981 \emph{Spatial processes}, Pion;
Ord, J. K. 1975 Estimation methods for models of spatial interaction,
\emph{Journal of the American Statistical Association}, 70, 120-126;
Anselin, L. 1988 \emph{Spatial econometrics: methods and models.}
(Dordrecht: Kluwer); Anselin, L. 1995 SpaceStat, a software program for
the analysis of spatial data, version 1.80. Regional Research Institute,
West Virginia University, Morgantown, WV (\url{www.spacestat.com});
Anselin L, Bera AK (1998) Spatial dependence in linear regression models
with an introduction to spatial econometrics. In: Ullah A, Giles DEA
(eds) Handbook of applied economic statistics. Marcel Dekker, New York,
pp. 237-289.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}, with thanks to Andrew 
Bernat for contributions to the asymptotic standard error code.}

\seealso{\code{\link{lm}}, \code{\link{errorsarlm}}, 
\code{\link{eigenw}}, %\code{\link{log.spwdet}}
}

\examples{
require(spweights)
require(sptests)
data(oldcol)
COL.lag.eig <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb), method="eigen", quiet=FALSE)
%COL.lag.sp <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
% nb2listw(COL.nb), method="sparse", quiet=FALSE)
summary(COL.lag.eig)
%summary(COL.lag.sp)
}
\keyword{spatial}
