% Copyright 2002 by Roger S. Bivand
\name{summary.sarlm}
\alias{summary.sarlm}
\alias{print.sarlm}
\alias{print.summary.sarlm}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{summary method for class sarlm}
\description{
  Methods used for presenting the results of estimating spatial SAR models.
}
\usage{
summary.sarlm(object, ...)
print.sarlm(x, ...)
print.summary.sarlm(x, digits = max(5, .Options$digits - 3),
	signif.stars = FALSE, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{object}{\code{sarlm} object from \code{lagsarlm} or \code{errorsarlm}}
  \item{x}{\code{sarlm} object from \code{lagsarlm} or \code{errorsarlm} in
\code{print.sarlm}, summary object from \code{summary.sarlm} for 
\code{print.summary.sarlm}}
  \item{digits}{the number of significant digits to use when printing}
  \item{signif.stars}{logical. If TRUE, "significance stars" are printed
          for each coefficient.}
  \item{\dots}{further arguments passed to or from other methods}
}
\value{
The summary function \code{summary.sarlm} returns the \code{sarlm}
object augmented with a coefficient matrix with probability values
for coefficient asymptotic standard errors for type="error" and for
type="lag" or "mixed" when object\$ase=TRUE, or a coefficient matrix
with probability values for likelihood ratio tests between the model as
reported and models with independent variables dropped in turn.
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
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{errorsarlm}}, \code{\link{lagsarlm}},
\code{\link{summary.lm}}}

\examples{
require(spweights)
require(sptests)
data(oldcol)
COL.mix.eig <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb), type="mixed", method="eigen", quiet=FALSE)
summary(COL.mix.eig)
}
\keyword{spatial}
