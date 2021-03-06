% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/function.R
\name{SPH}
\alias{SPH}
\title{Random Generation for General Spherically Symmetric Distributions}
\usage{
SPH(n, p, D, ...)
}
\arguments{
\item{n}{number of rows (observations).}

\item{p}{number of columns (variables).}

\item{D}{random generation functions for some distributions (e.g., \code{rgamma}, \code{rbeta}).}

\item{...}{optional arguments passed to \code{D}.}
}
\value{
Returns univariate (\eqn{p=1}) or multivariate (\eqn{p>1}) random sample matrix.
}
\description{
Generate univariate or multivariate random sample for general spherically symmetric distributions.
}
\examples{
set.seed(12345)

## Generate 5X2 random sample matrix from SPH(Beta(1,1)) ##
SPH(n=5, p=2, D=rbeta, shape1=1, shape2=1)


## Power calculation against bivariate (p=2) SPH(Beta(1,1)) distribution ##
## at sample size n=50 at one-sided alpha = 0.05 ##

# Zhou-Shao's test #
power.mvnTest(a=0.05, n=50, p=2, B=100, FUN=SPH, D=rbeta, shape1=1, shape2=1)

}
\references{
Chmielewski, M. A. (1981). Elliptically symmetric distributions: A review and bibliography. \emph{International Statistical Review/Revue Internationale de Statistique}, 67-74.

Henze, N., & Zirkler, B. (1990). A class of invariant consistent tests for multivariate normality. \emph{Communications in statistics-Theory and Methods}, 19(10), 3595-3617.
}
