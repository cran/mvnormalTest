% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/function.R
\name{power.usw}
\alias{power.usw}
\title{Power Calculation using the Univariate Shapiro-Wilk Test Statistic}
\usage{
power.usw(a, n, p = 1, B = 1000, FUN, ...)
}
\arguments{
\item{a}{significance level (\eqn{\alpha}).}

\item{n}{number of rows (observations).}

\item{p}{p=1 for univariate.}

\item{B}{number of Monte Carlo simulations, default is 1000 (can increase B to increase the precision).}

\item{FUN}{self-defined function for generate multivariate distribution. See example.}

\item{...}{optional arguments passed to \code{FUN}.}
}
\value{
Returns a numeric value of the estimated empirical power (value between 0 and 1).
}
\description{
Empirical power calculation using univariate Shapiro-Wilk test statistic.
}
\examples{
set.seed(12345)

## Power calculation against univariate (p=1) independent Beta(1, 1) distribution ##
## at sample size n=50 at one-sided alpha = 0.05 ##

power.usw(a = 0.05, n = 50, p = 1,  B = 100, FUN=IMMV, D1=runif)

}
\references{
Shapiro, S. S., & Wilk, M. B. (1965). An analysis of variance test for normality (complete samples). \emph{Biometrika}, 52(3/4), 591-611.
}
