% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/function.R
\name{msk}
\alias{msk}
\title{Bowman and Shenton Test for Multivariate Normality}
\usage{
msk(X, B = 1000)
}
\arguments{
\item{X}{an \eqn{n*p} numeric matrix or data frame.}

\item{B}{number of Monte Carlo simulations for null distribution, default is 1000 (increase B to increase the precision of p-value).}
}
\value{
Returns a list with two objects:
\describe{
\item{\code{mv.test}}{results of the Bowman and Shenton test, i.e., test statistic, \emph{p}-value, and multivariate normality summary (YES, if \emph{p}-value>0.05).}
\item{\code{uv.shapiro}}{a dataframe with \eqn{p} rows detailing univariate Shapiro-Wilk tests. Columns in the dataframe contain test statistics \emph{W}, \emph{p}-value,and univariate normality summary (YES, if \emph{p}-value>0.05).}
}
}
\description{
It computes Bowman and Shenton (1975)'s test statistic (MSK) and its corresponding
 p-value for multivariate normality. The statistic is calculated based on a combination of
 multivariate skewness (MS) and kurtosis (MK) such that \eqn{MSK=MS+|MK|^2}. For formulas of MS and MK,
 please refer to Mardia (1970). The corresponding p-value of the statistic is computed based on a
 simulated null distribution of MSK. The skewness statistic (MS) will be adjusted for sample size \eqn{n < 20}.
}
\examples{
set.seed(12345)

## Data from gamma distribution
X = matrix(rgamma(50*4,shape =  2),50)
msk(X, B=100)

## load the ubiquitous multivariate iris data ##
## (first 50 observations of columns 1:4) ##
iris.df = iris[1:50, 1:4]
msk(iris.df, B=100)

}
\references{
Bowman, K. O., & Shenton, L. R. (1975). Omnibus test contours for departures from normality based on \eqn{\sqrt b_1} and \eqn{b_2}. \emph{Biometrika}, 62(2), 243-250.

Shapiro, S. S., & Wilk, M. B. (1965). An analysis of variance test for normality (complete samples). \emph{Biometrika}, 52(3/4), 591-611.

Mardia, K. V. (1970). Measures of multivariate skewness and kurtosis with applications. \emph{Biometrika}, 57(3), 519-530.

Doornik, J. A., & Hansen, H. (2008). An omnibus test for univariate and multivariate normality. \emph{Oxford Bulletin of Economics and Statistics}, 70, 927-939.

Zhou, M., & Shao, Y. (2014). A powerful test for multivariate normality. \emph{Journal of applied statistics}, 41(2), 351-363.
}
\seealso{
\code{\link{power.msk}}, \code{\link{mvnTest}}, \code{\link{faTest}}, \code{\link{msw}}, \code{\link{mardia}}, \code{\link{mhz}}, \code{\link{mvn}}
}
