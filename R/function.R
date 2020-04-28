#' A Powerful Test for Multivariate Normality (Zhou-Shao's Test)
#'
#' @description A simple and powerful test for multivariate normality with a combination of multivariate
#' kurtosis (MK) and Shapiro-Wilk which was proposed by Zhou and Shao (2014). The \emph{p}-value of the test
#' statistic (\eqn{T_n}) is computed based on a simulated null distribution of \eqn{T_n}. Details see Zhou and Shao (2014).
#' @param X an \eqn{n*p} data matrix or data frame, where \eqn{n} is number of rows (observations) and \eqn{p} is number of columns (variables) and \eqn{n>p}.
#' @param B number of Monte Carlo simulations for null distribution, default is 1000 (increase B to increase the precision of \emph{p}-value).
#' @param pct percentiles of MK to get \eqn{c_1} and \eqn{c_2} described in the reference paper, default is (0.01, 0.99).
#' @return Returns a list with two objects:
#' \describe{
#' \item{\code{mv.test}}{results of the Zhou-Shao's test for multivariate normality , i.e., test statistic \eqn{T_n},
#' \emph{p}-value (under H0, i.e. multivariate normal, that \eqn{T_n} is at least as extreme as the observed value), and multivariate normality summary (YES, if \emph{p}-value>0.05).}
#' \item{\code{uv.shapiro}}{a dataframe with \eqn{p} rows detailing univariate Shapiro-Wilk tests. Columns in the dataframe contain test statistics \emph{W}, \emph{p}-value,and univariate normality summary (YES, if \emph{p}-value>0.05).}
#' }
#' @references Zhou, M., & Shao, Y. (2014). A powerful test for multivariate normality. \emph{Journal of applied statistics}, 41(2), 351-363.
#' @references Shapiro, S. S., & Wilk, M. B. (1965). An analysis of variance test for normality (complete samples). \emph{Biometrika}, 52(3/4), 591-611.
#' @import stats
#' @seealso \code{\link{power.mvnTest}}, \code{\link{msk}}, \code{\link{mardia}}, \code{\link{msw}}, \code{\link{faTest}}, \code{\link{mhz}}
#' @export
#' @examples
#' set.seed(12345)
#'
#' ## Data from gamma distribution ##
#' X = matrix(rgamma(50*4,shape =  2),50)
#' mvnTest(X, B=100)
#'
#' ## load the ubiquitous multivariate iris data ##
#' ## (first 50 observations of columns 1:4) ##
#' iris.df = iris[1:50, 1:4]
#' mvnTest(iris.df, B=100)
#'
mvnTest <- function(X, B=1000, pct=c(0.01,0.99)) {
  X <- as.matrix(X)
  X <- X[complete.cases(X), ]
  X <- as.matrix(X)
  n  = NROW(X)
  p  = NCOL(X)
  mk = Dsn_MK(n, p, B)

  ##obtain c1 and c2 from the null distribution of MK
  cri_mk = quantile(mk, pct)

  ## Component of the Tn test statistic
  SWM2 <- function(X){
    SWS <- function(x) shapiro.test(x)$statistic
    X <- as.matrix(X)
    p = NCOL(X)
    StdX = STD(X)
    return((mean(apply(StdX, 2, SWS)) + mean(sort(apply(StdX%*%t(StdX), 2, SWS))[1:p]))/2)
  }

  ## Tn statistic function
  SWM2MK <- function(X, cri_mk) {
    Mk      <- MK(X)
    Index   <- (Mk > cri_mk[2]|Mk < cri_mk[1])
    return( (1 - Index)*(1 - SWM2(X)) + Index)
  }

  ## Simulate the null distribution of Tn
  Dsn_SWM2MK <- function(n, p, cri_mk, B) {
    s <- rep(0,B)
    for(i in 1:B) {
      s[i] <- SWM2MK(matrix(rnorm(n*p), n), cri_mk)
    }
    return(s)
  }

  ## obtain Tn statistic
  Tn = SWM2MK(X, cri_mk)

  ##obtain the null distribution sample from Tn
  ##(note that cri_mk only depends on the null distribution of MK (i.e., data independent), and can be calculated in advance
  swm2mk = Dsn_SWM2MK(n, p, cri_mk, B)

  ##calculate the p-value, i.e., the chance (under H0, i.e. multivariate normal) that Tn is at least as extreme as the observed value
  p.value <- mean(swm2mk > Tn)
  result <- ifelse(p.value>0.05,"YES","NO")
  output <- noquote(c(round(Tn,4),round(p.value,4),result))
  names(output) <- c("Tn","p-value","Result")

  final <- list(output,univ(X))
  names(final) <- c("mv.test","uv.shapiro")

  return(final)
}



#' Rotational Robust Shapiro-Wilk Type (SWT) Test for Multivariate Normality (FA Test of Fattorini)
#'
#' @description It computes FA Test proposed by Fattorini (1986). This test would be more rotationally robust than other
#' SWT tests such as Royston (1982) H test and the test proposed by Villasenor-Alva and Gonzalez-Estrada (2009).
#' The \emph{p}-value of the test statistic is computed based on a
#' simulated null distribution of the statistic.
#' @param X an \eqn{n*p} data matrix or data frame, where \eqn{n} is number of rows (observations) and \eqn{p} is number of columns (variables) and \eqn{n>p}.
#' @param B number of Monte Carlo simulations for null distribution, default is 1000 (increase B to increase the precision of \emph{p}-value).
#' @return Returns a list with two objects:
#' \describe{
#' \item{\code{mv.test}}{results of the FA test for multivariate normality, i.e., test statistic,
#' \emph{p}-value, and multivariate normality summary (YES, if \emph{p}-value>0.05).}
#' \item{\code{uv.shapiro}}{a dataframe with \eqn{p} rows detailing univariate Shapiro-Wilk tests. Columns in the dataframe contain test statistics \emph{W}, \emph{p}-value,and univariate normality summary (YES, if \emph{p}-value>0.05).}
#' }
#' @references Fattorini, L. (1986). Remarks on the use of Shapiro-Wilk statistic for testing multivariate normality. \emph{Statistica}, 46(2), 209-217.
#' @references Lee, R., Qian, M., & Shao, Y. (2014). On rotational robustness of Shapiro-Wilk type tests for multivariate normality. \emph{Open Journal of Statistics}, 4(11), 964.
#' @references Shapiro, S. S., & Wilk, M. B. (1965). An analysis of variance test for normality (complete samples). \emph{Biometrika}, 52(3/4), 591-611.
#' @references Royston, J. P. (1982). An extension of Shapiro and Wilk's W test for normality to large samples. \emph{Journal of the Royal Statistical Society: Series C (Applied Statistics)}, 31(2), 115-124.
#' @references Villasenor Alva, J. A., & Estrada, E. G. (2009). A generalization of Shapiro–Wilk's test for multivariate normality. \emph{Communications in Statistics—Theory and Methods}, 38(11), 1870-1883.
#' @references Zhou, M., & Shao, Y. (2014). A powerful test for multivariate normality. \emph{Journal of applied statistics}, 41(2), 351-363.
#' @import stats
#' @seealso \code{\link{power.faTest}}, \code{\link{mvnTest}}, \code{\link{msk}}, \code{\link{mardia}}, \code{\link{msw}}, \code{\link{mhz}}
#' @export
#' @examples
#' set.seed(12345)
#'
#' ## Data from gamma distribution ##
#' X = matrix(rgamma(50*4,shape =  2),50)
#' faTest(X, B=100)
#'
#' ## load the ubiquitous multivariate iris data ##
#' ## (first 50 observations of columns 1:4) ##
#' iris.df = iris[1:50, 1:4]
#' faTest(iris.df, B=100)
#'
faTest <- function(X, B=1000) {

  X <- as.matrix(X)
  X <- X[complete.cases(X), ]
  X <- as.matrix(X)
  n  = NROW(X)
  p  = NCOL(X)

  SW <- function(x){shapiro.test(x)$statistic}

  FA <- function(X) {
    X <- as.matrix(X)
    n <- NROW(X)
    p <- NCOL(X)
    mu <- apply(X,2,mean)
    nSinver <- solve((n-1)*cov(X))
    Y <- X%*%t((X-matrix(rep(mu,n),ncol=p,byrow=T))%*%nSinver)
    return(1 - min(apply(Y,2,SW)))
  }

  Dsn_FA <- function(n, p, B){
    s <- rep(0, B)

    for(i in 1:B) {
      s[i] <- FA(matrix(rnorm(n*p), n))
    }
    return(s)

  }

  fa = FA(X)
  Dsn_fa = Dsn_FA(n, p, B)

  p.value <- mean(Dsn_fa > fa)
  result <- ifelse(p.value>0.05,"YES","NO")
  output <- noquote(c(round(fa,4),round(p.value,4),result))
  names(output) <- c("Statistic","p-value","Result")

  final <- list(output,univ(X))
  names(final) <- c("mv.test","uv.shapiro")

  return(final)
}



#' Shapiro-Wilk Type (SWT) Tests for Multivariate Normality
#'
#' @description The SWT-based tests for multivariate normality including Royston's H test and the test proposed
#' by Villasenor-Alva and Gonzalez-Estrada (2009).
#' @param X an \eqn{n*p} numeric matrix or data frame, the number of \eqn{n} must be between 3 and 5000, n>p.
#' @return Returns a list with two objects:
#' \describe{
#' \item{\code{mv.test}}{a result table of multivariate normality tests,
#' including the name of the test, test statistic,
#' \emph{p}-value, and multivariate normality summary (Yes, if \emph{p}-value>0.05). Note that the test results
#' of \code{Royston} will not be reported if \eqn{n > 2000} or \eqn{n < 3} and the test results
#' of Villasenor-Alva and Gonzalez-Estrada (\code{VAGE}) will not be reported if  \eqn{n > 5000} or \eqn{n < 12}.}
#' \item{\code{uv.shapiro}}{a dataframe with \eqn{p} rows detailing univariate Shapiro-Wilk tests.
#' Columns in the dataframe contain test statistics \emph{W}, \emph{p}-value,and univariate normality summary (YES, if \emph{p}-value>0.05).}
#' If the number of variable is \eqn{p=1}, only univariate Shapiro-wilk's test result will be produced.}
#' @references Shapiro, S. S., & Wilk, M. B. (1965). An analysis of variance test for normality (complete samples). \emph{Biometrika}, 52(3/4), 591-611.
#' @references Royston, J. P. (1982). An extension of Shapiro and Wilk's W test for normality to large samples. \emph{Journal of the Royal Statistical Society: Series C (Applied Statistics)}, 31(2), 115-124.
#' @references Villasenor Alva, J. A., & Estrada, E. G. (2009). A generalization of Shapiro–Wilk's test for multivariate normality. \emph{Communications in Statistics—Theory and Methods}, 38(11), 1870-1883.
#' @references Lee, R., Qian, M., & Shao, Y. (2014). On rotational robustness of Shapiro-Wilk type tests for multivariate normality. \emph{Open Journal of Statistics}, 4(11), 964.
#' @import stats nortest moments
#' @seealso \code{\link{power.mswR}}, \code{\link{power.mswV}}, \code{\link{mvnTest}}, \code{\link{faTest}}, \code{\link{msk}}, \code{\link{mardia}}, \code{\link{mhz}}, \code{\link{mvn}}, \code{\link{shapiro.test}}
#' @export
#' @examples
#' set.seed(12345)
#'
#' ## Data from gamma distribution
#' X = matrix(rgamma(50*4,shape =  2),50)
#' msw(X)
#'
#' ## Data from normal distribution
#' X = matrix(rnorm(50*4,mean = 2 , sd = 1),50)
#' msw(X)
#'
#' ## load the ubiquitous multivariate iris data ##
#' ## (first 50 observations of columns 1:4) ##
#' iris.df = iris[1:50, 1:4]
#' msw(iris.df)
#'
msw <- function (X) {

  X <- as.matrix(X)
  X <- X[complete.cases(X), ]
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)

  VAGE<- function (X)
  {
  dname <- deparse(substitute(X))
  if (is.vector(X) == TRUE)
    X = cbind(X)
  stopifnot(is.matrix(X))
  n <- nrow(X)
  if (n < 12 || n > 5000)
    stop("Sample size must be between 12 and 5000.")
  p <- ncol(X)
  if (n <= p)
    stop("Sample size must be larger than vector dimension.")
  if (n > p) {
    x <- scale(X, scale = FALSE)
    eigenv <- eigen(var(X), symmetric = TRUE)
    e.vec <- as.matrix(eigenv$vectors)
    sqrS <- e.vec %*% diag(1/sqrt(eigenv$values), ncol = p) %*%
      t(e.vec)
    z <- t(sqrS %*% t(x))
    w <- rep(NA, p)
    for (k in 1:p) {
      w[k] <- shapiro.test(z[, k])$statistic
    }
    wast <- mean(w)
    y <- log(n)
    w1 <- log(1 - wast)
    m <- -1.5861 - 0.31082 * y - 0.083751 * y^2 + 0.0038915 *
      y^3
    s <- exp(-0.4803 - 0.082676 * y + 0.0030302 * y^2)
    s2 <- s^2
    sigma2 <- log((p - 1 + exp(s2))/p)
    mu1 <- m + s2/2 - sigma2/2
    p.value <- pnorm(w1, mean = mu1, sd = sqrt(sigma2), lower.tail = FALSE)
    results <- list(statistic = c(MVW = wast), p.value = p.value,
                    method = "Generalized Shapiro-Wilk test for Multivariate Normality by Villasenor-Alva and Gonzalez-Estrada",
                    data.name = dname)
    class(results) = "htest"
    return(results)
  }
  }

  royston <- function (a)
  {
    p = dim(a)[2]
    n = dim(a)[1]
    z <- matrix(nrow = p, ncol = 1)
    z <- as.data.frame(z)
    w <- matrix(nrow = p, ncol = 1)
    w <- as.data.frame(w)
    if (n <= 3) {
      stop("n must be greater than 3")
    }
    else if (n >= 4 || n <= 11) {
      x = n
      g = -2.273 + 0.459 * x
      m = 0.544 - 0.39978 * x + 0.025054 * x^2 - 0.0006714 *
        x^3
      s = exp(1.3822 - 0.77857 * x + 0.062767 * x^2 - 0.0020322 *
                x^3)
      for (i in 1:p) {
        a2 = a[, i]
        {
          if (kurtosis(a2) > 3) {
            w = sf.test(a2)$statistic
          }
          else {
            w = shapiro.test(a2)$statistic
          }
        }
        z[i, 1] = (-log(g - (log(1 - w))) - m)/s
      }
    }
    if (n > 2000) {
      stop("n must be less than 2000")
    }
    else if (n >= 12 || n <= 2000) {
      x = log(n)
      g = 0
      m = -1.5861 - 0.31082 * x - 0.083751 * x^2 + 0.0038915 *
        x^3
      s = exp(-0.4803 - 0.082676 * x + 0.0030302 * x^2)
      for (i in 1:p) {
        a2 = a[, i]
        {
          if (kurtosis(a2) > 3) {
            w = sf.test(a2)$statistic
          }
          else {
            w = shapiro.test(a2)$statistic
          }
        }
        z[i, 1] = ((log(1 - w)) + g - m)/s
      }
    }
    else {
      stop("n is not in the proper range")
    }
    u = 0.715
    v = 0.21364 + 0.015124 * (log(n))^2 - 0.0018034 * (log(n))^3
    l = 5
    C = cor(a)
    NC = (C^l) * (1 - (u * (1 - C)^u)/v)
    T = sum(sum(NC)) - p
    mC = T/(p^2 - p)
    edf = p/(1 + (p - 1) * mC)
    Res <- matrix(nrow = p, ncol = 1)
    Res <- as.data.frame(Res)
    for (i in 1:p) {
      Res = (qnorm((pnorm(-z[, ]))/2))^2
    }
    RH = (edf * (sum(Res)))/p
    pv = pchisq(RH, edf, lower.tail = FALSE)
    {
      if (pv < 0.05) {
        dist.check = ("Data analyzed have a non-normal distribution.")
      }
      else {
        dist.check = ("Data analyzed have a normal distribution.")
      }
    }
    res = structure(list(H = RH, p.value = pv, distribution = dist.check))
    res
  }

  if (p == 1){
    final <- list(univ(X))
    names(final) <- c("uv.shapiro")
    return(final)
  } else if (p > 1){
    if (n <= p){
      stop("Sample size must be larger than vector dimension.")}
    if (n > p){
      if (n <= 3) {
        stop("n must be greater than 3")
      }else if(n > 3 & n < 12){
        rtest <-  royston(X)
        r.result <- ifelse(rtest$p.value > 0.05, "YES", "NO")
        r.output <- noquote(c("Royston",round(rtest$H,4),round(rtest$p.value,4),r.result))
        names(r.output) <- c("Test","Statistic","p-value","Result")

        final <- list(r.output,univ(X))
        names(final) <- c("mv.test","uv.shapiro")
        return(final)
      }else if(n >= 12 & n <= 2000){
        rtest <-  royston(X)
        r.result <- ifelse(rtest$p.value > 0.05, "YES", "NO")
        r.output <- noquote(c("Royston",round(rtest$H,4),round(rtest$p.value,4),r.result))
        names(r.output) <- c("Test","Statistic","p-value","Result")

        vtest <- VAGE(X)
        v.result <- ifelse(vtest$p.value > 0.05, "YES", "NO")
        v.output <- noquote(c("VA-GE",round(vtest$statistic,4),round(vtest$p.value,4),v.result))
        names(v.output) <- c("Test","Statistic","p-value","Result")

        combind <- noquote(rbind(r.output,v.output))
        rownames(combind)<-NULL
        final <- list(combind,univ(X))
        names(final) <- c("mv.test","uv.shapiro")
        return(final)
      }else if(n >2000 & n <=5000){
        vtest <- VAGE(X)
        v.result <- ifelse(vtest$p.value > 0.05, "YES", "NO")
        v.output <- noquote(c("VA-GE",round(vtest$statistic,4),round(vtest$p.value,4),v.result))
        names(v.output) <- c("Test","Statistic","p-value","Result")

        final <- list(v.output,univ(X))
        names(final) <- c("mv.test","uv.shapiro")
        return(final)
      }
      if(n >5000)
        stop("Sample size must be between 3 and 5000.")
    }
  }
}



#' Mardia Test (Skewness and Kurtosis) for Multivariate Normality
#'
#' @description It computes Mardia (1970)'s multivariate skewness and kurtosis statistics and their corresponding
#'  p-value. Both p-values of skewness and kurtosis statistics should be greater than 0.05 to conclude
#'  multivariate normality. The skewness statistic will be adjusted for sample size \eqn{n < 20}.
#' @param X an \eqn{n*p} numeric matrix or data frame.
#' @param std if \code{TRUE}, the data matrix or data frame will be standardized via normalizing
#' the covariance matrix by \eqn{n}.
#' @return Returns a list with two objects:
#' \describe{
#' \item{\code{mv.test}}{results of the Mardia test, i.e., test statistic, \emph{p}-value, and multivariate normality summary (YES, if both skewness and kurtosis \emph{p}-value>0.05).}
#' \item{\code{uv.shapiro}}{a dataframe with \eqn{p} rows detailing univariate Shapiro-Wilk tests. Columns in the dataframe contain test statistics \emph{W}, \emph{p}-value,and univariate normality summary (YES, if \emph{p}-value>0.05).}
#' }
#' @references Mardia, K. V. (1970). Measures of multivariate skewness and kurtosis with applications. \emph{Biometrika}, 57(3), 519-530.
#' @references Shapiro, S. S., & Wilk, M. B. (1965). An analysis of variance test for normality (complete samples). \emph{Biometrika}, 52(3/4), 591-611.
#' @references Doornik, J. A., & Hansen, H. (2008). An omnibus test for univariate and multivariate normality. \emph{Oxford Bulletin of Economics and Statistics}, 70, 927-939.
#' @references Zhou, M., & Shao, Y. (2014). A powerful test for multivariate normality. \emph{Journal of applied statistics}, 41(2), 351-363.
#' @import stats
#' @seealso \code{\link{mvnTest}}, \code{\link{faTest}}, \code{\link{msw}}, \code{\link{msk}}, \code{\link{mhz}}, \code{\link{mvn}}
#' @export
#' @examples
#' set.seed(12345)
#'
#' ## Data from gamma distribution
#' X = matrix(rgamma(50*4,shape =  2),50)
#' mardia(X)
#'
#' ## Data from normal distribution
#' X = matrix(rnorm(50*4,mean = 2 , sd = 1),50)
#' mardia(X)
#'
#' ## load the ubiquitous multivariate iris data ##
#' ## (first 50 observations of columns 1:4) ##
#' iris.df = iris[1:50, 1:4]
#' mardia(iris.df)
#'
mardia <- function (X, std = TRUE){
  X = as.data.frame(X)
  dname <- deparse(substitute(X))
  data <- X[complete.cases(X), ]
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  data.org <- data
  data <- scale(data, scale = FALSE)
  if (std) {
    S <- ((n - 1)/n) * cov(data)
  }
  else {
    S <- cov(data)
  }
  D <- data %*% solve(S, tol = 1e-25) %*% t(data)
  g1p <- sum(D^3)/n^2
  g2p <- sum(diag((D^2)))/n
  df <- p * (p + 1) * (p + 2)/6
  k = ((p + 1) * (n + 1) * (n + 3))/(n * ((n + 1) * (p + 1) -
                                            6))
  if (n < 20) {
    skew <- n * k * g1p/6
    p.skew <- pchisq(skew, df, lower.tail = FALSE)
  }
  else {
    skew <- n * g1p/6
    p.skew <- pchisq(skew, df, lower.tail = FALSE)
  }
  kurt <- (g2p - p * (p + 2)*(n-1)/(n+1)) * sqrt(n/(8 * p * (p + 2)))
  p.kurt <- 2 * (1 - pnorm(abs(kurt)))
  skewMVN = ifelse(p.skew > 0.05, "YES", "NO")
  kurtoMVN = ifelse(p.kurt > 0.05, "YES", "NO")
  MVN = ifelse(p.kurt > 0.05 && p.skew > 0.05, "YES",
               "NO")
  result <- cbind.data.frame(test = "Mardia", g1p = g1p,
                             chi.skew = skew, p.value.skew = p.skew, skewnewss = skewMVN,
                             g2p = g2p, z.kurtosis = kurt, p.value.kurt = p.kurt,
                             kurtosis = kurtoMVN, MVN = MVN)
  resultSkewness = cbind.data.frame(Test = "Skewness",
                                    Statistic = as.factor(round(skew,4)), `p-value` = as.factor(round(p.skew,4)),
                                    Result = skewMVN)
  resultKurtosis = cbind.data.frame(Test = "Kurtosis",
                                    Statistic = as.factor(round(kurt,4)), `p-value` = as.factor(round(p.kurt,4)),
                                    Result = kurtoMVN)
  MVNresult = cbind.data.frame(Test = "MV Normality", Statistic = NA,
                               `p-value` = NA, Result = MVN)
  result = rbind.data.frame(resultSkewness, resultKurtosis,
                            MVNresult)

  final <- list(result,univ(X))
  names(final) <- c("mv.test","uv.shapiro")

  return(final)
}



#' Bowman and Shenton Test for Multivariate Normality
#'
#' @description It computes Bowman and Shenton (1975)'s test statistic (MSK) and its corresponding
#'  p-value for multivariate normality. The statistic is calculated based on a combination of
#'  multivariate skewness (MS) and kurtosis (MK) such that \eqn{MSK=MS+|MK|^2}. For formulas of MS and MK,
#'  please refer to Mardia (1970). The corresponding p-value of the statistic is computed based on a
#'  simulated null distribution of MSK. The skewness statistic (MS) will be adjusted for sample size \eqn{n < 20}.
#' @param X an \eqn{n*p} numeric matrix or data frame.
#' @param B number of Monte Carlo simulations for null distribution, default is 1000 (increase B to increase the precision of p-value).
#' @return Returns a list with two objects:
#' \describe{
#' \item{\code{mv.test}}{results of the Bowman and Shenton test, i.e., test statistic, \emph{p}-value, and multivariate normality summary (YES, if \emph{p}-value>0.05).}
#' \item{\code{uv.shapiro}}{a dataframe with \eqn{p} rows detailing univariate Shapiro-Wilk tests. Columns in the dataframe contain test statistics \emph{W}, \emph{p}-value,and univariate normality summary (YES, if \emph{p}-value>0.05).}
#' }
#' @references Bowman, K. O., & Shenton, L. R. (1975). Omnibus test contours for departures from normality based on \eqn{\sqrt b_1} and \eqn{b_2}. \emph{Biometrika}, 62(2), 243-250.
#' @references Shapiro, S. S., & Wilk, M. B. (1965). An analysis of variance test for normality (complete samples). \emph{Biometrika}, 52(3/4), 591-611.
#' @references Mardia, K. V. (1970). Measures of multivariate skewness and kurtosis with applications. \emph{Biometrika}, 57(3), 519-530.
#' @references Doornik, J. A., & Hansen, H. (2008). An omnibus test for univariate and multivariate normality. \emph{Oxford Bulletin of Economics and Statistics}, 70, 927-939.
#' @references Zhou, M., & Shao, Y. (2014). A powerful test for multivariate normality. \emph{Journal of applied statistics}, 41(2), 351-363.
#' @import stats
#' @seealso \code{\link{power.msk}}, \code{\link{mvnTest}}, \code{\link{faTest}}, \code{\link{msw}}, \code{\link{mardia}}, \code{\link{mhz}}, \code{\link{mvn}}
#' @export
#' @examples
#' set.seed(12345)
#'
#' ## Data from gamma distribution
#' X = matrix(rgamma(50*4,shape =  2),50)
#' msk(X, B=100)
#'
#' ## load the ubiquitous multivariate iris data ##
#' ## (first 50 observations of columns 1:4) ##
#' iris.df = iris[1:50, 1:4]
#' msk(iris.df, B=100)
#'
msk <- function (X, B=1000)
{
  X = as.data.frame(X)
  dname <- deparse(substitute(X))
  data <- X[complete.cases(X), ]
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)

  msk_1 <- function (X)
  {
    X = as.data.frame(X)
    dname <- deparse(substitute(X))
    data <- X[complete.cases(X), ]
    data <- as.matrix(data)
    n <- nrow(data)
    p <- ncol(data)
    data.org <- data
    data <- scale(data, scale = FALSE)
    S <- ((n - 1)/n) * cov(data)
    D <- data %*% solve(S, tol = 1e-25) %*% t(data)
    g1p <- sum(D^3)/n^2
    g2p <- sum(diag((D^2)))/n
    df <- p * (p + 1) * (p + 2)/6
    k = ((p + 1) * (n + 1) * (n + 3))/(n * ((n + 1) * (p + 1) -
                                              6))
    if (n < 20) {
      skew <- n * k * g1p/6
    }
    else {
      skew <- n * g1p/6
    }
    kurt <- (g2p - p * (p + 2)*(n-1)/(n+1)) * sqrt(n/(8 * p * (p + 2)))
    msk <- skew + kurt^2

    return(msk)
  }


  dsn_msk_1 <- function(n, p, B) {
    s <- rep(0,B)
    for(i in 1:B) {
      s[i] <- msk_1(matrix(rnorm(n*p), n))
    }
    return(s)
  }

  msk <- msk_1(X)
  dsn_msk_1 <-  dsn_msk_1(n, p, B)

  p.value <- mean(dsn_msk_1 > msk)
  result <- ifelse(p.value>0.05,"YES","NO")
  output <- noquote(c(round(msk,4),round(p.value,4),result))
  names(output) <- c("Statistic","p-value","Result")

  final <- list(output,univ(X))
  names(final) <- c("mv.test","uv.shapiro")

  return(final)
}



#' Henze-Zirkler Test for Multivariate Normality
#'
#' @description It computes a multiviariate normality test based on a non-negative functional distance which
#' was proposed by Henze and Zirkler (1990). Under the null hypothesis the test statistic is approximately
#' log-normally distributed.
#' @param X an \eqn{n*p} numeric matrix or data frame.
#' @return Returns a list with two objects:
#' \describe{
#' \item{\code{mv.test}}{results of the Henze-Zirkler test, i.e., test statistic, \emph{p}-value, and multivariate normality summary (YES, if \emph{p}-value>0.05).}
#' \item{\code{uv.shapiro}}{a dataframe with \eqn{p} rows detailing univariate Shapiro-Wilk tests. Columns in the dataframe contain test statistics \emph{W}, \emph{p}-value,and univariate normality summary (YES, if \emph{p}-value>0.05).}
#' }
#' @references Henze, N., & Zirkler, B. (1990). A class of invariant consistent tests for multivariate normality. \emph{Communications in statistics-Theory and Methods}, 19(10), 3595-3617.
#' @references Shapiro, S. S., & Wilk, M. B. (1965). An analysis of variance test for normality (complete samples). \emph{Biometrika}, 52(3/4), 591-611.
#' @references Zhou, M., & Shao, Y. (2014). A powerful test for multivariate normality. \emph{Journal of applied statistics}, 41(2), 351-363.
#' @import stats
#' @seealso \code{\link{power.mhz}}, \code{\link{mvnTest}}, \code{\link{faTest}}, \code{\link{msw}}, \code{\link{msk}}, \code{\link{mardia}}, \code{\link{mvn}}
#' @export
#' @examples
#' set.seed(12345)
#'
#' ## Data from gamma distribution
#' X = matrix(rgamma(50*4,shape =  2),50)
#' mhz(X)
#'
#' ## Data from normal distribution
#' X = matrix(rnorm(50*4,mean = 2 , sd = 1),50)
#' mhz(X)
#'
#' ## load the ubiquitous multivariate iris data ##
#' ## (first 50 observations of columns 1:4) ##
#' iris.df = iris[1:50, 1:4]
#' mhz(iris.df)
mhz <- function(X){

  X   <- as.matrix(X)
  X   <- X[complete.cases(X),]
  X   <- as.matrix(X)
  n   <- NROW(X)
  p   <- NCOL(X)

  HZ <- function(X) {
    HZoptb <- function(n, p)  ((2*p + 1)*n/4)^(1/(p + 4))/sqrt(2)

    X   <- as.matrix(X)
    n   <- NROW(X)
    p   <- NCOL(X)
    S   <- var(X)*(n-1)/n
    b   <- HZoptb(n,p)
    if(n == 1)  return(NA)
    if( min(eigen(S)$values) < 1e-3)  return(4*n)
    else {
      Y <- STD(X)
      DY    <- as.matrix(dist(Y, diag = TRUE, upper = TRUE))
      part1 <- sum(exp(- b^2 * DY^2/2))/n
      part2 <- (-2*(1+b^2)^{-p/2}*sum(exp(-b^2/(2*(1+b^2))*apply(Y^2,1,sum))))
      part3 <- (1+2*b^2)^{-p/2}*n

      statistic <- part1+part2+part3

      mu <- 1 - ((1 + 2 * b^2)^(-p/2)) * (1 + (p * b^2)/(1 +
                                                           2 * b^2) + (p * (p + 2) * b^4)/(2 * (1 + 2 * b^2)^2))
      w.b <- (1 + b^2) * (1 + 3 * b^2)
      sigma.sq <- 2 * (1 + 4 * b^2)^(-p/2) + 2 * (1 + 2 * b^2)^(-p) *
        (1 + (2 * p * b^4)/(1 + (2 * b^2))^2 + (3 * p *
                                                  (p + 2) * b^8)/(4 * (1 + 2 * b^2)^4)) - 4 *
        (w.b^(-p/2)) * (1 + (3 * p * b^4)/(2 * w.b) +
                          (p * (p + 2) * b^8)/(2 * w.b^2))
      p.mu <- log(sqrt(mu^4/(sigma.sq + mu^2)))
      p.sigma <- sqrt(log((sigma.sq + mu^2)/mu^2))
      p.value <- 1 - (plnorm(statistic, p.mu, p.sigma))

      result <- ifelse(p.value>0.05,"YES","NO")
      output <- noquote(c(round(statistic,4),round(p.value,4),result))
      names(output) <- c("Statistic","p-value","Result")
      return(output)
    }
  }

  hz <- HZ(X)

  final <- list(hz,univ(X))
  names(final) <- c("mv.test","uv.shapiro")
  return(final)
}



#' Power Calculation using the Zhou-Shao's Multivariate Normality Test Statistic (\eqn{T_n})
#'
#' @description Empirical power calculation using the Zhou-Shao's multivariate normality test Statistic \eqn{T_n}.
#' @param a significance level (\eqn{\alpha}).
#' @param n number of rows (observations).
#' @param p number of columns (variables), \eqn{n>p}.
#' @param B number of Monte Carlo simulations, default is 1000 (can increase B to increase the precision).
#' @param pct percentiles of MK to get c1 and c2 described in the reference paper,default is (0.01, 0.99).
#' @param FUN self-defined function for generate multivariate distribution. See example.
#' @param ... optional arguments passed to \code{FUN}.
#' @return Returns a numeric value of the estimated empirical power (value between 0 and 1).
#' @references Zhou, M., & Shao, Y. (2014). A powerful test for multivariate normality. \emph{Journal of applied statistics}, 41(2), 351-363.
#' @import stats
#' @export
#' @examples
#' set.seed(12345)
#'
#' ## Power calculation against bivariate (p=2) independent Beta(1, 1) distribution ##
#' ## at sample size n=50 for Tn at one-sided alpha = 0.05 ##
#'
#' power.mvnTest(a = 0.05, n = 50, p = 2,  B = 100, pct = c(0.01, 0.99), FUN=IMMV, D1=runif)
#'
power.mvnTest <- function(a, n, p, B=1000, pct=c(0.01,0.99), FUN,...){

  mk = Dsn_MK(n, p, B)

  ##obtain c1 and c2 from the null distribution of MK
  cri_mk = quantile(mk, pct)

  ## Component of the Tn test statistic
  SWM2 <- function(X){
    SWS <- function(x) shapiro.test(x)$statistic
    X <- as.matrix(X)
    p = NCOL(X)
    StdX = STD(X)
    return((mean(apply(StdX, 2, SWS)) + mean(sort(apply(StdX%*%t(StdX), 2, SWS))[1:p]))/2)
  }

  ## Tn statistic function
  SWM2MK <- function(X, cri_mk) {
    Mk      <- MK(X)
    Index   <- (Mk > cri_mk[2]|Mk < cri_mk[1])
    return( (1 - Index)*(1 - SWM2(X)) + Index)
  }

  FUN  <- match.fun(FUN)
  s    <- rep(0,B)
  for(i in 1:B){
    X    <- FUN(n, p, ...)
    s[i] <- SWM2MK(X, cri_mk)
  }

  tn  <- Dsn_SWM2MK(n, p, B, pct)
  cri <- quantile(tn, (1-a))

  return(sum(s > cri)/B)
}



#' Power Calculation using the Fattorini's FA Test Statistic
#'
#' @description Empirical power calculation using the Fattorini's FA Test Statistic.
#' @param a significance level (\eqn{\alpha}).
#' @param n number of rows (observations).
#' @param p number of columns (variables), \eqn{n>p}.
#' @param B number of Monte Carlo simulations, default is 1000 (can increase B to increase the precision).
#' @param FUN self-defined function for generate multivariate distribution. See example.
#' @param ... optional arguments passed to \code{FUN}.
#' @return Returns a numeric value of the estimated empirical power (value between 0 and 1).
#' @references Fattorini, L. (1986). Remarks on the use of Shapiro-Wilk statistic for testing multivariate normality. \emph{Statistica}, 46(2), 209-217.
#' @import stats
#' @export
#' @examples
#' set.seed(12345)
#'
#' ## Power calculation against bivariate (p=2) independent Beta(1, 1) distribution ##
#' ## at sample size n=50 at one-sided alpha = 0.05 ##
#'
#' power.faTest(a = 0.05, n = 50, p = 2,  B = 100, FUN=IMMV, D1=runif)
#'
power.faTest <- function(a, n, p, B=1000, FUN,...) {

  ## SW statistic ##
  SW <- function(x)  shapiro.test(x)$statistic

  ## FA statistic ##
  FA <- function(X) {
    X <- as.matrix(X)
    n <- NROW(X)
    p <- NCOL(X)
    mu <- apply(X,2,mean)
    nSinver <- solve((n-1)*cov(X))
    Y <- X%*%t((X-matrix(rep(mu,n),ncol=p,byrow=T))%*%nSinver)
    return(1 - min(apply(Y,2,SW)))
  }

  FUN <- match.fun(FUN)
  s   <- rep(0,B)
  for(i in 1:B) {
    s[i] <- FA(FUN(n, p, ...))
  }

  tn  <- Dsn_FA(n, p, B)
  cri <- quantile(tn, (1-a))

  return(sum(s > cri)/B)
}



#' Power Calculation using the Henze-Zirkler Test Statistic
#'
#' @description Empirical power calculation using the Henze-Zirkler Test Statistic.
#' @param a significance level (\eqn{\alpha}).
#' @param n number of rows (observations).
#' @param p number of columns (variables), \eqn{n>p}.
#' @param B number of Monte Carlo simulations, default is 1000 (can increase B to increase the precision).
#' @param FUN self-defined function for generate multivariate distribution. See example.
#' @param ... optional arguments passed to \code{FUN}.
#' @return Returns a numeric value of the estimated empirical power (value between 0 and 1).
#' @references Henze, N., & Zirkler, B. (1990). A class of invariant consistent tests for multivariate normality. \emph{Communications in statistics-Theory and Methods}, 19(10), 3595-3617.
#' @import stats
#' @export
#' @examples
#' set.seed(12345)
#'
#' ## Power calculation against bivariate (p=2) independent Beta(1, 1) distribution ##
#' ## at sample size n=50 at one-sided alpha = 0.05 ##
#'
#' power.mhz(a = 0.05, n = 50, p = 2,  B = 100, FUN=IMMV, D1=runif)
#'
power.mhz <- function(a, n, p, B=1000, FUN,...) {
  FUN <- match.fun(FUN)
  s   <-  rep(0,B)
  for(i in 1:B) {
    s[i] <- as.numeric(mhz(FUN(n, p, ...))$mv.test[1])
  }

  tn  <- Dsn_HZ(n, p, B)
  cri <- quantile(tn, (1-a))

  return(sum(s > cri)/B)
}



#' Power Calculation using the Bowman and Shenton Test Statistic
#'
#' @description Empirical power calculation using Bowman and Shenton Test Statistic.
#' @param a significance level (\eqn{\alpha}).
#' @param n number of rows (observations).
#' @param p number of columns (variables), \eqn{n>p}.
#' @param B number of Monte Carlo simulations, default is 1000 (can increase B to increase the precision).
#' @param FUN self-defined function for generate multivariate distribution. See example.
#' @param ... optional arguments passed to \code{FUN}.
#' @return Returns a numeric value of the estimated empirical power (value between 0 and 1).
#' @references Bowman, K. O., & Shenton, L. R. (1975). Omnibus test contours for departures from normality based on \eqn{\sqrt b_1} and \eqn{b_2}. \emph{Biometrika}, 62(2), 243-250.
#' @import stats
#' @export
#' @examples
#' set.seed(12345)
#'
#' ## Power calculation against bivariate (p=2) independent Beta(1, 1) distribution ##
#' ## at sample size n=50 at one-sided alpha = 0.05 ##
#'
#' power.msk(a = 0.05, n = 50, p = 2,  B = 100, FUN=IMMV, D1=runif)
#'
power.msk <- function(a, n, p, B=1000, FUN,...){
  FUN  <- match.fun(FUN)
  s    <- rep(0,B)
  for(i in 1:B){
    X    <- FUN(n, p, ...)
    s[i] <- MSMK(X)
  }

  tn  <- Dsn_MSMK(n, p, B)
  cri <- quantile(tn, (1-a))

  return(sum(s > cri)/B)
}



#' Power Calculation using the SWT-based Royston Test Statistic
#'
#' @description Empirical power calculation using Royston test statistic.
#' @param a significance level (\eqn{\alpha}).
#' @param n number of rows (observations).
#' @param p number of columns (variables), \eqn{n>p}.
#' @param B number of Monte Carlo simulations, default is 1000 (can increase B to increase the precision).
#' @param FUN self-defined function for generate multivariate distribution. See example.
#' @param ... optional arguments passed to \code{FUN}.
#' @return Returns a numeric value of the estimated empirical power (value between 0 and 1).
#' @references Royston, J. P. (1982). An extension of Shapiro and Wilk's W test for normality to large samples. \emph{Journal of the Royal Statistical Society: Series C (Applied Statistics)}, 31(2), 115-124.
#' @import stats
#' @export
#' @examples
#' set.seed(12345)
#'
#' ## Power calculation against bivariate (p=2) independent Beta(1, 1) distribution ##
#' ## at sample size n=50 at one-sided alpha = 0.05 ##
#'
#' power.mswR(a = 0.05, n = 50, p = 2,  B = 100, FUN=IMMV, D1=runif)
#'
power.mswR <- function(a, n, p, B=1000, FUN,...){
  FUN  <- match.fun(FUN)
  s    <- rep(0,B)
  for(i in 1:B){
    X    <- FUN(n, p, ...)
    s[i] <- as.numeric(msw(X)$mv.test[1,2])
  }

  tn  <- Dsn_royston(n, p, B)
  cri <- quantile(tn, (1-a))

  return(sum(s > cri)/B)
}



#' Power Calculation using the SWT-based Villasenor-Alva and Gonzalez-Estrada (VAGE) Test Statistic
#'
#' @description Empirical power calculation using VAGE test statistic.
#' @param a significance level (\eqn{\alpha}).
#' @param n number of rows (observations).
#' @param p number of columns (variables), \eqn{n>p}.
#' @param B number of Monte Carlo simulations, default is 1000 (can increase B to increase the precision).
#' @param FUN self-defined function for generate multivariate distribution. See example.
#' @param ... optional arguments passed to \code{FUN}.
#' @return Returns a numeric value of the estimated empirical power (value between 0 and 1).
#' @references Villasenor Alva, J. A., & Estrada, E. G. (2009). A generalization of Shapiro–Wilk's test for multivariate normality. \emph{Communications in Statistics—Theory and Methods}, 38(11), 1870-1883.
#' @import stats
#' @export
#' @examples
#' set.seed(12345)
#'
#' ## Power calculation against bivariate (p=2) independent Beta(1, 1) distribution ##
#' ## at sample size n=50 at one-sided alpha = 0.05 ##
#'
#' power.mswV(a = 0.05, n = 50, p = 2,  B = 100, FUN=IMMV, D1=runif)
#'
power.mswV <- function(a, n, p, B=1000, FUN,...){
  FUN  <- match.fun(FUN)
  s    <- rep(0,B)
  for(i in 1:B){
    X    <- FUN(n, p, ...)
    s[i] <- as.numeric(msw(X)$mv.test[2,2])
  }

  tn  <- Dsn_VAGE(n, p, B)
  cri <- quantile(tn, a)

  return(sum(s < cri)/B)
}
# power.mswV <- function(a, n, p, B, FUN,...){
#   FUN  <- match.fun(FUN)
#   s    <- rep(0,B)
#   for(i in 1:B){
#     X    <- FUN(n, p, ...)
#     s[i] <- as.numeric(mvShapiro.Test(X)$statistic)
#   }
#
#   tn  <- Dsn_VAGE(n, p, B)
#   cri <- quantile(tn, a)
#
#   return(sum(s < cri)/B)
# }


#' Power Calculation using the Univariate Shapiro-Wilk Test Statistic
#'
#' @description Empirical power calculation using univariate Shapiro-Wilk test statistic.
#' @param a significance level (\eqn{\alpha}).
#' @param n number of rows (observations).
#' @param p p=1 for univariate.
#' @param B number of Monte Carlo simulations, default is 1000 (can increase B to increase the precision).
#' @param FUN self-defined function for generate multivariate distribution. See example.
#' @param ... optional arguments passed to \code{FUN}.
#' @return Returns a numeric value of the estimated empirical power (value between 0 and 1).
#' @references Shapiro, S. S., & Wilk, M. B. (1965). An analysis of variance test for normality (complete samples). \emph{Biometrika}, 52(3/4), 591-611.
#' @import stats
#' @export
#' @examples
#' set.seed(12345)
#'
#' ## Power calculation against univariate (p=1) independent Beta(1, 1) distribution ##
#' ## at sample size n=50 at one-sided alpha = 0.05 ##
#'
#' power.usw(a = 0.05, n = 50, p = 1,  B = 100, FUN=IMMV, D1=runif)
#'
power.usw <- function(a, n, p=1, B=1000, FUN,...){
  FUN  <- match.fun(FUN)
  s    <- rep(0,B)
  for(i in 1:B){
    X    <- FUN(n, p, ...)
    s[i] <- as.numeric(msw(X)$uv.shapiro[1,1])
  }

  tn  <- Dsn_usw(n, p, B)
  cri <- quantile(tn, a)

  return(sum(s < cri)/B)
}



#' Random Generation for the Spherically Symmetric Pearson Type II Distribution
#'
#' @description Generate univariate or multivariate random sample for the spherically symmetric Pearson type II distribution.
#' @param n number of rows (observations).
#' @param p number of columns (variables).
#' @param s shape parameter, \eqn{s>-1}.
#' @return Returns univariate (\eqn{p=1}) or multivariate (\eqn{p>1}) random sample matrix.
#' @references Kotz, S. (1975). Multivariate distributions at a cross road. In A Modern Course on Statistical Distributions in Scientific Work (pp. 247-270). Springer, Dordrecht.
#' @references Henze, N., & Zirkler, B. (1990). A class of invariant consistent tests for multivariate normality. \emph{Communications in statistics-Theory and Methods}, 19(10), 3595-3617.
#' @import stats
#' @export
#' @examples
#' set.seed(12345)
#'
#' ## Generate 5X2 random sample matrix from PSII(s=1) ##
#' PSII(n=5, p=2, s=1)
#'
#'
#' ## Power calculation against bivariate (p=2) PSII(s=1) distribution ##
#' ## at sample size n=50 at one-sided alpha = 0.05 ##
#'
#' # Zhou-Shao's test #
#' power.mvnTest(a = 0.05, n = 50, p = 2,  B = 100, FUN = PSII, s = 1)
#'
PSII <- function(n, p, s) {

  V      <- matrix(rnorm(n*p), n)
  Vnew   <- V/sqrt(apply(V^2, 1, sum))
  repara <- 1/(s + 1)
  r    <- sqrt(rbeta(n, shape1 = p/2, shape2 = s + 1))

  return(r*Vnew)
}



#' Random Generation for the Spherically Symmetric Pearson Type VII Distribution
#'
#' @description Generate univariate or multivariate random sample for the spherically symmetric Pearson type VII distribution.
#' @param n number of rows (observations).
#' @param p number of columns (variables).
#' @param s shape parameter, \eqn{s > p/2}.
#' @return Returns univariate (\eqn{p=1}) or multivariate (\eqn{p>1}) random sample matrix.
#' @references Kotz, S. (1975). Multivariate distributions at a cross road. In A Modern Course on Statistical Distributions in Scientific Work (pp. 247-270). Springer, Dordrecht.
#' @references Henze, N., & Zirkler, B. (1990). A class of invariant consistent tests for multivariate normality. \emph{Communications in statistics-Theory and Methods}, 19(10), 3595-3617.
#' @import stats
#' @export
#' @examples
#' set.seed(12345)
#'
#' ## Generate 5X2 random sample matrix from PSVII(s=3) ##
#' PSVII(n=5, p=2, s=3)
#'
#'
#' ## Power calculation against bivariate (p=2) PSVII(s=3) distribution ##
#' ## at sample size n=50 at one-sided alpha = 0.05 ##
#'
#' # Zhou-Shao's test #
#' power.mvnTest(a = 0.05, n = 50, p = 2,  B = 100, FUN = PSVII, s = 3)
#'
PSVII <- function(n, p, s) {

  V      <- matrix(rnorm(n*p), n)
  Vnew   <- V/sqrt(apply(V^2, 1, sum))
  repara <- - 1/(s - 1)
  r    <- sqrt(1/rbeta(n, shape1 = s - p/2, shape2 = p/2) - 1)

  return(r*Vnew)
}



#' Random Generation for General Spherically Symmetric Distributions
#'
#' @description Generate univariate or multivariate random sample for general spherically symmetric distributions.
#' @param n number of rows (observations).
#' @param p number of columns (variables).
#' @param D random generation functions for some distributions (e.g., \code{rgamma}, \code{rbeta}).
#' @param ... optional arguments passed to \code{D}.
#' @return Returns univariate (\eqn{p=1}) or multivariate (\eqn{p>1}) random sample matrix.
#' @references Chmielewski, M. A. (1981). Elliptically symmetric distributions: A review and bibliography. \emph{International Statistical Review/Revue Internationale de Statistique}, 67-74.
#' @references Henze, N., & Zirkler, B. (1990). A class of invariant consistent tests for multivariate normality. \emph{Communications in statistics-Theory and Methods}, 19(10), 3595-3617.
#' @import stats
#' @export
#' @examples
#' set.seed(12345)
#'
#' ## Generate 5X2 random sample matrix from SPH(Beta(1,1)) ##
#' SPH(n=5, p=2, D=rbeta, shape1=1, shape2=1)
#'
#'
#' ## Power calculation against bivariate (p=2) SPH(Beta(1,1)) distribution ##
#' ## at sample size n=50 at one-sided alpha = 0.05 ##
#'
#' # Zhou-Shao's test #
#' power.mvnTest(a=0.05, n=50, p=2, B=100, FUN=SPH, D=rbeta, shape1=1, shape2=1)
#'
SPH <- function(n, p, D,...){

  FUN  <- match.fun(D)
  r    <- FUN(n, ...)
  V    <- matrix(rnorm(n*p), n)
  Vnew <- V/sqrt(apply(V^2, 1, sum))

  return(r*Vnew)
}



#' Random Generation for Distribution with Independent Marginals
#'
#' @description Generate univariate or multivariate random sample for distribution with independent marginals such that \eqn{D_1 \otimes D_2}.
#' \eqn{D_1 \otimes D_2} denotes the distribution having independent marginal distributions \eqn{D_1} and \eqn{D_2}. This function can generate
#' multivariate random samples only from distribution \eqn{D_1} or from both \eqn{D_1} and \eqn{D_2}.
#' @param n number of rows (observations).
#' @param p total number of columns (variables).
#' @param q number of columns from distribution \code{D1} if generate multivariate samples from independent marginal distribution \eqn{D_1} and \eqn{D_2}.
#' Default is \code{NULL}, i.e., generating samples only from one distribution.
#' @param D1 random generation function for 1st distribution (e.g., \code{rnorm}, \code{rbeta}).
#' @param D2 random generation function for 2nd distribution (e.g., \code{rnorm}, \code{rbeta}).
#' @param D1.args a list of optional arguments passed to \code{D1}.
#' @param D2.args a list of optional arguments passed to \code{D2}.
#' @return Returns univariate (\eqn{p=1}) or multivariate (\eqn{p>1}) random sample matrix.
#' @references Zhou, M., & Shao, Y. (2014). A powerful test for multivariate normality. \emph{Journal of applied statistics}, 41(2), 351-363.
#' @references Henze, N., & Zirkler, B. (1990). A class of invariant consistent tests for multivariate normality. \emph{Communications in statistics-Theory and Methods}, 19(10), 3595-3617.
#' @import stats
#' @export
#' @examples
#' set.seed(12345)
#'
#' ## Generate 5X2 random sample matrix from IMMV(N(0,1),Beta(1,2)) ##
#' IMMV(n=5, p=2, q=1, D1=rbeta, D1.args=list(shape1=1,shape2=2), D2=rnorm)
#'
#'
#' ## Power calculation against bivariate (p=2) IMMV(Gamma(5,1)) distribution ##
#' ## at sample size n=50 at one-sided alpha = 0.05 ##
#'
#' # Zhou-Shao's test #
#' power.mvnTest(a=0.05, n=50, p=2, B=100, FUN=IMMV, D1=rgamma, D1.args=list(shape=5, rate=1))
#'
#' ## Power calculation against bivariate (p=2) IMMV(N(0,1),Beta(1,2)) distribution ##
#' ## at sample size n=50 at one-sided alpha = 0.05 ##
#'
#' # Zhou-Shao's test #
#' power.mvnTest(a=0.05, n=50, p=2, B=100, FUN=IMMV, q=1, D1=rbeta, D1.args=list(shape1=1,shape2=2),
#' D2=rnorm)
#'
IMMV <- function(n,p,q=NULL,D1,D2=NULL,D1.args=list(),D2.args=list()){
  if(is.null(q)){
    FUN1 <- match.fun(D1)
    ds1 <- matrix(do.call(FUN1, c(n*p, D1.args)),n)
    return(ds1)
  }else{
    FUN1 <- match.fun(D1)
    ds1 <- matrix(do.call(FUN1, c(n*q, D1.args)),n)
    FUN2 <- match.fun(D2)
    ds2 <- matrix(do.call(FUN2, c(n*(p-q), D2.args)),n)
    return(cbind(ds1,ds2))
  }
}



#' Random Generation for the Normal Mixture Distribution
#'
#' @description Generate univariate or multivariate random sample for the normal mixture distribution with density
#' \eqn{\lambda N(0,\sum_1)+(1-\lambda)N(bl, \sum_2)}, where \eqn{l} is the column vector with all elements being 1,
#' \eqn{\sum_i=(1-\rho_i)I+\rho_ill^T} for \eqn{i=1,2}. \eqn{\rho} has to satisfy \eqn{\rho > -1/(p-1)} in order to make the
#' covariance matrix meaningful.
#' @param n number of rows (observations).
#' @param p total number of columns (variables).
#' @param lambda weight parameter to allocate the proportions of the mixture, \eqn{0<\lambda<1}.
#' @param mu2 is \eqn{bl} of \eqn{N(bl, \sum_2)}.
#' @param rho1 parameter in \eqn{\sum_1}.
#' @param rho2 parameter in \eqn{\sum_2}.
#' @return Returns univariate (\eqn{p=1}) or multivariate (\eqn{p>1}) random sample matrix.
#' @references Zhou, M., & Shao, Y. (2014). A powerful test for multivariate normality. \emph{Journal of applied statistics}, 41(2), 351-363.
#' @import stats
#' @export
#' @examples
#' set.seed(12345)
#'
#' ## Generate 5X2 random sample matrix from MVNMIX(0.5,4,0,0) ##
#' MVNMIX(n=5, p=2, lambda=0.5, mu2=4, rho1=0, rho2=0)
#'
#'
#' ## Power calculation against bivariate (p=2) MVNMIX(0.5,4,0,0) distribution ##
#' ## at sample size n=50 at one-sided alpha = 0.05 ##
#'
#' # Zhou-Shao's test #
#' power.mvnTest(a=0.05, n=50, p=2, B=100, FUN=MVNMIX, lambda=0.5, mu2=4, rho1=0, rho2=0)
#'
MVNMIX <- function(n, p, lambda, mu2, rho1 = 0, rho2 = 0) {
  A1 <- diag(sqrt(1 - rho1), p) +
    (sqrt(rho1*p + 1 - rho1) - sqrt(1 - rho1))/p * matrix(rep(1, p*p), p)
  A2 <- diag(sqrt(1 - rho2), p) +
    (sqrt(rho2*p + 1 - rho2) - sqrt(1 - rho2))/p * matrix(rep(1, p*p), p)
  u <- runif(n)
  S <- NULL
  for( i in 1:n) {
    if(u[i] < lambda)  X <- t(A1%*%rnorm(p))
    else               X <- mu2 + t(A2%*%rnorm(p))
    S <- rbind(S, X)
  }
  return(S)
}



#' Random Generation for the Copula Generated Distributions
#'
#' @description Generate univariate or multivariate random sample for the Copula Generated Distributions.
#' @param n number of rows (observations).
#' @param p total number of columns (variables).
#' @param c name of an Archimedean copula, choosing from "\code{clayton}" (default), "\code{frank}", or "\code{gumbel}".
#' @param param number (numeric) specifying the copula parameter.
#' @param invF inverse function (quantile function, e.g. \code{qnorm}).
#' @param ... optional arguments passed to \code{invF}.
#' @return univariate (\eqn{p=1}) or multivariate (\eqn{p>1}) random sample.
#' @references Yan, J. (2007). Enjoy the joy of copulas: with a package copula. \emph{Journal of Statistical Software}, 21(4), 1-21.
#' @importFrom copula claytonCopula frankCopula gumbelCopula rCopula
#' @export
#' @examples
#' set.seed(12345)
#'
#' ## Generate 5X2 random sample matrix from Clayton(0.5, qnorm) ##
#' copulas(n=50, p=2, c="clayton", param=0.5, invF=qnorm)
#'
#'
#' ## Power calculation against bivariate (p=2) Clayton(0.5, qnorm) distribution ##
#' ## at sample size n=50 at one-sided alpha = 0.05 ##
#'
#' # Zhou-Shao's test #
#' power.mvnTest(a=0.05, n=50, p=2, B=100, FUN=copulas, c="clayton", param=0.5, invF=qnorm)
#'
copulas <- function(n, p, c="clayton", param, invF, ...) {
  if(c=="clayton"){
    obj  <- claytonCopula(param, dim = p)
  }else if(c=="frank"){
    obj  <- frankCopula(param, dim = p)
  }else if(c=="gumbel"){
    obj  <- gumbelCopula(param, dim = p)
  }

  u    <- rCopula(n=n,copula=obj)
  invF <- match.fun(invF)
  s    <- apply(u, 2, invF, ...)
  return(s)
}








##### some inner functions for above functions ########
## standardize data ##
STD <- function(X){
  n     <- NROW(X)
  p     <- NCOL(X)
  if(p == 1) {
    X <- as.numeric(X)
    return(as.matrix((X-mean(X))/sd(X)*sqrt(n/(n-1))))}
  else {
    X1    <- (X-matrix(rep(apply(X,2,mean),n),ncol=p,byrow=T))
    EIGEN <- eigen(cov(X)*(n-1)/n)
    S1    <- EIGEN$vectors%*%diag((sqrt(EIGEN$values))^{-1})%*%t(EIGEN$vectors)
    return(X1%*%t(S1))
  }
}


## MK Test Statistic ##
MK <- function(X) {
  n <- NROW(X)
  p <- NCOL(X)
  if(p==1) {
    mu<-mean(X)
    b21<-n*sum((X-mu)^4)/(sum((X-mu)^2))^2
    return(sqrt(n/24)*(b21-3*(n-1)/(n+1)))
  }
  else
  {
    stdX <- STD(X)
    b    <- apply(stdX^2,1,sum)
    return(sqrt(n/(8*p*(p+2)))*(mean(b^2)-p*(p+2)*(n-1)/(n+1)))
  }
}


## Simulate Null Distribution of the MK Statistic ##
Dsn_MK <- function(n, p, B) {
  s <- rep(0,B)
  for(i in 1:B) {
    s[i] <- MK(matrix(rnorm(n*p), n))
  }
  return(s)
}


## Simulate Null Distribution of the Zhou-Shao's Test Statistic (Tn) ##
Dsn_SWM2MK <- function(n, p, B, pct) {
  mk = Dsn_MK(n, p, B)
  ##obtain c1 and c2 from the null distribution of MK
  cri_mk = quantile(mk, pct)
  ## Component of the Tn test statistic
  SWM2 <- function(X){
    SWS <- function(x) shapiro.test(x)$statistic
    X <- as.matrix(X)
    p = NCOL(X)
    StdX = STD(X)
    return((mean(apply(StdX, 2, SWS)) + mean(sort(apply(StdX%*%t(StdX), 2, SWS))[1:p]))/2)
  }
  ## Tn statistic function
  SWM2MK <- function(X, cri_mk) {
    Mk      <- MK(X)
    Index   <- (Mk > cri_mk[2]|Mk < cri_mk[1])
    return( (1 - Index)*(1 - SWM2(X)) + Index)
  }
  s <- rep(0,B)
  for(i in 1:B) {
    s[i] <- SWM2MK(matrix(rnorm(n*p), n), cri_mk)
  }
  return(s)
}


## Simulate Null Distribution of FA Test Statistic ##
Dsn_FA <- function(n, p, B){

  ## SW statistic ##
  SW <- function(x)  shapiro.test(x)$statistic

  ## FA statistic ##
  FA <- function(X) {
    X <- as.matrix(X)
    n <- NROW(X)
    p <- NCOL(X)
    mu <- apply(X,2,mean)
    nSinver <- solve((n-1)*cov(X))
    Y <- X%*%t((X-matrix(rep(mu,n),ncol=p,byrow=T))%*%nSinver)
    return(1 - min(apply(Y,2,SW)))
  }

  s <- rep(0, B)

  for(i in 1:B) {
    s[i] <- FA(matrix(rnorm(n*p), n))
  }
  return(s)

}


## Simulate Null Distribution of HZ Test Statistic ##
Dsn_HZ <- function(n, p, B) {
  s  <- rep(0, B)
  for( i in 1:B) {
    s[i] <- as.numeric(mhz(matrix(rnorm(n*p), n))$mv.test[1])
  }
  return(s)
}


## MSK statistic function ##
MSMK <- function (X)
{
  dataframe = as.data.frame(X)
  dname <- deparse(substitute(X))
  data <- X[complete.cases(X), ]
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  data.org <- data
  data <- scale(data, scale = FALSE)
  S <- ((n - 1)/n) * cov(data)
  D <- data %*% solve(S, tol = 1e-25) %*% t(data)
  g1p <- sum(D^3)/n^2
  g2p <- sum(diag((D^2)))/n
  df <- p * (p + 1) * (p + 2)/6
  k = ((p + 1) * (n + 1) * (n + 3))/(n * ((n + 1) * (p + 1) -
                                            6))
  if (n < 20) {
    skew <- n * k * g1p/6
  }
  else {
    skew <- n * g1p/6
  }
  kurt <- (g2p - p * (p + 2)*(n-1)/(n+1)) * sqrt(n/(8 * p * (p + 2)))
  msk <- skew + kurt^2

  return(msk)
}


## Simulate Null Distribution of MSK Test Statistic ##
Dsn_MSMK <- function(n, p, B) {
  s <- rep(0,B)
  for(i in 1:B) {
    s[i] <- MSMK(matrix(rnorm(n*p), n))
  }
  return(s)
}


## Simulate Null Distribution of MSW-Royston Test Statistic ##
Dsn_royston <- function(n, p, B) {
  s <- rep(0,B)
  for(i in 1:B) {
    s[i] <- as.numeric(msw(matrix(rnorm(n*p), n))$mv.test[1,2])
  }
  return(s)
}


## Simulate Null Distribution of MSW-VAGE Test Statistic ##
Dsn_VAGE <- function(n, p, B) {
  s <- rep(0,B)
  for(i in 1:B) {
    s[i] <- as.numeric(msw(matrix(rnorm(n*p), n))$mv.test[2,2])
  }
  return(s)
}
# Dsn_VAGE <- function(n, p, B) {
#   s <- rep(0,B)
#   for(i in 1:B) {
#     s[i] <- as.numeric(mvShapiro.Test(matrix(rnorm(n*p), n))$statistic)
#   }
#   return(s)
# }


## Univariate Shapiro-Wilk test ##
univ <- function(X){
  univariate_1 <- t(sapply(as.data.frame(X), function(x) shapiro.test(x))[1:2,])
  uvresult <- ifelse(univariate_1[,2]>0.05, "Yes", "No")
  univariate_2 <- rbind(round(Reduce(cbind,univariate_1[,1]),4),round(Reduce(cbind,univariate_1[,2]),4),uvresult)
  colnames(univariate_2)<-names(univariate_1[,1])
  univariate <- noquote(t(univariate_2))
  colnames(univariate) <- c("W","p-value","UV.Normality")
  return(univariate)
}

## Simulate Null Distribution of univariate SW Test Statistic ##
Dsn_usw <- function(n, p, B) {
  s <- rep(0,B)
  for(i in 1:B) {
    s[i] <- as.numeric(msw(matrix(rnorm(n*p), n))$uv.shapiro[1])
  }
  return(s)
}
