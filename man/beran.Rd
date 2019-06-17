\name{beran}

\alias{beran}

\title{Compute Beran's Estimator of the Conditional Survival
}

\description{This function computes the Beran nonparametric estimator of
  the conditional survival function.
}

\usage{beran(x, t, d, dataset, x0, h, local = TRUE, testimate = NULL,
conflevel = 0L, cvbootpars = if (conflevel == 0 && !missing(h)) NULL
else controlpars())
}

\arguments{
  \item{x}{If \code{dataset} is missing, a numeric object giving the
    covariate values. If \code{dataset} is a data frame, it is
    interpreted as the name of the variable corresponding to the
    covariate in the data frame.
  }
  \item{t}{If \code{dataset} is missing, a numeric object giving the
    observed times. If \code{dataset} is a data frame, it is interpreted
    as the name of the variable corresponding to the observed times in
    the data frame.
  }
  \item{d}{If \code{dataset} is missing, an integer object giving the
    values of the uncensoring indicator. Censored observations must be
    coded as 0, uncensored ones as 1. If \code{dataset} is a data frame,
    it is interpreted as the name of the variable corresponding to the
    uncensoring indicator. 
  }
  \item{dataset}{An optional data frame in which the variables named in
    \code{x}, \code{t} and \code{d} are interpreted. If it is missing,
    \code{x}, \code{t} and \code{d} must be objects of the workspace.
  }
  \item{x0}{A numeric vector of covariate values where the survival
    estimates will be computed.
  }
  \item{h}{A numeric vector of bandwidths. If it is missing the default
    is to use the cross-validation bandwidth computed by the
    \code{berancv} function.
  }
  \item{local}{A logical value, \code{TRUE} by default, specifying
    whether local or global bandwidths are used.
  }
  \item{testimate}{A numeric vector specifying the times at which the
    survival is estimated. By default it is \code{NULL}, and then the
    survival is estimated at the times given by \code{t}.
  }
  \item{conflevel}{A value controlling whether bootstrap confidence
    intervals (CI) of the survival are to be computed. With the default
    value, 0L, the CIs are not computed. If a numeric value 
    between 0 and 1 is passed, it specifies the confidence level of the
    CIs.
  }
  \item{cvbootpars}{A list of parameters controlling the bootstrap when
    computing the CIs of the survival: \code{B}, the number of bootstrap
    resamples, and \code{nnfrac}, the fraction of the sample size that
    determines the order of the nearest neighbor used for choosing a
    pilot bandwidth. If \code{h} is missing the list of parameters is
    extended to be the same used for computing the cross-validation
    bandwidth (see the help of \code{berancv} for details). The default
    is the value returned by the \code{controlpars} function called
    without arguments. In case the CIs are not computed and \code{h} is
    not missing the default is \code{NULL}.
  }
}

\details{This function computes the kernel type product-limit estimator
  of the conditional survival function \eqn{S(t | x) = P(Y > t | X = x)}
  under censoring, using the Nadaraya-Watson weights. The kernel used is
  the Epanechnikov. If the smoothing parameter \eqn{h} is not provided,
  then the cross-validation bandwidth selector in Geerdens et al. (2018)
  is used. The function is available only for one continuous covariate
  \eqn{X}.}

\value{An object of S3 class 'npcure'. Formally, a list of components:
  \item{type}{The constant string "survival".
  }
  \item{local}{The value of the \code{local} argument.
  }
  \item{h}{The value of the \code{h} argument, unless this is missing,
    in which case its value is that of the cross-validation bandwidth.
  }
  \item{x0}{The value of the \code{x0} argument.
  }
  \item{testim}{The numeric vector of time values where the survival function
    is estimated.
  }
  \item{S}{A list whose components are the estimates of the survival function
    for each one of the covariate values, i.e., those specified by the
    \code{x0} argument. The survival estimates are given at the times
    determined by the \code{testimate} argument.
  }
}

\references{Beran, R. (1981). Nonparametric regression with randomly
  censored survival data. Technical report, University of California,
  Berkeley.

  Geerdens, C., Acar, E. F., Janssen, P. (2018). Conditional copula
  models for right-censored clustered event time data. Biostatistics,
  19(2): 247-262. \url{https://doi.org/10.1093/biostatistics/kxx034}.
}

\author{Ignacio López-de-Ullibarri [aut, cre],
  Ana López-Cheda [aut],
  Maria Amalia Jácome [aut]
}


\seealso{\code{\link{controlpars}}, \code{\link{berancv}}}

\examples{
## Some artificial data
set.seed(123)
n <- 50
x <- runif(n, -2, 2) ## Covariate values
y <- rweibull(n, shape = .5*(x + 4)) ## True lifetimes
c <- rexp(n) ## Censoring values
p <- exp(2*x)/(1 + exp(2*x)) ## Probability of being susceptible
u <- runif(n)
t <- ifelse(u < p, pmin(y, c), c) ## Observed times
d <- ifelse(u < p, ifelse(y < c, 1, 0), 0) ## Uncensoring indicator
data <- data.frame(x = x, t = t, d = d)

## Survival estimates for covariate values 0, 0.5 using...
## ... (a) global bandwidths 0.3, 0.5, 1.
## By default, the estimates are computed at the observed times
x0 <- c(0, .5)
S1 <- beran(x, t, d, data, x0 = x0, h = c(.3, .5, 1), local = FALSE) 

## Plot predicted survival curves for covariate value 0.5
plot(S1$testim, S1$S$h0.3$x0.5, type = "s", xlab = "Time", ylab =
"Survival", ylim = c(0, 1)) 
lines(S1$testim, S1$S$h0.5$x0.5, type = "s", lty = 2)
lines(S1$testim, S1$S$h1$x0.5, type = "s", lty = 3)
## The true survival curve is plotted for reference
p0 <- exp(2*x0[2])/(1 + exp(2*x0[2]))
lines(S1$testim, 1 - p0 + p0*pweibull(S1$testim, shape = .5*(x0[2] + 4),
lower.tail = FALSE), col = 2)
legend("topright", c("Estimate, h = 0.3", "Estimate, h = 0.5",
"Estimate, h = 1", "True"), lty = c(1:3, 1), col = c(rep(1, 3), 2))

## As before, but with estimates computed at fixed times 0.1, 0.2,...,1
S2 <- beran(x, t, d, data, x0 = x0, h = c(.3, .5, 1), local = FALSE,
testimate = .1*(1:10))

## ... (b) local bandwidths 0.3, 0.5.
## Note that the length of the covariate vector x0 and the bandwidth h
## must be the same.
S3 <- beran(x, t, d, data, x0 = x0, h = c(.3, .5), local = TRUE)

## ... (c) the cross-validation (CV) bandwidth selector (the default
## when the bandwidth argument is not provided). 
## The CV bandwidth is searched in a grid of 150 bandwidths (hl = 150)
## between 0.2 and 2 times the standardized interquartile range
## of the covariate values (hbound = c(.2, 2)).
## 95% confidence intervals are also given.
S4 <- beran(x, t, d, data, x0 = x0, conflevel = .95, cvbootpars =
controlpars(hl = 150, hbound = c(.2, 2))) 
     
## Plot of predicted survival curve and confidence intervals for
## covariate value 0.5 
plot(S4$testim, S4$S$x0.5, type = "s", xlab = "Time", ylab = "Survival",
ylim = c(0, 1))
lines(S4$testim, S4$conf$x0.5$lower, type = "s", lty = 2)
lines(S4$testim, S4$conf$x0.5$upper, type = "s", lty = 2)
lines(S4$testim, 1 - p0 + p0 * pweibull(S4$testim, shape = .5*(x0[2] +
4), lower.tail = FALSE), col = 2) 
legend("topright", c("Estimate with CV bandwidth", "95\% CI limits",
"True"), lty = c(1, 2, 1), col = c(1, 1, 2))

\donttest{
## Example with the dataset 'bmt' in the 'KMsurv' package
## to study the survival of patients aged 25 and 40.
data("bmt", package = "KMsurv")
x0 <- c(25, 40)
S <- beran(z1, t2, d3, bmt, x0 = x0, conflevel = .95)
## Plot of predicted survival curves and confidence intervals
plot(S$testim, S$S$x25, type = "s", xlab = "Time", ylab = "Survival",
ylim = c(0, 1))
lines(S$testim, S$conf$x25$lower, type = "s", lty = 2)
lines(S$testim, S$conf$x25$upper, type = "s", lty = 2)
lines(S$testim, S$S$x40, type = "s", lty = 1, col = 2)
lines(S$testim, S$conf$x40$lower, type = "s", lty = 2, col = 2)
lines(S$testim, S$conf$x40$upper, type = "s", lty = 2, col = 2)
legend("topright", c("Age 25: Estimate", "Age 25: 95\% CI limits",
"Age 40: Estimate", "Age 40: 95\% CI limits"), lty = 1:2,
col = c(1, 1, 2, 2)) 
}
}

\keyword{survival}

\keyword{nonparametric}