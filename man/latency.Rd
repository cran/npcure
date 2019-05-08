\name{latency}

\alias{latency}

\title{Compute Nonparametric Estimator of the Conditional Latency}

\description{This function computes the nonparametric estimator of the
  conditional latency function proposed by López-Cheda et al. (2017).
}

\usage{latency(x, t, d, dataset, x0, h, local = TRUE, testimate = NULL,
conflevel = 0L, bootpars = if (conflevel == 0 && !missing(h)) NULL else
npcure::controlpars())
}

\arguments{
  \item{x}{If \code{dataset} is missing, a numeric object giving the
    covariate values. If \code{dataset} is a data frame, it is interpreted
    as the name of the variable corresponding to the covariate in the data
    frame.
  }
  \item{t}{If \code{dataset} is missing, a numeric object giving the
    observed times. If \code{dataset} is a data frame, it is interpreted
    as the name of the variable corresponding to the observed times in the
    data frame.
  }
  \item{d}{If \code{dataset} is missing, an integer object giving the
  values of the uncensoring indicator. Censored observations must be
  coded as 0, uncensored ones as 1. If\code{dataset} is a data frame, it
  is interpreted as the name of the variable corresponding to the
  uncensoring indicator.
}
  \item{dataset}{An optional data frame in which the variables named in
    \code{x}, \code{t} and \code{d} are interpreted. If it is missing,
    \code{x}, \code{t} and \code{d} must be objects of the workspace.
  }
  \item{x0}{A numeric vector of covariate values where the latency
    estimates will be computed.
  }
  \item{h}{A numeric vector of bandwidths. If it is missing the default
    is to use the local bootstrap bandwidth computed by the
    \code{latencyhboot} function.
  }
  \item{local}{A logical value, \code{TRUE} by default, specifying
    whether local or global bandwidths are used.
  }
  \item{testimate}{A numeric vector specifying the times at which the
    latency is estimated. By default it is \code{NULL}, and then the
    latency is estimated at the times given by \code{t}.
  }
  \item{conflevel}{A value controlling whether bootstrap confidence
    intervals (CI) of the latency are to be computed. With the default
    value, 0L, the CIs are not computed. If a numeric value 
    between 0 and 1 is passed, it specifies the confidence level of the
    CIs.
  }
  \item{bootpars}{A list of parameters controlling the bootstrap when
    computing the CIs of the latency: \code{B}, the number of bootstrap
    resamples, and \code{nnfrac}, the fraction of the sample size that
    determines the order of the nearest neighbor used for choosing a
    pilot bandwidth. If \code{h} is missing the list of parameters is
    extended to be the same used for computing the bootstrap bandwidth
    (see the help of \code{latencyhboot} for details). The default is
    the value returned by the \code{controlpars} function called without
    arguments. In case the CIs are not computed and \code{h} is not
    missing the default is \code{NULL}. 
  }
}

\details{The function computes the nonparametric estimator of the
  conditional latency \eqn{S_0(t | X = x_0) = P(Y>t | Y<\infty,
  X=x_0)}{S_0(t | X = x_0) = P(Y > t | Y < infinity, X = x_0)}
  proposed by López-Cheda et al. (2017). It is only available for a
  continuous covariate \eqn{X}.
}

\value{An object of S3 class 'npcure'. Formally, a list of components:
  \item{type}{The constant string "latency".
  }
  \item{local}{The value of the \code{local} argument.
  }
  \item{h}{The value of the \code{h} argument, unless this is missing,
    in which case its value is that of the cross-validation bandwidth.
  }
  \item{x0}{The value of the \code{x0} argument.
  }
  \item{testim}{The numeric vector of time values where the latency
    function is estimated.
  }
  \item{S}{A list whose components are the estimates of the latency
    function for each one of the covariate values, i.e., those specified
    by the \code{x0} argument. The latency estimates are given at the
    times determined by the \code{testimate} argument.
  }
  \item{conf}{A list of components \code{lower} and \code{upper} giving
    the lower and the upper limits of the confidence intervals,
    respectively.}
}

\references{López-Cheda, A., Jácome, M. A., Cao,
  R. (2017). Nonparametric latency estimation for mixture cure
  models. Test, 26:
  353–376. \url{https://doi.org/10.1007/s11749-016-0515-1}.
}

\author{Ignacio López-de-Ullibarri [aut, cre],
  Ana López-Cheda [aut],
  Maria Amalia Jácome [aut]
}

\seealso{\code{\link{controlpars}}, \code{\link{latencyhboot}}
}

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

## Latency estimates for covariate value 0.5...
x0 <- .5

## ... (a) with global bandwidths 0.5, 1, 2.
## By default, estimates are computed at the time values of 't'
S1 <- latency(x, t, d, data, x0 = x0, h = c(.5, 1, 2), local = FALSE)
plot(S1$testim, S1$S$h0.5$x0.5, type = "S", xlab = "Time", ylab =
"Latency", ylim = c(0, 1))
lines(S1$testim, S1$S$h1$x0.5, type = "S", lty = 2)
lines(S1$testim, S1$S$h2$x0.5, type = "S", lty = 3)
## The true latency curve is plotted as reference
lines(S1$testim, pweibull(S1$testim, shape = .5*(x0 + 4), lower.tail =
FALSE), lwd = 2)
legend("topright", c(paste("Estimate, ", c("h = 0.5", "h = 1",
"h = 2")), "True"), lty = c(1:3, 1), lwd = c(rep(1, 3), 2))

## As before, but with estimates computed at times 0.1, 0.2,..., 1
S2 <- latency(x, t, d, data, x0 = x0, h = c(.5, 1, 2), local = FALSE,
testimate = .1*(1:10))

## ... (b) with local bandwidth 2.
S3 <- latency(x, t, d, data, x0 = x0, h = 2, local = TRUE)
#### Note that with only one covariate value the results with
#### 'local = FALSE' and 'local = TRUE' coincide, but the output formats
#### differ slightly. Compare with
S3 <- latency(x, t, d, data, x0 = x0, h = 2, local = FALSE)

## ... (c) with local bootstrap bandwidth
b  <- latencyhboot(x, t, d, data, x0 = x0)
S4 <- latency(x, t, d, data, x0 = x0, h = b$h)

## ... (d) when the bandwidth is not specified, the bootstrap bandwidth
#### selector given by the 'latencyhboot' function is used by default.
#### The computation of 95% confidence bands based on 1999 bootstrap
#### resamples is also illustrated
S5 <- latency(x, t, d, data, x0 = x0, conflevel = .95, bootpars =
controlpars(B = 1999))
plot(S5$testim, S5$S$x0, type = "S", xlab = "Time", ylab = "Latency",
ylim = c(0, 1))
lines(S5$testim, S5$conf$x0$lower, type = "S", lty = 2)
lines(S5$testim, S5$conf$x0$upper, type = "S", lty = 2)
lines(S5$testim, pweibull(S5$testim, shape = .5*(x0 + 4), lower.tail =
FALSE), lwd = 2)
legend("topright", c("Estimate", "95\% CI limits", "True"), lty = c(1,
2, 1), lwd = c(1, 1, 2))

\donttest{
## Example with the dataset 'bmt' of the 'KMsurv' package
## to study the survival of the uncured patients aged 25 and 40
data("bmt", package = "KMsurv")
x0 <- c(25, 40)
S <- latency(z1, t2, d3, bmt, x0 = x0, conflevel = .95)
## Plot of predicted latency curves and confidence bands
plot(S$testim, S$S$x25, type = "S", xlab = "Time", ylab = "Latency",
ylim = c(0,1))
lines(S$testim, S$conf$x25$lower, type = "S", lty = 2)
lines(S$testim, S$conf$x25$upper, type = "S", lty = 2)
lines(S$testim, S$S$x40, type = "S", lty = 1, lwd = 2)
lines(S$testim, S$conf$x40$lower, type = "S", lty = 2, lwd = 2)
lines(S$testim, S$conf$x40$upper, type = "S", lty = 2, lwd = 2)
legend("topright", c("25 years", "95\% CI limits", "40 years",
"95\% CI limits"), lty = c(1, 2, 1, 2), lwd = c(1, 1, 2, 2))
}
}

\keyword{survival}

\keyword{nonparametric}