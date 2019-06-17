\name{probcure}

\alias{probcure}

\title{Compute Nonparametric Estimator of the Conditional Probability of
  Cure
}

\description{This function computes the nonparametric estimator of the
  conditional probability of cure proposed by Xu and Peng (2014).
}

\usage{probcure(x, t, d, dataset, x0, h, local = TRUE, conflevel = 0L,
bootpars = if (conflevel == 0 && !missing(h)) NULL else controlpars())
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
    uncensoring indicator in the data frame.
  }
  \item{dataset}{An optional data frame in which the variables named in
    \code{x}, \code{t} and \code{d} are interpreted. If it is missing,
    \code{x}, \code{t} and \code{d} must be objects of the workspace.
  }
  \item{x0}{A numeric vector of covariate values where the estimates
    of cure probability will be computed.
  }
  \item{h}{A numeric vector of bandwidths. If it is missing the default
  is to use the local bootstrap bandwidth computed by the
  \code{probcurehboot} function.
}
  \item{local}{A logical value, \code{TRUE} by default, specifying
    whether local or global bandwidths are used.
  }
  \item{conflevel}{A value controlling whether bootstrap confidence
    intervals (CI) of the cure probability are to be computed. With the
    default value, 0L, the CIs are not computed. If a numeric value
    between 0 and 1 is passed, it specifies the confidence level of the
    CIs.
  }
  \item{bootpars}{A list of parameters controlling the bootstrap when
    computing either the CIs of the cure probability or the bootstrap
    bandwidth (if \code{h} is missing): \code{B}, the number of bootstrap
    resamples, and \code{nnfrac}, the fraction of the sample size that
    determines the order of the nearest neighbor used for choosing a pilot
    bandwidth. The default is the value returned by the \code{controlpars}
    function called without arguments. If the CIs are not computed and
    \code{h} is not missing the default is \code{NULL}.
  } 
}

\details{The function computes the nonparametric estimator of the
  conditional cure probability \eqn{q(X = x_0) \equiv 1 - p(X =
  x_0)= P(Y=\infty | X=x_0)}{q(X = x_0) := 1 - p(X = x_0)= P(Y =
  infinity | X = x_0)} proposed by Xu and Peng (2014), and also studied
  by López-Cheda et al (2017). It is only available for a continuous
  covariate \eqn{X}.
}

\value{An object of S3 class 'npcure'. Formally, a list of components:
  \item{estimate}{The constant string "cure".
  }
  \item{local}{The value of the \code{local} argument.
  }
  \item{h}{The value of the \code{h} argument, unless this is missing,
    in which case its value is that of the bootstrap bandwidth.
  }
  \item{x0}{The value of the \code{x0} argument.
  }
  \item{q}{A list with the estimates of the probability of cure.
  }
  \item{conf}{A list of components \code{lower} and \code{upper} giving
    the lower and the upper limits of the confidence intervals,
    respectively.
  }
}

\references{López-Cheda, A., Cao, R., Jácome, M. A., Van Keilegom,
  I. (2017). Nonparametric incidence estimation and bootstrap bandwidth
  selection in mixture cure models. Computational Statistics & Data
  Analysis, 105:
  144–165. \url{https://doi.org/10.1016/j.csda.2016.08.002}.
  
  Xu, J., Peng, Y. (2014). Nonparametric cure rate estimation with
  covariates. The Canadian Journal of Statistics 42:
  1-17. \url{https://doi.org/10.1002/cjs.11197}.
}

\author{Ignacio López-de-Ullibarri [aut, cre],
  Ana López-Cheda [aut],
  Maria Amalia Jácome [aut]
}

\seealso{\code{\link{controlpars}}, \code{\link{probcurehboot}}
}

\examples{
## Some artificial data
set.seed(123)
n <- 50
x <- runif(n, -2, 2) ## Covariate values
y <- rweibull(n, shape = 0.5 * (x + 4)) ## True lifetimes
c <- rexp(n) ## Censoring values
p <- exp(2*x)/(1 + exp(2*x)) ## Probability of being susceptible
u <- runif(n)
t  <- ifelse(u < p, pmin(y, c), c) ## Observed times
d  <- ifelse(u < p, ifelse(y < c, 1, 0), 0) ## Uncensoring indicator
data <- data.frame(x = x, t = t, d = d)

## Covariate values where cure probability is estimated 
x0 <- seq(-1.5, 1.5, by = 0.1) 

## Nonparametric estimates of cure probability at 'x0'...

## ... (a) with global bandwidths 1, 1.5, 2  
q1 <- probcure(x, t, d, data, x0 = x0, h = c(1, 1.5, 2), local = FALSE)

#### Plot predicted cure probabilities at 'x0' for each bandwidth in 'h'
#### (the true cure probability is displayed for reference)
plot(q1$x0, q1$q$h1, type = "l", xlab = "Covariate", ylab =
"Cure probability", ylim = c(0, 1))
lines(q1$x0, q1$q$h1.5, lty = 2)
lines(q1$x0, q1$q$h2, lty = 3)
lines(q1$x0, 1 - exp(2*q1$x0)/(1 + exp(2*q1$x0)), col = 2)
legend("topright", c(paste("Estimate, ", c("h = 1", "h = 1.5",
"h = 2")), "True"), lty = c(1, 2, 3, 1), col = c(1, 1, 1, 2))

## ... (b) with local bandwidths (default)
#### (the vectors passed to 'x0' and 'h' must have the same length)
q2 <- probcure(x, t, d, data, x0 = x0, h = seq(1, 2.5, along = x0)) 

## ... (c) with local bootstrap bandwidths (based on 1999 booostrap
#### resamples). Besides, 95% confidence intervals are computed  and
#### smoothed (with a 15-th order moving average)
set.seed(1) ## Not needed, just for reproducibility.
hb <- probcurehboot(x, t, d, data, x0 = x0, bootpars = controlpars(B =
1999, hsmooth = 15))
q3 <- probcure(x, t, d, data, x0 = x0, h = hb$hsmooth, conflevel = .95,
bootpars = controlpars(B = 1999))

\donttest{
## ... (d) If the bandwidth is not specified, the local bootstrap
#### bandwidth is used (same results as in (c))
set.seed(1) ## Not needed, just for reproducibility.
q4 <- probcure(x, t, d, data, x0 = x0, conflevel = .95, bootpars =
controlpars(B = 1999, hsmooth = 15))

#### Plot of the estimated cure probabilities evaluated at 'x0'
#### (true cure rate displayed as reference)
plot (q4$x0, q4$q, type = "l", ylim = c(0, 1), xlab = "Covariate X",
ylab = "Cure probability")
lines(q4$x0, q4$conf$lower, lty = 2)
lines(q4$x0, q4$conf$upper, lty = 2)
lines(q4$x0, 1-exp(2 * q4$x0)/(1 + exp(2 * q4$x0)), col = 2)
legend("topright", c("Estimate", "95\% CI limits", "True"),
lty = c(1, 2, 1), col = c(1, 1, 2))

## Example with the dataset 'bmt' in the 'KMsurv' package
## to study the probability of cure as a function of the age (z1).
data("bmt", package = "KMsurv")
x0 <- seq(quantile(bmt$z1, .05), quantile(bmt$z1, .95), length.out = 100)
q.age <- probcure(z1, t2, d3, bmt, x0 = x0, conflevel = .95, bootpars =
controlpars(B = 1999, hsmooth = 10)) 

## Plot of estimated cure probability and confidence intervals
par(mar = c(5, 4, 4, 5) + .1)
plot(q.age$x0, q.age$q, type = "l", ylim = c(0, 1), xlab =
"Patient age (years)", ylab = "Cure probability")
lines(q.age$x0, q.age$conf$lower, lty = 2)
lines(q.age$x0, q.age$conf$upper, lty = 2)
## The estimated density of age (z1) is added for reference
par(new = TRUE)
d.age <- density(bmt$z1)
plot(d.age, xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = 2,
main = "", zero.line = FALSE)
mtext("Density", side = 4, col = 2, line = 3) 
axis(4, ylim = c(0, max(d.age$y)), col = 2, col.axis = 2)
legend("topright", c("Estimate", "95\% CI limits", "Covariate density"),
lty = c(1, 2, 1), col = c(1, 1, 2))
}
}

\keyword{survival}

\keyword{nonparametric}