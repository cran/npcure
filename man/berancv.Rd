\name{berancv}

\alias{berancv}

\title{Compute the Cross-Validation Bandwidth for Beran's Estimator
  of the Conditional Survival
}

\description{This function computes the cross-validation bandwidth
  for Beran's estimator of the conditional survival function.
}

\usage{berancv(x, t, d, dataset, x0, cvpars = controlpars())
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
    \code{x}, \code{t} and \code{indicator} are interpreted. If it is
    missing, \code{x}, \code{t} and \code{indicator} must be objects
    of the workspace.
  }
  \item{x0}{A numeric vector of covariate values where the local
    cross-validation bandwidth will be computed.
  }
  \item{cvpars}{A list of parameters controlling the process of bandwidth
    selection. The default is the value returned by the \code{controlpars}
    function called without arguments. See the help for
    \code{controlpars} for details.
  }
}

\details{The cross-validation (CV) bandwidth is taken as the largest
  local minimizer of the leave-one-out cross-validated criterion in
  Geerdens et al. (2018). Let \eqn{F^{(-i)}(t | x_i)}, \eqn{i = 1,
  \ldots, n}{i = 1, \dots, n} be the Beran estimator obtained using the
  data points \eqn{(x_j, t_j, d_j)}, \eqn{j = 1, \ldots, i-1, i+1,
  \ldots, n}{j = 1, \dots, i-1, i+1, \dots, n}. For the CV criterion,
  the differences \eqn{I(t_i \leq t_j)-F^{(-i)}(t_j|x_i)}{I(t_i <=
  t_j)-F^{(-i)}(t_j | x_i)} are computed only for the so-called 'useful
  pairs' of observed times \eqn{(t_i, t_j)}. A pair \eqn{(T_i, T_j)} is
  useful if the value of the indicator \eqn{I(T_i \leq T_j)}{I(T_i <=
  T_j)} gives an unambiguous correct value for the indicator \eqn{I(Y_i
  \leq Y_j)}{I(Y_i <= Y_j)} which contains the corresponding true
  (possibly unknown) event times, see Geerdens et al. (2018) for
  details. Gannoun et al. (2007) apply a similar criterion to perform
  bandwidth selection for the Beran estimator, but they consider only
  the pairs of true (uncensored) event times. Note that the inclusion of
  useful pairs of observed times would be especially advantageous if the
  censoring rate is high.
}

\value{An object of S3 class 'npcure'. Formally, a list of components:
  \item{type}{The constant character string
  c("Cross-validation bandwidth", "survival").
  }
  \item{x0}{Grid of covariate values.
  }
  \item{h}{Selected local cross-validation bandwidths.
  }
  \item{hgrid}{Grid of bandwidths used (optional).
  }
}

\references{Gannoun, A., Saracco, J., Yu, K. (2007). Comparison of
  kernel estimators of conditional distribution function and quantile
  regression under censoring. Statistical Modeling, 7:
  329-344. \url{https://doi.org/10.1177/1471082X0700700404}.

  Geerdens, C., Acar, E. F., Janssen, P. (2018). Conditional copula
  models for right-censored clustered event time data. Biostatistics,
  19(2): 247-262. \url{https://doi.org/10.1093/biostatistics/kxx034}.
}

\author{Ignacio López-de-Ullibarri [aut, cre],
  Ana López-Cheda [aut],
  Maria Amalia Jácome [aut]
}

\seealso{\code{\link{beran}}, \code{\link{controlpars}},
  \code{\link{hpilot}}
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

## Computation of cross-validation (CV) local bandwidth for Beran's
## estimator of survival for covariate values 0, 1, ...
#### ... with the default control parameters (passed through 'cvpars')
x0 <- c(0, 1)
hcv <- berancv(x, t, d, data, x0 = x0)

#### ... changing the default 'cvpars' by calling the 'controlpars()'
#### function: 
#### (a) the CV local bandwidth is searched in a grid of 150 bandwidths
#### ('hl = 150') between 0.2 and 4 times the standardized interquartile
#### range of the covariate values of x ('hbound = c(.2, 4'))
#### (b) all the grid bandwidths are saved ('hsave = TRUE')
hcv <- berancv(x, t, d, data, x0 = x0, cvpars = controlpars(hbound =
c(.2, 4), hl = 150, hsave = TRUE))

## Survival estimates for covariate values 0, 1, with CV local bandwidth
S1 <- beran(x, t, d, data, x0 = x0,  h = hcv$h)
## Plot predicted survival curves for covariate values 0, 1
plot (S1$testim, S1$S$x0, type = "s", xlab = "Time", ylab = "Survival",
ylim = c(0, 1))
lines(S1$testim, S1$S$x1, type = "s", lty = 2)
## The survival curves are displayed for reference
p0 <- exp(2*x0)/(1 + exp(2*x0))
lines(S1$testim, 1 - p0[1] + p0[1]*pweibull(S1$testim, shape = .5*(x0[1]
+ 4), lower.tail = FALSE), col = 2)
lines(S1$testim, 1 - p0[2] + p0[2]*pweibull(S1$testim, shape = .5*(x0[2]
+ 4), lower.tail = FALSE), lty = 2, col = 2)
legend("topright", c("Estimate, x = 0", "True, x = 0",
"Estimate, x = 1", "True, x = 1"), lty = c(1, 1, 2, 2), col = 1:2)

\donttest{
## Example with the dataset 'bmt' of the 'KMsurv' package to study the
## survival of patients aged 25 and 40.
data("bmt", package = "KMsurv")
x0 <- c(25, 40)
hcv <- berancv(z1, t2, d3, bmt, x0 = x0, cvpars = controlpars(hbound =
c(.2, 4), hl = 150, hsave = TRUE))
S <- beran(z1, t2, d3, bmt, x0 = x0, h = hcv$h, conflevel = .95)
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