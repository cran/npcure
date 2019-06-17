\name{testcov}

\alias{testcov}

\title{Covariate Significance Test of the Cure Probability}

\description{This function carries out a significance test of covariate
  effect on the probability of cure.
}

\usage{testcov(x, t, d, dataset, bootpars = controlpars())
}

\arguments{
  \item{x}{If \code{dataset} is missing, an object giving the covariate
    values, whose type can be numeric, integer, factor or character.
    If\code{dataset} is a data frame, it is interpreted as the name of
    the variable corresponding to the covariate in the data frame.
  }
  \item{t}{If \code{dataset} is missing, a numeric object giving the
    observed times. If \code{dataset} is a data frame, it is
    interpreted as the name of the variable corresponding to the
    observed times in the data frame.
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
  \item{bootpars}{A list of parameters controlling the test. Currently,
    the only accepted component is \code{B}, the number of bootstrap
    resamples. The default is the value returned by the
    \code{controlpars} function called without arguments.
  }
}

\details{The function computes a statistic, based on the method by
  Delgado and González-Manteiga (2004), to test whether a covariate X
  has an effect on the cure probability (\eqn{H_1}: cure probability =
  q(x)) or not (\eqn{H_0}: cure probability = q). Since the cure rate
  can be written as the regression function \eqn{E(\nu | X = x) =
  q(x)}{E(nu | X = x) = q(x)}, where \eqn{\nu}{nu} is the cure
  indicator, the procedure is carried out as a significance test in
  nonparametric regression. The challenge of the test is that the cure
  indicator \eqn{\nu}{nu} is only partially observed due to censoring:
  for the uncensored observations \eqn{\nu = 0}{nu = 0}, but it is
  unknown if a censored individual will be eventually cured or not
  (\eqn{\nu}{nu} unknown). The approach consists in expressing the cure
  rate as a regression function with another response, not observable
  but estimable, say \eqn{\eta}{eta}. The estimated values of
  \eqn{\eta}{eta} depend on suitable estimates of the conditional 
  distribution of the censoring variable and \eqn{\tau(x)}{tau(x)}, an
  unknown time beyond which a subject could be considered as cured; see
  López-Cheda (2018). For the computation of the values of
  \eqn{\eta}{eta}, the censoring distribution is estimated
  unconditionally using the Kaplan-Meier product-limit estimator. The
  time \eqn{\tau}{tau} is estimated by the largest uncensored time.

  The test statistic is a weighted mean of the difference between the
  observations of \eqn{\eta_i}{eta_i} and the values of the conditional
  mean of \eqn{\eta}{eta} under the null hypothesis:
  
  \deqn{T_{n(x)} = 1/n \sum (\eta_i - \bar{\eta}) I(x_i \leq x)}{T_{n(x)
  = 1/n sum(eta_i - mean(eta)) I(x_i <= x)}}.

  For a qualitative covariate there is no natural way to order the values
  \eqn{x_i}. In principle, this makes impossible to compute the indicator
  function in the test statistic. The problem is solved by considering all
  the possible combinations of the covariate values, and by computing the
  test statistic for each 'ordered' combination. \eqn{T_{n(x)}} is taken
  as the largest value of these test statistics.
  
  The distribution of the test statistic under the null hypothesis is
  approximated by bootstrap, using independent naive resampling.
}

\value{An object of S3 class 'npcure'. Formally, a list of components:
  \item{type}{The constant character string c("test", "Covariate").
  }
  \item{x}{The name of the covariate.
  }
  \item{CM}{The result of the Cramer-von Mises test: a list with
    components \code{statistic}, the test statistic, and \code{pvalue},
    the  p-value.
  }
  \item{KS}{The result of the Kolmogorov-Smirnov test: a list with
    components \code{statistic}, the test statistic, and \code{pvalue},
    the p-value.
  }
}

\references{Delgado M. A., González-Manteiga W. (2001). Significance
  testing in nonparametric regression based on the bootstrap. Annals of
  Statistics, 29: 1469-1507.
  
  López-Cheda A. (2018). \emph{Nonparametric Inference in Mixture Cure
    Models}. PhD dissertation, Universidade da Coruña. Spain.
}

\author{Ignacio López-de-Ullibarri [aut, cre],
  Ana López-Cheda [aut],
  Maria Amalia Jácome [aut]
}

\seealso{\code{\link{controlpars}}}

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

## Test of the significance of the covariate 'x'
testcov(x, t, d, data)

## Test carried out with 1999 bootstrap resamples (the default is 999)
testcov(x, t, d, data, bootpars = controlpars(B = 1999))

\donttest{
## How to apply the test repeatedly when there is more than one
## covariate... 
## ... 'y' is another continuous covariate of the data frame 'data'
data$y <- runif(n, -1, 1)
namecovar <- c("x", "y")
## ... testcov is called from a 'for' loop
for (i in 1:length(namecovar)) {
   result <- testcov(data[, namecovar[i]], data$t, data$d)
   print(result)
}

## In the previous example, testcov() was called without using the
## argument 'dataset'. To use it, the 'for' must be avoided
testcov(x, t, d, data)
testcov(y, t, d, data)

## Non-numeric covariates can also be tested...
## ... 'z' is a nominal covariate of the data frame 'data'
data$z <- rep(factor(letters[1:5]), each = n/5)
testcov(z, t, d, data)

## Example with the dataset 'bmt' in the 'KMsurv' package
## to study the effect on the probability of cure of...
## ... (a) a continuous covariate (z1 = age of the patient)
data("bmt", package = "KMsurv")
set.seed(1) ## Not needed, just for reproducibility.
testcov(z1, t2, d3, bmt, bootpars = controlpars(B = 4999))

## ... (b) a qualitative covariate (z3 = patient gender)
set.seed(1) ## Not needed, just for reproducibility.
testcov(z3, t2, d3, bmt, bootpars = controlpars(B = 4999))
}
}

\keyword{survival}

\keyword{nonparametric}