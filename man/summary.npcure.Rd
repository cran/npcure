\name{summary.npcure}

\alias{summary.npcure}

\title{Summary Method for Objects of Class 'npcure'
}

\description{This function prints a summary of a 'npcure' object.
}

\usage{\method{summary}{npcure}(object, \dots)
}

\arguments{
  \item{object}{An object of class 'npcure'.
  }
  \item{\dots}{Further optional arguments for the default method
    (i.e., \code{summary.default}) of the \code{summary} generic
    function. 
  }
}

\value{A compact summary showing the components of the object.
}

\author{Ignacio López-de-Ullibarri [aut, cre],
  Ana López-Cheda [aut],
  Maria Amalia Jácome [aut]
}

\seealso{\code{\link{print.npcure}}}

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

## Calling 'summary()' with an object of class 'npcure' created by
## 'latency()' 
S1 <- latency(x, t, d, data, x0 = c(0, .5), h = c(1, 1.5))
summary(S1)

## If needed, the number of significant digits of the output can be set
summary(S1, digits = 5)

## Calling 'summary()' with an object created by 'probcure()'
q1 <- probcure(x, t, d, data, x0 = c(0, .5), h = c(.5, 1, 1.5), local =
FALSE, conflevel = .95)
summary(q1)
}
