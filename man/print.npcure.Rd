\name{print.npcure}

\alias{print.npcure}

\title{Print Method for Objects of Class 'npcure'
}

\description{This function implements a print method for 'npcure'
  objects.
} 

\usage{\method{print}{npcure}(x, how, head = FALSE, \dots)
}

\arguments{
  \item{x}{An object of class 'npcure'.
  }
  \item{how}{A character string with values "wide" or "long". If
    missing, the function itself chooses a convenient default.
  }
  \item{head}{A logical value that controls whether the function's
    output must be abbreviated (\code{TRUE}) or not (\code{FALSE},
    the default). 
  }
  \item{\dots}{Further optional arguments. Excepting for \code{n}, which
  controls how many lines are printed when \code{head = TRUE}, these are
  the arguments for the default method (i.e., \code{print.default}) of
  the \code{print} generic function.
  }
}
    
\value{A formatted output.}

\author{Ignacio López-de-Ullibarri [aut, cre],
  Ana López-Cheda [aut],
  Maria Amalia Jácome [aut]
}

\seealso{\code{\link{summary.npcure}}}

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

## Calling 'print()' with an object of class 'npcure' created by
## 'latency()' 
S1 <- latency(x, t, d, data, x0 = c(0, .5), h = c(1, 1.5))

## In this case (latency estimation with local bandwidths and without
## confidence bands), the 'wide' format is used by default
S1
print(S1, how = "wide")
print(S1, how = "long")

## How to control the number of significant digits of the output, and
## how to abbreviate the output 
print(S1, digits = 5, head = TRUE, n = 4)

## Calling 'print()' with a 'npcure' object created by 'probcure()'
q1 <- probcure(x, t, d, data, x0 = c(0, .5), h = c(.5, 1, 1.5), local =
FALSE, conflevel = .95)

## Only the 'long' format is available when confidence bands are
## computed
q1
print(q1, how = "long")
print(q1, how = "wide")
}
