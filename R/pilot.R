## pilot bandwidth for bootstrap bandwidth selector
hpilot <- function(x,
                   x0,
                   nnfrac = 0.25) {
    k <- floor(length(x)*nnfrac)
    x <- sort(x)
    n <- length(x)
    pilot <- numeric(length = length(x0))
    for (i in 1:length(x0)) {
        xi <- x0[i]
        nNNleft <- sum(x < xi) 
        nNNright <- sum(x > xi)
        if (nNNleft >= k && nNNright >= k) pilot[i] <- (x[x > xi][k] - rev(x[x < xi])[k])*0.5*(100/n)^(1/9)
        else {
            if (nNNleft >= k) pilot[i] <- (x[n] - rev(x[x < xi])[k])*0.5*(100/n)^(1/9)
            else {
                if (nNNright >= k) pilot[i] <- (x[x > xi][k] - x[1])*0.5*(100/n)^(1/9)
                else pilot[i] <- (x[n] - x[1])*(100/n)^(1/9)
            }
        }
    }
    pilot
}
