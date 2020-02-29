## summary method
summary.npcure <- function(object, ...) {
    dots <- list(...)
    cat(paste("An object of class '", class(object), "'\n\n", sep = ""))
    if (is.null(dots$give.attr))
        str(object, give.attr = FALSE, ...)
    else
        str(object, ...)
}

## print method
print.npcure <- function(x, how, head = FALSE, ...) {
    dots <- list(...)
    if (is.null(dots$digits))
        dots$digits <- getOption("digits")
    if (head && is.null(dots$n))
        dots$n <- 6
    if (length(x$type) == 1) {
        if (!missing(how)) {
            if (how == "wide" && !is.null(x$conf))
                warning("The option 'how = 'long'' is the only available for 'npcure' objects with 'type' component equal to 'cure', 'latency' or 'survival' and a 'conf' component")
            if (how == "long" && x$local && x$type == "cure")
                warning("The option 'how = 'wide'' is the only available for 'npcure' objets with components 'local' and 'type' equal to 'TRUE' and 'cure', respectively")
        }
        cat("\nBandwidth type: ")
        if (x$local)
            cat("local\n\n")
        else
            cat("global\n\n")
        if (x$type == "cure") {
            cat("Conditional", x$type, "estimate:\n")
            if (head && length(x$x0) < dots$n)
                head <- FALSE
            if (x$local) { ## local
                if (is.null(x$conf)) { ## without CI
                    df <- as.data.frame(x[c("h", "x0", "q")])
                    dimnames(df)[[2]][3] <- x$type
                }
                else { ## with CI
                    df <- as.data.frame(x[c("h", "x0", "q", "conf")])
                    dimnames(df)[[2]][3:5] <- c(x$type, sprintf("lower %s%% CI", formatC(x$conflevel*100)), sprintf("upper %s%% CI", formatC(x$conflevel*100)))                
                }
                if (head) {
                    print(df[1:dots$n, ], row.names = FALSE, digits = dots$digits)
                    cat("...\n")
                }
                else {
                    print(df, row.names = FALSE, digits = dots$digits)
                }
            }
            else { ## cure, global
                if (is.null(x$conf)) { ## without CI
                    if (!missing(how) && how == "long") {## long
                        for (i in 1:length(x$h)) {
                            cat("\nh =", format(x$h[i], digits = dots$digits), "\n")
                            df <- data.frame(x["x0"], x$q[i])
                            dimnames(df)[[2]][2] <- x$type
                            if (head) {
                                print(df[1:dots$n, ], row.names = FALSE, digits = dots$digits)
                                cat("...\n")
                            }
                            else {
                                print(df, row.names = FALSE, digits = dots$digits)
                            }
                        }
                    }
                    else { ## wide
                        df <- data.frame(x["x0"]) 
                        for (i in 1:length(x$h))
                            df <- cbind(df, x$q[i])
                        dimnames(df)[[2]] <- c("x0", formatC(paste("h = ", x$h, sep = ""), getOption("digits") + 1))
                        if (head) {
                            print(df[1:dots$n, ], row.names = FALSE, digits = dots$digits)
                            cat("...\n")
                        }
                        else {
                            print(df, row.names = FALSE, digits = dots$digits)
                        }
                    }
                }
                else { ## long by default
                    for (i in 1:length(x$h)) {
                        cat("\nh =", format(x$h[i], digits = dots$digits), "\n")
                        df <- data.frame(x["x0"], x$q[i], x$conf[[i]])
                        dimnames(df)[[2]] <- c("x0", x$type, sprintf("lower %s%% CI", formatC(x$conflevel*100)), sprintf("upper %s%% CI", formatC(x$conflevel*100)))
                        if (head) {
                            print(df[1:dots$n, ], row.names = FALSE, digits = dots$digits)
                            cat("...\n")
                        }
                        else {
                            print(df, row.names = FALSE, digits = dots$digits)
                        }
                    }
                }
            }
        }
        else { ## latency or beran
            cat("Covariate (x0):", format(x$x0, digits = dots$digits))
            cat("\nBandwidth (h): ", format(x$h, digits = dots$digits), "\n\n")
            if (x$type == "latency")
                cat("Conditional ")
            else
                cat("Beran's conditional ")
            cat(x$type, "estimate:\n")
            if (head && length(x$testim) < dots$n)
                head <- FALSE
            if (x$local) { ## latency or beran, local 
                if (is.null(x$conf)) {
                    if (!missing(how) && how == "long") { ## long
                        for (i in 1:length(x$h)) {
                            cat("\nx0 =", format(x$x0[i], digits = dots$digits))
                            cat("\nh =", format(x$h[i], digits = dots$digits), "\n")
                            if (length(x$h) == 1)
                                df <- data.frame(x["testim"], x$S)
                            else
                                df <- data.frame(x["testim"], x$S[i])
                            dimnames(df)[[2]] <- c("time", x$type)
                            if (head) {
                                print(df[1:dots$n, ], row.names = FALSE, digits = dots$digits)
                                cat("...\n")
                            }
                            else {
                                print(df, row.names = FALSE, digits = dots$digits)
                            }
                        }
                    }
                    else { ## wide
                        df <- data.frame(x["testim"])
                        for (i in 1:length(x$h))
                            df <- cbind(df, x$S[i])
                        dimnames(df)[[2]] <- c("time", formatC(paste("x0 = ", x$x0, sep = ""), getOption("digits") + 1))
                        cat("\n")
                        if (head) {
                            print(df[1:dots$n, ], row.names = FALSE, digits = dots$digits)
                            cat("...\n")
                        }
                        else {
                            print(df, row.names = FALSE, digits = dots$digits)
                        }
                    }
                }
                else { ## latency or beran, local, with CI, only long by default 
                    for (i in 1:length(x$h)) {
                        cat("\nx0 =", format(x$x0[i], digits = dots$digits), "\n")
                        if (length(x$h) == 1)
                            df <- data.frame(x["testim"], x$S, x$conf)
                        else
                            df <- data.frame(x["testim"], x$S[i], x$conf[[i]])
                        dimnames(df)[[2]] <- c("time", x$type, sprintf("lower %s%% CI", formatC(x$conflevel*100)), sprintf("upper %s%% CI", formatC(x$conflevel*100)))
                        if (head) {
                            print(df[1:dots$n, ], row.names = FALSE, digits = dots$digits)
                            cat("...\n")
                        }
                        else {
                            print(df, row.names = FALSE, digits = dots$digits)
                        }
                    }
                }
            }
            else { ## latency or beran, global
                if (is.null(x$conf)) { ## without CI
                    if (!missing(how) && how == "long") { ## long
                        for (i in 1:length(x$x0)) {
                            for (j in 1:length(x$h)) {
                                cat("\nx0 =", format(x$x0[i], digits = dots$digits))
                                cat("\nh =", format(x$h[j], digits = dots$digits), "\n")
                                df <- data.frame(x["testim"], x$S[[j]][i])
                                dimnames(df)[[2]] <- c("time", x$type)
                                if (head) {
                                    print(df[1:dots$n, ], row.names = FALSE, digits = dots$digits)
                                    cat("...\n")
                                }
                                else {
                                    print(df, row.names = FALSE, digits = dots$digits)
                                }
                            }
                        }  
                    }
                    else { ## wide
                        for (i in 1:length(x$x0)) {
                            cat("\nx0 =", format(x$x0[i], digits = dots$digits), "\n")
                            df <- data.frame(x["testim"]) 
                            for (j in 1:length(x$h))
                                df <- cbind(df, x$S[[j]][i])
                            dimnames(df)[[2]] <- c("time", formatC(paste("h = ", x$h, sep = ""), dots$digits + 1))
                            if (head) {
                                print(df[1:dots$n, ], row.names = FALSE, digits = dots$digits)
                                cat("...\n")
                            }
                            else {
                                print(df, row.names = FALSE, digits = dots$digits)
                            }
                        }
                    }
                }
                else { ## latency or beran, global, with CI
                    for (i in 1:length(x$x0)) { ## only long (default)
                        for (j in 1:length(x$h)) {
                            cat("\nx0 =", format(x$x0[i], digits = dots$digits), "\n")
                            cat("h =", format(x$h[j], digits = dots$digits), "\n")
                            df <- data.frame(x["testim"])
                            df <- cbind(df, x$S[[j]][i], x$conf[[j]][i])
                            dimnames(df)[[2]] <- c("time", x$type, sprintf("lower %s%% CI", formatC(x$conflevel*100)), sprintf("upper %s%% CI", formatC(x$conflevel*100)))
                            if (head) {
                                print(df[1:dots$n, ], row.names = FALSE, digits = dots$digits)
                                cat("...\n")
                            }
                            else {
                                print(df, row.names = FALSE, digits = dots$digits)
                            }
                        }
                    }
                }
            }
        }
    }
    else { ## of if (length(x$type) == 1)
        if (!missing(how))
            warning("The option 'how' is ignored for 'npcure' objects with 'type' component of length 2")
        if (x$type[1] == "test") { ## test
            cat(x$type[2], x$type[1], "\n\n")
            if (x$type[2] == "Covariate") {
                cat("Covariate: ", x$x, "\n")
                df <- data.frame(test = c("Cramer-von Mises", "Kolmogorov-Smirnov"), statistic = c(x$CM$statistic, x$KS$statistic), p.value =  c(x$CM$pvalue, x$KS$pvalue))
            }
            else {
                df <- data.frame(statistic = x$aux$statistic, n = x$aux$n,  p.value = x$pvalue)
            }
            print(df, row.names = FALSE, digits = dots$digits)
        }
        else { ## bandwidth
            cat(x$type[1], "(h) for", if (x$type[2] == "survival") "Beran's", x$type[2], "estimator conditional on covariate x0:\n")
            if (head && length(x$x0) < dots$n)
                head <- FALSE
            if (is.null(x$hsmooth)) {
                df <- data.frame(x[c("x0", "h")])
                dimnames(df)[[2]] <- c("x0", "h")
            }
            else {
                df <- data.frame(x[c("x0", "h", "hsmooth")])
                dimnames(df)[[2]] <- c("x0", "h", "h.smooth")
            }
            cat("\n")
            if (head) {
                print(df[1:dots$n, ], row.names = FALSE, digits = dots$digits)
                cat("...\n")
            }
            else {
                print(df, row.names = FALSE, digits = dots$digits)
            }
        }
    }
}
