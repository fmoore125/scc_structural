## Overview:
##
## The core function defined here is `generate.pdf`, which takes draws
## from a distribution estimated from its mean and quantile values.
##
## Throughout, the `mu` parameter is assumed to be the mean, and may
## be NA, `qs` is a vector of quantiles between 0 and 1, and `as` is
## the corresponding vector of locations of those quantiles.
##
## The distribution estimation process is as follows:
##
## 1. If only a central value is given, all draws match it.
##
## 2. If only min & max, or only min, max, and central are given, it
##    is assumed to be a triangular distribution with the given
##    extreme values.
##
## 3. If a min or max is given, the distribution is bottom- or
##    top-coded with that value.
##
## 4. If a .001 or .999 quantile is given, the distribution is
##    truncated to those values.
##
## 5. If a central value (mean or median) and one other
##    (non-truncating) quantile is given, the distribution is assumed
##    to be Gaussian.
##
## 6. If a central value (mean or median) and two other
##    (non-truncating) quantiles are given, the distribution is
##    assumed to be either a Skew normal or a exponentially modified
##    normal, whichever produces a better fit.
##
## 7. Otherwise, if `allow.mixed` is true, the distribution is assumed
##    to be a mixture of up to k - 2 Gaussians, where k is the number
##    of fitting values. If `allow.mixed` is false, a skew-normal or
##    exponentially modified normal fit is attempted.
##
## 8. If cases 3 - 7 are used, an alternnative model consisting of a
##    piecewise uniform distribution with weights from the spans
##    between quantiles, combined with a fitted tail distribution, is
##    tried as an alternative, and the best-scoring distribution is
##    returned.

library(sn)
library(EnvStats)

allow.mixed <- F

##### Checks for special cases

## Discrete distributions

is.one.value <- function(mu, qs, as)
(!is.na(mu) && length(qs) == 0) || (length(qs) == 1 && (is.na(mu) || mu == as)) ||
    ((is.na(mu) || mu == as[1]) && min(qs) == 0 && max(qs) == 1 && diff(as) == 0)



is.two.values <- function(mu, qs, as) {
    if (length(qs) == 2 && any(qs == 0) && any(qs == 1))
        return(is.na(mu) || mu == mean(as))
    return(F)
}

is.three.values <- function(mu, qs, as) {
    if (length(qs) == 2 && any(qs == 0) && any(qs == 1))
        return(!is.na(mu) && mu != mean(as))
    return(length(qs) == 3 && all(qs %in% c(0, .5, 1)) && is.na(mu))
}

## Characteristics of continuous distributions

is.bounded <- function(qs)
    any(qs %in% c(0, .001, .999, 1))

is.normal <- function(mu, qs, as) {
    return((!is.na(mu) && length(qs) == 1 && qs != .5 && !is.bounded(qs)) ||
           (length(qs) == 2 && any(qs == .5) && !is.bounded(qs)
               && (is.na(mu) || mu == as[qs == .5])) ||
           (is.na(mu) && length(qs) == 2 && !any(qs == .5) && !is.bounded(qs)))
}

is.skewnormal <- function(mu, qs, as) {
    return((!is.na(mu) && length(qs) == 2 && !is.bounded(qs)) ||
           (length(qs) == 3 && any(qs == .5) && !is.bounded(qs)
               && (is.na(mu) || mu == as[qs == .5])) ||
           (!is.na(mu) && length(qs) == 1 && qs == .5 && as != mu))
}

## Other chekcs

check.valid <- function(mu, qs, as) {
    if (length(qs) > 0) {
        if (length(qs) != length(as))
            return(F)
        if (!all(order(qs) == 1:length(qs)))
            return(F)
        if (!all(order(as) == 1:length(as)))
            return(F)
        if (!is.na(mu)) {
            if (qs[1] == 0 && mu < as[1])
                return(F)
            if (qs[length(qs)] == 1 && mu > as[length(as)])
                return(F)
        }
    }
    return(T)
}

##### Fitting functions

## Generate a PDF with a mean mu, and the value at quantiles qs equal
## to as, Then take N draws from that distribution
last.solution <- NA
generate.pdf <- function(mu, qs, as, N) {
    last.solution <<- NA

    if (!check.valid(mu, qs, as)) {
        last.solution <<- "invalid"
        return(NULL)
    }

    ## Simple distributions

    if (is.one.value(mu, qs, as)) {
        last.solution <<- "delta"
        return(rep(get.central(mu, qs, as), N))
    }

    if (is.two.values(mu, qs, as) || is.three.values(mu, qs, as)) {
        last.solution <<- "triangle"
        central <- get.central(mu, qs, as)
        if (central == as[qs == 0] || central == as[qs == 1]) {
            if (central == as[qs == 0])
                return(rtri(N, as[qs == 0], as[qs == 1], central + 1e-6))
            else if (central == as[qs == 1])
                return(rtri(N, as[qs == 0], as[qs == 1], central - 1e-6))
        } else
            return(rtri(N, as[qs == 0], as[qs == 1], central))
    }

    values <- generate.general.pdf(mu, qs, as, N)

    ## Try piecewise uniform as a backup
    values.alt <- generate.piecewise.pdf(mu, qs, as, N)
    if (!is.null(values.alt)) {
        score <- score.dist.draws(mu, qs, as, values)
        score.alt <- score.dist.draws(mu, qs, as, values.alt)
        if (score.alt < score) {
            if (is.bounded(qs))
                last.solution <<- "piecewise bounded"
            else {
                ## Estimate best tail
                if (min(qs) > 0) {
                    left.fit <- fit.left.tail(mu, qs, as)
                    values.alt[values.alt < min(as)] <- generate.left.tail.pdf(left.fit, sum(values.alt < min(as)), min(as))
                } else {
                    values.alt[values.alt < min(as)] <- min(as)
                    left.fit <- list(type='bound')
                }

                if (min(qs) < 1) {
                    right.fit <- fit.right.tail(mu, qs, as)
                    values.alt[values.alt > max(as)] <- generate.right.tail.pdf(right.fit, sum(values.alt > max(as)), max(as))
                } else {
                    values.alt[values.alt > max(as)] <- max(as)
                    right.fit <- list(type='bound')
                }

                last.solution <<- paste0("piecewise ", left.fit$type, "-", right.fit$type)
            }
            values <- values.alt
        }
    }

    values
}

## Always produces bounded tails
generate.piecewise.pdf <- function(mu, qs, as, N) {
    qs2 <- qs
    as2 <- as
    if (!is.na(mu) && !any(qs == .5)) {
        qs2 <- c(qs2, .5)
        as2 <- c(as2, mu)
    }

    if (length(qs2) <= 1)
        return(NULL)

    as2 <- as2[order(qs2)]
    qs2 <- sort(qs2)

    if (any(as2 != sort(as2)))
        return(NULL)

    if (is.bounded(qs2)) {
        if (qs2[1] != 0)
            qs2[1] <- 0
        if (qs2[length(qs2)] != 1)
            qs2[length(qs2)] <- 1
    }

    quants <- runif(N)
    values <- rep(NA, N)
    below <- quants < qs2[1]
    values[below] <- as2[1]
    for (ii in 2:length(qs2)) {
        within <- quants >= qs2[ii-1] & quants < qs2[ii]
        values[within] <- runif(sum(within), as2[ii-1], as2[ii])
    }
    above <- quants >= qs2[length(qs2)]
    values[above] <- as2[length(as2)]

    values
}

## Handle skew-normal and exp-modified normal cases
generate.skewnormal.pdf <- function(mu, qs, as, N) {
    ## xi = 0, omega = 1, alpha = 1
    ## mu = 0.5627814
    ## qs = c(.1, .75)
    ## as = c(-0.4772229, 1.1075878)
    if (!is.na(mu) && any(qs == .5) && as[qs == .5] != mu) {
        xi0 <- as[qs == .5]
        omega0 <- abs(mu - xi0)
        alpha0 <- sign(mu - xi0)
    } else {
        xi0 <- ifelse(!is.na(mu), mu, as[qs == .5])
        omega0 <- max(abs(as[qs != .5] - mu))
        alpha0 <- 0
    }
    ## Try skew normal
    result1 <- optim(c(xi0, omega0, alpha0), function(params) {
        delta <- params[3] / sqrt(1 + params[3]^2)
        mu.pred <- params[1] + abs(params[2]) * delta * sqrt(2 / pi)
        as.pred <- tryCatch({
            qsn(qs, params[1], abs(params[2]), params[3])
        }, error=function(e) {
            NA * qs
        })
        score.dist(mu, mu.pred, as, as.pred)
    })
    ## Try exponentially modified Gaussian distribution
    result2 <- optim(c(xi0, omega0, 1), function(params) {
        mu.pred <- params[1] + 1 / abs(params[3])
        as1 <- qnorm(qs, params[1], abs(params[2])) + qexp(.5, abs(params[3]))
        as2 <- qnorm(.5, params[1], abs(params[2])) + qexp(qs, abs(params[3]))
        as.pred <- (as1 + as2) / 2
        score.dist(mu, mu.pred, as, as.pred)
    })
    ## Re-evaluate result2's score
    result2.draws <- rnorm(1e6, result2$par[1], abs(result2$par[2])) + rexp(1e6, abs(result2$par[3]))
    result2$value <- score.dist.draws(mu, qs, as, result2.draws)

    if (result1$value < result2$value) {
        last.solution <<- "skew-normal"
        return(rsn(N, result1$par[1], result1$par[2], result1$par[3]))
    } else {
        last.solution <<- "exp-normal"
        return(rnorm(N, result2$par[1], abs(result2$par[2])) + rexp(N, abs(result2$par[3])))
    }
}

## Handle general cases
generate.general.pdf <- function(mu, qs, as, N) {
    ## Truncated distributions

    if (is.bounded(qs)) {
        ## Top/bottom-coding
        if (any(qs == 0 | qs == 1)) {
            vals <- generate.pdf(mu, qs[!(qs %in% c(0, 1))],
                                 as[!(qs %in% c(0, 1))], N)
            ## Apply shift to maintain results
            result <- optimize(function(shift) {
                newvals <- vals + shift
                if (0 %in% qs)
                    newvals[newvals < as[qs == 0]] <- as[qs == 0]
                if (1 %in% qs)
                    newvals[newvals > as[qs == 1]] <- as[qs == 1]
                score.dist.draws(mu, qs, as, newvals)
            }, c(-1, 1) * max(abs(vals)))

            newvals <- vals + result$minimum
            if (0 %in% qs)
                newvals[newvals < as[qs == 0]] <- as[qs == 0]
            if (1 %in% qs)
                newvals[newvals > as[qs == 1]] <- as[qs == 1]

            last.solution <<- paste(last.solution, "ext-coded")
            return(newvals)
        }

        ## Truncate distribution
        if (any(qs == .001 | qs == .999)) {
            vals <- generate.pdf(mu, qs[!(qs %in% c(.001, .999))],
                                 as[!(qs %in% c(.001, .999))], N)
            ## Apply shift to maintain results
            result <- optimize(function(shift) {
                newvals <- vals + shift
                valid <- rep(T, N)
                if (.001 %in% qs)
                    valid <- newvals >= as[qs == .001]
                if (.999 %in% qs)
                    valid <- valid & (newvals <= as[qs == .999])
                score.dist.draws(mu, qs, as, newvals[valid])
            }, c(-1, 1) * max(abs(vals)))

            newvals <- vals + result$minimum
            valid <- rep(T, N)
            if (.001 %in% qs)
                valid <- newvals >= as[qs == .001]
            if (.999 %in% qs)
                valid <- valid & (newvals <= as[qs == .999])

            last.solution <<- paste(last.solution, "trucated")
            return(sample(vals[valid], N, replace=T))
        }
    }

    ## None of the following have trunctions

    if (is.normal(mu, qs, as)) {
        if (!is.na(mu) || any(qs == .5)) {
            ## We have a central value: easy case
            mu <- ifelse(!is.na(mu), mu, as[qs == .5])
            sigma <- (as[qs != .5] - mu) / qnorm(qs[qs != .5])
        } else {
            mu0 <- get.central(mu, qs, as)
            sigma0 <- mean((as - mu0) / qnorm(qs))
            result <- optim(c(mu0, sigma0), function(params) {
                mu.pred <- NA
                as.pred <- qnorm(qs, params[1], params[2])
                score.dist(mu, mu.pred, as, as.pred)
            })
            mu <- result$par[1]
            sigma <- result$par[2]
        }

        last.solution <<- "normal"
        return(rnorm(N, mu, sigma))
    }

    if (is.skewnormal(mu, qs, as) || !allow.mixed) {
        return(generate.skewnormal.pdf(mu, qs, as, N))
    }

    ## Estimate general mixture of Gaussians

    ## Take an average sd
    mus <- get.central(mu, qs, as)
    if (length(qs) == 1) {
        sigmas <- abs(mus[1] - as)
    } else {
        sigmas <- mean(abs((as[qs != .5] - mus[1]) / qnorm(qs[qs != .5])))
    }
    scales <- 1

    Nobs <- length(qs) - is.na(mu)

    ## Add on up to Nobs - 2 additional entries
    for (ii in 1:min(1, (Nobs - 2))) {
        ## Check which quantile fails most
        pred.qs <- get.pred.qs(scales, mus, sigmas, as)
        errors <- pred.qs - qs
        sqrerrs <- errors^2
        if (sum(sqrerrs) < 1e-6)
            break

        if (errors[which.max(sqrerrs)] < 0) {
            ## Want more mass below this point
            if (qs[which.max(sqrerrs)] == min(qs)) {
                new.mu <- as[which.max(sqrerrs)] / 2
                mu.q <- qs[which.max(sqrerrs)] / 2
            } else {
                new.mu <- mean(as[qs < qs[which.max(sqrerrs)]])
                mu.q <- mean(qs[qs < qs[which.max(sqrerrs)]])
            }
        } else {
            ## Want more mass above this point
            if (qs[which.max(sqrerrs)] == max(qs)) {
                new.mu <- 2 * as[which.max(sqrerrs)]
                mu.q <- 1 - (1 - qs[which.max(sqrerrs)]) / 2
            } else {
                new.mu <- mean(as[qs > qs[which.max(sqrerrs)]])
                mu.q <- mean(qs[qs > qs[which.max(sqrerrs)]])
            }
        }
        new.sigma <- abs(as[which.max(sqrerrs)] - new.mu) / (abs(qnorm(abs(qs[which.max(sqrerrs)] - mu.q))) + .5)

        ## Fit the distribution
        result <- optim(.5, function(rescale) {
            if (!is.na(mu)) {
                ## Require that rescale shift mu + (1 - rescale) new.mu = mu
                shift <- (mu - (1 - rescale) * new.mu) / (rescale * mu)
                new.pred.qs <- get.pred.qs(c(rescale * scales, 1 - rescale), c(shift * mus, new.mu), c(sigmas, new.sigma), as)
            } else
                new.pred.qs <- rescale * pred.qs + (1 - rescale) * pnorm(as, new.mu, new.sigma)
            new.sqrerrs <- (new.pred.qs - qs)^2
            sum(new.sqrerrs)
        }, method='L-BFGS-B', lower=0, upper=1)

        if (!is.na(mu)) {
            shift <- (mu - (1 - result$par) * new.mu) / (result$par * mu)
            mus <- c(shift * mus, new.mu)
        } else {
            mus <- c(mus, new.mu)
        }
        sigmas <- c(sigmas, new.sigma)
        scales <- c(result$par * scales, (1 - result$par))
    }

    ## Take draws from distributions
    values <- c()
    for (kk in 2:length(mus))
        values <- c(values, rnorm(floor(N * scales[kk]), mus[kk], sigmas[kk]))

    last.solution <<- "mixed gaussian"
    c(values, rnorm(N - length(values), mus[1], sigmas[1]))
}

fit.distribution <- function(init, qdist, mudist, mu, qs, as) {
    if (!is.null(mudist)) {
        optim(init, function(params) {
            as.fit <- qdist(qs, params)
            if (is.null(as.fit))
                return(Inf)
            score.dist(mu, mudist(params), as, as.fit)
        })
    } else {
        optim(init, function(params) {
            as.fit <- qdist(qs, params)
            if (is.null(as.fit))
                return(Inf)
            score.dist(mu, NA, as, as.fit)
        })
    }
}

fit.gaussian <- function(mu, qs, as) {
    if (!is.na(mu)) {
        result <- fit.distribution(c(mu, abs(mu)), function(qs, params) {
            qnorm(qs, params[1], abs(params[2]))
        }, function(params) {
            params[1]
        }, mu, qs, as)
    } else {
        result <- fit.distribution(c(mean(as), abs(mean(as))), function(qs, params) {
            qnorm(qs, params[1], abs(params[2]))
        }, NULL, mu, qs, as)
    }

    if (result$par[2] < 0)
        result$par[2] <- 0

    result
}

## Fit tail distribution(s)
fit.right.tail <- function(mu, qs, as) {
    aboves <- as > get.central(mu, qs, as) | qs > .5
    if (sum(aboves) == 0) {
        ## No tail information: Fit Gaussian to whole
        result <- fit.gaussian(mu, qs, as)
        return(list('type'='gaussian', 'params'=result$par, 'score'=result$val))
    }

    if (sum(aboves) > 2)
        mu.include <- NA
    else
        mu.include <- get.central(mu, qs, as)

    ## Gaussian
    result <- fit.gaussian(mu.include, qs[aboves], as[aboves])
    best.fit <- list('type'='gaussian', 'params'=result$par, 'score'=result$val)

    if (sum(aboves) >= 1) {
        init <- c(ifelse(is.na(mu.include), mean(as), mu.include), max(as))
        if (diff(init) > 0) {
            result <- fit.distribution(init, function(qs, params) {
                if (params[1] > params[2])
                    return(NULL)
                qtri(qs, params[1] - diff(params), params[2], params[1])
            }, function(params) {
                params[1]
            }, mu.include, qs[aboves], as[aboves])

            if (result$val < best.fit$score && result$par[2] > max(as)) # Ensure that can draw above max
                best.fit <- list('type'='triangle', 'params'=result$par, 'score'=result$val)
        }

        ## Limit to above median (drop above means)
        as2 <- as[qs > .5]
        qs2 <- qs[qs > .5]
        if ((length(as2) + !is.na(mu.include)) > 2) {
            log.flip.cdf.y <- log(1 - qs2)
            mod <- lm(log.flip.cdf.y ~ as2)
            lambda <- abs(mod$coeff[2])
            if (is.na(lambda))
                lambda <- 1 # Why not?

            result <- fit.distribution(c(ifelse(is.na(mu.include), mean(as2), mu.include), lambda), function(qs2, params) {
                params[1] + qexp(1 - 2*(1 - qs2), abs(params[2]))
            }, function(params) {
                params[1]
            }, mu.include, qs2, as2)

            if (result$val < best.fit$score)
                best.fit <- list('type'='exponential', 'params'=result$par, 'score'=result$val)
        }
    }

    best.fit
}

fit.left.tail <- function(mu, qs, as) {
    best.fit <- fit.right.tail(-mu, 1 - qs, -as)
    if (best.fit$type == 'gaussian')
        best.fit$params[1] <- -best.fit$params[1]
    else if (best.fit$type == 'triangle')
        best.fit$params <- -best.fit$params
    else if (best.fit$type == 'exponential')
        best.fit$params[1] <- -best.fit$params[1]

    best.fit
}

## Generate values from tail distributions
generate.right.tail.pdf <- function(right.fit, N, maxas) {
    if (is.null(right.fit))
        return(rep(maxas, N))
    if (right.fit$type == 'gaussian')
        values <- rnorm(1e4, right.fit$params[1], right.fit$params[2])
    else if (right.fit$type == 'triangle')
        values <- rtri(1e4, right.fit$params[1] - diff(right.fit$params), right.fit$params[2], right.fit$params[1])
    else if (right.fit$type == 'exponential')
        values <- right.fit$params[1] + rexp(1e4, abs(right.fit$params[2]))

    values <- values[values > maxas]
    if (length(values) < N)
        values <- c(values, generate.right.tail.pdf(right.fit, N - length(values), maxas))

    values[1:N]
}

generate.left.tail.pdf <- function(left.fit, N, minas) {
    if (is.null(left.fit))
        return(rep(minas, N))
    if (left.fit$type == 'gaussian')
        values <- rnorm(1e4, left.fit$params[1], left.fit$params[2])
    else if (left.fit$type == 'triangle')
        values <- rtri(1e4, left.fit$params[2], left.fit$params[1] - diff(left.fit$params), left.fit$params[1])
    else if (left.fit$type == 'exponential')
        values <- left.fit$params[1] - rexp(1e4, abs(left.fit$params[2]))

    values <- values[values < minas]
    if (length(values) < N)
        values <- c(values, generate.left.tail.pdf(left.fit, N - length(values), minas))

    values[1:N]
}

## Score a proposed distribution against the data
score.dist <- function(mu.true, mu.pred, as.true, as.pred) {
    if (is.na(mu.true))
        return(sqrt(mean((as.true - as.pred)^2)))
    else
        return(sqrt(((mu.true - mu.pred)^2 + sum((as.true - as.pred)^2)) / (1 + length(as.true))))
}

## Score a set of draws from a distribution
score.dist.draws <- function(mu.true, qs, as.true, draws) {
    mu.pred <- mean(draws)
    as.pred <- quantile(draws, qs)
    score.dist(mu.true, mu.pred, as.true, as.pred)
}

## Get the predicted quantiles given a mixed Gaussian distribution
get.pred.qs <- function(scales, mus, sigmas, as) {
    ## scales are the weights for k Gaussians
    ## mus and sigmas describe those k Gaussians
    ## as are values at which the quantiles should be estimated
    stopifnot(sum(scales) == 1)

    pred.qs <- rep(0, length(as))
    for (kk in 1:length(mus))
        pred.qs <- pred.qs + scales[kk] * pnorm(as, mus[kk], sigmas[kk])

    pred.qs
}

## Get the best central value we can
get.central <- function(mu, qs, as) {
    ifelse(!is.na(mu), mu,
           ifelse(any(qs == .5), as[qs == .5], mean(as)))
}

## Return a data frame with validation checks
validate.pdf <- function(mu, qs, as, generate.pdf.func=generate.pdf) {
    draws <- generate.pdf(mu, qs, as, 1e6)
    data.frame(metric=c('mean', qs), desired=c(mu, as), observed=c(mean(draws), quantile(draws, qs)))
}

if (F) {
    ## Tests
    qs <- .5
    as <- 50
    mu <- 100

    validate.pdf(mu, qs, as)

    qs <- c(0, .5, 1)
    as <- c(10, 50, 200)
    mu <- 100

    validate.pdf(mu, qs, as)

    qs <- c(.25, .75)
    as <- c(0, 100)
    mu <- 75

    validate.pdf(mu, qs, as)

    qs <- c(.25, .5, .75)
    as <- c(10, 50, 100)
    mu <- 100

    validate.pdf(mu, qs, as)

    qs <- c(.1, .25, .5, .75)
    as <- c(0, 10, 50, 100)
    mu <- 100

    validate.pdf(mu, qs, as)

    ## Test tail fits

    ## Fit Gaussian tails
    fit.right.tail(3, c(.1, .25, .75, .9), qnorm(c(.1, .25, .75, .9), 3, 2))
    fit.left.tail(3, c(.1, .25, .75, .9), qnorm(c(.1, .25, .75, .9), 3, 2))

    ## Fit exponential tails
    qq <- seq(0, 1, length.out=100)
    aa <- c(3 - qexp(1 - 2*qq[1:50], .5), 3 + qexp(1 - 2*(1 - qq[51:100]), .5))
    plot(aa, qq)

    fit.right.tail(3, c(.1, .25, .75, .9), c(3 - qexp(c(.8, .5), .5), 3 + qexp(c(.5, .8), .5)))
    fit.left.tail(3, c(.1, .25, .75, .9), c(3 - qexp(c(.8, .5), .5), 3 + qexp(c(.5, .8), .5)))

    ## Fit triangular tails
    fit.right.tail(3, c(.1, .25, .75, .9), qtri(c(.1, .25, .75, .9), 0, 6, 3))
    fit.left.tail(3, c(.1, .25, .75, .9), qtri(c(.1, .25, .75, .9), 0, 6, 3))

    ## Test skewnormal
    validate.pdf(10, c(.025, .975), c(0, 100), generate.skewnormal.pdf)
    validate.pdf(50, c(.025, .975), c(0, 100), generate.skewnormal.pdf)
    validate.pdf(90, c(.025, .975), c(0, 100), generate.skewnormal.pdf)
}
