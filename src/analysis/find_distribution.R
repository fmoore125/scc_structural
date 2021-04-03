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
##    is assumed to be a uniform distribution over those 2 or 3
##    values.
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
## 7. Otherwise, the distribution is assumed to be a mixture of up to
##    k - 2 Gaussians, where k is the number of fitting values.
##
## 8. If cases 3 - 7 are used, an alternnative model consisting of a
##    piecewise uniform distribution with weights from the spans
##    between quantiles is tried as an alternative, and the
##    best-scoring distribution is returned.

library(sn)

##### Checks for special cases

## Discrete distributions

is.one.value <- function(mu, qs, as)
    (!is.na(mu) && length(qs) == 0) || (length(qs) == 1 && qs == .5 && (is.na(mu) || mu == as))

is.two.values <- function(mu, qs, as) {
    if (length(qs) == 2 && any(qs == 0) && any(qs == 1))
        return(is.na(mu) || mu == mean(as))
    return(F)
}

is.three.values <- function(mu, qs, as) {
    if (length(qs) == 2 && any(qs == 0) && any(qs == 1))
        return(!is.na(mu) && mu != mean(as))
    return(F)
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

##### Fitting functions

## Generate a PDF with a mean mu, and the value at quantiles qs equal
## to as, Then take N draws from that distribution
last.solution <- NA
generate.pdf <- function(mu, qs, as, N) {
    last.solution <<- NA

    ## Discrete distributions

    if (is.one.value(mu, qs, as)) {
        last.solution <<- "delta"
        return(rep(get.central(mu, qs, as), N))
    }

    if (is.two.values(mu, qs, as)) {
        last.solution <<- "discrete"
        return(c(rep(as[1], N/2), rep(as[2], N/2)))
    }

    if (is.three.values(mu, qs, as)) {
        last.solution <<- "discrete"
        return(c(rep(as[1], floor(N/3)), rep(as[2], floor(N/3)), rep(mu, N - 2*floor(N/3))))
    }

    values <- generate.general.pdf(mu, qs, as, N)

    if (!is.bounded(qs)) {
        ## Try piecewise uniform as a backup
        qs2 <- qs
        as2 <- as
        if (!is.na(mu) && !any(qs == .5)) {
            qs2 <- c(qs2, .5)
            as2 <- c(as2, mu)
        }

        values.alt <- generate.piecewise.pdf(qs2, as2, N)
        if (!is.null(values.alt)) {
            score <- score.dist.draws(mu, qs, as, values)
            score.alt <- score.dist.draws(mu, qs, as, values.alt)
            if (score.alt < score) {
                last.solution <<- "piecewise"
                values <- values.alt
            }
        }
    }

    values
}

generate.piecewise.pdf <- function(qs, as, N) {
    as <- as[order(qs)]
    qs <- sort(qs)

    if (any(as != sort(as)))
        return(NULL)

    if (qs[1] != 0)
        qs[1] <- 0
    if (qs[length(qs)] != 1)
        qs[length(qs)] <- 1

    quants <- runif(N)
    values <- rep(NA, N)
    for (ii in 2:length(qs)) {
        within <- quants >= qs[ii-1] & quants < qs[ii]
        values[within] <- runif(sum(within), as[ii-1], as[ii])
    }
    values[is.na(values)] <- as[length(as)]

    values
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
                    valid <- valid & (newvals <= as[qs == .001])
                score.dist.draws(mu, qs, as, newvals[valid])
            }, c(-1, 1) * max(abs(vals)))

            newvals <- vals + result$minimum
            valid <- rep(T, N)
            if (.001 %in% qs)
                valid <- newvals >= as[qs == .001]
            if (.999 %in% qs)
                valid <- valid & (newvals <= as[qs == .001])

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

    if (is.skewnormal(mu, qs, as)) {
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
        # Try exponentially modified Gaussian distribution
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
        new.sigma <- abs(as[which.max(sqrerrs)] - new.mu) / (qnorm(abs(qs[which.max(sqrerrs)] - mu.q)) + .5)

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

## Score a proposed distribution against the data
score.dist <- function(mu.true, mu.pred, as.true, as.pred) {
    if (is.na(mu.true))
        return(sum((as.true - as.pred)^2))
    else
        return((mu.true - mu.pred)^2 + sum((as.true - as.pred)^2))
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
validate.pdf <- function(mu, qs, as) {
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
}
