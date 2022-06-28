setwd("~/research/scciams/scc_structural")

source("src/analysis/find_distribution.R")
source("src/data_cleaining_scripts/cleaning_master.R")

get.hill <- function(qs, as) {
    aboves <- (as > get.central(mu, qs, as) | qs > .5) & qs < 1

    if (sum(aboves) < 2)
        return(NULL)

    soln <- optimize(function(xi) {
        qjk <- (max(1 - qs[aboves]) / (1 - qs[aboves]))^(1 / xi) * min(as[aboves])
        max(abs(as[aboves] - qjk))
    }, c(0, 10))

    if (sum(aboves) == 2)
        return(data.frame(kk=sum(aboves), xi=soln$minimum, objective=soln$objective, se=NA))

    xis <- c()
    for (ii in which(aboves)) {
        aboves2 <- aboves
        aboves2[ii] <- F
        soln2 <- optimize(function(xi) {
            qjk <- (max(1 - qs[aboves2]) / (1 - qs[aboves2]))^(1 / xi) * min(as[aboves2])
            max(abs(as[aboves2] - qjk))
        }, c(0, 10))
        xis <- c(xis, soln2$minimum)
    }

    data.frame(kk=sum(aboves), xi=soln$minimum, objective=soln$objective, se=sd(xis))
}

qs <- c(.75, .9, .95)
get.hill(qs, qpareto(qs, 10, .1))
get.hill(qs, qpareto(qs, 10, .5))
get.hill(qs, qpareto(qs, 10, 1))
get.hill(qs, qpareto(qs, 10, 10))
get.hill(qs, qunif(qs))
get.hill(qs, qunif(qs, -1, 10))
get.hill(qs, qnorm(qs))
get.hill(qs, qnorm(qs, 10, 100))
get.hill(qs, qnorm(qs, 10, .1))
get.hill(qs, qexp(qs))
get.hill(qs, qexp(qs, .1))
get.hill(qs, qexp(qs, 10))
qs <- c(.75, .9, .95, .99)
get.hill(qs, qnorm(qs))
get.hill(qs, qexp(qs))

all.qs <- c(0, 0.001, 0.01, .025, .05, .1, .17, .25, .5, .75, .83, .9, .95, .975, .99, 0.999, 1)
all.as.cols <- which(names(dat) == 'Min'):which(names(dat) == 'Max')

results <- data.frame()
for (ii in 1:nrow(dat)) {
    print(ii)
    all.as <- t(dat[ii, all.as.cols])
    qs <- all.qs[!is.na(all.as)]
    as <- all.as[!is.na(all.as)]
    mu <- dat$`Central Value ($ per ton CO2)`[ii]

    soln <- get.hill(qs, as)
    if (is.null(soln)) {
        ## results <- rbind(results, data.frame(ii, hill=NA, note="Insufficient data"))
        next
    }

    results <- rbind(results, cbind(ii=ii, soln))
}

