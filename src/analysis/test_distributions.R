setwd("~/research/scciams/scc_structural")

source("src/analysis/find_distribution.R")
source("src/data_cleaining_scripts/cleaning_standardizedollaryears.R")

all.qs <- c(0, .025, .05, .1, .17, .25, .5, .75, .83, .9, .95, .975, .99, 1)
all.as.cols <- which(names(dat) == 'Min'):which(names(dat) == 'Max')

results <- data.frame()
for (ii in 1:nrow(dat)) {
    print(ii)
    all.as <- t(dat[ii, all.as.cols])
    qs <- all.qs[!is.na(all.as)]
    as <- all.as[!is.na(all.as)]
    mu <- dat$`Central Value ($ per ton CO2)`[ii]
    if (is.na(mu) && length(qs) == 0) {
        results <- rbind(results, data.frame(ii, central=NA, score=NA, solution="no data"))
        next
    }

    values <- generate.pdf(mu, qs, as, 1e6)
    score <- score.dist.draws(mu, qs, as, values)

    results <- rbind(results, data.frame(ii, central=get.central(mu, qs, as), score, solution=last.solution))
}

library(ggplot2)

ggplot(results, aes(pmin(2.5, sqrt(score) / central), fill=solution)) +
    geom_histogram() + scale_y_continuous(expand=c(0, 0)) +
    theme_bw() + xlab("Root-Sum-Squared Errors / Central Estimate")

ggplot(results, aes(ifelse(sqrt(score) / central < 1, sqrt(score) / central, NA), fill=solution)) +
    geom_histogram() + coord_cartesian(ylim=c(0, 100)) + scale_y_continuous(expand=c(0, 0)) +
    scale_x_continuous(expand=c(0, 0)) +
    theme_bw() + xlab("Root-Sum-Squared Errors / Central Estimate")
