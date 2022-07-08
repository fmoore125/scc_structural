setwd("~/research/scciams/scc_structural")

source("src/analysis/find_distribution.R")
source("src/data_cleaining_scripts/cleaning_master.R")

all.qs <- c(0, 0.001, 0.01, .025, .05, .1, .17, .25, .5, .75, .83, .9, .95, .975, .99, 0.999, 1)
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
    if (is.null(values)) {
        results <- rbind(results, data.frame(ii, central=NA, score=NA, solution=last.solution))
        next
    }

    score <- score.dist.draws(mu, qs, as, values)
    results <- rbind(results, data.frame(ii, central=get.central(mu, qs, as), score, solution=last.solution))
}

library(ggplot2)

ggplot(results, aes(pmin(2.5, score / central), fill=solution)) +
    geom_histogram() + scale_y_continuous(expand=c(0, 0)) +
    theme_bw() + xlab("RMSE / Central Estimate")

ggplot(results, aes(ifelse(score / central < 1, score / central, NA), fill=solution)) +
    geom_histogram() + coord_cartesian(ylim=c(0, 100)) + scale_y_continuous(expand=c(0, 0)) +
    scale_x_continuous(expand=c(0, 0)) +
    theme_bw() + xlab("RMSE / Central Estimate")

ggplot(subset(results, !(solution %in% c('delta', 'delta ext-coded', 'discrete'))), aes(ifelse(score / central < 1, score / central, NA), fill=solution)) +
    geom_histogram() + scale_y_continuous(expand=c(0, 0)) +
    scale_x_continuous(expand=c(0, 0)) + scale_fill_discrete("Distribution:") +
    theme_bw() + xlab("RMSE / Central Estimate")
ggsave("figures/distsoln.pdf", width=6.5, height=3.5)

results$short <- sapply(results$solution, function(xx) strsplit(xx, ' ')[[1]][1])
ggplot(subset(results, !(short %in% c('delta'))), aes(ifelse(score / central > 0 & score / central < 1, score / central, NA), fill=short)) +
    geom_histogram() + scale_y_continuous(expand=c(0, 0)) +
    scale_x_continuous(expand=c(0, 0)) + scale_fill_discrete("Distribution:") +
    theme_bw() + xlab("RMSE / Central Estimate")

table(results$solution)

sum(results$score / results$central > 1, na.rm=T)
sum(results$score / results$central < .01, na.rm=T)
sum(rowSums(!is.na(dat[, which(names(dat) == 'Min'):which(names(dat) == 'Max')])) == 0 & !is.na(dat$`Central Value ($ per ton CO2)`))

results[results$score > abs(results$central),]
