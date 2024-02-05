## setwd("~/research/scciams/scc_structural")

source("src/analysis/randomforest_dists_load.R")
library(ggplot2)

if (F) {
    ad.test(draw ~ branch, data=data.frame(draw=c(rnorm(1e3, 0, 1), rnorm(1e3, 0, 1)), branch=rep(c('A', 'B'), each=1e3)), method='asymptotic')
    ad.test(draw ~ branch, data=data.frame(draw=c(rnorm(1e3, 0, 1), rep(0, 1e3)), branch=rep(c('A', 'B'), each=1e3)), method='asymptotic')
    ad.test(draw ~ branch, data=data.frame(draw=c(rnorm(1e3, 1, 1), rnorm(1e3, 0, 1)), branch=rep(c('A', 'B'), each=1e3)), method='asymptotic')
    ad.test(draw ~ branch, data=data.frame(draw=c(rnorm(1e3, 0, 1), rnorm(1e3, 0, 2)), branch=rep(c('A', 'B'), each=1e3)), method='asymptotic')
}

if (.Platform$OS.type != "unix") {
    pdffunc <- grDevices::cairo_pdf
} else
    pdffunc <- pdf

## Make tree with all values
tree <- make.tree(dat, incrows, dist, cols, mtry=length(cols), minleaf=5)
save(tree, file="outputs/rfdistsmodel.RData")

sum(map.tree(tree, function(tree) 1))

pdf.test <- data.frame()
for (ii in which(incrows)) {
    print(ii)
    draws <- predict.tree(tree, dat[ii,], dist)
    pdf.test <- rbind(pdf.test, data.frame(row=ii, quant=all.qs, pred=quantile(draws, all.qs),
                                 true=quantile(dist$draw[dist$row == ii], all.qs)))
}

ggplot(pdf.test, aes(true, pred, group=row)) +
    geom_line(linewidth=.1) + geom_abline(slope=1, col='red') +
    scale_x_log10("True SCC quantile") + scale_y_log10("Predicted SCC quantile") +
    theme_bw()
ggsave("figures/rfdists-test.pdf", width=6.5, height=5)

calc.mean <- function(tree) {
    mean(dist$draw[dist$row %in% tree$rows])
}

calc.quantile <- function(tree, qq) {
    quantile(dist$draw[dist$row %in% tree$rows], qq)
}

pdf.tree <- apply.tree(tree, function(init, branch, branches, tree) {
    if (!is.na(tree$split) && tree$split == 'terminal')
        splittext <- NA
    else if (is.na(tree$split))
        splittext <- paste(tree$col, "is NA")
    else if (is.numeric(tree$split))
        splittext <- paste(tree$col, "<", round(tree$split, 2))
    else
        splittext <- tree$col
    if (length(branches) == 1) { # Initial split
        data.frame(branch, column=tree$col, split=tree$split, splittext, allowwidth=init$allowwidth, count=length(tree$rows), mu=calc.mean(tree), min=calc.quantile(tree, 0), max=calc.quantile(tree, 1), parentx=init$x, parenty=init$y, x=init$x, y=init$y)
    } else if (!is.na(tree$split) && tree$split == 'terminal') {
        data.frame(branch, column=NA, split=tree$split, splittext, allowwidth=0, count=length(tree$rows), mu=calc.mean(tree), min=calc.quantile(tree, 0), max=calc.quantile(tree, 1), parentx=init$x, parenty=init$y, x=init$x - init$allowwidth * ((2*which(branches == branch)-1) / (2*length(branches)) - .5), y=init$y - 1)
    } else {
        data.frame(branch, column=tree$col, split=tree$split, splittext, allowwidth=init$allowwidth / length(branches), count=length(tree$rows), mu=calc.mean(tree), min=calc.quantile(tree, 0), max=calc.quantile(tree, 1), parentx=init$x, parenty=init$y, x=init$x - init$allowwidth * ((2*which(branches == branch)-1) / (2*length(branches)) - .5), y=init$y - 1)
    }
}, data.frame(x=0, y=0, allowwidth=2), "", "")
pdf.tree$slanty <- ifelse(pdf.tree$y >= -3, pdf.tree$y, ifelse(pdf.tree$y == -4, pdf.tree$y - pdf.tree$x - 1, -(pdf.tree$y + 4)**2 + pdf.tree$y - pdf.tree$x - 1))
pdf.tree$slantparenty <- ifelse(pdf.tree$parenty >= -3, pdf.tree$parenty, ifelse(pdf.tree$parenty == -4, pdf.tree$parenty - pdf.tree$parentx - 1, -(pdf.tree$parenty + 4)**2 + pdf.tree$parenty - pdf.tree$parentx - 1))
pdf.tree$branchtext <- ifelse(pdf.tree$branch %in% c('TRUE', 'FALSE', 'NA'), pdf.tree$branch, ifelse(pdf.tree$branch %in% c('1', '1.0'), "Yes", "No"))

ggplot(pdf.tree, aes(x, slanty)) +
    geom_segment(aes(xend=x, yend=slantparenty)) +
    geom_segment(aes(y=slantparenty, xend=parentx, yend=slantparenty)) +
    geom_label(aes(label=paste0("$", round(mu), " (", count, ")"))) +
    geom_label(aes(label=branchtext, y=slanty+.25)) +
    geom_text(data=subset(pdf.tree, slanty > -6), aes(label=splittext, y=slanty-.25)) +
    ylim(-8, 0) + theme_bw() +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
          axis.text.x=element_blank(), axis.text.y=element_blank(),
          axis.ticks.x=element_blank(), axis.ticks.y=element_blank())
ggsave("figures/rfdists-tree.pdf", width=20, height=10)

plot.split(tree, dist, incrows) + xlim(-1000, 2000)

forest <- make.forest(dat, incrows, dist, cols, mtry=10)
save(tree, forest, file="outputs/rfdistsmodel.RData")

## Prune trees with fewer than 7 nodes and large unclassified dists
tbl <- data.frame()
for (ii in 1:length(forest)) {
    print(ii)
    mins <- map.tree(forest[[ii]], function(tree) if (!is.na(tree$split) && tree$split == 'terminal') min(dist$draw[dist$row %in% tree$rows]) else NA)
    maxs <- map.tree(forest[[ii]], function(tree) if (!is.na(tree$split) && tree$split == 'terminal') max(dist$draw[dist$row %in% tree$rows]) else NA)
    cnts <- map.tree(forest[[ii]], function(tree) if (!is.na(tree$split) && tree$split == 'terminal') length(tree$rows) else NA)
    nodes <- sum(map.tree(forest[[ii]], function(tree) 1))
    tbl <- rbind(tbl, data.frame(ii, nodes=nodes, drop=nodes < 7 || any(maxs > 10000 & mins < 0 & cnts > 50))) # NA means keep
}
sum(is.na(tbl$drop))

for (ii in 1:length(forest)) {
    if (!is.na(tbl$drop[ii]) && tbl$drop[ii] == T)
        forest[[ii]] <- NA
}

save(forest, file="outputs/rfdistsmodel-final.RData")

pdf <- data.frame()
for (ii in sample(which(incrows))) {
    if (ii %in% pdf$row)
        next
    print(ii)
    quants.pred <- predict.forest(forest, dat[ii,], dist)
    quants.true <- quantile(dist$draw[dist$row == ii], all.qs)
    pdf <- rbind(pdf, data.frame(row=ii, quant=all.qs, pred=quants.pred, true=quants.true))
}

ggplot(pdf, aes(true, pred, group=row)) +
    geom_line(linewidth=.1, alpha=.5) + geom_abline(slope=1, col='red') +
    scale_x_log10("True SCC quantile") + scale_y_log10("Predicted SCC quantile") +
    theme_bw()
ggsave("figures/rfdists-test-forest.pdf", width=6.5, height=5)

## Variable importance
## load("rfdistsmodel.RData")

calc.varimport.tree(tree) %>% group_by(column) %>% summarize(import=max(import)) %>% arrange(desc(import))

varimport <- calc.varimport.forest(forest)
label.chg <- c('Other Market Failure?'='Feature: Other Market Failure',
               'Income Elasticity'='Uncertainty: Income Elasticity',
               'Equilibrium Climate Sensitivity'='Uncertainty: Equilibrium Climate Sensitivity',
               'Population Growth'='Uncertainty: Population Growth',
               'Tipping Point Magnitude'='Uncertainty: Tipping Point Magnitude',
               'Inequality Aversion'='Structural: Inequality Aversion',
               'Earth system'='Structural: Earth System',
               'Tipping Points'='Structural: Climate Tipping Points',
               'TFP Growth'='Uncertainty: TFP Growth',
               'Declining Discounting?'='Feature: Declining Discounting',
               'Ambiguity/Model Uncertainty'='Structural: Ambiguity/Model Uncertainty',
               'Carbon Cycle2'='Uncertainty: Carbon Cycle',
               'Adaptation Rates'='Uncertainty: Adaptation Rates',
               'Backstop Price?'='Feature: Backstop Price',
               'Tipping Points2'='Structural: Damage Tipping Points',
               'EMUC2'='Uncertainty: EMUC',
               'PRTP2'='Uncertainty: PRTP',
               'Limitedly-Substitutable Goods'='Structural: Limitedly-Substitutable Goods',
               'Damage Function'='Uncertainty: Damage Function',
               'log.scc.synth'='Quantity: Damage Function-based SCC',
               'Year'='Quantity: Publication Year',
               'discountrate'='Quantity: Standardized Discount Rate',
               'Learning'='Structural: Learning',
               'Persistent / Growth Damages'='Structural: Persistent/Growth Damages',
               'sccyearformerge'='Quantity: SCC Year',
               'Market Only Damages'='Feature: Market Only Damages',
               'Constant Discount Rate'='Uncertainty: Constant Discount Rate',
               'Emissions Growth'='Uncertainty: Emissions Growth',
               'Epstein-Zin'='Structural: Epstein-Zin Preferences',
               'Risk Aversion (EZ Utility)'='Uncertainty: EZ Risk Aversion',
               'Transient Climate Response'='Uncertainty: Transient Climate Response')
varimport$label <- label.chg[varimport$column]
varimport$label <- factor(varimport$label, levels=varimport$label[order(varimport$import)])

ggplot(varimport, aes(label, import)) +
    coord_flip() +
    geom_col() + theme_bw() +
    scale_y_continuous("Variable Importance", labels=scales::percent) + xlab(NULL)
ggsave("figures/rfdists-varimport.pdf", width=6.5, height=5)
write.csv(varimport, "outputs/rfdists-varimport.csv", row.names=F)

## Pre-compute quantiles

forest <- pre.compute.forest(forest, dist)
save(forest, file="outputs/rfdistsmodel-final-precompute.RData")
