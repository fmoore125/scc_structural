## setwd("~/research/scciams/scc_structural")

source("src/data_cleaining_scripts/cleaning_master.R")
source("src/analysis/all_scc_lib.R")
library(ggplot2)
library(kSamples)

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

dat$`SCC Year` <- as.numeric(dat$`SCC Year`)

dat$log.scc.2020usd <- log(dat$`Central Value ($ per ton CO2)`)
dat$log.scc.2020usd[!is.finite(dat$log.scc.2020usd)] <- NA

source("src/analysis/damage_funcs_lib.R")
dat <- multivar.prep(dat)

dist <- read.csv(file="outputs/distribution_v2.csv")

#drop outlier Nordhaus row
todrop=which(dat$`Central Value ($ per ton CO2)`>70000)
if(is.finite(todrop)) dist=dist[-which(dist$row==todrop),]

dat$row <- 1:nrow(dat)
dat$`Earth system` <- ifelse(dat$`Carbon Cycle` == "1.0" | dat$`Carbon Cycle` == "1", "1",
                      ifelse(dat$`Climate Model` == "1.0" | dat$`Climate Model` == "1", "1", "0"))
dat$`Inequality Aversion`[dat$`Inequality Aversion` == "Calibrated"] <- "1.0"
dat$`Inequality Aversion`[dat$`Inequality Aversion` == "1.0"] <- "1"
dat$`Persistent / Growth Damages`[dat$`Persistent / Growth Damages` == "Calibrated"] <- "1.0"
dat$`Persistent / Growth Damages`[dat$`Persistent / Growth Damages` == "1.0"] <- "1"
dat$`Tipping Points2`[dat$`Tipping Points2` == "-1.0" | dat$`Tipping Points2` == "-1"] <- "0"
dat$`TFP Growth`[is.na(dat$`TFP Growth`)] <- "0"
dat$`Population Growth`[is.na(dat$`Population Growth`)] <- "0"
dat$`Emissions Growth`[is.na(dat$`Emissions Growth`)] <- "0"
dat$`Transient Climate Response`[is.na(dat$`Transient Climate Response`)] <- "0"
dat$`Carbon Cycle2`[is.na(dat$`Carbon Cycle2`)] <- "0"
dat$`Equilibrium Climate Sensitivity`[is.na(dat$`Equilibrium Climate Sensitivity`)] <- "0"
dat$`Tipping Point Magnitude`[is.na(dat$`Tipping Point Magnitude`)] <- "0"
dat$`Damage Function`[is.na(dat$`Damage Function`)] <- "0"
dat$`Adaptation Rates`[is.na(dat$`Adaptation Rates`)] <- "0"
dat$`Income Elasticity`[is.na(dat$`Income Elasticity`)] <- "0"
dat$`Constant Discount Rate`[is.na(dat$`Constant Discount Rate`)] <- "0"
dat$`EMUC2`[is.na(dat$`EMUC2`)] <- "0"
dat$`PRTP2`[is.na(dat$`PRTP2`)] <- "0"
dat$`Risk Aversion (EZ Utility)`[is.na(dat$`Risk Aversion (EZ Utility)`)] <- "0"
dat$`Declining Discounting?`[is.na(dat$`Declining Discounting?`)] <- "0"
dat$log.scc.synth[dat$missing.scc.synth] <- NA # We can handle this

incrows <- !is.na(dat$sccyearformerge) & dat$sccyearformerge <= 2100 &
    dat$row %in% dist$row # Some have no distribution

score.split <- function(dist2) {
    if (nrow(dist2) > 1e4) {
        ## Max of three checks
        max.score <- 0
        for (ii in 1:3) {
            score <- score.split(dist2[sample(nrow(dist2), 1e4),])
            max.score <- max(score, max.score)
            if (max.score > 0.1)
                return(max.score)
        }
        return(max.score)
    }
    ns <- table(dist2$branch)
    if (min(ns) > 1e3)
        dist2 <- dist2[sample(nrow(dist2), 1e3),]
    else if (max(ns) > 1e3) {
        included <- rep(F, nrow(dist2))
        rands <- runif(nrow(dist2))
        for (br in unique(dist2$branch)) {
            currrows <- dist2$branch == br
            currn <- sum(currrows, na.rm=T)
            if (currn <= 1e3)
                included[currrows] <- T
            else
                included[currrows & rands < 1e3 / currn] <- T
        }
        dist2 <- dist2[included,]
    }

    ad.res <- ad.test(draw ~ branch, data=dist2, method='asymptotic')
    prod(sqrt(ad.res$ad[, 3])) ** (min(ns) / sum(ns)) # Penalize if uneven split
}

## dist2 <- dist
## dist2$branch <- dist2$row %in% dat$row[dat$Year < 2015]
## ## dist2$branch <- runif(nrow(dist2)) > .5
## ## dist2$branch <- dist2$row %in% dat$row[dat$Year < 2021]
## score.split(dist2)

get.numeric.branches <- function(datcol, splitpt) {
    if (is.na(splitpt))
        as.character(is.na(datcol)) # TRUE or FALSE
    else
        ifelse(is.na(datcol), "NA", as.character(datcol < splitpt)) # TRUE, FALSE, NA
}

## dist2 <- dist
## dist2$branch <- get.numeric.branches(dat$discountrate, 2.419026)[dist2$row]
## dist2$branch <- get.numeric.branches(dat$sccyearformerge, 2015)[dist2$row]
## score.split(dist2)

score.column.split <- function(dat, incrows, dist2, col, minleaf, verbose=T) {
    values <- dat[incrows, col]
    if (is.numeric(values)) {
        allquants <- quantile(values, c(0, .1, .25, .5, .75, .9, runif(1), 1), na.rm=T)
        splitpts <- unique(allquants[2:7])
        splitpts <- splitpts[splitpts != allquants[1] & splitpts != allquants[8]]
        best.score <- 1
        best.split <- NULL
        for (splitpt in splitpts) {
            if (verbose)
                print(c(col, splitpt))
            datbranch <- get.numeric.branches(dat[, col], splitpt)
            if (length(unique(datbranch)) == 1)
                next
            if (min(table(datbranch)) < minleaf)
                next
            dist2$branch <- datbranch[dist2$row]
            score <- score.split(dist2)
            if (score < best.score) {
                best.score <- score
                best.split <- splitpt
            }
        }
        if (any(is.na(values))) {
            if (verbose)
                print(c(col, NA))
            datbranch <- get.numeric.branches(dat[, col], NA)
            if (length(unique(datbranch)) == 1)
                next
            if (min(table(datbranch)) < minleaf)
                next
            dist2$branch <- datbranch[dist2$row]
            score <- score.split(dist2)
            if (score < best.score) {
                best.score <- score
                best.split <- NA
            }
        }
        if (!is.null(best.split))
            list(split=best.split, score=best.score)
        else
            NULL
    } else {
        if (verbose)
            print(col)
        if (length(unique(values)) == 1)
            return(list(split='categorical', score=1))
        if (min(table(dat[dist2$row, col])) < minleaf)
            NULL
        dist2$branch <- dat[dist2$row, col]
        list(split='categorical', score=score.split(dist2))
    }
}

## score.column.split(dat, T, dist, "Limitedly-Substitutable Goods", 5)
## score.column.split(dat, T, dist, "discountrate", 5)

make.tree <- function(dat, incrows, dist2, cols, ancestors=0, mtry=floor(sqrt(length(cols))), minleaf=5, minsplit=ceiling(0.01*nrow(dat)), verbose=T) {
    if (sum(incrows) <= minsplit)
        return(list(split="terminal", reason='minsplit', rows=which(incrows)))
    best.split <- list(score=1)
    for (col in sample(cols, mtry)) {
        this.split <- tryCatch({
            score.column.split(dat, incrows, dist2, col, minleaf, verbose=verbose)
        }, error=function(e) {
            NULL
        })
        if (!is.null(this.split) && this.split$score < best.split$score) {
            this.split$column <- col
            best.split <- this.split
        }
    }
    if (best.split$score > 0.05)
        return(list(split="terminal", reason='nosplit', rows=which(incrows)))

    if (verbose) {
        print(ancestors)
        print(best.split)
        print(sum(incrows))
    }

    ## Label every row, even if not included, so branches[dist2$row] works
    if (is.numeric(dat[, best.split$column])) {
        branches <- get.numeric.branches(dat[, best.split$column], best.split$split)
    } else {
        branches <- dat[, best.split$column]
    }
    dist2$branch <- branches[dist2$row]

    children <- list()
    for (branch in unique(branches[incrows])) {
        if (is.na(branch))
            next
        if (sum(incrows & !is.na(branches) & branches == branch) == 0) {
            print(best.split)
            print(branch)
            FAIL
        }
        child <- make.tree(dat, incrows & !is.na(branches) & branches == branch,
                           dist2[!is.na(dist2$branch) & dist2$branch == branch,], cols, ancestors=ancestors+1,
                           mtry=mtry, minsplit=minsplit, verbose=verbose)
        children[[as.character(branch)]] <- child
    }

    best.split$children <- children
    best.split$rows <- which(incrows)
    best.split
}

cols <- c("Tipping Points", "Tipping Points2", "Persistent / Growth Damages", "Epstein-Zin",
          "Ambiguity/Model Uncertainty", "Limitedly-Substitutable Goods", "Inequality Aversion",
          "Learning", "Earth system", "TFP Growth", "Population Growth", "Emissions Growth",
          "Transient Climate Response", "Carbon Cycle2", "Equilibrium Climate Sensitivity",
          "Tipping Point Magnitude", "Damage Function", "Adaptation Rates", "Income Elasticity",
          "Constant Discount Rate", "EMUC2", "PRTP2", "Risk Aversion (EZ Utility)",
          "Backstop Price?", "Declining Discounting?", "Market Only Damages", "Other Market Failure?",
          "sccyearformerge", "discountrate", "log.scc.synth", "Year")
if (F) {
    for (col in cols)
        print(c(col, unique(dat[, col])[1:min(3, length(unique(dat[, col])))], "NAs:", sum(is.na(dat[, col]))))
}
## Make tree with all values
tree <- make.tree(dat, incrows, dist, cols, mtry=length(cols), minleaf=5)
save(tree, file="outputs/rfdistsmodel.RData")

map.tree <- function(tree, func) {
    results <- func(tree)
    for (branch in unique(names(tree$children)))
        results <- c(results, map.tree(tree$children[[branch]], func))

    results
}

sum(map.tree(tree, function(tree) 1))

predict.tree <- function(tree, datpred, dist, ndraw=1e3) {
    if (nrow(datpred) > 1) {
        draws <- c()
        for (ii in 1:nrow(datpred))
            draws <- c(draws, predict.tree(tree, datpred[ii,], dist, ndraw=ndraw))
        return(draws)
    }

    if (!is.na(tree$split) && tree$split == 'terminal')
        return(sample(dist$draw[dist$row %in% tree$rows], ndraw, replace=T))

    if (!is.na(tree$split) && tree$split == 'categorical')
        branch <- as.character(datpred[, tree$col])
    else if (is.numeric(tree$split) || is.na(tree$split))
        branch <- get.numeric.branches(datpred[, tree$col], tree$split)
    else
        return(NULL) # unknown split

    if (is.null(tree$children[[branch]])) # branch not available in training data
        return(sample(dist$draw[dist$row %in% tree$rows], ndraw, replace=T))
    return(predict.tree(tree$children[[branch]], datpred, dist, ndraw=ndraw))
}

all.qs <- c(0, 0.001, 0.01, .025, .05, .1, .17, .25, .5, .75, .83, .9, .95, .975, .99, 0.999, 1)

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

apply.tree <- function(tree, func, init, branch, branches) {
    results <- func(init, branch, branches, tree)
    if (!is.na(tree$split) && tree$split != 'terminal') {
        for (branch in unique(names(tree$children)))
            results <- rbind(results, apply.tree(tree$children[[branch]], func, results[1,], branch, unique(names(tree$children))))
    }

    results
}

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

## Show the distributional split
plot.split <- function(tree, dist2, incrows) {
    ## Label every row, even if not included, so branches[dist2$row] works
    if (is.numeric(dat[, tree$column])) {
        branches <- get.numeric.branches(dat[, tree$column], tree$split)
    } else {
        branches <- dat[, tree$column]
    }
    dist2$branch <- branches[dist2$row]

    print(score.split(dist2))

    ggplot(dist2[incrows,], aes(draw)) +
        geom_density(aes(colour=branch))
}

plot.split(tree, dist, incrows) + xlim(-1000, 2000)

make.forest <- function(dat, incrows, dist2, cols, trees=500, mtry=floor(sqrt(length(cols))), minleaf=5, minsplit=ceiling(0.01*nrow(dat))) {
    alltrees <- list()
    for (tt in 1:trees) {
        print(tt)
        bagrows <- sample(which(incrows), sum(incrows), replace=T)
        incrows2 <- incrows & (dat$row %in% bagrows) # Don't actually up-weight multiply sampled rows
        tree <- make.tree(dat, incrows2, dist2, cols, mtry=mtry, minleaf=minleaf, minsplit=minsplit, verbose=F)
        alltrees[[tt]] <- tree
    }
    alltrees
}

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

predict.forest <- function(forest, datpred, dist, ndraw=1e5, quants=all.qs) {
    qmat <- matrix(NA, length(forest), length(quants))
    for (ii in 1:length(forest)) {
        if (!is.list(forest[[ii]]))
            next
        draws <- predict.tree(forest[[ii]], datpred, dist, ndraw=ndraw)
        qvals <- quantile(draws, quants)
        qmat[ii,] <- qvals
    }

    colMeans(qmat, na.rm=T)
}

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

calc.varimport.tree <- function(tree, baseimport=1) {
    results <- data.frame(column=tree$column, import=baseimport)
    for (child in names(tree$children)) {
        childtree <- tree$children[[child]]
        if (is.na(childtree$split) || childtree$split != 'terminal') {
            rows <- calc.varimport.tree(childtree, baseimport=baseimport * length(childtree$rows) / length(tree$rows))
            results <- rbind(results, rows)
        }
    }
    results
}

calc.varimport.tree(tree) %>% group_by(column) %>% summarize(import=max(import)) %>% arrange(desc(import))

calc.varimport.forest <- function(forest) {
    alltrees <- data.frame()
    count <- 0
    for (ii in 1:length(forest)) {
        if (!is.list(forest[[ii]]))
            next
        results <- calc.varimport.tree(forest[[ii]])
        alltrees <- rbind(alltrees, calc.varimport.tree(tree) %>% group_by(column) %>% summarize(import=max(import)))
        count <- count + 1
    }

    alltrees %>% group_by(column) %>% summarize(import=sum(import) / count)
}

varimport <- calc.varimport.forest(forest)
varimport$column <- factor(varimport$column, levels=varimport$column[order(varimport$import)])

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
               'sccyearformerge'='Quantity: SCC Year')
varimport$label <- label.chg[varimport$column]
varimport$label <- factor(varimport$label, levels=varimport$label[order(varimport$import)])

ggplot(varimport, aes(label, import)) +
    coord_flip() +
    geom_col() + theme_bw() +
    scale_y_continuous("Variable Importance", labels=scales::percent) + xlab(NULL)
ggsave("figures/rfdists-varimport.pdf", width=6.5, height=5)

## Construct synthetic SCC

idealdat <- dat # include all rows

## set obvious ones
for (cc in which(names(dat) == 'TFP Growth'):which(names(dat) == 'Risk Aversion (EZ Utility)'))
    idealdat[, cc] <- "1"
idealdat$`Backstop Price?` <- "0"
idealdat$`Declining Discounting?` <- "1"
idealdat$`Market Only Damages` <- "0"
idealdat$`Other Market Failure?` <- "0"
## for damage-function-based scc - draw from literature values
idealdat$log.scc.synth[idealdat$missing.scc.synth] <- sample(idealdat$log.scc.synth[!idealdat$missing.scc.synth],
                                                            sum(idealdat$missing.scc.synth),replace=TRUE)
idealdat$Year <- 2020
idealdat$sccyearformerge <- 2020

bayespost <- read.csv("data/expert_survey/meta-analysis-distribution.csv")
bayespost2 <- dcast(bayespost, iterations ~ question, value.var='prob')
bayespost3 <- bayespost2[sample(nrow(bayespost2), nrow(dat), replace=T),]
bayespost3[, -1] <- ifelse(matrix(runif((ncol(bayespost3)-1)*nrow(bayespost3)), nrow(bayespost3), ncol(bayespost3)-1) < bayespost3[, -1], "1", "0")

idealdat <- cbind(idealdat[-c(which(names(dat) == 'Carbon Cycle'):which(names(dat) == 'Learning'), which(names(dat) == 'Earth system'))], bayespost3)
names(idealdat)[names(idealdat) == 'Tipping Points: Climate'] <- 'Tipping Points'
names(idealdat)[names(idealdat) == 'Tipping Points: Damages'] <- 'Tipping Points2'
names(idealdat)[names(idealdat) == 'Limited Substitutability'] <- 'Limitedly-Substitutable Goods'
names(idealdat)[names(idealdat) == 'Earth System'] <- 'Earth system'

discountsurvey <- read.csv("data/Drupp_et_al_2018_AEJ_Constant_SDR.csv")
idealdat$discountrate <- sample(discountsurvey$SDR[!is.na(discountsurvey$SDR)], nrow(idealdat), replace=TRUE)

predict.forest.all <- function(expdat, incrows, outfile) {
    pdf <- data.frame()
    for (ii in sample(which(incrows))) {
        if (ii %in% pdf$row)
            next
        print(ii)
        quants.pred <- predict.forest(forest, expdat[ii,], dist)
        pdf <- rbind(pdf, data.frame(row=ii, quant=all.qs, pred=quants.pred))
    }

    ## Construct Monte Carlo draws
    allsamp <- c()
    for (row in unique(pdf$row)) {
        inv.cdf <- approx(pdf$quant[pdf$row == row], pdf$pred[pdf$row == row])
        uu <- runif(1000)
        samples <- inv.cdf$y[findInterval(uu, inv.cdf$x)]
        allsamp <- c(allsamp, samples)
    }

    print(quantile(allsamp, c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1.)))
    print(mean(allsamp))

    save(pdf, allsamp, file=outfile)
}

predict.forest.all(idealdat, incrows, "outputs/rf_experiments/RFD_best.RData")

##----A. No structural changes (classic DICE assumptions)-----

## set obvious ones
dicedat <- idealdat
for (cc in which(names(dicedat) == 'TFP Growth'):which(names(dicedat) == 'Risk Aversion (EZ Utility)'))
    dicedat[, cc] <- "0"
dicedat$`Backstop Price?` <- "0"
dicedat$`Declining Discounting?` <- "0"
dicedat$`Market Only Damages` <- "0"
dicedat$`Other Market Failure?` <- "0"

## for damage-function-based scc - draw from literature values
rel <- dat$`Damage Function Info: Model, Commonly-Used Function, or Function`%in%c("DICE-2007","DICE-2013R","DICE-2016R2","DICE 2007","DICE 2010","DICE 2013","DICE 2013R","DICE 2016","DICE2007","DICE2010","DICE2013","DICE2016R")
dicedat$log.scc.synth <- sample(dat$log.scc.synth[rel & !dat$missing.scc.synth],
                                nrow(dat),replace=TRUE)

for (cc in which(names(dicedat) == 'Ambiguity/Model Uncertainty'):which(names(dicedat) == 'Tipping Points2'))
    dicedat[, cc] <- "0"

dicedat$discountrate <- 4.6

predict.forest.all(dicedat, 1:nrow(dat) %in% sample(which(incrows), 100), "outputs/rf_experiments/RFD_A_dice.RData")

##----B. EPA assumptions -----
epadat <- dicedat

## structural changes to Earth System
epadat$`Earth system` <- "1"
## parametric uncertainty in tfp growth, pop growth, earth system and damage functions
epadat$`TFP Growth` <- "1"
epadat$`Population Growth` <- "1"
epadat$`Emissions Growth` <- "1"
epadat$`Transient Climate Response` <- "1"
epadat$`Carbon Cycle2` <- "1"
epadat$`Equilibrium Climate Sensitivity` <- "1"
epadat$`Damage Function` <- "1"
## central discount rate of 2%
epadat$discountrate <- 2
## damage function from Howard and Sterner
rel <- dat$`Damage Function Info: Model, Commonly-Used Function, or Function`%in%c("HowardSterner","HowardSterner (0.007438*T^2)")
epadat$log.scc.synth <- sample(dat$log.scc.synth[rel & !dat$missing.scc.synth],
                               nrow(dat),replace=TRUE)

predict.forest.all(epadat, 1:nrow(dat) %in% sample(which(incrows), 100), "outputs/rf_experiments/RFD_B_epa.RData")

##---- C. All structural changes and no structural changes
alldat <- idealdat
for (cc in which(names(alldat) == 'Ambiguity/Model Uncertainty'):which(names(alldat) == 'Tipping Points2'))
    alldat[, cc] <- "1"

predict.forest.all(alldat, 1:nrow(dat) %in% sample(which(incrows), 100), "outputs/rf_experiments/RFD_C_all.RData")

