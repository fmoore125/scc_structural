library(kSamples)

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

map.tree <- function(tree, func) {
    results <- func(tree)
    for (branch in unique(names(tree$children)))
        results <- c(results, map.tree(tree$children[[branch]], func))

    results
}

## dist may be NULL to use pre-computed quantiles
predict.tree <- function(tree, datpred, dist, ndraw=1e3) {
    if (nrow(datpred) > 1) {
        draws <- c()
        for (ii in 1:nrow(datpred))
            draws <- c(draws, predict.tree(tree, datpred[ii,], dist, ndraw=ndraw))
        return(draws)
    }

    if (!is.na(tree$split) && tree$split == 'terminal') {
        if (is.null(dist))
            return(tree$quants)
        else
            return(sample(dist$draw[dist$row %in% tree$rows], ndraw, replace=T))
    }

    if (!is.na(tree$split) && tree$split == 'categorical')
        branch <- as.character(datpred[, tree$col])
    else if (is.numeric(tree$split) || is.na(tree$split))
        branch <- get.numeric.branches(datpred[, tree$col], tree$split)
    else
        return(NULL) # unknown split

    if (is.null(tree$children[[branch]])) { # branch not available in training data
        if (is.null(dist))
            return(tree$quants)
        else
            return(sample(dist$draw[dist$row %in% tree$rows], ndraw, replace=T))
    }
    return(predict.tree(tree$children[[branch]], datpred, dist, ndraw=ndraw))
}

apply.tree <- function(tree, func, init, branch, branches) {
    results <- func(init, branch, branches, tree)
    if (!is.na(tree$split) && tree$split != 'terminal') {
        for (branch in unique(names(tree$children)))
            results <- rbind(results, apply.tree(tree$children[[branch]], func, results[1,], branch, unique(names(tree$children))))
    }

    results
}

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

all.qs <- c(0, 0.001, 0.01, .025, .05, .1, .17, .25, .5, .75, .83, .9, .95, .975, .99, 0.999, 1)

## dist may be NULL to use pre-computed quantiles
predict.forest <- function(forest, datpred, dist, ndraw=1e5, quants=all.qs) {
    qmat <- matrix(NA, length(forest), length(quants))
    for (ii in 1:length(forest)) {
        if (!is.list(forest[[ii]]))
            next
        if (is.null(dist))
            qvals <- predict.tree(forest[[ii]], datpred, NULL, ndraw=ndraw)
        else {
            draws <- predict.tree(forest[[ii]], datpred, dist, ndraw=ndraw)
            qvals <- quantile(draws, quants)
        }
        qmat[ii,] <- qvals
    }

    colMeans(qmat, na.rm=T)
}

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

pre.compute.tree <- function(tree, dist, ndraw=1e5, quants=all.qs) {
    draws <- sample(dist$draw[dist$row %in% tree$rows], ndraw, replace=T)
    tree$quants <- quantile(draws, quants)
    for (branch in unique(names(tree$children)))
        tree$children[[branch]] <- pre.compute.tree(tree$children[[branch]], dist, ndraw=ndraw, quants=quants)
    tree
}

pre.compute.forest <- function(forest, dist, ndraw=1e5, quants=all.qs) {
    new.forest <- list()
    for (ii in 1:length(forest)) {
        print(ii)
        if (!is.list(forest[[ii]]))
            new.forest[[ii]] <- forest[[ii]]
        else
            new.forest[[ii]] <- pre.compute.tree(forest[[ii]], dist, ndraw=ndraw, quants=quants)
    }
    new.forest
}
