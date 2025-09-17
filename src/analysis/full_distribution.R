## Translates quantile information into Monte Carlo distribution draws for all literature review observations.

## setwd("~/research/scciams/scc_structural")

library(readxl)
library(tidyverse)
library(viridisLite)
library(data.table)

source("src/analysis/find_distribution.R")
source("src/data_cleaining_scripts/cleaning_master.R")

if (F) {
    ## Make table of available quantiles
    tbl <- data.frame()
    for (col in which(names(dat) == 'Min'):which(names(dat) == 'Max')) {
        tbl <- rbind(tbl, data.frame(quantile=names(dat)[col], count=sum(!is.na(dat[, col])), percent=paste0(round(mean(!is.na(dat[, col])) * 100, 1), '%')))
    }
    library(xtable)
    print(xtable(tbl), include.rownames=F)
}

set.seed(12345)

all.qs <- c(0,0.001,0.01, .025, .05, .1, .17, .25, .5, .75, .83, .9, .95, .975, .99,0.999, 1)
all.as.cols <- which(names(dat) == 'Min'):which(names(dat) == 'Max')

#start by generating distributions for each row
sumstats <- data.frame()
dists=list()
for (ii in 1:nrow(dat)) {
    print(ii)
    all.as <- t(dat[ii, all.as.cols])
    ## add a minimum quantile value for some PAGE values to prevent fitting algorithm from producing spurious negative values
    if(dat$`Damage Function Info: Model, Commonly-Used Function, or Function`[ii]%in%c("PAGE2009","PAGE2002")&sum(is.na(all.as[1:2]))==2){ all.as[2]=0}
    qs <- all.qs[!is.na(all.as)]
    as <- all.as[!is.na(all.as)]
    mu <- dat$`Central Value ($ per ton CO2)`[ii]

    if (is.na(mu) && length(qs) == 0) {
        sumstats <- rbind(sumstats, data.frame(ii, central=NA, score=NA, solution="no data"))
        next
    }

    dists[[ii]] <- generate.pdf(mu, qs, as, 1e6)
    if (is.null(dists[[ii]])) {
        sumstats <- rbind(sumstats, data.frame(ii, central=NA, score=NA, solution=last.solution))
        next
    }

    score <- score.dist.draws(mu, qs, as, dists[[ii]])
    sumstats <- rbind(sumstats, data.frame(ii, central=get.central(mu, qs, as), score, solution=last.solution))
}

sumstats$solshort <- sapply(sumstats$solution, function(ss) strsplit(ss, " ")[[1]][1])
ggplot(subset(sumstats, !(solshort %in% c('delta', 'discrete'))), aes(ifelse(score / abs(central) < 1, score / abs(central), NA), fill=solshort)) +
    geom_histogram() + scale_y_continuous(expand=c(0, 0)) +
    scale_x_continuous(expand=c(0, 0)) + scale_fill_discrete("Distribution:") +
    theme_bw() + xlab("RMSE / Central Estimate")
ggsave("figures/distsoln.pdf", width=6.5, height=3.5)

if (F) {
    ## Other checks on distributions
    table(results$solution)

    sum(results$score / results$central > 1, na.rm=T)
    sum(results$score / results$central < .01, na.rm=T)
    sum(rowSums(!is.na(dat[, which(names(dat) == 'Min'):which(names(dat) == 'Max')])) == 0 & !is.na(dat$`Central Value ($ per ton CO2)`))
}


coauthorweights=read.csv(file="src/analysis/paper_covariance/paperweightings.csv")
coauthorweights$prob=coauthorweights$weight/sum(coauthorweights$weight)
citationweights=read.csv(file="src/analysis/paper_covariance/citations.csv")
citationweights$prob=citationweights$normalized_peryear/sum(citationweights$normalized_peryear)

dat$ID_number=as.integer(dat$ID_number)
papers=unique(dat$ID_number)

## set both to false for unweighted distribution, set one to false and the other to true for either citation or coauthor weighting
nsamp=1e7
for (weighting_option in c('unweighted', 'coauthors', 'citations')) {
    weighting_coauthors <- weighting_option == 'coauthors'
    weighting_citations <- weighting_option == 'citations'

    dist=matrix(nrow=nsamp,ncol=2)

    for(i in 1:nsamp){
        if(i%%10000==0) print(i)

        if(weighting_coauthors==FALSE&weighting_citations==FALSE) paper=sample(papers,1) #if no independence weighting, sample papers with equal probability
        if(weighting_coauthors==TRUE&weighting_citations==FALSE) paper=sample(coauthorweights$ID_number,1,prob=coauthorweights$prob) #weigthing is inversely proportional to degree of shared authorship
        if(weighting_coauthors==FALSE&weighting_citations==TRUE) paper=sample(citationweights$ID_number,1,prob=citationweights$prob) #weigthing is proportional to citations

        ## draw from rows for each paper
        rows=which(dat$ID_number==paper)
        row=ifelse(length(rows)==1,rows,sample(rows,1))

        draw=sample(dists[[row]],1)
        dist[i,]=c(draw,row)
    }

    colnames(dist)=c("draw","row")
    if(weighting_coauthors==FALSE&weighting_citations==FALSE) fwrite(dist,file="outputs/distribution_v2.csv")
    if(weighting_coauthors==TRUE&weighting_citations==FALSE) fwrite(dist,file="outputs/distribution_coauthorweighted_v2.csv")
    if(weighting_coauthors==FALSE&weighting_citations==TRUE) fwrite(dist,file="outputs/distribution_citationweighted_v2.csv")
}

