## setwd("~/research/scciams/scc_structural/")

library(lfe)
library(glmnet)

source("src/data_cleaining_scripts/cleaning_master.R")
source("src/analysis/all_scc_lib.R")
source("src/analysis/damage_funcs_lib.R")

do.forcepreds <- T

df <- multivar.prep(dat)
df <- subset(df, !is.na(df$`Central Value ($ per ton CO2)`) & !is.na(df$discountrate) & df$`Central Value ($ per ton CO2)` > 0)
df$`Persistent / Growth Damages` <- df$`Persistent / Growth Damages` != '0'
df$`Tipping Points2` <- df$`Tipping Points2` == '1'
df$`Inequality Aversion` <- df$`Inequality Aversion` != '0'
df$`Earth System` <- df$`Carbon Cycle` != '0' | df$`Climate Model` != '0'
df$`Declining Discounting?` <- !is.na(df$`Declining Discounting?`)
df$`TFP Growth` <- !is.na(df$`TFP Growth`)
df$`Population Growth` <- !is.na(df$`Population Growth`)
df$`Emissions Growth` <- !is.na(df$`Emissions Growth`)
df$`Transient Climate Response` <- !is.na(df$`Transient Climate Response`)
df$`Carbon Cycle2` <- !is.na(df$`Carbon Cycle2`)
df$`Equilibrium Climate Sensitivity` <- !is.na(df$`Equilibrium Climate Sensitivity`)
df$`Tipping Point Magnitude` <- !is.na(df$`Tipping Point Magnitude`)
df$`Damage Function` <- !is.na(df$`Damage Function`)
df$`Adaptation Rates` <- !is.na(df$`Adaptation Rates`)
df$`Income Elasticity` <- !is.na(df$`Income Elasticity`)
df$`Constant Discount Rate` <- !is.na(df$`Constant Discount Rate`)
df$`EMUC2` <- !is.na(df$`EMUC2`)
df$`PRTP2` <- !is.na(df$`PRTP2`)
df$`Risk Aversion (EZ Utility)` <- !is.na(df$`Risk Aversion (EZ Utility)`)

df$dicemodel <- grepl("DICE", df$`Base IAM (if applicable)`) | grepl("DICE", df$`IAM Calibrated To (if applicable)`)
df$fundmodel <- grepl("FUND", df$`Base IAM (if applicable)`) | grepl("FUND", df$`IAM Calibrated To (if applicable)`)
df$pagemodel <- grepl("PAGE", df$`Base IAM (if applicable)`) | grepl("PAGE", df$`IAM Calibrated To (if applicable)`)

df$sccyear_from2020 <- as.numeric(df$`SCC Year`)-2020
df$logscc <- log(df$`Central Value ($ per ton CO2)`)

## Make sure we have just either numerical or 2-factor options
sapply(names(df)[-1], function(col) length(unique(df[, col])))
sapply(names(df)[-1], function(col) is.numeric(df[, col]))
sapply(names(df)[-1], function(col) any(is.na(df[, col])))

## Construct all interactions for LASSO
preds <- c('Earth System', 'Tipping Points', 'Tipping Points2', 'Persistent / Growth Damages',
           'Epstein-Zin', 'Ambiguity/Model Uncertainty', 'Limitedly-Substitutable Goods', 'Inequality Aversion',
           'Learning', 'TFP Growth', 'Population Growth', 'Emissions Growth', 'Transient Climate Response',
           'Carbon Cycle2', 'Equilibrium Climate Sensitivity', 'Tipping Point Magnitude', 'Damage Function',
           'Adaptation Rates', 'Income Elasticity', 'Constant Discount Rate', 'EMUC2', 'PRTP2',
           'Risk Aversion (EZ Utility)', 'sccyear_from2020', 'discountrate', 'Declining Discounting?',
           'dicemodel', 'fundmodel', 'pagemodel', 'Backstop Price?', 'Other Market Failure?', 'log.scc.synth', 'missing.scc.synth')

allpreds <- c()
todemean <- c('logscc')
for (pp in 1:length(preds)) {
    print(pp / length(preds))

    pred <- preds[pp]
    if (is.numeric(df[, pred])) {
        df[, paste0(pred, '2')] <- df[, pred]^2
        if (!do.forcepreds) {
            allpreds <- c(allpreds, paste0(c(pred, paste0(pred, '2')), '.dm'))
            todemean <- c(todemean, pred, paste0(pred, '2'))
        } else {
            allpreds <- c(allpreds, paste0(pred, '2.dm'))
            todemean <- c(todemean, paste0(pred, '2'))
        }
    } else if (!do.forcepreds) {
        allpreds <- c(allpreds, pred)
    }

    for (pp2 in pp:length(preds)) {
        pred2 <- preds[pp2]
        if (pred == pred2)
            next # keep this logic, so don't need to handle pp+1 > length case

        if (is.numeric(df[, pred]) && is.numeric(df[, pred2])) {
            df[, paste0(pred, ':', pred2)] <- df[, pred] * df[, pred2]
            allpreds <- c(allpreds, paste0(pred, ':', pred2, '.dm'))
            todemean <- c(todemean, paste0(pred, ':', pred2))
        } else if ((is.numeric(df[, pred]) && !is.numeric(df[, pred2])) ||
                   (!is.numeric(df[, pred]) && is.numeric(df[, pred2]))) {
            if (is.numeric(df[, pred])) {
                numpred <- pred
                facpred <- pred2
            } else {
                numpred <- pred2
                facpred <- pred
            }
            df[, paste0(numpred, ':', facpred, '=', unique(df[, facpred])[1])] <- df[, numpred] * (df[, facpred] == unique(df[, facpred])[1])
            allpreds <- c(allpreds, paste0(numpred, ':', facpred, '=', unique(df[, facpred])[1], '.dm'))
            todemean <- c(todemean, paste0(numpred, ':', facpred, '=', unique(df[, facpred])[1]))
        } else {
            ## Both factors
            df[, paste0(pred, ':', pred2)] <- factor(paste(df[, pred], df[, pred2]))
            allpreds <- c(allpreds, paste0(pred, ':', pred2))
        }
    }
}

## Removed paper fixed effects
df$ID_number <- factor(df$ID_number)
if (!do.forcepreds) {
    dfdm <- demeanlist(df[, todemean], list(df$ID_number))
    for (pred in todemean)
        df[, paste0(pred, '.dm')] <- dfdm[, pred]
} else {
    for (pred in todemean) {
        mod <- felm(as.formula(paste0('`', pred, '` ~ `', paste(preds, collapse='` + `'), '` | ID_number')), data=df)
        df[, paste0(pred, '.dm')] <- mod$resid
    }
}

XX <- data.matrix(df[, allpreds])
cv_model <- cv.glmnet(XX, df$logscc.dm, alpha=1)
plot(cv_model)
best_lambda <- cv_model$lambda.min

best_model <- glmnet(XX, df$logscc.dm, alpha=1, lambda=best_lambda)

included <- as.logical(!(abs(coef(best_model)) < 1e-10))
results <- data.frame(pred=rownames(coef(best_model))[included], coef=coef(best_model)[included])

if (do.forcepreds) {
    usepreds <- c(preds, gsub("\\.dm$", "", results$pred))
    usepreds <- usepreds[-which(usepreds %in% c('Transient Climate Response', 'Constant Discount Rate'))]

    ## All of my final factors have only 2 levels, so turn to 0/1
    numdf <- df
    for (pred in usepreds) {
        if (is.numeric(numdf[, pred]))
            next
        else if (is.logical(numdf[, pred]))
            numdf[, pred] <- as.numeric(numdf[, pred])
        else if (is.factor(numdf[, pred])) {
            stopifnot(all(levels(numdf[, pred]) == c('0', '1')))
            numdf[, pred] <- as.numeric(numdf[, pred] == '1')
        } else if (is.character(numdf[, pred])) {
            if (all(unique(numdf[, pred]) == c('0', '1')))
                numdf[, pred] <- as.numeric(numdf[, pred] == '1')
            else if (all(unique(numdf[, pred]) == c('0', '1.0')))
                numdf[, pred] <- as.numeric(numdf[, pred] == '1.0')
            else
                OTHERERROR
        } else {
            ERROR
        }
    }
    lassomod <- felm(as.formula(paste0('logscc ~ `', paste(usepreds, collapse='` + `'), '` | ID_number')), data=numdf, keepCX=T)
    summary(lassomod)
    newres <- data.frame(pred=usepreds, coef=coef(lassomod))

    ## Now try to use it
    source("src/analysis/multivariate_prep.R")

    discountsurvey=read.csv("data/Drupp_et_al_2018_AEJ_Constant_SDR.csv")
    bayespost=read.csv("data/expert_survey/meta-analysis-distribution.csv")
    bayespost.wide <- dcast(bayespost, iterations ~ question, value.var='prob')

    bayespost.rows <- sample(1:nrow(bayespost.wide), nrow(distreg), replace=T)

    preddf1 <- tibble(`Earth System`=as.numeric(bayespost.wide$`Earth System`[bayespost.rows] > runif(nrow(distreg))),
                      `Tipping Points`=as.numeric(bayespost.wide$`Tipping Points: Climate`[bayespost.rows] > runif(nrow(distreg))),
                      `Tipping Points2`=as.numeric(bayespost.wide$`Tipping Points: Damages`[bayespost.rows] > runif(nrow(distreg))),
                      `Persistent / Growth Damages`=as.numeric(bayespost.wide$`Persistent / Growth Damages`[bayespost.rows] > runif(nrow(distreg))),
                      `Epstein-Zin`=as.numeric(bayespost.wide$`Epstein-Zin`[bayespost.rows] > runif(nrow(distreg))),
                      `Ambiguity/Model Uncertainty`=as.numeric(bayespost.wide$`Ambiguity/Model Uncertainty`[bayespost.rows] > runif(nrow(distreg))),
                      `Limitedly-Substitutable Goods`=as.numeric(bayespost.wide$`Ambiguity/Model Uncertainty`[bayespost.rows] > runif(nrow(distreg))),
                      `Inequality Aversion`=as.numeric(bayespost.wide$`Inequality Aversion`[bayespost.rows] > runif(nrow(distreg))),
                      `Learning`=as.numeric(bayespost.wide$`Learning`[bayespost.rows] > runif(nrow(distreg))))
    preddf2 <- tibble(`TFP Growth`=1, `Population Growth`=1, `Emissions Growth`=1, `Transient Climate Response`=1,
                      `Carbon Cycle2`=1, `Equilibrium Climate Sensitivity`=1,
                      `Tipping Point Magnitude`=as.numeric(preddf1$`Tipping Points` == 1 | preddf1$`Tipping Points2` == 1),
                      `Damage Function`=1, `Adaptation Rates`=1, `Income Elasticity`=1,
                      `Constant Discount Rate`=as.numeric(!is.na(dat$`Constant Discount Rate (%)`[distreg$row])),
                      EMUC2=as.numeric(!is.na(dat$EMUC[distreg$row])), PRTP2=as.numeric(!is.na(dat$PRTP[distreg$row])),
                      `Risk Aversion (EZ Utility)`=preddf1$`Epstein-Zin`, sccyear_from2020=0,
                      discountrate=sample(discountsurvey$SDR[-which(is.na(discountsurvey$SDR))], nrow(distreg), replace=T),
                      `Declining Discounting?`=distreg$declining,
                      dicemodel=distreg$dicemodel, fundmodel=distreg$fundmodel, pagemodel=distreg$pagemodel,
                      `Backstop Price?`=distreg$backstop, `Other Market Failure?`=distreg$failure,
                      log.scc.synth=distreg$log.scc.synth, missing.scc.synth=as.numeric(distreg$missing.scc.synth))
    preddf <- cbind(preddf1, preddf2)

    for (col in rownames(coef(best_model))[included]) {
        if (col == "sccyear_from20202.dm")
            preddf$`sccyear_from20202` <- preddf$sccyear_from2020^2
        else if (col == "discountrate2.dm")
            preddf$`discountrate2` <- preddf$discountrate^2
        else if (col == "log.scc.synth2.dm")
            preddf$`log.scc.synth2` <- preddf$log.scc.synth^2
        else if (col == 'sccyear_from2020:discountrate.dm')
            preddf$`sccyear_from2020:discountrate` <- 0
        else if (col == 'discountrate:log.scc.synth.dm')
            preddf$`discountrate:log.scc.synth` <- preddf$discountrate * preddf$log.scc.synth
        else {
            parts <- str_split(gsub("\\.dm$", "", col), ":|=")[[1]]
            stopifnot(parts[1] %in% c('discountrate', 'sccyear_from2020', 'log.scc.synth'))
            match <- list('TRUE'=T, 'FALSE'=F, '0'=0, '1'=1)[[parts[3]]]
            preddf[, gsub("\\.dm$", "", col)] <- preddf[, parts[1]] * (preddf[, parts[2]] == match)
        }
    }

    source("~/projects/research-common/R/felm-tools.R")
    fes <- getfe(lassomod)
    dat2 <- dat %>% left_join(fes, by=c('ID_number'='idx'))

    distreg$logdraw.post <- NA
    for (row1 in seq(1, nrow(preddf), by=1e4)) {
        print(row1)
        resls <- predict.felm(lassomod, preddf[row1 - 1 + (1:1e4),], se.fit=T)
        distreg$logdraw.post[row1 - 1 + (1:1e4)] <- rnorm(nrow(resls$fit), resls$fit$fit, resls$se.fit)
    }

    distreg$draw.post <- exp(dat2$effect[distreg$row] + distreg$logdraw.post + var(lassomod$resid)[1, 1] / 2)

    ggplot(distreg) +
        geom_density(aes(x=draw, colour='Original')) +
        geom_density(aes(x=draw.post, colour='LASSO')) +
        scale_x_continuous("SCC (% per ton CO2)", limits=c(-100, 2000), expand=c(0, 0)) +
        scale_y_continuous(expand=c(0, 0)) +
        scale_colour_discrete(NULL) +
        theme_bw() + theme(legend.justification=c(1,1), legend.position=c(1,1))

    distreg$draw.post[distreg$draw.post < quantile(distreg$draw.post, .005, na.rm=T) | distreg$draw.post > quantile(distreg$draw.post, .995, na.rm=T)] <- NA

    quantile(distreg$draw.post, na.rm=T)
    mean(distreg$draw.post, na.rm=T)
}

