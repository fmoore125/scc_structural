## setwd("~/research/scciams/scc_structural/")

source("src/analysis/multivariate_prep.R")

## Drop invalid points, convert to DF
distreg <- as.data.frame(subset(distreg, draw > 0))
## Subsample, to accommodate new columns in memory
distreg <- distreg[sample(nrow(distreg), floor(nrow(distreg) / 10)),]

## Make sure we have just either numerical or 2-factor options
sapply(names(distreg)[-1], function(col) length(unique(distreg[, col])))

## Make sure 0-1 columns are factors
for (pred in c('declining', 'dicemodel', 'fundmodel', 'pagemodel', 'backstop', 'failure'))
    distreg[, pred] <- factor(distreg[, pred])

## Construct all interactions for LASSO
preds <- c('Earth_system_struc', 'Tipping.Points_struc', 'Tipping.Points2_struc', 'Persistent...Growth.Damages_struc',
           'Epstein.Zin_struc', 'Ambiguity.Model.Uncertainty_struc', 'Limitedly.Substitutable.Goods_struc', 'Inequality.Aversion_struc',
           'Learning_struc', 'TFP.Growth_param', 'Population.Growth_param', 'Emissions.Growth_param', 'Transient.Climate.Response_param',
           'Carbon.Cycle2_param', 'Equilibrium.Climate.Sensitivity_param', 'Tipping.Point.Magnitude_param', 'Damage.Function_param',
           'Adaptation.Rates_param', 'Income.Elasticity_param', 'Constant.Discount.Rate_param', 'EMUC2_param', 'PRTP2_param',
           'Risk.Aversion..EZ.Utility._param', 'sccyear_from2020', 'discountrate', 'declining', 'dicemodel', 'fundmodel', 'pagemodel',
           'backstop', 'failure', 'log.scc.synth', 'missing.scc.synth')

allpreds <- c()
todemean <- c('draw')
for (pp in 1:length(preds)) {
    print(pp / length(preds))

    pred <- preds[pp]
    if (is.numeric(distreg[, pred])) {
        distreg[, paste0(pred, '2')] <- distreg[, pred]^2
        allpreds <- c(allpreds, paste0(c(pred, paste0(pred, '2')), '.dm'))
        todemean <- c(todemean, pred, paste0(pred, '2'))
    } else {
        allpreds <- c(allpreds, pred)
    }

    for (pp2 in pp:length(preds)) {
        pred2 <- preds[pp2]
        if (pred == pred2)
            next # keep this logic, so don't need to handle pp+1 > length case

        if (is.numeric(distreg[, pred]) && is.numeric(distreg[, pred2])) {
            distreg[, paste0(pred, ':', pred2)] <- distreg[, pred] * distreg[, pred2]
            allpreds <- c(allpreds, paste0(pred, ':', pred2, '.dm'))
            todemean <- c(todemean, paste0(pred, ':', pred2))
        } else if ((is.numeric(distreg[, pred]) && !is.numeric(distreg[, pred2])) ||
                   (!is.numeric(distreg[, pred]) && is.numeric(distreg[, pred2]))) {
            if (is.numeric(distreg[, pred])) {
                numpred <- pred
                facpred <- pred2
            } else {
                numpred <- pred2
                facpred <- pred
            }
            distreg[, paste0(numpred, ':', facpred, '=', unique(distreg[, facpred])[1])] <- distreg[, numpred] * (distreg[, facpred] == unique(distreg[, facpred])[1])
            distreg[, paste0(numpred, ':', facpred, '=', unique(distreg[, facpred])[2])] <- distreg[, numpred] * (distreg[, facpred] == unique(distreg[, facpred])[2])
            allpreds <- c(allpreds, paste0(numpred, ':', facpred, '=', unique(distreg[, facpred])[1], '.dm'), paste0(numpred, ':', facpred, '=', unique(distreg[, facpred])[2], '.dm'))
            todemean <- c(todemean, paste0(numpred, ':', facpred, '=', unique(distreg[, facpred])[1]), paste0(numpred, ':', facpred, '=', unique(distreg[, facpred])[2]))
        } else {
            ## Both factors
            distreg[, paste0(pred, ':', pred2)] <- factor(paste(distreg[, pred], distreg[, pred2]))
            allpreds <- c(allpreds, paste0(pred, ':', pred2))
        }
    }
}

## Removed paper fixed effects
library(lfe)

distregdm <- demeanlist(distreg[, todemean], list(distreg$paper))
for (pred in todemean)
    distreg[, paste0(pred, '.dm')] <- distregdm[, pred]

library(glmnet)

XX <- data.matrix(distreg[, allpreds])
cv_model <- cv.glmnet(XX[1:ceiling(877261 / 2), ], distreg$draw[1:ceiling(877261 / 2)], alpha=1)
## cv_model <- cv.glmnet(XX[, 1:327], distreg$draw, alpha=1)
plot(cv_model)
best_lambda <- cv_model$lambda.min

best_model <- glmnet(XX, distreg$draw, alpha=1, lambda=best_lambda)

