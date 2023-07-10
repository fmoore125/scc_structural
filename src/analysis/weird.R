## setwd("~/research/scciams/scc_structural/")

library(lfe)
library(MASS)

do.bagging <- 'none'

source("src/analysis/multivariate_prep.R")

mod_paperfe=felm(log(draw)~Earth_system_struc+Tipping.Points_struc+Tipping.Points2_struc+Persistent...Growth.Damages_struc+
            Epstein.Zin_struc+Ambiguity.Model.Uncertainty_struc+Limitedly.Substitutable.Goods_struc+Inequality.Aversion_struc+
            Learning_struc+TFP.Growth_param+Population.Growth_param+Emissions.Growth_param+Transient.Climate.Response_param+Carbon.Cycle2_param+
            Equilibrium.Climate.Sensitivity_param+Tipping.Point.Magnitude_param+Damage.Function_param+Adaptation.Rates_param+Income.Elasticity_param+
            Constant.Discount.Rate_param+EMUC2_param+PRTP2_param+Risk.Aversion..EZ.Utility._param+sccyear_from2020+I(sccyear_from2020^2)+discountrate+I(discountrate^2)+declining+dicemodel+fundmodel+pagemodel+backstop+failure+log.scc.synth + missing.scc.synth | paper, data=distreg[which(distreg$draw>0),])

valid <- !is.nan(coef(mod_paperfe))
coefmu <- coef(mod_paperfe)[valid]
coefvv <- vcov(mod_paperfe)[valid, valid]

discountsurvey=read.csv("data/Drupp_et_al_2018_AEJ_Constant_SDR.csv")
bayespost=read.csv("data/expert_survey/meta-analysis-distribution.csv")
bayespost.wide <- dcast(bayespost, iterations ~ question, value.var='prob')

## Estimate how large the difference will be:
sum(coef(mod_paperfe)[grepl("_struc", names(coef(mod_paperfe))) & valid]) # 1.386852
sum(coef(mod_paperfe)[grepl("_param", names(coef(mod_paperfe))) & valid]) # 2.437915

bases <- coef(mod_paperfe)[names(coef(mod_paperfe)) == 'discountrate'] * dat$discountrate +
    coef(mod_paperfe)[names(coef(mod_paperfe)) == 'I(discountrate^2)'] * dat$discountrate^2

newdiscs <- sample(discountsurvey$SDR[-which(is.na(discountsurvey$SDR))], nrow(dat), replace=T)
news <- coef(mod_paperfe)[names(coef(mod_paperfe)) == 'discountrate'] * newdiscs +
    coef(mod_paperfe)[names(coef(mod_paperfe)) == 'I(discountrate^2)'] * newdiscs

mean(news - bases, na.rm=T) # 0.4920397


num.struc <- (distreg$Tipping.Points_struc == 'Yes') + (distreg$Tipping.Points2_struc == 'Yes') +
    (distreg$Persistent...Growth.Damages_struc == 'Yes') + (distreg$Epstein.Zin_struc == 'Yes') +
    (distreg$Ambiguity.Model.Uncertainty_struc == 'Yes') + (distreg$Limitedly.Substitutable.Goods_struc == 'Yes') +
    (distreg$Inequality.Aversion_struc == 'Yes') + (distreg$Learning_struc == 'Yes') +
    (distreg$Earth_system_struc == 'Yes')

num.param <- (distreg$TFP.Growth_param == '1') + (distreg$Population.Growth_param == '1') +
    (distreg$Emissions.Growth_param == '1') + (distreg$Transient.Climate.Response_param == '1') +
    (distreg$Carbon.Cycle2_param == '1') + (distreg$Equilibrium.Climate.Sensitivity_param == '1') +
    (distreg$Tipping.Point.Magnitude_param == '1') + (distreg$Damage.Function_param == '1') +
    (distreg$Constant.Discount.Rate_param == '1') + (distreg$EMUC2_param == '1') +
    (distreg$PRTP2_param == '1') + (distreg$Risk.Aversion..EZ.Utility._param == '1')

quantile(num.struc)
quantile(num.param)

if (do.bagging == 'none') {
    distreg$draw.post <- NA # All entries
} else if (do.bagging == 'half') {
    distreg$draw.half <- NA # Maximum observed in the data, across struc and param
} else if (do.bagging == 'ones') {
    distreg$draw.ones <- NA # 75th percentile is 1 from each struc and param
}

for (ii in 1:nrow(dat)) {
    print(ii / nrow(dat))

    jjs <- which(distreg$row == ii)
    if (length(jjs) == 0)
        next

    discrates <- sample(discountsurvey$SDR[-which(is.na(discountsurvey$SDR))], length(jjs), replace=T)
    bayesset <- as.matrix(bayespost.wide[sample(nrow(bayespost.wide), length(jjs), replace=T),])
    bayesbol <- as.data.frame(matrix(runif(nrow(bayesset) * (ncol(bayesset)-1)), nrow(bayesset), ncol(bayesset)-1) < bayesset[, -1])

    coeffs <- mvrnorm(length(jjs), coefmu, coefvv)

    preddf <- data.frame(Earth_system_struc=c(!is.na(dat$`Carbon Cycle`[ii]) || !is.na(dat$`Climate Model`[ii]),
                                              bayesbol$`Earth System`),
                         Tipping.Points_struc=c(!is.na(dat$`Tipping Points`[ii]),
                                                bayesbol$`Tipping Points: Climate`),
                         Tipping.Points2_struc=c(!is.na(dat$`Tipping Points2`[ii]),
                                                 bayesbol$`Tipping Points: Damages`),
                         Persistent...Growth.Damages_struc=c(!is.na(dat$`Persistent / Growth Damages`[ii]),
                                                             bayesbol$`Persistent / Growth Damages`),
                         Epstein.Zin_struc=c(!is.na(dat$`Epstein-Zin`[ii]),
                                             bayesbol$`Epstein-Zin`),
                         Ambiguity.Model.Uncertainty_struc=c(!is.na(dat$`Ambiguity/Model Uncertainty`[ii]),
                                                             bayesbol$`Ambiguity/Model Uncertainty`),
                         Limitedly.Substitutable.Goods_struc=c(!is.na(dat$`Limitedly-Substitutable Goods`[ii]),
                                                               bayesbol$`Limited Substitutability`),
                         Inequality.Aversion_struc=c(!is.na(dat$`Inequality Aversion`[ii]),
                                                     bayesbol$`Inequality Aversion`),
                         Learning_struc=c(!is.na(dat$`Learning`[ii]),
                                          bayesbol$`Learning`),
                         TFP.Growth_param=c(!is.na(dat$`TFP Growth`[ii]), rep(T, length(jjs))),
                         Population.Growth_param=c(!is.na(dat$`Population Growth`[ii]), rep(T, length(jjs))),
                         Emissions.Growth_param=c(!is.na(dat$`Emissions Growth`[ii]), rep(T, length(jjs))),
                         ## Transient.Climate.Response_param=c(!is.na(dat$`Transient Climate Response`[ii]), rep(T, length(jjs))),
                         Carbon.Cycle2_param=c(!is.na(dat$`Carbon Cycle2`[ii]), rep(T, length(jjs))),
                         Equilibrium.Climate.Sensitivity_param=c(!is.na(dat$`Equilibrium Climate Sensitivity`[ii]), rep(T, length(jjs))),
                         Tipping.Point.Magnitude_param=c(!is.na(dat$`Tipping Point Magnitude`[ii]),
                                                         bayesbol$`Tipping Points: Climate` | bayesbol$`Tipping Points: Damages`),
                         Damage.Function_param=c(!is.na(dat$`Damage Function`[ii]), rep(T, length(jjs))),
                         Adaptation.Rates_param=c(!is.na(dat$`Adaptation Rates`[ii]), rep(T, length(jjs))),
                         Income.Elasticity_param=c(!is.na(dat$`Income Elasticity`[ii]), rep(T, length(jjs))),
                         ## Constant.Discount.Rate_param=c(!is.na(dat$`Constant Discount Rate`[ii]), rep(!is.na(dat$`Constant Discount Rate (%)`[ii]), length(jjs))),
                         MUC2_param=c(!is.na(dat$`EMUC2`[ii]), rep(!is.na(dat$`PRTP`[ii]), length(jjs))),
                         PRTP2_param=c(!is.na(dat$`PRTP2`[ii]), rep(!is.na(dat$`EMUC`[ii]), length(jjs))),
                         Risk.Aversion..EZ.Utility._param=c(!is.na(dat$`Risk Aversion (EZ Utility)`[ii]),
                                                            bayesbol$`Epstein-Zin`),
                         sccyear_from2020=c(as.numeric(dat$`SCC Year`[ii])-2020, rep(0, length(jjs))),
                         sccyear_from2020=c(as.numeric(dat$`SCC Year`[ii])-2020, rep(0, length(jjs)))^2,
                         discountrate=c(dat$discountrate[ii], discrates),
                         discountrate2=c(dat$discountrate[ii], discrates)^2,
                         declining=rep(!is.na(dat$`Declining Discounting?`[ii]), length(jjs)+1),
                         dicemodel=rep(grepl('DICE', dat$`Base IAM (if applicable)`[ii]) || grepl('DICE', dat$`IAM Calibrated To (if applicable)`[ii]), length(jjs)+1),
                         fundmodel=rep(grepl('FUND', dat$`Base IAM (if applicable)`[ii]) || grepl('FUND', dat$`IAM Calibrated To (if applicable)`[ii]), length(jjs)+1),
                         pagemodel=rep(grepl('PAGE', dat$`Base IAM (if applicable)`[ii]) || grepl('PAGE', dat$`IAM Calibrated To (if applicable)`[ii]), length(jjs)+1),
                         backstop=rep(!is.na(dat$`Backstop Price?`[ii]), length(jjs)+1),
                         failure=rep(!is.na(dat$`Other Market Failure?`[ii]), length(jjs)+1),
                         log.scc.synth=rep(dat$log.scc.synth[ii], length(jjs)+1),
                         missing.scc.synth=rep(dat$missing.scc.synth[ii], length(jjs)+1))

    ## XXX: Try randomly reducing the number of differences we allow (like bagging)
    if (do.bagging == 'half') {
        bagged <- matrix(runif(ncol(preddf) * length(jjs)), length(jjs), ncol(preddf)) > .5
        for (kk in 1:ncol(preddf))
            preddf[c(F, !bagged[, kk]), kk] <- preddf[1, kk]
    } else if (do.bagging == 'ones') {
        ## bagged <- t(matrix(!grepl("_struc|_param", names(preddf)), ncol(preddf), length(jjs)))
        ## bagged[sample(grep("_struc", names(preddf)), 1)] <- T
        ## bagged[sample(grep("_param", names(preddf)), 1)] <- T
        keep.struc <- sample(grep("_struc", names(preddf)), length(jjs), replace=T)
        for (kk in grep("_struc", names(preddf)))
            preddf[c(F, keep.struc != kk), kk] <- preddf[1, kk]
        keep.param <- sample(grep("_param", names(preddf)), length(jjs), replace=T)
        for (kk in grep("_param", names(preddf)))
            preddf[c(F, keep.param != kk), kk] <- preddf[1, kk]
    }

    preds <- as.matrix(preddf) %*% coeffs[1, ] # Just use first row
    logdiffs <- preds[-1] - preds[1]
    if (do.bagging == 'none') {
        distreg$draw.post[jjs] <- distreg$draw[jjs] * exp(logdiffs)
    } else if (do.bagging == 'half') {
        distreg$draw.half[jjs] <- distreg$draw[jjs] * exp(logdiffs)
    } else if (do.bagging == 'ones') {
        distreg$draw.ones[jjs] <- distreg$draw[jjs] * exp(logdiffs)
    }

    ## for (kk in 1:length(jjs)) {
    ##     ## XXX: Try randomly reducing the number of differences we allow (like bagging)
    ##     bagged <- runif(ncol(preddf)) > .5
    ##     mypreddf <- preddf[c(1, kk+1),]
    ##     mypreddf[2, !bagged] <- mypreddf[1, !bagged]

    ##     logdiff <- diff(as.matrix(mypreddf) %*% coeffs[kk, ])
    ##     distreg$draw.half[jjs[kk]] <- distreg$draw[jjs[kk]] * exp(logdiff[1, 1])
    ## }
}

quantile(distreg$draw, na.rm=T)
mean(distreg$draw, na.rm=T)


if (do.bagging == 'none') {
    distreg$draw.post[distreg$draw.post < quantile(distreg$draw.post, .005) | distreg$draw.post > quantile(distreg$draw.post, .995)] <- NA
    quantile(distreg$draw.post, na.rm=T)
    mean(distreg$draw.post, na.rm=T)
} else if (do.bagging == 'half') {
    distreg$draw.half[distreg$draw.half < quantile(distreg$draw.half, .005) | distreg$draw.half > quantile(distreg$draw.half, .995)] <- NA

    quantile(distreg$draw.half, na.rm=T)
    mean(distreg$draw.half, na.rm=T)
} else {
    distreg$draw.ones[distreg$draw.ones < quantile(distreg$draw.ones, .005) | distreg$draw.ones > quantile(distreg$draw.ones, .995)] <- NA

    quantile(distreg$draw.ones, na.rm=T)
    mean(distreg$draw.ones, na.rm=T)
}

ggplot(distreg) +
    ## geom_histogram(aes(x=draw, after_stat(density), fill='Original'), position='dodge', alpha=0.5) +
    ## geom_histogram(aes(x=draw.half, after_stat(density), fill='Half Bags'), position='dodge', alpha=0.5) +
    ## geom_histogram(aes(x=draw.ones, after_stat(density), fill='One + One'), position='dodge', alpha=0.5) +
    geom_density(aes(x=draw, colour='Original')) +
    geom_density(aes(x=draw.half, colour='Half Bags')) +
    geom_density(aes(x=draw.ones, colour='One + One')) +
    scale_x_continuous("SCC (% per ton CO2)", limits=c(-100, 2000), expand=c(0, 0)) +
    scale_y_continuous(expand=c(0, 0)) +
    ## scale_fill_discrete(NULL) +
    scale_colour_discrete(NULL) +
    theme_bw() + theme(legend.justification=c(1,1), legend.position=c(1,1))
