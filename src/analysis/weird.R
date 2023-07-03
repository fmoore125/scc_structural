## setwd("~/research/scciams/scc_structural/")

library(lfe)
library(MASS)

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

distreg$draw.post <- NA

for (jj in 1:nrow(distreg)) {
    if (jj %% 1000 == 0)
        print(jj / nrow(distreg))

    ii <- distreg$row[jj]

    discrate <- sample(discountsurvey$SDR[-which(is.na(discountsurvey$SDR))], 1)
    bayesset <- subset(bayespost, iterations == sample(bayespost$iterations, 1))
    bayesset$bool <- runif(nrow(bayesset)) < bayesset$prob

    coeffs <- mvrnorm(1, coefmu, coefvv)

    preddf <- data.frame(Earth_system_struc=c(!is.na(dat$`Carbon Cycle`[ii]) || !is.na(dat$`Climate Model`[ii]),
                                              bayesset$bool[bayesset$question == 'Earth System']),
                         Tipping.Points_struc=c(!is.na(dat$`Tipping Points`[ii]),
                                                bayesset$bool[bayesset$question == 'Tipping Points: Climate']),
                         Tipping.Points2_struc=c(!is.na(dat$`Tipping Points2`[ii]),
                                                 bayesset$bool[bayesset$question == 'Tipping Points: Damages']),
                         Persistent...Growth.Damages_struc=c(!is.na(dat$`Persistent / Growth Damages`[ii]),
                                                             bayesset$bool[bayesset$question == 'Persistent / Growth Damages']),
                         Epstein.Zin_struc=c(!is.na(dat$`Epstein-Zin`[ii]),
                                             bayesset$bool[bayesset$question == 'Epstein-Zin']),
                         Ambiguity.Model.Uncertainty_struc=c(!is.na(dat$`Ambiguity/Model Uncertainty`[ii]),
                                                             bayesset$bool[bayesset$question == 'Ambiguity/Model Uncertainty']),
                         Limitedly.Substitutable.Goods_struc=c(!is.na(dat$`Limitedly-Substitutable Goods`[ii]),
                                                               bayesset$bool[bayesset$question == 'Limited Substitutability']),
                         Inequality.Aversion_struc=c(!is.na(dat$`Inequality Aversion`[ii]),
                                                     bayesset$bool[bayesset$question == 'Inequality Aversion']),
                         Learning_struc=c(!is.na(dat$`Learning`[ii]),
                                          bayesset$bool[bayesset$question == 'Learning']),
                         TFP.Growth_param=c(!is.na(dat$`TFP Growth`[ii]), T),
                         Population.Growth_param=c(!is.na(dat$`Population Growth`[ii]), T),
                         Emissions.Growth_param=c(!is.na(dat$`Emissions Growth`[ii]), T),
                         ## Transient.Climate.Response_param=c(!is.na(dat$`Transient Climate Response`[ii]), T),
                         Carbon.Cycle2_param=c(!is.na(dat$`Carbon Cycle2`[ii]), T),
                         Equilibrium.Climate.Sensitivity_param=c(!is.na(dat$`Equilibrium Climate Sensitivity`[ii]), T),
                         Tipping.Point.Magnitude_param=c(!is.na(dat$`Tipping Point Magnitude`[ii]),
                                                         bayesset$bool[bayesset$question == 'Tipping Points: Climate'] ||
                                                         bayesset$bool[bayesset$question == 'Tipping Points: Damages']),
                         Damage.Function_param=c(!is.na(dat$`Damage Function`[ii]), T),
                         Adaptation.Rates_param=c(!is.na(dat$`Adaptation Rates`[ii]), T),
                         Income.Elasticity_param=c(!is.na(dat$`Income Elasticity`[ii]), T),
                         ## Constant.Discount.Rate_param=c(!is.na(dat$`Constant Discount Rate`[ii]), !is.na(dat$`Constant Discount Rate (%)`[ii])),
                         EMUC2_param=c(!is.na(dat$`EMUC2`[ii]), !is.na(dat$`PRTP`[ii])),
                         PRTP2_param=c(!is.na(dat$`PRTP2`[ii]), !is.na(dat$`EMUC`[ii])),
                         Risk.Aversion..EZ.Utility._param=c(!is.na(dat$`Risk Aversion (EZ Utility)`[ii]),
                                                            bayesset$bool[bayesset$question == 'Epstein-Zin']),
                         sccyear_from2020=c(as.numeric(dat$`SCC Year`[ii])-2020, 0),
                         sccyear_from2020=c(as.numeric(dat$`SCC Year`[ii])-2020, 0)^2,
                         discountrate=c(dat$discountrate[ii], discrate),
                         discountrate2=c(dat$discountrate[ii], discrate)^2,
                         declining=rep(!is.na(dat$`Declining Discounting?`[ii]), 2),
                         dicemodel=rep(grepl('DICE', dat$`Base IAM (if applicable)`[ii]) || grepl('DICE', dat$`IAM Calibrated To (if applicable)`[ii]), 2),
                         fundmodel=rep(grepl('FUND', dat$`Base IAM (if applicable)`[ii]) || grepl('FUND', dat$`IAM Calibrated To (if applicable)`[ii]), 2),
                         pagemodel=rep(grepl('PAGE', dat$`Base IAM (if applicable)`[ii]) || grepl('PAGE', dat$`IAM Calibrated To (if applicable)`[ii]), 2),
                         backstop=rep(!is.na(dat$`Backstop Price?`[ii]), 2),
                         failure=rep(!is.na(dat$`Other Market Failure?`[ii]), 2),
                         log.scc.synth=rep(dat$log.scc.synth[ii], 2),
                         missing.scc.synth=rep(dat$missing.scc.synth[ii], 2))

    logdiff <- diff(as.matrix(preddf) %*% coeffs)
    distreg$draw.post[jj] <- distreg$draw[jj] * exp(logdiff[1, 1])
}
