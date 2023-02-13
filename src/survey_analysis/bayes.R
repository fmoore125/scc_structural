## setwd("~/research/scciams/scc_structural")

library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

load("data/expert_survey/fig2surveydata.rdat")

df <- data.frame()
for (cc in 2:ncol(fig2dat_qual)) {
    values <- fig2dat_qual[!is.na(fig2dat_qual[, cc]), cc]

    level = ifelse(values == "Strongly Disagree", 1,
              ifelse(values == "Disagree", 2,
              ifelse(values == "Neither Agree nor Disagree", 3,
              ifelse(values == "Agree", 4,
              ifelse(values == "Strongly Agree", 5, NA)))))

    df <- rbind(df, data.frame(expert=fig2dat_qual[!is.na(fig2dat_qual[, cc]), 1],
                               question=names(fig2dat_qual)[cc], level))
}

df$expert <- factor(df$expert)
df$question <- factor(df$question)

logit <- function(p) log(p / (1 - p))
inv.logit <- function(x) exp(x)/(1+exp(x))

library(reshape2)
library(dplyr)
library(ggplot2)

stan.code <- "
data {
  int<lower=0> N; // # answers
  int<lower=0> J; // # experts
  int<lower=0> K; // # questions;
  int<lower=1,upper=J> expert[N];
  int<lower=1,upper=K> question[N];
  int<lower=1,upper=5> level[N]; // 1 - 5
  ordered[4] hyperc;
  real<lower=0> hcsprior;
}
parameters {
  // hyper assessment of agreement
  real mu[K];
  real<lower=0> tau[K];
  // disagreement about scale
  real<lower=1> hypercscale;
  real<lower=0, upper=1> sigmac;
  // expert understanding of scale
  ordered[4] c[J];
  // expert-question opinion shifter
  real eta[N];
}
transformed parameters {
  vector[N] theta;
  for (ii in 1:N)
    theta[ii] = mu[question[ii]] + tau[question[ii]] * eta[ii];
}
model {
  target += normal_lpdf(eta | 0, 1);
  target += exponential_lpdf(hypercscale | hcsprior);
  tau ~ cauchy(0, 1);
  for (ii in 1:N) {
    level[ii] ~ ordered_logistic(theta[ii] * hypercscale, c[expert[ii]] * hypercscale);
  }
  for (jj in 1:J)
    c[jj] ~ normal(hyperc, sigmac);
}"

fit <- stan(model_code=stan.code, data=list(N=nrow(df), J=length(levels(df$expert)),
                                            K=length(levels(df$question)), expert=as.numeric(df$expert),
                                            question=as.numeric(df$question), level=df$level,
                                            hyperc=c(logit(.2), logit(.4), logit(.6), logit(.8)), hcsprior=.001))

la <- extract(fit, permute=T)

pdf <- melt(la$mu)
names(pdf)[2] <- 'question'
pdf$question <- levels(df$question)[pdf$question]
pdf$prob <- inv.logit(pdf$value)

df$answer <- c('Strongly Disagree', 'Disagree', 'Neither...', 'Agree', 'Strongly Agree')[df$level]
df$answer <- factor(df$answer, levels=c('Strongly Disagree', 'Disagree', 'Neither...', 'Agree', 'Strongly Agree'))
df$prob <- inv.logit(colMeans(la$theta))

df$naive <- (df$level - .5) / 5
df2 <- df %>% group_by(question) %>% summarize(naive1=mean(naive), naive2=inv.logit(mean(logit(naive))))

## pdf <- read.csv("outputs/meta-analysis-distribution.csv")

ggplot(pdf, aes(prob)) +
    facet_wrap(~ question) +
    geom_vline(data=df, aes(xintercept=prob, colour=answer), alpha=.5) +
    geom_vline(data=df2, aes(xintercept=naive2, colour='Simple average'), size=1) +
    geom_density() + scale_x_continuous(expand=c(0, 0)) + theme_bw() +
    scale_colour_manual("Expert response", breaks=c('Strongly Disagree', 'Disagree', 'Neither...', 'Agree', 'Strongly Agree', 'Simple average'), values=c('#e41a1c', '#ff7f00', '#4daf4a', '#377eb8', '#984ea3', '#808080')) +
    xlab("Probability of inclusion")
ggsave("figures/meta-analysis.pdf", width=9, height=4)

write.csv(pdf, "outputs/meta-analysis-distribution.csv", row.names=F)

library(xtable)
pp <- summary(fit)
xtable(pp$summary)

## Look at covariance

wpdf <- dcast(pdf, iterations ~ question)
