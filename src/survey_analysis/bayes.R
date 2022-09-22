## setwd("~/research/scciams/scc_structural")

library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

load("data/expert_survey/fig2surveydata.rdat")

## stan.code <- "
## data {
##   int<lower=0> J; // # experts
##   real plo[J]; // low prob. true
##   real phi[J]; // high prob. true
## }
## parameters {
##   real<lower=0> alpha; // hyper beta alpha
##   real<lower=0> beta; // hyper beta beta
##   real eta[J]; // true expert opinion
## }
## model {
##   eta ~ uniform(plo, phi);
##   eta ~ beta(alpha, beta);
## }"

stan.code <- "
data {
  int<lower=0> J; // # experts
  int<lower=1,upper=5> level[J]; // 1 - 5
}
parameters {
  real mu;
  real<lower=0> sigma;
  ordered[4] c;
  real eta[J]; // expert opinion shifter
}
model {
  for (jj in 1:J)
    level[jj] ~ ordered_logistic(eta[jj], c);
  eta ~ normal(mu, sigma);
}"


for (cc in 2:ncol(fig2dat_qual)) {
    values <- fig2dat_qual[!is.na(fig2dat_qual[, cc]), cc]
    ## plotrue = ifelse(values == "Strongly Disagree", 0,
    ##           ifelse(values == "Disagree", 0.05,
    ##           ifelse(values == "Neither Agree nor Disagree", 0.25,
    ##           ifelse(values == "Agree", 0.5,
    ##           ifelse(values == "Strongly Agree", 0.75, NA)))))
    ## phitrue = ifelse(values == "Strongly Disagree", 0.25,
    ##           ifelse(values == "Disagree", 0.5,
    ##           ifelse(values == "Neither Agree nor Disagree", 0.75,
    ##           ifelse(values == "Agree", 0.95,
    ##           ifelse(values == "Strongly Agree", 1, NA)))))

    level = ifelse(values == "Strongly Disagree", 1,
              ifelse(values == "Disagree", 2,
              ifelse(values == "Neither Agree nor Disagree", 3,
              ifelse(values == "Agree", 4,
              ifelse(values == "Strongly Agree", 5, NA)))))

    fit <- stan(model_code=stan.code, data=list(J=length(ptrue), level=level))
}

stan.code <- "
data {
  int<lower=0> N; // # answers
  int<lower=0> J; // # experts
  int<lower=0> K; // # questions;
  int<lower=1,upper=J> expert[N];
  int<lower=1,upper=K> question[N];
  int<lower=1,upper=5> level[N]; // 1 - 5
}
parameters {
  real mu[K];
  real<lower=0> sigma[K];
  ordered[4] c;
  real eta[N]; // expert-question opinion shifter
}
model {
  for (ii in 1:N) {
    level[ii] ~ ordered_logistic(eta[ii], c);
    eta[ii] ~ normal(mu[question[ii]], sigma[question[ii]]);
  }
}"

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

fit <- stan(model_code=stan.code, data=list(N=nrow(df), J=length(levels(df$expert)),
                                            K=length(levels(df$question)), expert=as.numeric(df$expert),
                                            question=as.numeric(df$question), level=df$level))

stan.code <- "
data {
  int<lower=0> N; // # answers
  int<lower=0> J; // # experts
  int<lower=0> K; // # questions;
  int<lower=1,upper=J> expert[N];
  int<lower=1,upper=K> question[N];
  int<lower=1,upper=5> level[N]; // 1 - 5
}
parameters {
  // hyper assessment of agreement
  real mu[K];
  real<lower=0> sigma[K];
  // hyper understanding of scale
  ordered[4] hyperc;
  real<lower=0> sigmac;
  // expert understanding of scale
  ordered[4] c[J];
  // expert-question opinion shifter
  real eta[N];
}
model {
  for (ii in 1:N) {
    level[ii] ~ ordered_logistic(eta[ii], c[expert[ii]]);
    eta[ii] ~ normal(mu[question[ii]], sigma[question[ii]]);
  }
  for (jj in 1:J)
    c[jj] ~ normal(hyperc, sigmac);

}"

fit <- stan(model_code=stan.code, data=list(N=nrow(df), J=length(levels(df$expert)),
                                            K=length(levels(df$question)), expert=as.numeric(df$expert),
                                            question=as.numeric(df$question), level=df$level))

logit <- function(p) log(p / (1 - p))

stan.code <- "
data {
  int<lower=0> N; // # answers
  int<lower=0> J; // # experts
  int<lower=0> K; // # questions;
  int<lower=1,upper=J> expert[N];
  int<lower=1,upper=K> question[N];
  int<lower=1,upper=5> level[N]; // 1 - 5
  ordered[4] hyperc;
}
parameters {
  // hyper assessment of agreement
  real mu[K];
  real<lower=0> tau[K];
  // disagreement about scale
  real<lower=0, upper=1> sigmac; // force c to roughy follow hyperc
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
  for (ii in 1:N) {
    level[ii] ~ ordered_logistic(theta[ii], c[expert[ii]]);
  }
  for (jj in 1:J)
    c[jj] ~ normal(hyperc, sigmac);
}"

fit <- stan(model_code=stan.code, data=list(N=nrow(df), J=length(levels(df$expert)),
                                            K=length(levels(df$question)), expert=as.numeric(df$expert),
                                            question=as.numeric(df$question), level=df$level,
                                            hyperc=c(logit(.2), logit(.4), logit(.6), logit(.8))))

la <- extract(fit, permute=T)

inv.logit <- function(x) exp(x)/(1+exp(x))

library(reshape2)
pdf <- melt(la$mu)
names(pdf)[2] <- 'question'
pdf$question <- levels(df$question)[pdf$question]
pdf$prob <- inv.logit(pdf$value)

df$answer <- c('Strongly Disagree', 'Disagree', 'Neither...', 'Agree', 'Strongly Agree')[df$level]
df$answer <- factor(df$answer, levels=c('Strongly Disagree', 'Disagree', 'Neither...', 'Agree', 'Strongly Agree'))
df$prob <- inv.logit(colMeans(la$theta))

ggplot(pdf, aes(prob)) +
    facet_wrap(~ question) +
    geom_vline(data=df, aes(xintercept=prob, colour=answer), alpha=.5) +
    geom_density() + scale_x_continuous(expand=c(0, 0)) + theme_bw() +
    scale_colour_discrete("Expert response") + xlab("Probability of inclusion")

df[df$question == "Earth System" & df$answer == "Strongly Disagree",]
which(levels(df$expert) == "570")
which(levels(df$question) == "Earth System")

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
  for (ii in 1:N) {
    level[ii] ~ ordered_logistic(theta[ii] * hypercscale, c[expert[ii]] * hypercscale);
  }
  for (jj in 1:J)
    c[jj] ~ normal(hyperc, sigmac);
}"

fit <- stan(model_code=stan.code, data=list(N=nrow(df), J=length(levels(df$expert)),
                                            K=length(levels(df$question)), expert=as.numeric(df$expert),
                                            question=as.numeric(df$question), level=df$level,
                                            hyperc=c(logit(.2), logit(.4), logit(.6), logit(.8)), hcsprior=.01))

la <- extract(fit, permute=T)

pdf <- melt(la$mu)
names(pdf)[2] <- 'question'
pdf$question <- levels(df$question)[pdf$question]
pdf$prob <- inv.logit(pdf$value)

df$answer <- c('Strongly Disagree', 'Disagree', 'Neither...', 'Agree', 'Strongly Agree')[df$level]
df$answer <- factor(df$answer, levels=c('Strongly Disagree', 'Disagree', 'Neither...', 'Agree', 'Strongly Agree'))
df$prob <- inv.logit(colMeans(la$theta))

df$naive <- (df$level - .5) / 5

library(dplyr)
df2 <- df %>% group_by(question) %>% summarize(naive=mean(naive))

ggplot(pdf, aes(prob)) +
    facet_wrap(~ question) +
    geom_vline(data=df, aes(xintercept=prob, colour=answer), alpha=.5) +
    geom_vline(data=df2, aes(xintercept=naive, colour='Simple average')) +
    geom_density() + scale_x_continuous(expand=c(0, 0)) + theme_bw() +
    scale_colour_manual("Expert response", breaks=c('Strongly Disagree', 'Disagree', 'Neither...', 'Agree', 'Strongly Agree', 'Simple average'), values=c('#e41a1c', '#ff7f00', '#4daf4a', '#377eb8', '#984ea3', '#000000')) +
    xlab("Probability of inclusion")

## Don't change expert's opinions

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
  for (ii in 1:N) {
    level[ii] ~ ordered_logistic(theta[ii] * hypercscale, c[expert[ii]] * hypercscale);
  }
  for (jj in 1:J)
    c[jj] ~ normal(hyperc, sigmac);
}"

fit <- stan(model_code=stan.code, data=list(N=nrow(df), J=length(levels(df$expert)),
                                            K=length(levels(df$question)), expert=as.numeric(df$expert),
                                            question=as.numeric(df$question), level=df$level,
                                            hyperc=c(logit(.2), logit(.4), logit(.6), logit(.8)), hcsprior=.01))
