# Analysis of synthetic data
# From Qian et al 2021 PsychMehtods

rm(list = ls())

library(tidyverse)
library(xtable)
library(geepack)
library(ggplot2)
source("xgeepack.R")
source("estimate.R")
source("estimators.R")


generate_synthetic_data <- function(N, T = 210, p = 0.4) {
  data <- data.frame(
    userid = rep(1:N, each = T),
    decision.index.nogap = rep(1:T, times = N),
    study.day.nogap = rep(rep(0:((T/5)-1), rep(5, (T/5))), times = N),
    study.day.square = rep(rep((0:((T/5)-1))**2, rep(5, (T/5))), times = N),
    avail = rbinom(N * T, 1, 1)
  )
  
  data$send <- with(data, ifelse(avail == 1, rbinom(N * T, 1, p), 0))
  # Model parameters
  alpha <- c(2, 0.01, 0.01) # Example values for B_t coefficients
  alpha <- c(2, 0, 0)
  beta <- c(0.2, 0, 0) # Example values for Z_t coefficients
  B_t <- function(t) c(1, floor((t-1)/5), floor((t-1)/5)^2)
  Z_t <- function(t) c(1, floor((t-1)/5), floor((t-1)/5)^2)
  
  # Generate steps
  data$steps <- numeric(N*T)
  for (i in 1:N) {
    for (t in 1:T) {
      idx <- (i - 1) * T + t
      B <- B_t(t)
      Z <- Z_t(t)
      data$steps[idx] <- sum(B * alpha) + sum(Z * beta) * (data$send[idx] - p) + rnorm(1, 0, 1)
    }
  }
  
  return(data)
}

# continuous outcome estimator
p_values <- numeric(0)
for (i in 1:200) {
  print(i)
  set.seed(i+1001)
  synthetic_data <- generate_synthetic_data(N = 43, T = 210)
  head(synthetic_data)
  
  synthetic_data$"(Intercept)" <- 1
  synthetic_data$"I(send - 0.4)" <- synthetic_data$send - 0.4
  
  xmat <- synthetic_data %>%
    transmute("(Intercept)" = .$"(Intercept)",
              "days" = .$"study.day.nogap",
              "days_sq" = .$"study.day.square",
              "I(send - 0.4)" = .$"I(send - 0.4)")
  
  fit_model1 <- geese.glm(x = as.matrix(xmat),
                          y = synthetic_data$steps, w = synthetic_data$avail, id = as.factor(synthetic_data$userid),
                          family = gaussian(), corstr = "independence")
  p_values <- c(p_values, estimate(fit_model1)[4,8])
}


# binary outcome estimator
p_values_2 <- numeric(0)
for (i in 1:500) {
  print(i)
  set.seed(i+1001)
  
  sample_size <- 34
  synthetic_data <- generate_synthetic_data(N = sample_size, T = 50)
  synthetic_data$Y <- as.integer(synthetic_data$steps > 2)
  synthetic_data$prob_A <- 0.4
  
  fit_wcls <- weighted_centered_least_square(
    dta = synthetic_data,
    id_varname = "userid",
    decision_time_varname = "decision.index.nogap",
    treatment_varname = "send",
    outcome_varname = "Y",
    control_varname = NULL,
    moderator_varname = NULL,
    avail_varname = "avail",
    rand_prob_varname = "prob_A",
    rand_prob_tilde_varname = NULL,
    rand_prob_tilde = 0.4,
    estimator_initial_value = NULL
  )
  p_value <- 2 * pt(abs(fit_wcls$beta_hat) / fit_wcls$beta_se_adjusted, sample_size - 1, lower.tail = FALSE)
  p_values_2 <- c(p_values_2, p_value)
}
