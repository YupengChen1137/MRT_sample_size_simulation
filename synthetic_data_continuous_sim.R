# Analysis of synthetic data
# From Qian et al 2021 PsychMehtods

rm(list = ls())

library(tidyverse)
library(xtable)
library(geepack)
library(ggplot2)
library(MASS)
source("xgeepack.R")
source("estimate.R")
source("estimators.R")


generate_synthetic_data <- function(N, Time = 210, p = 0.4) {
  data <- data.frame(
    userid = rep(1:N, each = Time),
    decision.index.nogap = rep(1:Time, times = N),
    study.day.nogap = rep(rep(0:((Time/5)-1), rep(5, (Time/5))), times = N),
    study.day.square = rep(rep((0:((Time/5)-1))**2, rep(5, (Time/5))), times = N),
    avail = rbinom(N * Time, 1, 1)
  )
  
  data$send <- with(data, ifelse(avail == 1, rbinom(N * Time, 1, p), 0))
  # Model parameters
  alpha <- c(2, 0.01, 0.01) # Example values for B_t coefficients
  #alpha <- c(2, 0, 0)
  beta <- c(0, 0.00964, -0.000172)
  beta <- c(0, 0.00476, 0) # Example values for Z_t coefficients
  B_t <- function(t) c(1, floor((t-1)/5), floor((t-1)/5)^2)
  Z_t <- function(t) c(1, floor((t-1)/5), floor((t-1)/5)^2)
  
  # Generate steps
  data$steps <- numeric(N*Time)
  for (i in 1:N) {
    for (t in 1:Time) {
      idx <- (i - 1) * Time + t
      B <- B_t(t)
      Z <- Z_t(t)
      data$steps[idx] <- sum(B * alpha) + sum(Z * beta) * (data$send[idx] - p) + rnorm(1, 0, 1)
    }
  }
  
  return(data)
}

# continuous outcome estimator

Time <- 210  # Total number of decision points (42 days * 5 times per day)
rho <- 0.4  # Randomization probability
availability <- 1  # Constant availability

N <- 19
alpha.0 = 0.05; p = 3; q = 3;
multiple = p*(N-q-1)/(N-p-q)
crit_value = multiple*qf((1-alpha.0), df1 = p, df2 = N - p - q)
crit_value_2 = qf((1-alpha.0)/multiple, df1 = p, df2 = N - p - q)

Z <- matrix(0, nrow = Time, ncol = 3)
for (t in 1:Time) {
  Z[t,] <- c(1, floor((t-1)/5), floor((t-1)/5)^2)
}
Q <- t(Z) %*% diag(rep(availability,Time) * rep(rho, Time) * rep((1-rho), Time)) %*% Z

trial <- 200
test_stats <- numeric(0)
for (i in 1:trial) {
  print(i)
  set.seed(i+1001)
  synthetic_data <- generate_synthetic_data(N = N, T = Time)
  
  
  synthetic_data$"(Intercept)" <- 1
  synthetic_data$"I(send - 0.4)" <- synthetic_data$send - 0.4
  
  xmat <- synthetic_data %>%
    transmute("(Intercept)" = .$"(Intercept)",
              "days" = .$"study.day.nogap",
              "days.square" = .$"study.day.square",
              "I(send - 0.4)" = .$"I(send - 0.4)",
              "I(send - 0.4):study.day.nogap" = .$"I(send - 0.4)" * .$"study.day.nogap",
              "I(send - 0.4):study.day.square" = .$"I(send - 0.4)" * .$"study.day.square")
  
  fit_model1 <- geese.glm(x = as.matrix(xmat),
                          y = synthetic_data$steps, w = synthetic_data$avail, id = as.factor(synthetic_data$userid),
                          family = gaussian(), corstr = "independence")
  # alpha_hat <- estimate(fit_model1)[1:3,1]
  beta_hat <- fit_model1$coefficients[4:6]
  
  
  W_sum <- matrix(0,3,3)
  for (i in 1:N) {
    summation_term <- numeric(3)
    for (t in 1:Time) {
      idx <- (i - 1) * Time + t
      summation_term <- summation_term + fit_model1$residuals[idx] * synthetic_data$avail[idx] * (synthetic_data$send[idx] - rho) * Z[t,]
    }
    W_sum <- W_sum + summation_term %*% t(summation_term)
  }
  W <- W_sum / N
  test_stats <- c(test_stats, N * t(beta_hat) %*% (Q %*% solve(W) %*% Q) %*% beta_hat)
  
}
