# Analysis of synthetic data
# From Qian et al 2021 PsychMehtods

rm(list = ls())

library(tidyverse)
library(xtable)
library(geepack)
library(ggplot2)
source("xgeepack.R")
source("estimate.R")
# # Read the data into a data frame
# s_data <- data.frame(read_csv("synthetic_data_37subject_210time.csv"))
# 
# # Calculate the average availability for each time point
# average_avail <- aggregate(avail ~ decision.index.nogap, data = s_data, mean)
# 
# # Plot the average availability for each time point
# ggplot(average_avail, aes(x = decision.index.nogap, y = avail)) +
#   geom_line() +
#   geom_point() +
#   labs(title = "Average Availability Over Time",
#        x = "Time Point (decision.index.nogap)",
#        y = "Average Availability") +
#   theme_minimal()
# 
# var_sample <- var(average_avail$avail)
# var_bernoulli <- mean(average_avail$avail) * (1 - mean(average_avail$avail)) / 37
generate_synthetic_header <- function(N, T = 210, p = 0.4) {
  data <- data.frame(
    userid = rep(1:N, each = T),
    decision.index.nogap = rep(1:T, times = N),
    study.day.nogap = rep(rep(0:((T/5)-1), rep(5, (T/5))), times = N),
    study.day.square = rep(rep((0:((T/5)-1))**2, rep(5, (T/5))), times = N),
    avail = rbinom(N * T, 1, 0.5)
  )
  
  data$send <- with(data, ifelse(avail == 1, rbinom(N * T, 1, p), 0))
  # Model parameters
  alpha <- c(2.5, 0.01, 0.005) # Example values for B_t coefficients
  beta <- c(0, 0.00964, -0.000172) # Example values for Z_t coefficients
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

# Example usage
p_values_2 <- numeric(0)
p_values_1 <- numeric(0)
p_values <- numeric(0)
for (i in 1:200) {
  print(i)
  set.seed(i+1001)
  synthetic_data <- generate_synthetic_header(N = 43, T = 210)
  head(synthetic_data)
  
  synthetic_data$"(Intercept)" <- 1
  synthetic_data$"I(send - 0.4)" <- synthetic_data$send - 0.4
  
  xmat <- synthetic_data %>%
    transmute("(Intercept)" = .$"(Intercept)",
              "days" = .$"study.day.nogap",
              "days_sq" = .$"study.day.square",
              "I(send - 0.4)" = .$"I(send - 0.4)")
  
  fit_model2 <- geese.glm(x = as.matrix(xmat),
                          y = synthetic_data$steps, w = synthetic_data$avail, id = as.factor(synthetic_data$userid),
                          family = gaussian(), corstr = "independence")
  p_values_2 <- c(p_values_2, estimate(fit_model2)[4,8])
}


