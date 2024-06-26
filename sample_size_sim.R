# Continuous Outcome MRT sample size calculator
# From Liao et al 2016

rm(list = ls())

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(MASS)
library(pwr)

# Set parameters
Time <- 210  # Total number of decision points (42 days * 5 times per day)
rho <- 0.4  # Randomization probability
availability <- 1  # Constant availability

create_z <- function(t) {
  Z <- matrix(0, nrow = 3, ncol = 1)
  Z <- c(1, floor((t-1)/5), floor((t-1)/5)^2)
  return(Z)
}

#calculate matrix Z and Q
Z <- matrix(0, nrow = Time, ncol = 3)
for (t in 1:Time) {
  Z[t,] <- c(1, floor((t-1)/5), floor((t-1)/5)^2)
}
Q <- t(Z) %*% diag(rep(availability,Time) * rep(rho, Time) * rep((1-rho), Time)) %*% Z



# Standardized treatment effect (d)
d <- c(0, 0.00964, -0.000172)
d <- c(0, 0.00476, 0)
# d <- c(0, 0.01, -0.0002)
# beta_t <- matrix(rep(0.1, 210), ncol = 1)
# d <- solve(t(Z) %*% diag(rep(availability,Time)) %*% Z) %*% t(Z)  %*% diag(rep(availability,Time)) %*% beta_t

# Verify the max effect day
effect_by_day <- c()
for (i in 1:(Time/5)) {
  effect_by_day <- c(effect_by_day, t(create_z(i*5)) %*% d)
}
print(paste0("Maximum effect: ", which.max(effect_by_day)))

# Verify the d_bar
d_t <- numeric(0)
for (i in 1:Time) {
  d_t <- c(d_t, t(create_z(i)) %*% d)
}
print(paste0("d bar: ", sum(d_t)/ Time))

plot(1:length(d_t), d_t)

# Function to calculate sample size for given power and alpha
calculate_sample_size <- function(beta, alpha, d, Q) {
  p <- 3
  q <- 3
  # Calculate power function for given N
  power_function <- function(N) {
    
    #calculate the noncentraility parameter
    c_N <- N * t(d) %*% Q %*% d
    df1 <- p
    df2 <- N - q - p
    
    # Calculate the critical value for the F-distribution
    inv.f <- qf(1-alpha,df1,df2)
    #inv.f <- qf((N-p-q) *(1-alpha) / p / (N-q-1),df1,df2)
    # Evaluate the left-hand side of the formula
    calc_power = 1-pf(inv.f,df1,df2,ncp = c_N)
    #calc_power = 1 - p * (N-q-1) * pf(inv.f,df1,df2,ncp = c_N) / (N-p-q)
    #calc_power = 1 - pf( p * (N-q-1) * inv.f / (N-p-q),df1,df2,ncp = c_N)
    #print(paste0("Sample size: ", N, ", power: ", calc_power))
    return(calc_power - 0.8)
  }
  
  # Debugging: Print values at interval endpoints
  print(paste("Value at N=10:", power_function(10)))
  print(paste("Value at N=100:", power_function(1000)))
  
  # Calculate sample size for given power
  result <- uniroot(power_function, interval = c(10, 1000))
  

  return(ceiling(result$root))
}

# Desired power and significance level
beta <- 0.8
alpha <- 0.05

# Calculate sample size
sample_size <- calculate_sample_size(beta, alpha, d, Q)
print(paste("Required sample size: ", sample_size))
 

library(ggplot2)
sample_size <- 19
# Define parameters
df1 <- 3   # degrees of freedom 1
df2 <- sample_size - 3 -3   # degrees of freedom 2
c_N <- sample_size * t(d) %*% Q %*% d    # non-centrality parameter

# Define the range for x
x <- seq(0, 80, length.out = 1000)

# Calculate the pdf values for the F-distribution
pdf_values <- df(x, df1, df2, c_N) * 3 * (sample_size - 3 -1)  / (sample_size -3-3)

# Create a data frame for plotting
data <- data.frame(x = x, pdf = pdf_values)

# Plot the pdf
ggplot(data, aes(x = x, y = pdf)) +
  geom_line(color = "blue") +
  labs(title = "PDF of the F-Distribution",
       x = "x",
       y = "Density") +
  theme_minimal()
