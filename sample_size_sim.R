# Load necessary libraries
library(dplyr)
library(ggplot2)
library(MASS)
library(pwr)

# Set parameters
T <- 210  # Total number of decision points (42 days * 5 times per day)
rho <- 0.4  # Randomization probability
availability <- 0.5  # Constant availability

# Standardized Proximal Main Effect (assuming a quadratic form)
Z <- matrix(0, nrow = T, ncol = 3)
for (t in 1:T) {
  Z[t,] <- c(1, floor((t-1)/5), floor((t-1)/5)^2)
}


create_z <- function(t) {
  Z <- matrix(0, nrow = 3, ncol = 1)
  Z <- c(1, floor((t-1)/5), floor((t-1)/5)^2)
  return(Z)
}

#calculate matrix Q
Q <- matrix(0, nrow = 3, ncol = 3)
for (i in 1:T) {
  Z_t <- create_z(i)
  Q <- Q + (Z_t) %*% t(Z_t) * availability * rho * (1 - rho)
}

# Standardized treatment effect (d)
d <- c(0, 0.00964, -0.000172)

# Verify the d_bar
d_bar <- 0
for (i in 1:T) {
  d_bar <- d_bar + t(create_z(i)) %*% d / T
}
print(d_bar)

# Function to calculate sample size for given power and alpha
calculate_sample_size <- function(power, alpha, d_bar, Q) {
  p <- 3
  q <- 3
  # Calculate power function for given N
  power_function <- function(N) {
    
    #calculate the noncentraility parameter
    c_N <- N * (t(d) %*% Q %*% d)
    print(N)
    df1 <- p
    df2 <- N - q - p
    
    # Calculate the critical value for the F-distribution
    F_critical <- qf((N - q - p) * (1 - alpha) / (p * (N - q - 1)), df1 = df1, df2 = df2, ncp = 0)
    
    # Evaluate the left-hand side of the formula
    projected_power <- (p * (N - q - 1) / (N - q - p)) * pf(F_critical, df1 = df1, df2 = df2, ncp = c_N)
    return(projected_power)
  }
  
  beta.0 <- power
  max.iters = 10000
  N = 100
  i = 1
  
  while (i < max.iters ) {
    if ( power_function(N) < 1- beta.0) {
      if (power_function(N-1) > 1 - beta.0) {
        break
      } else {N = N - 1}
    } else if (power_function(N) > 1 - beta.0) {
      N = N+1
    } else if (power_function(N) == 1 - beta.0) {break}
    i = i+1
  }
  
  return(N)
  
  # # Debugging: Print values at interval endpoints
  # print(paste("Value at N=10:", power_function(10)))
  # print(paste("Value at N=1000:", power_function(50)))
  # 
  # # Calculate sample size for given power
  # result <- uniroot(power_function, interval = c(10, 50))
  # 
  # return(ceiling(result$root))
}

# Desired power and significance level
desired_power <- 0.2
alpha <- 0.05

# Calculate sample size
sample_size <- calculate_sample_size(desired_power, alpha, d, Q)
 print(paste("Required sample size: ", sample_size))
