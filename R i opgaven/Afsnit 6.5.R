library(quantreg)
library(dplyr)
library(stats)
library(ggplot2)
library(tidyr)
library(stats)
library(Rlab)
library(fGarch)
library(VaRES)
library(dplyr)
library(gridExtra)
library(knitr)
library(patchwork)
#Estimating VaR using a linear QR
# Simulate a time series of returns
set.seed(123)
n <- 250
returns <- rnorm(n, mean = 0, sd = 1)  
data <- data.frame(
  returns = returns,
  lag_1 = dplyr::lag(returns, 1),
  lag_2 = dplyr::lag(returns, 2)
) %>% na.omit()
#View(data) 
# Define quantile level for VaR
alpha <- 0.025  # 5% quantile
 

hqr <- function(params) {
  beta0 <- params[1]
  beta1 <- params[2]
  beta2 <- params[3]
  
  residuals <- data$returns - (beta0 + beta1 * data$lag_1 + beta2 * data$lag_2)
  sum(residuals * (alpha - ifelse(residuals < 0, 1, 0)))
}

# Initial guess & optimization
initial_guess <- c(0, 0, 0)
result <- optim(initial_guess, hqr, method = "BFGS")

# Extract estimated parameters
beta0 <- result$par[1]
beta1 <- result$par[2]
beta2 <- result$par[3]

# Compute estimated Value at Risk (VaR)
VaR_e <- beta0 + beta1 * data$lag_1[nrow(data)] + beta2 * data$lag_2[nrow(data)]


#______________________________________________________________________

#ES estimator (5)
ES_e5 <- function(t,X,tau,k0,VaR_e){
  r_k <- (X-VaR_e)*(tau-ifelse(X<VaR_e,1,0))
  1/(t-k0)*sum(X[k0:(t-1)])-1/(tau*(t-k0))*sum(r_k[k0:(t-1)])
}

-ES_e5(n,returns,alpha,n/2,VaR_e)


#______________________________________________________________________
#ES estimator (8)
expectile_loss <- function(m, returns, theta) {
  u <- returns - m
  loss <- sum(abs(theta - (u < 0)) * u^2)  # Expectile loss function
  return(loss)
}


mu <- function(theta, X) {
  result <- optim(par = median(X),  # Initial guess (median)
                  fn = expectile_loss,
                  returns = X,
                  theta = theta,
                  method = "BFGS")  
  return(result$par)  
}


find_theta_e <- function(X, VaR_e, tol = 0.0001, lower = 0.002, upper = 0.495) {
  while ((upper - lower) > tol) {
    theta_mid <- (lower + upper) / 2  
    mu_mid <- mu(theta_mid, X)  
    

    if (mu_mid > VaR_e) {
      upper <- theta_mid
    } else {
      lower <- theta_mid
    }
  }
  return((lower + upper) / 2)  
} 

ES_e8 <- function(t,X,tau,theta_e,VaR_e){
  -(1+theta_e/((1-2*theta_e)*tau))*VaR_e+theta_e/((1-2*theta_e)*tau)*1/(t-1)*sum(X[1:(t-1)])
}

theta_e <- find_theta_e(returns,VaR_e)

ES_e8(250,returns,alpha,theta_e,VaR_e)





#______________________________________________________________________
#alt ovenfor i en funktion
estimate_VaR_ES <- function(gamma,df, alpha, n, ES_t, M) {
  VaR_values <- numeric(M)
  ES5_values <- numeric(M)
  ES8_values <- numeric(M)
  rel_errors_ES5 <- numeric(M)
  rel_errors_ES8 <- numeric(M)
  
  for (i in 1:M) {
    # Generate new returns for each simulation
    simulated_returns <- gamma*rt(n,df)
    
    # Ensure correct length of returns
    if (length(simulated_returns) < n) stop("Not enough return data.")
    
    # Create lagged returns
    data <- data.frame(
      returns = simulated_returns,
      lag_1 = dplyr::lag(simulated_returns, 1),
      lag_2 = dplyr::lag(simulated_returns, 2)
    ) %>% na.omit()
    
    # Quantile regression loss function
    hqr <- function(params) {
      residuals <- data$returns - (params[1] + params[2] * data$lag_1 + params[3] * data$lag_2)
      sum(residuals * (alpha - ifelse(residuals < 0, 1, 0)))
    }
    
    # Optimization for beta parameters
    result <- optim(c(0, 0, 0), hqr, method = "BFGS")
    beta0 <- result$par[1]
    beta1 <- result$par[2]
    beta2 <- result$par[3]
    
    # Estimate Value at Risk (VaR)
    VaR_e <- beta0 + beta1 * data$lag_1[nrow(data)] + beta2 * data$lag_2[nrow(data)]
    VaR_values[i] <- VaR_e
    
    # Expected Shortfall Estimation - Method (5)
    ES_e5 <- function(t, X, tau, k0, VaR_e) {
      r_k <- (X - VaR_e) * (tau - ifelse(X < VaR_e, 1, 0))
      mean(X[k0:(t - 1)]) - sum(r_k[k0:(t - 1)]) / (tau * (t - k0))
    }
    
    ES_5 <- -ES_e5(n, simulated_returns, alpha, n/2, VaR_e)
    ES5_values[i] <- ES_5
    rel_errors_ES5[i] <- abs(ES_5 - ES_t) / ES_t
    
    # Expectile loss function for Method (8)
    expectile_loss <- function(m, returns, theta) {
      u <- returns - m
      sum(abs(theta - (u < 0)) * u^2)
    }
    
    mu <- function(theta, X) {
      optim(par = median(X), fn = expectile_loss, returns = X, theta = theta, method = "BFGS")$par
    }
    
    find_theta_e <- function(X, VaR_e, tol = 1e-4, lower = 0.002, upper = 0.495) {
      while ((upper - lower) > tol) {
        theta_mid <- (lower + upper) / 2
        if (mu(theta_mid, X) > VaR_e) {
          upper <- theta_mid
        } else {
          lower <- theta_mid
        }
      }
      (lower + upper) / 2
    }
    
    theta_e <- find_theta_e(simulated_returns, VaR_e)
    
    ES_e8 <- function(t, X, tau, theta_e, VaR_e) {
      -(1 + theta_e / ((1 - 2 * theta_e) * tau)) * VaR_e +
        theta_e / ((1 - 2 * theta_e) * tau) * mean(X[1:(t - 1)])
    }
    
    ES_8 <- ES_e8(250, simulated_returns, alpha, theta_e, VaR_e)
    ES8_values[i] <- ES_8
    rel_errors_ES8[i] <- abs(ES_8 - ES_t) / ES_t
  }
  
  mean_VaR <- mean(VaR_values)
  mean_ES5 <- mean(ES5_values)
  mean_ES8 <- mean(ES8_values)
  mean_rel_error_ES5 <- mean(rel_errors_ES5)
  mean_rel_error_ES8 <- mean(rel_errors_ES8)
  
  cat("Mean VaR across", M, "simulations:", mean_VaR, "\n")
  cat("Mean ES (Method 5) across", M, "simulations:", mean_ES5, "\n")
  cat("Mean ES (Method 8) across", M, "simulations:", mean_ES8, "\n")
  cat("Mean Relative Error (ES5) across", M, "simulations:", mean_rel_error_ES5, "\n")
  cat("Mean Relative Error (ES8) across", M, "simulations:", mean_rel_error_ES8, "\n")
  
  return(list(
    Mean_VaR = mean_VaR,
    Mean_ES_method_5 = mean_ES5,
    Mean_ES_method_8 = mean_ES8,
    Mean_Relative_Error_ES5 = mean_rel_error_ES5,
    Mean_Relative_Error_ES8 = mean_rel_error_ES8
  ))
}

#______________________________________________________________________
#Plot for ændringer i gamma 
set.seed(1234)
df <- 100
alpha <- 0.025  
n <- 250   
m <- 4000    


gamma_vals <- seq(0.8, 1.5, by = 0.05)

compute_ES <- function(gamma) {
  -integrate(function(q) gamma * qt(q, df), 0, alpha)$value / alpha
}


rel_err_es5 <- numeric(length(gamma_vals))
rel_err_es8 <- numeric(length(gamma_vals))

# Main loop
for (i in seq_along(gamma_vals)) {
  gamma <- gamma_vals[i]
  ES_t <- compute_ES(gamma)  
  
  result <- estimate_VaR_ES(gamma, df, alpha, n, ES_t, m)
  
  rel_err_es5[i] <- result$Mean_Relative_Error_ES5
  rel_err_es8[i] <- result$Mean_Relative_Error_ES8
}

# Plot
par(mfrow = c(1,1))
plot(gamma_vals, rel_err_es5, type = "l", col = "blue", lwd = 2,
     ylim = range(c(0, 0.5)),
     xlab = expression(gamma), ylab = "Relativ fejl",
     main = "")
lines(gamma_vals, rel_err_es8, col = "red", lwd = 2)
legend("topright", legend = c("ES_q", "ES_qe"),
       col = c("blue", "red"), lwd = 2)



#______________________________________________________________________
#Plot for ændringer i df
set.seed(123)

estimate_VaR_ES <- function(gamma, df, alpha, n, ES_t, M) {
  VaR_values <- numeric(M)
  ES5_values <- numeric(M)
  ES8_values <- numeric(M)
  rel_errors_ES5 <- numeric(M)
  rel_errors_ES8 <- numeric(M)
  rel_err_ES5_cond <- numeric()
  rel_err_ES8_cond <- numeric()
  condition_met_ES5 <- 0
  condition_met_ES8 <- 0
  
  for (i in 1:M) {
    simulated_returns <- gamma * rt(n, df)
    
    if (length(simulated_returns) < n) stop("Not enough return data.")
    
    data <- data.frame(
      returns = simulated_returns,
      lag_1 = dplyr::lag(simulated_returns, 1),
      lag_2 = dplyr::lag(simulated_returns, 2)
    ) %>% na.omit()
    
    hqr <- function(params) {
      residuals <- data$returns - (params[1] + params[2] * data$lag_1 + params[3] * data$lag_2)
      sum(residuals * (alpha - ifelse(residuals < 0, 1, 0)))
    }
    
    result <- optim(c(0, 0, 0), hqr, method = "BFGS")
    beta0 <- result$par[1]
    beta1 <- result$par[2]
    beta2 <- result$par[3]
    
    VaR_e <- beta0 + beta1 * data$lag_1[nrow(data)] + beta2 * data$lag_2[nrow(data)]
    VaR_values[i] <- VaR_e
    
    ES_e5 <- function(t, X, tau, k0, VaR_e) {
      r_k <- (X - VaR_e) * (tau - ifelse(X < VaR_e, 1, 0))
      mean(X[k0:(t - 1)]) - sum(r_k[k0:(t - 1)]) / (tau * (t - k0))
    }
    
    ES_5 <- -ES_e5(n, simulated_returns, alpha, n / 2, VaR_e)
    ES5_values[i] <- ES_5
    rel_errors_ES5[i] <- abs(ES_5 - ES_t) / ES_t
    
    expectile_loss <- function(m, returns, theta) {
      u <- returns - m
      sum(abs(theta - (u < 0)) * u^2)
    }
    
    mu <- function(theta, X) {
      optim(par = median(X), fn = expectile_loss, returns = X, theta = theta, method = "BFGS")$par
    }
    
    find_theta_e <- function(X, VaR_e, tol = 1e-4, lower = 0.002, upper = 0.495) {
      while ((upper - lower) > tol) {
        theta_mid <- (lower + upper) / 2
        if (mu(theta_mid, X) > VaR_e) {
          upper <- theta_mid
        } else {
          lower <- theta_mid
        }
      }
      (lower + upper) / 2
    }
    
    theta_e <- find_theta_e(simulated_returns, VaR_e)
    
    ES_e8 <- function(t, X, tau, theta_e, VaR_e) {
      -(1 + theta_e / ((1 - 2 * theta_e) * tau)) * VaR_e +
        theta_e / ((1 - 2 * theta_e) * tau) * mean(X[1:(t - 1)])
    }
    
    ES_8 <- ES_e8(250, simulated_returns, alpha, theta_e, VaR_e)
    ES8_values[i] <- ES_8
    rel_errors_ES8[i] <- abs(ES_8 - ES_t) / ES_t
    
    if (ES_5 < ES_t) {
      rel_err_ES5_cond <- c(rel_err_ES5_cond, rel_errors_ES5[i])
      condition_met_ES5 <- condition_met_ES5 + 1
    }
    if (ES_8 < ES_t) {
      rel_err_ES8_cond <- c(rel_err_ES8_cond, rel_errors_ES8[i])
      condition_met_ES8 <- condition_met_ES8 + 1
    }
  }
  
  condition_percent_ES5 <- 100 * condition_met_ES5 / M
  condition_percent_ES8 <- 100 * condition_met_ES8 / M
  
  return(list(
    Mean_VaR = mean(VaR_values),
    Mean_ES_method_5 = mean(ES5_values),
    Mean_ES_method_8 = mean(ES8_values),
    Mean_Relative_Error_ES5 = mean(rel_errors_ES5),
    Mean_Relative_Error_ES8 = mean(rel_errors_ES8),
    Rel_Error_ES5_Condition = rel_err_ES5_cond,
    Rel_Error_ES8_Condition = rel_err_ES8_cond,
    Condition_Percent_ES5 = condition_percent_ES5,
    Condition_Percent_ES8 = condition_percent_ES8
  ))
}

# Simulation parameters
df_vals <- c(3, 5, 10, 15, 20, 30, 50, 100)
alpha <- 0.025
n <- 250
m <- 4000
gamma <- 1

# Function to compute ES_t
compute_ES <- function(df) {
  -integrate(function(q) gamma * qt(q, df), 0, alpha)$value / alpha
}

# Storage
rel_err_es5 <- numeric(length(df_vals))
rel_err_es8 <- numeric(length(df_vals))
rel_err_ES5_all <- list()
rel_err_ES8_all <- list()
condition_percents_ES5 <- numeric(length(df_vals))
condition_percents_ES8 <- numeric(length(df_vals))

# Loop over df
for (i in seq_along(df_vals)) {
  df <- df_vals[i]
  ES_t <- compute_ES(df)
  
  result <- estimate_VaR_ES(gamma, df, alpha, n, ES_t, m)
  
  rel_err_es5[i] <- result$Mean_Relative_Error_ES5
  rel_err_es8[i] <- result$Mean_Relative_Error_ES8
  
  rel_err_ES5_all[[i]] <- result$Rel_Error_ES5_Condition
  rel_err_ES8_all[[i]] <- result$Rel_Error_ES8_Condition
  
  condition_percents_ES5[i] <- result$Condition_Percent_ES5
  condition_percents_ES8[i] <- result$Condition_Percent_ES8
}

# --- Plot 1: Line plot of average relative errors
x_pos <- seq_along(df_vals)

plot(x_pos, rel_err_es5, type = "l", col = "blue", lwd = 2, pch = 19,
     ylim = range(c(0, 0.5)),
     xlab = "df", ylab = "Relative Error",
     xaxt = "n", main = "")
axis(1, at = x_pos, labels = df_vals)
lines(x_pos, rel_err_es8, type = "l", col = "red", lwd = 2, pch = 17)
legend("topright", legend = c("ES_q", "ES_qe"), col = c("blue", "red"), lwd = 2)

# --- Plot 2: Boxplots of conditional relative errors
df_labels_ES5 <- rep(df_vals, times = sapply(rel_err_ES5_all, length))
df_labels_ES8 <- rep(df_vals, times = sapply(rel_err_ES8_all, length))
es5_errors <- unlist(rel_err_ES5_all)
es8_errors <- unlist(rel_err_ES8_all)

plot_data <- data.frame(
  df = as.factor(c(df_labels_ES5, df_labels_ES8)),
  rel_error = c(es5_errors, es8_errors),
  Method = c(rep("ES5", length(es5_errors)), rep("ES8", length(es8_errors)))
)

# Remove rows where rel_error > 5
plot_data <- plot_data %>%
  filter(rel_error <= 5)


# --- Boxplot (with ES_q and ES_qe)
plot_data$Method <- factor(plot_data$Method, levels = c("ES5", "ES8"), labels = c("ES_q", "ES_qe"))

p1 <- ggplot(plot_data, aes(x = df, y = rel_error, fill = Method)) +
  geom_boxplot(position = position_dodge(0.8), outlier.alpha = 0.3) +
  scale_fill_manual(values = c("blue", "red")) +
  labs(
    title = "",
    x = "",
    y = "Relativ fejl",
    fill = "Metode"
  ) +
  theme_minimal()+
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    #legend.position = "none"  # Hide legend here
  )

# --- Barplot (also with ES_q and ES_qe)
bar_df <- data.frame(
  df = rep(df_vals, 2),
  method = rep(c("ES_q", "ES_qe"), each = length(df_vals)),
  percent = c(condition_percents_ES5, condition_percents_ES8)
)

p2 <- ggplot(bar_df, aes(x = factor(df), y = percent, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("blue", "red")) +
  labs(
    title = "",
    x = "df",
    y = "Procent",
    fill = "Metode"
  ) +
  theme_minimal()

# --- Combine vertically
p1 / p2 + plot_layout(heights = c(2, 1))

#______________________________________________________________________
#For forskellige GARCH fordelinger 
estimate_VaR_ES_matrix_error <- function(X_matrix, alpha, n, ES_H) {
  M <- nrow(X_matrix)  # Number of simulations
  
  # Storage for results
  VaR_estimates <- numeric(M)
  ES_5_estimates <- numeric(M)
  ES_8_estimates <- numeric(M)
  relative_errors_5 <- numeric(M)
  relative_errors_8 <- numeric(M)
  
  for (m in 1:M) {
    # Extract one simulation's returns
    returns <- X_matrix[m, ]
    
    # Ensure correct length of returns
    if (length(returns) < n) next
    
    # Create lagged returns
    data <- data.frame(
      returns = returns,
      lag_1 = dplyr::lag(returns, 1),
      lag_2 = dplyr::lag(returns, 2)
    ) %>% na.omit()
    
    # Quantile regression loss function
    hqr <- function(params) {
      residuals <- data$returns - (params[1] + params[2] * data$lag_1 + params[3] * data$lag_2)
      sum(residuals * (alpha - ifelse(residuals < 0, 1, 0)))
    }
    
    # Optimization for beta parameters
    result <- optim(c(0, 0, 0), hqr, method = "BFGS")
    beta0 <- result$par[1]
    beta1 <- result$par[2]
    beta2 <- result$par[3]
    
    # Estimate Value at Risk (VaR)
    VaR_e <- beta0 + beta1 * data$lag_1[nrow(data)] + beta2 * data$lag_2[nrow(data)]
    VaR_estimates[m] <- VaR_e
    
    # Expected Shortfall Estimation - Method (5)
    ES_e5 <- function(t, X, tau, k0, VaR_e) {
      r_k <- (X - VaR_e) * (tau - ifelse(X < VaR_e, 1, 0))
      mean(X[k0:(t - 1)]) - sum(r_k[k0:(t - 1)]) / (tau * (t - k0))
    }
    
    ES_5 <- -ES_e5(n, returns, alpha, n / 2, VaR_e)
    ES_5_estimates[m] <- ES_5
    
    # Expectile loss function for Method (8)
    expectile_loss <- function(m, returns, theta) {
      u <- returns - m
      sum(abs(theta - (u < 0)) * u^2)
    }
    
    mu <- function(theta, X) {
      optim(par = median(X), fn = expectile_loss, returns = X, theta = theta, method = "BFGS")$par
    }
    
    find_theta_e <- function(X, VaR_e, tol = 1e-4, lower = 0.002, upper = 0.495) {
      while ((upper - lower) > tol) {
        theta_mid <- (lower + upper) / 2
        if (mu(theta_mid, X) > VaR_e) {
          upper <- theta_mid
        } else {
          lower <- theta_mid
        }
      }
      (lower + upper) / 2
    }
    
    theta_e <- find_theta_e(returns, VaR_e)
    
    ES_e8 <- function(t, X, tau, theta_e, VaR_e) {
      -(1 + theta_e / ((1 - 2 * theta_e) * tau)) * VaR_e +
        theta_e / ((1 - 2 * theta_e) * tau) * mean(X[1:(t - 1)])
    }
    
    ES_8 <- ES_e8(250, returns, alpha, theta_e, VaR_e)
    ES_8_estimates[m] <- ES_8
    
    # Compute relative errors using provided ES_H instead of integral
    ES_t <- ES_H[m,250]  # Use precomputed ES_H values for corresponding simulation
    relative_errors_5[m] <- abs(ES_5 - ES_t) / ES_t
    relative_errors_8[m] <- abs(ES_8 - ES_t) / ES_t
  }
  
  # Compute means across simulations
  mean_VaR <- mean(VaR_estimates, na.rm = TRUE)
  mean_ES5 <- mean(ES_5_estimates, na.rm = TRUE)
  mean_ES8 <- mean(ES_8_estimates, na.rm = TRUE)
  mean_relative_error_5 <- mean(relative_errors_5, na.rm = TRUE)
  mean_relative_error_8 <- mean(relative_errors_8, na.rm = TRUE)
  
  return(list(
    Mean_VaR = mean_VaR,
    Mean_ES5 = mean_ES5,
    Mean_ES8 = mean_ES8,
    Mean_Relative_Error_ES5 = mean_relative_error_5,
    Mean_Relative_Error_ES8 = mean_relative_error_8
  ))
}

# --- PARAMETERS ---
T <- 250
M <- 250
burn_in <- 250
T_total <- T + burn_in
alpha <- 0.025
w <- 0.01
a <- 0.1
b <- 0.85
df_values <- c(3, 5, 10, 15, 20, 30, 50, 100)

# --- STORAGE ---
X_list <- list()
Var_list <- list()
ES_list <- list()
result_list <- list()

# --- SIMULATION FUNCTION ---
simulate_garch_paths <- function(M, T_total, df) {
  epsilons <- matrix(rstd(M * T_total, mean = 0, sd = 1, nu = df), nrow = M, ncol = T_total)
  X <- matrix(0, nrow = M, ncol = T_total)
  sigma2 <- matrix(0, nrow = M, ncol = T_total)
  sigma2[,1] <- 0.05
  
  for (t in 2:T_total) {
    sigma2[,t] <- w + a * X[,t-1]^2 + b * sigma2[,t-1]
    X[,t] <- sqrt(sigma2[,t]) * epsilons[,t]
  }
  
  list(X = X[, (burn_in + 1):T_total],
       sigma2 = sigma2[, (burn_in + 1):T_total])
}

# --- MAIN LOOP OVER DF ---
for (df in df_values) {
  # Simulate paths
  sim <- simulate_garch_paths(M, T_total, df)
  X_post <- sim$X
  sigma2_post <- sim$sigma2
  
  # Compute VaR and ES for each simulation
  VaR_post <- -qstd(alpha, mean = 0, sd = sqrt(sigma2_post), nu = df)
  
  ES_post_mat <- matrix(NA, nrow = M, ncol = T)
  for (i in 1:M) {
    for (j in 1:T) {
      ES_post_mat[i, j] <- -1 / alpha * integrate(
        function(q) qstd(q, mean = 0, sd = sqrt(sigma2_post[i, j]), nu = df),
        lower = 0, upper = alpha
      )$value
    }
  }
  
  # Store each simulation result
  X_list[[as.character(df)]] <- X_post
  Var_list[[as.character(df)]] <- VaR_post
  ES_list[[as.character(df)]] <- ES_post_mat
  
  # Run your estimate function
  result <- estimate_VaR_ES_matrix_error(X_post, alpha, n = T, ES_H = ES_post_mat)
  result_list[[as.character(df)]] <- result
}

# Extract relative errors from result_list
rel_err_es5 <- sapply(result_list, function(res) res$Mean_Relative_Error_ES5)
rel_err_es8 <- sapply(result_list, function(res) res$Mean_Relative_Error_ES8)


# Create evenly spaced x-axis positions
x_pos <- seq_along(df_values)

# Plot
plot(x_pos, rel_err_es5, type = "l", col = "blue", lwd = 2, pch = 19,
     xaxt = "n",  # turn off default x-axis
     ylim = range(c(0,0.5), na.rm = TRUE),
     xlab = "df", ylab = "Relativ fejl",
     main = "GARCH(1,1)")
axis(1, at = x_pos, labels = df_values)

# Add ES8 line
lines(x_pos, rel_err_es8, col = "red", lwd = 2, pch = 17)

# Add legend
legend("topright", legend = c("ES_q (ES5)", "ES_qe (ES8)"),
       col = c("blue", "red"), lwd = 2, pch = c(19, 17))


#______________________________________________________________________
  
