library(quantreg)
library(dplyr)
library(stats)
library(ggplot2)
library(tidyr)
library(stats)
library(Rlab)
library(fGarch)
library(VaRES)
library(gridExtra)
library(knitr)
library(patchwork)
library(reshape2)

#______________________________________________________________________
#Sammenlign QR og Z-test 
compute_risk_measures_rolling <- function(X_full, alpha = 0.025, window_size = 500) {
  M <- nrow(X_full)
  T_total <- ncol(X_full)
  T_eval <- 250  # only evaluate the last 250 days
  
  VaR_array <- matrix(NA, nrow = M, ncol = T_eval)
  ES5_array <- matrix(NA, nrow = M, ncol = T_eval)
  ES8_array <- matrix(NA, nrow = M, ncol = T_eval)
  
  # --- Helper functions ---
  hqr <- function(params, data, alpha) {
    beta0 <- params[1]
    beta1 <- params[2]
    beta2 <- params[3]
    fitted <- beta0 + beta1 * data$lag_1 + beta2 * data$lag_2
    residuals <- data$returns - fitted
    if (any(!is.finite(residuals))) return(Inf)
    sum(residuals * (alpha - ifelse(residuals < 0, 1, 0)))
  }
  
  ES_e5 <- function(X, tau, k0, VaR_e) {
    r_k <- (X - VaR_e) * (tau - ifelse(X < VaR_e, 1, 0))
    t <- length(X)
    mean(X[k0:t]) - sum(r_k[k0:t]) / (tau * (t - k0 + 1))
  }
  
  expectile_loss <- function(m, returns, theta) {
    u <- returns - m
    sum(abs(theta - (u < 0)) * u^2)
  }
  
  mu <- function(theta, X) {
    if (any(!is.finite(X))) return(NA)
    result <- optim(par = median(X), fn = expectile_loss, returns = X, theta = theta, method = "BFGS")
    result$par
  }
  
  find_theta_e <- function(X, VaR_e, tol = 0.0001, lower = 0.002, upper = 0.495) {
    if (any(!is.finite(X)) || !is.finite(VaR_e)) return(NA)
    while ((upper - lower) > tol) {
      theta_mid <- (lower + upper) / 2
      mu_mid <- mu(theta_mid, X)
      if (!is.finite(mu_mid)) return(NA)
      if (mu_mid > VaR_e) {
        upper <- theta_mid
      } else {
        lower <- theta_mid
      }
    }
    (lower + upper) / 2
  }
  
  ES_e8 <- function(X, tau, theta_e, VaR_e) {
    t <- length(X)
    if (!is.finite(theta_e) || !is.finite(VaR_e)) return(NA)
    -(1 + theta_e / ((1 - 2 * theta_e) * tau)) * VaR_e +
      theta_e / ((1 - 2 * theta_e) * tau) * mean(X[1:(t - 1)])
  }
  
  # --- Main loop: evaluate last 250 days only ---
  for (m in 1:M) {
    returns <- X_full[m, ]
    for (t_eval in 1:T_eval) {
      t_full <- (T_total - T_eval) + t_eval  
      
      window_start <- t_full - window_size + 1
      window_end <- t_full
      window_returns <- returns[window_start:window_end]
      
      if (length(window_returns) < 3 || any(!is.finite(window_returns))) next
      
      data <- data.frame(
        returns = window_returns[3:length(window_returns)],
        lag_1 = window_returns[2:(length(window_returns) - 1)],
        lag_2 = window_returns[1:(length(window_returns) - 2)]
      )
      if (any(!is.finite(unlist(data)))) next
      
      initial_guess <- c(0, 0, 0)
      opt_result <- tryCatch(
        optim(initial_guess, hqr, data = data, alpha = alpha, method = "BFGS"),
        error = function(e) return(NULL)
      )
      if (is.null(opt_result)) next
      
      beta0 <- opt_result$par[1]
      beta1 <- opt_result$par[2]
      beta2 <- opt_result$par[3]
      
      VaR_e <- beta0 + beta1 * window_returns[length(window_returns) - 1] +
        beta2 * window_returns[length(window_returns) - 2]
      k0 <- floor(length(window_returns) / 2)
      ES5 <- -ES_e5(window_returns, alpha, k0, VaR_e)
      theta_e <- find_theta_e(window_returns, VaR_e)
      ES8 <- ES_e8(window_returns, alpha, theta_e, VaR_e)
      
      VaR_array[m, t_eval] <- VaR_e
      ES5_array[m, t_eval] <- ES5
      ES8_array[m, t_eval] <- ES8
    }
  }
  
  return(list(
    VaR = VaR_array,
    ES5 = ES5_array,
    ES8 = ES8_array
  ))
}

# Parameters
T <- 250
M <- 250
burn_in <- 1000
T_total <- T + burn_in
alpha <- 0.025

# GARCH(1,1) parameters
w <- 0.01
a <- 0.1
b <- 0.85

# Degrees of freedom to test under H1
df_values <- c(3, 5, 10,20, 50, 100)

set.seed(4444)
# --- STORAGE ---
X_list <- list()   
Var_list <- list()
ES_list <- list()
var_sd_list <- list() 
result_list <- list()       # Store risk measure results

# --- GARCH SIMULATION FUNCTION ---
simulate_garch_paths <- function(M, T_total, df) {
  epsilons <- matrix(rstd(M * T_total, mean = 0, sd = 1, nu = df), nrow = M, ncol = T_total)
  X <- matrix(0, nrow = M, ncol = T_total)
  sigma2 <- matrix(0, nrow = M, ncol = T_total)
  sigma2[,1] <- 0.05
  
  for (t in 2:T_total) {
    sigma2[,t] <- w + a * X[,t-1]^2 + b * sigma2[,t-1]
    X[,t] <- sqrt(sigma2[,t]) * epsilons[,t]
  }
  
  list(
    X = X,              # full path, including burn-in
    sigma2 = sigma2     # full conditional variances
  )
}
# --- MAIN LOOP OVER df-values ---
for (df in df_values) {
  message(paste("Running simulation for df =", df))
  
  # Simulate GARCH paths
  sim <- simulate_garch_paths(M, T_total, df)
  X_post <- sim$X
  sigma2_post <- sim$sigma2
  
  var_sd_post <- sqrt(w + a * sim$X[, (burn_in):(T_total-1)]^2 + b * sim$sigma2[, (burn_in):(T_total-1)])
  Var_post <- -qstd(alpha,0,var_sd_post, nu = as.numeric(df))
  ES_post_mat <- matrix(0, nrow = M, ncol = T)
  for (i in 1:M) {
    for (j in 1:T) {
      sd_val <- var_sd_post[i, j]
      ES_post_mat[i, j] <- -1 / alpha * integrate(
        function(q) qstd(q, 0, sd_val, nu = df),
        lower = 0, upper = alpha
      )$value
    }
  }
  
  X_list[[as.character(df)]]  <- X_post
  Var_list[[as.character(df)]]  <- Var_post
  ES_list[[as.character(df)]]  <- ES_post_mat
  var_sd_list[[as.character(df)]]  <- var_sd_post   # Store var_sd

  # Compute rolling risk measures
  result <- compute_risk_measures_rolling(X_post, alpha = alpha, window_size = 500)
  result_list[[as.character(df)]] <- result
}
X_list <- lapply(X_list, function(X) X[, (ncol(X) - T + 1):ncol(X)])


Z2 <- function(T, ES_matrix, Var_matrix, X_matrix) {
  I_t <- ifelse(X_matrix + Var_matrix < 0, 1, 0)  # Indicator function (same shape as X_matrix)
  Z_values <- rowSums((X_matrix * I_t) / (T * 0.025 * ES_matrix)) + 1  # Compute Z2 for each row (simulation)
  return(Z_values)
}
Z4 <- function(T, ES_matrix, Var_matrix, X_matrix) {
  I_t <- ifelse(X_matrix + Var_matrix < 0, 1, 0) 
  Z_values <- rowMeans(1-(Var_matrix/ES_matrix))+rowSums((X_matrix+Var_matrix)*I_t/(0.025*T*ES_matrix))
  return(Z_values)
}
T2 <- function(T, ES_matrix, Var_matrix, X_matrix) {
  M <- nrow(ES_matrix)  # Number of simulations
  T_values <- matrix(0, nrow = M, ncol = 2) 
  
  for (j in 1:M) {
    # Extract the j-th simulation path (length T vectors)
    X_t <- X_matrix[j, ]
    Var_t <- Var_matrix[j, ]
    ES_t <- ES_matrix[j, ]
    
    I_t <- ifelse(X_t + Var_t < 0, 1, 0)
    v1 <- alpha - I_t
    v2 <- Var_t - ES_t + (1 / alpha) * (-Var_t - X_t) * I_t
    v <- cbind(-v1, v2)
    
    Omega_sum <- matrix(0, nrow = 2, ncol = 2)
    
    for (i in 1:T) {
      vi <- matrix(v[i, ], ncol = 1)
      Omega_sum <- Omega_sum + vi %*% t(vi)
    }
    
    Omega <- Omega_sum / T
    zsum <- colSums(v)
    
    T_values[j, ] <- T^(-0.5) * diag(Omega)^(-0.5) * zsum
  }
  
  return(T_values)
}



# === Function to compute power ===
compute_power <- function(stat_mis, stat_true) {
  crit_val <- quantile(stat_true, probs = 0.05, type = 1)
  mean(stat_mis <= crit_val)
}

# === Storage containers ===
power_results <- list()

X_list[[3]]
df_values <- c("3", "5", "10", "20", "50", "100")
# === Loop over each df value ===
for (df in df_values) {
  X_df <- X_list[[df]]
  ES_true <- ES_list[[df]]
  Var_true <- Var_list[[df]]
  ES5 <- result_list[[df]]$ES5
  ES8 <- result_list[[df]]$ES8
  VaR <- -result_list[[df]]$VaR
  
  # Compute Z2 stats
  Z2_true <- Z2(T, ES_true, Var_true, X_df)
  Z2_5 <- Z2(T, ES5, VaR, X_df)
  Z2_8 <- Z2(T, ES8, VaR, X_df)
  
  # Compute Z4 stats
  Z4_true <- Z4(T, ES_true, Var_true, X_df)
  Z4_5 <- Z4(T, ES5, VaR, X_df)
  Z4_8 <- Z4(T, ES8, VaR, X_df)
  
  # Compute T2 stats
  T2_5 <- T2(T, ES5, VaR, X_df)
  T2_8 <- T2(T, ES8, VaR, X_df)
  
  # Calculate power
  power_results[[df]] <- data.frame(
    df = df,
    Z2_ES5 = compute_power(Z2_5, Z2_true),
    Z2_ES8 = compute_power(Z2_8, Z2_true),
    Z4_ES5 = compute_power(Z4_5, Z4_true),
    Z4_ES8 = compute_power(Z4_8, Z4_true),
    T2_ES5 = mean(2*apply(1-pnorm(T2_5),1,min) <= 0.05),  
    T2_ES8 = mean(2*apply(1-pnorm(T2_8),1,min) <= 0.05)
  )
}


power_df <- do.call(rbind, power_results)
rownames(power_df) <- NULL
df_levels <- c("3", "5", "10", "20", "50", "100")
power_df$df <- factor(power_df$df, levels = df_levels)


print(power_df)


power_long <- power_df %>%
  pivot_longer(cols = -df, names_to = "Test_ES", values_to = "Power") %>%
  mutate(
    Test_ES = gsub("ES5", "ES_q", Test_ES),
    Test_ES = gsub("ES8", "ES_qe", Test_ES)
  )


ggplot(power_long, aes(x = factor(df, levels = df_values), y = Power, color = Test_ES, group = Test_ES)) +
  geom_line() + 
  geom_point() +
  scale_y_continuous(limits = c(0, 1)) + 
  labs(
    title = "",
    x = "df", y = "Power"
  ) +
  theme_minimal()



error_plots <- list()

# Loop over each df
for (df_example in df_values) {
  # Subset matrices to last T columns
  ES_true <- (ES_list[[df_example]])
  VaR_true <- (Var_list[[df_example]])
  ES5 <- (result_list[[df_example]]$ES5)
  ES8 <- (result_list[[df_example]]$ES8)
  VaR_mis <- -(result_list[[df_example]]$VaR)
  
  # Flatten matrices and compute errors
  error_ES5 <- as.vector(ES5 - ES_true)
  error_ES8 <- as.vector(ES8 - ES_true)
  error_VaR <- as.vector(VaR_mis - VaR_true)
  
  # Create combined data frame for ggplot
  df_error <- data.frame(
    Error = c(error_ES5, error_ES8, error_VaR),
    Type = rep(c("ES_q - Sand", "ES_qe - Sand", "VaR - Sand"), each = length(error_ES5))
  )
  
  # Plot for this df
  p <- ggplot(df_error, aes(x = Error, fill = Type)) +
    geom_density(alpha = 0.5) +
    xlim(-1, 1) +
    labs(
      title = paste(" df =", df_example),
      x = "(Estimat - Sand)", y = ""
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Store the plot
  error_plots[[df_example]] <- p
}
wrap_plots(error_plots, ncol = 2)


result_list[[]]

# Select df = 10
df_to_plot <- "10"
X_df <- X_list[[df_to_plot]]
ES5 <- -result_list[[df_to_plot]]$ES5
ES8 <- -result_list[[df_to_plot]]$ES8
VaR <- result_list[[df_to_plot]]$VaR

# Use the first simulated path (row 1)
X_path <- X_df[6, ]
ES5_path <- ES5[6, ]
ES8_path <- ES8[6, ]
VaR_path <- VaR[6, ]


plot_data <- data.frame(
  Time = 1:length(X_path),
  Afkast = X_path,
  VaR = VaR_path,
  ES_q = ES5_path,
  ES_qe = ES8_path
)

# Convert to long format
plot_data <- melt(plot_data, id.vars = "Time",
                       variable.name = "Series", value.name = "Value")

ggplot(plot_data, aes(x = Time, y = Value, color = Series, linetype = Series)) +
  geom_line(size = 1) +
  scale_color_manual(values = c(
    "Afkast" = "black",
    "VaR" = "red",
    "ES_q" = "blue",
    "ES_qe" = "green"
  )) +
  scale_linetype_manual(values = c(
    "Afkast" = "solid",
    "VaR" = "dashed",
    "ES_q" = "dotted",
    "ES_qe" = "dotdash"
  )) +
  labs(
    title = "",     
    y = "Værdi",
    x = "Tid",
    color = NULL,
    linetype = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "top")

