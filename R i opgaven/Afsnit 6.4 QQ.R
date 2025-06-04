library(quantreg)
library(dplyr)
library(stats)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(ggpubr)
library(fGarch)  
generate_T2 <- function(T, df, gamma, alpha, ES, M) {
  T_values <- matrix(0, nrow = M, ncol = 2)
  for (j in 1:M) {
    X_t <- gamma * rt(T, df)
    VaR <- -qt(alpha, df)  # use consistent naming
    
    I_t <- ifelse(X_t + VaR < 0, 1, 0)
    v1 <- alpha - I_t
    v2 <- VaR- ES + (1 / alpha) * (-VaR - X_t) * I_t
    v <- cbind(-v1, v2)
    
    Omega_sum <- matrix(0, nrow = 2, ncol = 2)
    
    for (i in 1:T) {
      vi <- matrix(v[i, ], ncol = 1)  # Make it a column vector
      Omega_sum <- Omega_sum + vi %*% t(vi)
    }
    
    Omega <- Omega_sum/T
    zsum <- colSums(v)
    
    T_values[j, ] <- -(T^(-0.5) * diag(Omega)^(-0.5) * zsum)   #Obs her indsat negativt fortegn
  }
  return(T_values)
}


set.seed(1234)
T_list <- list(
  "T = 250" = generate_T2(250, 100, 1, 0.025, 2.378, 1000),
  "T = 500" = generate_T2(500, 100, 1, 0.025, 2.378, 1000),
  "T = 1000" = generate_T2(1000, 100, 1, 0.025, 2.378, 1000),
  "T = 5000" = generate_T2(5000, 100, 1, 0.025, 2.378, 1000)
)

par(mfrow = c(2, 2))  # 2 rækker, 2 kolonner

for (navn in c("T = 250", "T = 5000")) {
  T_vals <- T_list[[navn]]  # Hent de simulerede værdier
  
  # Komponent 1
  qqnorm(T_vals[, 1],
         main = paste(navn, "- Komponent 1"),
         xlab = "N(0,1) fraktiler",
         ylab = "T2 fraktiler")
  qqline(T_vals[, 1], col = "red")
  
  # Komponent 2
  qqnorm(T_vals[, 2],
         main = paste(navn, "- Komponent 2"),
         xlab = "N(0,1) fraktiler",
         ylab = "T2 fraktiler")
  qqline(T_vals[, 2], col = "blue")
}

# Optional: Print 5% quantile comparisons for T = 250 and T = 5000
cat("Empirisk 5%-kvantil for T = 250:", quantile(T_list[["T = 250"]][,1], 0.05, type = 1), "\n")
cat("Teoretisk 5%-kvantil (standard normal):", qnorm(0.04), "\n\n")

# Optional: Print 5% quantile comparisons for T = 250 and T = 5000
cat("Empirisk 5%-kvantil for T = 250:", quantile(T_list[["T = 250"]][,2], 0.05, type = 1), "\n")
cat("Teoretisk 5%-kvantil (standard normal):", qnorm(0.04), "\n\n")


cat("Empirisk 5%-kvantil for T = 5000:", quantile(T_list[["T = 5000"]][,1], 0.05, type = 1), "\n")
cat("Teoretisk 5%-kvantil (standard normal):", qnorm(0.04), "\n")

cat("Empirisk 5%-kvantil for T = 5000:", quantile(T_list[["T = 5000"]][,2], 0.05, type = 1), "\n")
cat("Teoretisk 5%-kvantil (standard normal):", qnorm(0.04), "\n")



# --- Parameters ---
M <- 1000
burn_in <- 250
alpha <- 0.025
df_h0 <- 10

# GARCH(1,1) parameters
w <- 0.01
a <- 0.1
b <- 0.85
alpha <- 0.025
df_h0 <- 10

# Function to simulate GARCH paths
simulate_garch_paths <- function(M, T_total, epsilons) {
  X <- matrix(0, nrow = M, ncol = T_total)
  sigma2 <- matrix(0, nrow = M, ncol = T_total)
  sigma2[, 1] <- 0.05
  
  for (t in 2:T_total) {
    sigma2[, t] <- w + a * X[, t - 1]^2 + b * sigma2[, t - 1]
    X[, t] <- sqrt(sigma2[, t]) * epsilons[, t]
  }
  
  list(X = X, sigma2 = sigma2)
}

# T2 test statistic function
T2 <- function(T, ES_matrix, Var_matrix, X_matrix) {
  M <- nrow(ES_matrix)  # Number of simulations
  T_values <- matrix(0, nrow = M, ncol = 2)  # Store 2 values for each simulation
  
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
    
    T_values[j, ] <- -(T^(-0.5) * diag(Omega)^(-0.5) * zsum)
  }
  
  return(T_values)
}

# --- Main loop for different T values ---
T_vals <- c(250,5000)
T2_results <- list()
qq_plots <- list()

set.seed(1234)
for (T in T_vals) {
  cat("Simulerer GARCH(1,1) for T =", T, "\n")
  
  T_total <- T + burn_in
  epsilons <- matrix(rstd(M * T_total, mean = 0, sd = 1, nu = df_h0), nrow = M, ncol = T_total)
  sim <- simulate_garch_paths(M, T_total, epsilons)
  
  X_post <- sim$X[, (burn_in + 1):T_total]
  sigma2_post <- sim$sigma2[, (burn_in + 1):T_total]
  X_pre <- sim$X[, burn_in:(T_total - 1)]
  sigma2_pre <- sim$sigma2[, burn_in:(T_total - 1)]
  var_sd_post <- sqrt(w + a * X_pre^2 + b * sigma2_pre)
  
  Var_post <- -qstd(alpha, mean = 0, sd = var_sd_post, nu = df_h0)
  
  # Compute ES numerically
  ES_post <- matrix(0, nrow = M, ncol = T)
  for (i in 1:M) {
    for (j in 1:T) {
      sd_val <- var_sd_post[i, j]
      ES_post[i, j] <- -1 / alpha * integrate(
        function(q) qstd(q, 0, sd_val, nu = df_h0),
        lower = 0, upper = alpha
      )$value
    }
  }
  
  # Compute T2 statistics
  T2_vals <- T2(T, ES_post, Var_post, X_post)
  T2_results[[paste0("T = ", T)]] <- T2_vals
  
  # Generate Danish Q-Q plot
  qq_plots[[paste0("T = ", T)]] <- ggqqplot(T2_vals,
                                            title = paste0("T2(T = ", T, ")"),
                                            xlab = "N(0,1) fraktiler",
                                            ylab = "T2 fraktiler",
                                            color = "blue") + theme_minimal(base_size = 14)
  
}
par(mfrow = c(2, 2))  # 2 rækker, 2 kolonner

for (navn in c("T = 250", "T = 5000")) {
  T_vals <- T2_results[[navn]]  # Hent de simulerede værdier
  
  # Komponent 1
  qqnorm(T_vals[, 1],
         main = paste(navn, "- Komponent 1"),
         xlab = "N(0,1) fraktiler",
         ylab = "T2 fraktiler")
  qqline(T_vals[, 1], col = "red")
  
  # Komponent 2
  qqnorm(T_vals[, 2],
         main = paste(navn, "- Komponent 2"),
         xlab = "N(0,1) fraktiler",
         ylab = "T2 fraktiler")
  qqline(T_vals[, 2], col = "blue")
}

# Optional: Print 5% quantile comparison
for (name in names(T2_results)) {
  cat("\n", name, "\n")
  cat("Empirisk 5%-kvantil:", quantile(T2_results[[name]][,2], 0.05, type = 1), "\n")
  cat("Teoretisk 5%-kvantil:", qnorm(0.05), "\n")
}

