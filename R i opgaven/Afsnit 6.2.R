library(ggplot2)
library(stats)
library(Rlab)
library(fGarch)
library(VaRES)
library(dplyr)
library(gridExtra)
library(knitr)
#For student t-distribution df 100 vs 10 
generate_Z1 <- function(T, df, alpha, M) {
  ES <- -integrate(function(q) qt(q, 100, 0), 0, alpha)$value / alpha
  Z_values <- numeric(M)
  for (i in 1:M) {
    X_t <- rt(T, df)
    I_t <- ifelse(X_t - qt(0.025, 100, 0) < 0, 1, 0)
    N_t <- sum(I_t)
    Z_values[i] <- (sum(X_t * I_t) / ES) / (ifelse(N_t < 1, 1, N_t)) + 1
  }
  return(Z_values)
}
generate_Z1_fixed <- function(T, n, df, alpha, M) {
  ES_val <- -integrate(function(q) qt(q,100), 0, alpha)$value / alpha
  a <- pt(qt(alpha, 100), df)
  Z_values <- numeric(M)
  for (i in 1:M) {
    X_t <- c(qt(runif(n, 0, a), df, 0), qt(runif(T - n, a, 1), df, 0))-(qt(0.025, df) -qt(0.025,100))
    I_t <- ifelse(X_t - qt(alpha,100) < 0, 1, 0)
    N_t <- sum(I_t)
    Z_values[i] <- (sum(X_t * I_t) / ES_val) / (ifelse(N_t < 1, 1, N_t)) + 1
  }
  return(Z_values)
} #OBS nu med shift 

generate_Z2 <- function(T, df, alpha, M) {
  ES <- -integrate(function(q) qt(q, 100), 0, alpha)$value / alpha
  Z_values <- numeric(M)
  for (i in 1:M) {
    X_t <- rt(T, df)
    I_t <- ifelse(X_t - qt(alpha,100) < 0, 1, 0)
    Z_values[i] <- sum(X_t * I_t) / (T * 0.025 * ES) + 1
  }
  return(Z_values)
}
generate_Z2_fixed <- function(T, n, df, alpha, M) {
  ES <- -integrate(function(q) qt(q,100), 0, alpha)$value / alpha
  a <- pt(qt(alpha, 100), df)
  Z_values <- numeric(M)
  for (i in 1:M) {
    X_t <- c(qt(runif(n, 0, a), df, 0), qt(runif(T - n, a, 1), df, 0))
    I_t <- ifelse(X_t - qt(alpha, 100) < 0, 1, 0)
    Z_values[i] <- sum(X_t * I_t) / (T * 0.025 * ES) + 1
  }
  return(Z_values)
}

generate_Z3 <- function(T, df, alpha, M) {
  Z_values <- numeric(M)
  for (i in 1:M) {
    Xdata <- rt(T, df)
    EST <- -1 / floor(alpha * T) * sum(sort(Xdata)[1:floor(alpha * T)])
    EV <- -T / floor(alpha * T) * integrate(function(p) pbeta(1 - p, T - floor(T * alpha), floor(T * alpha)) *
                                              qt(p, 100, 0), lower = 0, upper = 1)$value
    Z_values[i] <- -EST / EV + 1
  }
  return(Z_values)
}
generate_Z3_fixed <- function(T, n, df, alpha, M) {
  a <- pt(qt(alpha,100), df)
  Z_values <- numeric(M)
  for (i in 1:M) {
    Xdata <- c(qt(runif(n, 0, a), df, 0), qt(runif(T - n, a, 1), df, 0))
    EST <- -1 / floor(alpha * T) * sum(sort(Xdata)[1:floor(alpha * T)])
    EV <- -T / floor(alpha * T) * integrate(function(p) pbeta(1 - p, T - floor(T * alpha), floor(T * alpha)) *
                                              qt(p,100, 0), lower = 0, upper = 1)$value
    Z_values[i] <- -EST / EV + 1
  }
  return(Z_values)
}

generate_Z4 <- function(T, df, alpha, M) {
  ES <- -integrate(function(q) qt(q,100), 0, alpha)$value / alpha
  Z_values <- numeric(M)
  for (i in 1:M) {
    X_t <- rt(T, df)
    I_t <- ifelse(X_t - qt(alpha, df, 0) < 0, 1, 0)
    Z_values[i] <- sum(alpha * (ES + qt(alpha, df, 0)) + (X_t - qt(alpha, df, 0)) * I_t) / (T * 0.025 * ES)
  }
  return(Z_values)
}
generate_Z4_fixed <- function(T, n, df, alpha, M) {
  ES <- -integrate(function(q) qt(q,100), 0, alpha)$value / alpha
  a <- pt(qt(alpha, 100), df)
  Z_values <- numeric(M)
  for (i in 1:M) {
    X_t <- c(qt(runif(n, 0, a), df, 0), qt(runif(T - n, a, 1), df, 0))
    I_t <- ifelse(X_t - qt(alpha, 100, 0) < 0, 1, 0)
    Z_values[i] <- sum(alpha * (ES + qt(alpha, 100)) + (X_t - qt(alpha,100)) * I_t) / (T * 0.025 * ES)
  }
  return(Z_values)
}

generate_T2 <- function(T, df,alpha, M) {
  ES <- -integrate(function(q) qt(q, 100), 0, alpha)$value / alpha
  T_values <- matrix(0, nrow = M, ncol = 2)
  for (j in 1:M) {
    X_t <-  rt(T, df)
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
    
    T_values[j, ] <- T^(-0.5) * diag(Omega)^(-0.5) * zsum   
  }
  return(T_values)
}
generate_T2_fixed <- function(T, n, df, alpha, M) {
  ES <- -integrate(function(q) qt(q, 100), 0, alpha)$value / alpha
  a <- pt(qt(alpha, 100), df)
  T_values <- matrix(0, nrow = M, ncol = 2)
  for (j in 1:M) {
    X_t <- c(qt(runif(n, 0, a), df, 0), qt(runif(T - n, a, 1), df, 0))
    VaR <- -qt(alpha, 100)  # use consistent naming
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
    
    T_values[j, ] <- T^(-0.5) * diag(Omega)^(-0.5) * zsum  
  }
  return(T_values)
}


#-----------------------------
# Generic simulation & plot function
#-----------------------------

simulate_plot <- function(generate_Z, generate_Z_fixed, label, M = 10000, T = 250, alpha = 0.025) {
  n_seq <- 1:25
  set.seed(1)
  
  Z_H0 <- generate_Z(T, 100, alpha, M)  # df=100 for H0
  crit <- sort(Z_H0)[floor(M * 0.05)]
  
  ECDF_vals <- sapply(n_seq, function(n) {
    Z_H1 <- generate_Z_fixed(T, n, 10, alpha, M)  # df=10 for H1
    mean(Z_H1 < crit)
  })
  
  plot(n_seq, ECDF_vals, type = "l", col = "black", ylim = c(0, 1),
       main = paste(label),
       xlab = "n", ylab = "Power")
}


par(mfrow = c(1, 1))  

simulate_plot(generate_Z1, generate_Z1_fixed, "Z1")
abline(v = 9, col = "red", lwd = 2, lty = 2)
simulate_plot(generate_Z2, generate_Z2_fixed, "Z2")
abline(v = 9, col = "red", lwd = 2, lty = 2)
simulate_plot(generate_Z3, generate_Z3_fixed, "Z3")
abline(v = 9, col = "red", lwd = 2, lty = 2)
simulate_plot(generate_Z4, generate_Z4_fixed, "Z4")
abline(v = 9, col = "red", lwd = 2, lty = 2)


simulate_plot <- function(generate_Z, generate_Z_fixed, label, M = 10000, T = 250, alpha = 0.025) {
  n_seq <- 1:25
  set.seed(1)
  
  crit <-  0.05 
  
  ECDF_vals <- sapply(n_seq, function(n) {
    Z_H1 <- generate_Z_fixed(T, n, 10, alpha, M)  # df=10 for H1
    mean(2*apply(1-pnorm(Z_H1),1,min) <= crit)
  })
  
  plot(n_seq, ECDF_vals, type = "l", col = "black", ylim = c(0, 1),
       main = paste(label),
       xlab = "n", ylab = "Power")
}
simulate_plot(generate_T2, generate_T2_fixed, "T2")
abline(v = 9, col = "red", lwd = 2, lty = 2)

