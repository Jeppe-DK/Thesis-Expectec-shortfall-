library(ggplot2)
library(stats)
library(Rlab)
library(fGarch)
library(VaRES)
library(dplyr)
library(gridExtra)
library(knitr)
library(tidyr)

#Tests 
generate_Z2 <- function(T, df, sigma, alpha, M) {
  Z_values <- numeric(M)
  ES_val <- -integrate(function(q) qstd(q, mean = 0, sd = 1, nu = 25), 0, alpha)$value / alpha
  VaR_threshold <- qstd(alpha, mean = 0, sd = 1, nu = 25)
  for (i in 1:M) {
    X_t <- rstd(T, 0,sigma,df)
    I_t <- ifelse(X_t - VaR_threshold < 0, 1, 0)
    Z_values[i] <- sum(X_t * I_t) / (T*0.025*ES_val) + 1
  }
  return(Z_values)
}
generate_Z3 <- function(T, df, sigma, alpha, M) {
  Z_values <- numeric(M)
  for (i in 1:M){
    Xdata <- rstd(T,0,sigma, df)
    EST <- -1/floor(alpha*T)*sum(sort(Xdata)[1:floor(alpha*T)])
    EV <- -T/floor(alpha*T)*integrate(function(p) pbeta(1-p,T-floor(T*alpha),floor(T*alpha))*
                                        qstd(p,0,1,25),
                                      lower = 0,
                                      upper = 1)$value
    
    Z_values[i] <- -EST/EV+1
  }
  return(Z_values)
}
generate_Z4 <- function(T, df, sigma, alpha, M) {
  ES_val <- -integrate(function(q) qstd(q, mean = 0, sd = 1, nu = 25), 0, alpha)$value / alpha
  Z_values <- numeric(M)
  for (i in 1:M) {
    X_t <- rstd(T,0,sigma, df)
    I_t <- ifelse(X_t - qstd(alpha,0,1,25) < 0, 1, 0)
    Z_values[i] <- sum(alpha*(ES_val+qstd(alpha,0,1,25))+(X_t-qstd(alpha,0,1,25))*I_t)/(T*0.025*ES_val)
  }
  return(Z_values)
}
generate_VaR1 <- function(T,df,sigma,alpha,M){
  VaR1_values <- numeric(M)
  for (i in 1:M){
    XData <-rstd(T,0,sigma,df)
    I_t <- ifelse(XData - qstd(0.01,0,1,25) < 0, 1, 0)*(-1)
    VaR1_values[i] <-sum(I_t)
  }
  return(VaR1_values)
}
generate_T2 <- function(T, df, sigma, alpha, M) {
  T_values <- matrix(0, nrow = M, ncol = 2)
  ES <- -integrate(function(q) qstd(q, mean = 0, sd = 1, nu = 25), 0, alpha)$value / alpha
  for (j in 1:M) {
    X_t <- rstd(T, 0,sigma,df)
    VaR <- -qstd(alpha, mean = 0, sd = 1, nu = 25)  
    
    I_t <- ifelse(X_t + VaR < 0, 1, 0)
    v1 <- alpha - I_t
    v2 <- VaR- ES + (1 / alpha) * (-VaR - X_t) * I_t
    v <- cbind(-v1, v2)
    
    Omega_sum <- matrix(0, nrow = 2, ncol = 2)
    
    for (i in 1:T) {
      vi <- matrix(v[i, ], ncol = 1)  
      Omega_sum <- Omega_sum + vi %*% t(vi)
    }
    
    Omega <- Omega_sum/T
    zsum <- colSums(v)
    
    T_values[j, ] <- T^(-0.5) * diag(Omega)^(-0.5) * zsum   
  }
  return(T_values)
}


#Kritisk værdi er 4.3% 
set.seed(202503)
quantile(generate_VaR1(250,25,1,0.01,10000),seq(0.035,0.05,0.001),type=1)

# Generalized function
plot_power <- function(generator, title, df_seq, sigma_seq, M_H0 = 10000, M_H1 = 2500, alpha = 0.025) {
  set.seed(123)
  
  # Step 1: Generate critical value under H0
  Z_H0 <- generator(250, 100, 1, alpha, M_H0)
  critical_value <- quantile(Z_H0, 0.043, type = 1)
  
  # Step 2: Prepare grid
  grid <- expand.grid(df = df_seq, sigma = sigma_seq)
  
  # Step 3: Compute Power
  grid$acceptance_rate <- mapply(function(df_val, sigma_val) {
    Z_H1 <- generator(250, df_val, sigma_val, alpha, M_H1)
    ecdf(Z_H1)(critical_value)  # One-sided acceptance rate
  }, grid$df, grid$sigma)
  
  # --- Step 4: Tile plot with continuous color scale ---
  p1 <- ggplot(grid, aes(x = df, y = sigma, fill = acceptance_rate)) +
    geom_tile() +
    scale_fill_viridis_c(option = "D", limits = c(0,1), name = "Power", direction = -1) +
    labs(title = paste(title),
         x = "Frihedsgrader",
         y = "Sigma") +
    theme_minimal() +
    theme(legend.position = "right")
  
  print(p1)
} 
set.seed(2025)


#Z2
plot_power(generate_Z2, "Power for Z2",
           df_seq = seq(3, 25, by = 2),
           sigma_seq = seq(0.8, 1.6, by = 0.05))

#Z3
plot_power(generate_Z3, "Power for Z3",
           df_seq = seq(3, 25, by = 2),
           sigma_seq = seq(0.8, 1.6, by = 0.05))

#Z4
plot_power(generate_Z4, "Power for Z4",
           df_seq = seq(3, 25, by = 2),
           sigma_seq = seq(0.8, 1.6, by = 0.05))

# VaR1
plot_power(generate_VaR1, "Power for VaR1",
           df_seq = seq(3, 25, by = 2),
           sigma_seq = seq(0.8, 1.6, by = 0.05))

# T2
plot_power <- function(generator, title, df_seq, sigma_seq, M_H0 = 10000, M_H1 = 2500, alpha = 0.025) {
  set.seed(123)
  
  # Step 1: Generate critical value under H0
  Z_H0 <- generator(250, 100, 1, alpha, M_H0)
  critical_value <-  0.043 
  
  # Step 2: Prepare grid
  grid <- expand.grid(df = df_seq, sigma = sigma_seq)
  
  # Step 3: Compute Power
  grid$acceptance_rate <- mapply(function(df_val, sigma_val) {
    Z_H1 <- generator(250, df_val, sigma_val, alpha, M_H1)
    mean(2*apply(1-pnorm(Z_H1),1,min) <= critical_value)  
  }, grid$df, grid$sigma)
  
  # --- Step 4: Tile plot with continuous color scale ---
  p1 <- ggplot(grid, aes(x = df, y = sigma, fill = acceptance_rate)) +
    geom_tile() +
    scale_fill_viridis_c(option = "D", limits = c(0,1), name = "Power", direction = -1) +
    labs(title = paste(title),
         x = "Frihedsgrader",
         y = "Sigma") +
    theme_minimal() +
    theme(legend.position = "right")
  
  print(p1)
} 
plot_power(generate_T2, "Power for T2",
           df_seq = seq(3, 25, by = 2),
           sigma_seq = seq(0.8, 1.6, by = 0.05))



