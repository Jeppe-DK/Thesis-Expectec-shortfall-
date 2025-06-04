library(ggplot2)
library(stats)
library(Rlab)
library(fGarch)
library(VaRES)
library(dplyr)
library(gridExtra)
library(knitr)

#Tests 

generate_Z2 <- function(T, df, gamma, alpha,ES, M) {
  Z_values <- numeric(M)
  for (i in 1:M) {
    X_t <- rt(T, df)
    I_t <- ifelse(X_t -qt(alpha,df)*gamma < 0, 1, 0)
    Z_values[i] <- sum(X_t*I_t)/(T*0.025*ES)+1
  }
  return(Z_values)
}
generate_Z3 <- function(T, df, gamma, alpha, M) {
  Z_values <- numeric(M)
  for (i in 1:M){
    Xdata <- rt(T, df)
    EST <- -1/floor(alpha*T)*sum(sort(Xdata)[1:floor(alpha*T)])
    EV <- -T/floor(alpha*T)*integrate(function(p) pbeta(1-p,T-floor(T*alpha),floor(T*alpha))*
                                        qt(p,df),
                                      lower = 0,
                                      upper = 1)$value
    
    Z_values[i] <- -EST/EV+1
  }
  return(Z_values)
}
generate_Z4 <- function(T, df, gamma, alpha, ES, M) {
  Z_values <- numeric(M)
  for (i in 1:M) {
    X_t <- rt(T, df)
    I_t <- ifelse(X_t -qt(alpha,df)*gamma < 0, 1, 0)
    Z_values[i] <- sum(alpha*(ES+qt(alpha,df)*gamma)+(X_t-qt(alpha,df)*gamma)*I_t)/(T*0.025*ES)
  }
  return(Z_values)
}
generate_VaR1 <- function(T,df,gamma,alpha,VaR,M){
  VaR1_values <- numeric(M)
  
  for (i in 1:M){
    X_t <- rt(T, df)
    I_t <- ifelse(X_t -qt(alpha,df)*gamma < 0, 1, 0)*(-1)
    VaR1_values[i] <-sum(I_t)
  }
  return(VaR1_values)
}
generate_T2 <- function(T, df, gamma, alpha, ES, M) {
  T_values <- matrix(0, nrow = M, ncol = 2)
  for (j in 1:M) {
    X_t <-  rt(T, df)
    VaR <- -qt(alpha, df)*gamma 
    
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



-integrate(function(q) qt(q,5), 0, alpha)$value / alpha

# Fixed colors for each test
test_colors <- c(
  "Z1" = "#1b9e77",
  "Z2" = "#d95f02",
  "Z3" = "#7570b3",
  "Z4" = "#e7298a",
  "VaR1" = "#66a61e",
  "T2"="#e1d921"
)

# Fixed linetypes for each test
test_linetypes <- c(
  "Z1" = "solid",
  "Z2" = "dashed",
  "Z3" = "dotted",
  "Z4" = "dotdash",
  "VaR1" = "twodash"
)

# Fixed shapes for each test
test_shapes <- c(
  "Z1" = 16,  # filled circle
  "Z2" = 17,  # filled triangle
  "Z3" = 15,  # filled square
  "Z4" = 18,  # filled diamond
  "VaR1" = 8 , # star
  "T2" = 5 #
)




set.seed(1234)
gamma_seq <- seq(0.5, 1.5, by = 0.05)
M0 <- 5000
M1 <- 1000
T <- 250
alpha <- 0.025



power_list <- list()

# Simulate under H0 once (gamma = 1)
Z2_H0 <- generate_Z2(T, 100, 1, 0.025, 2.378, M0)
Z3_H0 <- generate_Z3(T, 100, 1, 0.025, M0)
Z4_H0 <- generate_Z4(T, 100, 1, 0.025, 2.378, M0)
VaR1_H0 <- generate_VaR1(T, 100, 1, 0.01, 2.36422, M0)
T2_H0 <- generate_T2(T, 100, 1, 0.025, 2.378, M0)

# Critical values at 4.2% level
crit_Z2 <- quantile(Z2_H0, probs = 0.042, type = 1)
crit_Z3 <- quantile(Z3_H0, probs = 0.042, type = 1)
crit_Z4 <- quantile(Z4_H0, probs = 0.042, type = 1)
crit_VaR1 <- quantile(VaR1_H0, probs = 0.042, type = 1)
crit_T2 <- 0.042 


# Loop over gamma values
for (g in gamma_seq) {
  
  Z2_sim <- generate_Z2(T, 100, g, 0.025, 2.378, M1)
  Z3_sim <- generate_Z3(T, 100, g, 0.025, M1)
  Z4_sim <- generate_Z4(T, 100, g, 0.025, 2.378, M1)
  VaR1_sim <- generate_VaR1(T, 100, g, 0.01, 2.36422, M1)
  T2_sim <- generate_T2(T, 100, g, 0.025, 2.378, M1)
  
  power_list[[length(power_list) + 1]] <- data.frame(
    Gamma = g,
    Test = c("Z2", "Z3", "Z4", "VaR1","T2"),
    Power = c(
      mean(Z2_sim <= crit_Z2),
      mean(Z3_sim <= crit_Z3),
      mean(Z4_sim <= crit_Z4),
      mean(VaR1_sim <= crit_VaR1),
      mean(2*apply(1-pnorm(T2_sim),1,min) <= crit_T2)
    )
  )
}


power_df <- do.call(rbind, power_list)




ggplot(power_df, aes(x = Gamma, y = Power,  color = Test, shape = Test, group = Test)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = test_colors) +
  scale_linetype_manual(values = test_linetypes) +
  scale_shape_manual(values = test_shapes) +
  scale_x_continuous(breaks = seq(0.5, 1.5, by = 0.1)) +
  labs(
    title = "gamma*ES og VaR/gamma, korrekt for gamma=1",
    x = expression(gamma),
    y = "Power",
    color = "Test"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  )




