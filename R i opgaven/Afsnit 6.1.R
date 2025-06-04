library(ggplot2)
library(stats)
library(Rlab)
library(fGarch)
library(VaRES)
library(dplyr)
library(gridExtra)
library(knitr)

#Tests 
generate_Z1 <- function(T, df,alpha,ES, M) {
  shift <- qt(0.025, df) -qt(0.025,100)  
  Z_values <- numeric(M)
  for (i in 1:M) {
    X_t <- rt(T,df)-shift
    I_t <- ifelse(X_t -qt(alpha,100) < 0, 1, 0)
    N_t <- sum(I_t)
    Z_values[i] <- (sum(X_t*I_t)/(ES))/(ifelse(N_t<1,1,N_t))+1
  }
  return(Z_values)
}

generate_Z2 <- function(T, df, gamma, alpha,ES, M) {
  Z_values <- numeric(M)
  for (i in 1:M) {
    X_t <- gamma*rt(T, df)
    I_t <- ifelse(X_t -qt(alpha,df) < 0, 1, 0)
    Z_values[i] <- sum(X_t*I_t)/(T*0.025*ES)+1
  }
  return(Z_values)
}

generate_Z2_df <- function(T, df,df1, gamma, alpha,ES, M) {
  Z_values <- numeric(M)
  for (i in 1:M) {
    X_t <- gamma*rt(T, df1)
    I_t <- ifelse(X_t -gamma*qt(alpha,df) < 0, 1, 0)
    Z_values[i] <- sum(X_t*I_t)/(T*0.025*ES)+1
  }
  return(Z_values)
}

generate_Z2_c <- function(T, df,alpha, ES, M) {
  shift <- qt(0.025, df) -qt(0.025,100) 
  Z_values <- numeric(M)
  for (i in 1:M) {
    X_t <- rt(T,df)-shift  
    I_t <- ifelse(X_t - q_2.5_ref < 0, 1, 0)  
    Z_values[i] <- sum(X_t * I_t) / (T * 0.025 * ES) + 1
  }
  return(Z_values)
}

generate_Z3 <- function(T, df, gamma, alpha, M) {
  Z_values <- numeric(M)
  for (i in 1:M){
    Xdata <- rt(T, df)*gamma
    EST <- -1/floor(alpha*T)*sum(sort(Xdata)[1:floor(alpha*T)])
    EV <- -T/floor(alpha*T)*integrate(function(p) pbeta(1-p,T-floor(T*alpha),floor(T*alpha))*
                                        qt(p,df),
                                      lower = 0,
                                      upper = 1)$value
    
    Z_values[i] <- -EST/EV+1
  }
  return(Z_values)
}


generate_Z3_df <- function(T, df,df1, gamma, alpha, M) {
  Z_values <- numeric(M)
  for (i in 1:M){
    Xdata <- rt(T, df1)*gamma
    EST <- -1/floor(alpha*T)*sum(sort(Xdata)[1:floor(alpha*T)])
    EV <- -T/floor(alpha*T)*integrate(function(p) pbeta(1-p,T-floor(T*alpha),floor(T*alpha))*
                                        qt(p,df),
                                      lower = 0,
                                      upper = 1)$value
    
    Z_values[i] <- -EST/EV+1
  }
  return(Z_values)
}

generate_Z3_c <- function(T, df,alpha, M) {
  shift <- qt(0.025, df) -qt(0.025,100)
  Z_values <- numeric(M)
  for (i in 1:M){
    Xdata <- rt(T, df)-shift
    EST <- -1/floor(alpha*T)*sum(sort(Xdata)[1:floor(alpha*T)])
    EV <- -T/floor(alpha*T)*integrate(function(p) pbeta(1-p,T-floor(T*alpha),floor(T*alpha))*
                                        qt(p,100),
                                      lower = 0,
                                      upper = 1)$value
    
    Z_values[i] <- -EST/EV+1
  }
  return(Z_values)
}


generate_Z4 <- function(T, df, gamma, alpha, ES, M) {
  Z_values <- numeric(M)
  for (i in 1:M) {
    X_t <- gamma*rt(T, df)
    I_t <- ifelse(X_t - qt(alpha,df) < 0, 1, 0)
    Z_values[i] <- sum(alpha*(ES+qt(alpha,df))+(X_t-qt(alpha,df))*I_t)/(T*0.025*ES)
  }
  return(Z_values)
}

generate_Z4_df <- function(T, df,df1, gamma, alpha, ES, M) {
  Z_values <- numeric(M)
  for (i in 1:M) {
    X_t <- gamma*rt(T, df1)
    I_t <- ifelse(X_t - qt(alpha,df) < 0, 1, 0)
    Z_values[i] <- sum(alpha*(ES+qt(alpha,df))+(X_t-qt(alpha,df))*I_t)/(T*0.025*ES)
  }
  return(Z_values)
}

generate_Z4_c <- function(T, df, gamma, alpha, ES, M) {
  shift <- qt(0.025, df) -qt(0.025,100)
  Z_values <- numeric(M)
  for (i in 1:M) {
    X_t <- gamma*rt(T, df)-shift
    I_t <- ifelse(X_t - q_2.5_ref < 0, 1, 0)
    Z_values[i] <- sum(alpha*(ES+qt(alpha,df))+(X_t-q_2.5_ref)*I_t)/(T*0.025*ES)
  }
  return(Z_values)
}


generate_VaR1 <- function(T,df,gamma,alpha,VaR,M){
  VaR1_values <- numeric(M)
  
  for (i in 1:M){
    XData <-rt(T,df)*gamma
    I_t <- ifelse(XData +VaR < 0, 1, 0)*(-1)
    VaR1_values[i] <-sum(I_t)
  }
  return(VaR1_values)
}

generate_VaR1_c <- function(T,df,VaR,M){
  shift <- qt(0.025, df) -qt(0.025,100)
  VaR1_values <- numeric(M) 
  for (i in 1:M){
    XData <-rt(T,df)-shift
    I_t <- ifelse(XData +VaR < 0, 1, 0)*(-1)
    VaR1_values[i] <-sum(I_t)
  }
  return(VaR1_values)
}

generate_T2 <- function(T, df, gamma, alpha, ES, M) {
  T_values <- matrix(0, nrow = M, ncol = 2)
  for (j in 1:M) {
    X_t <- gamma * rt(T, df)
    VaR <- -qt(alpha, df)  
    
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


generate_T2_df <- function(T, df,df1, gamma, alpha, ES, M) {
    T_values <- matrix(0, nrow = M, ncol = 2)
    for (j in 1:M) {
      X_t <- gamma * rt(T, df1)
      VaR <- -qt(alpha, df)  
      
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
      
      T_values[j, ] <- (T^(-0.5) * diag(Omega)^(-0.5) * zsum)  
    }
    return(T_values)
  }

generate_T2_c <- function(T, df, alpha, ES, M) {
  shift <- qt(0.025, df) -qt(0.025,100)
  T_values <- matrix(0, nrow = M, ncol = 2)
  for (j in 1:M) {
    X_t <- gamma * rt(T, df)-shift 
    VaR <- -q_2.5_ref  
    
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
    
    T_values[j, ] <- T^(-0.5)*diag(Omega)^(-0.5)*zsum  
  }
  return(T_values)
}





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


#Figur 1  skaleret t-fordeling df=100
gamma_values <- c(1, 1.1, 1.2, 1.3)
x_vals <- seq(-6, 6, length.out = 500)
plot_data <- data.frame()
for (gamma in gamma_values) {
  density_vals <- dt(x_vals / gamma, 100) / gamma
  plot_data <- rbind(plot_data, data.frame(x = x_vals, density = density_vals, gamma = as.factor(gamma)))
}
ggplot(plot_data, aes(x = x, y = density, color = gamma)) +
  labs(title = "", 
       x = "", y = "", color = "Gamma") +
  geom_line(size = 1) +
  theme_minimal() +
  theme(legend.position = "top")


gamma_values <- c(1, 1.1, 1.2,1.3)
alpha <- 0.025

compute_ES <- function(gamma) 
  -integrate(function(q) gamma * qt(q, 100), 0, alpha)$value / alpha

results <- data.frame(
  Gamma = gamma_values,
  VaR1 = -gamma_values * qt(0.01, 100),
  ES = sapply(gamma_values, compute_ES)
)
print(results)

#So for 0.05 quantile we get -5, and the lowest quantile that still generates -5 is 0.042 
set.seed(202503)
y <- generate_VaR1(250,100,1,0.01,2.36422,100000)
quantile(y,seq(0.035,0.05,0.001),type=1)



# Define gamma values to loop over
gamma_seq <- seq(0.8, 1.5, by = 0.05)
M0 <- 10000
M1 <- 2500
alpha <- 0.025


# Storage for results
power_list <- list()

# Simulate under H0 once (gamma = 1)
set.seed(1234)
Z2_H0 <- generate_Z2(250, 100, 1, 0.025, 2.378497, M0)
Z3_H0 <- generate_Z3(250, 100, 1, 0.025, M0)
Z4_H0 <- generate_Z4(250, 100, 1, 0.025, 2.378497, M0)
VaR1_H0 <- generate_VaR1(250, 100, 1, 0.01, 2.36422, M0)
T2_H0 <- generate_T2(250, 100, 1, 0.025, 2.378497, M0)

# Critical values at 4.2% level
crit_Z2 <- quantile(Z2_H0, probs = 0.042, type = 1)
crit_Z3 <- quantile(Z3_H0, probs = 0.042, type = 1)
crit_Z4 <- quantile(Z4_H0, probs = 0.042, type = 1)
crit_VaR1 <- quantile(VaR1_H0, probs = 0.042, type = 1)
crit_T2 <- 0.042 


# Loop over gamma values
for (g in gamma_seq) {
  
  Z2_sim <- generate_Z2(250, 100, g, 0.025, 2.378497, M1)
  Z3_sim <- generate_Z3(250, 100, g, 0.025, M1)
  Z4_sim <- generate_Z4(250, 100, g, 0.025, 2.378497, M1)
  VaR1_sim <- generate_VaR1(250, 100, g, 0.01, 2.36422, M1)
  T2_sim <- generate_T2(250, 100, g, 0.025, 2.378497, M1)
  
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

# Combine into one data frame
power_df <- do.call(rbind, power_list)




ggplot(power_df, aes(x = Gamma, y = Power,  color = Test, shape = Test, group = Test)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = test_colors) +
  scale_linetype_manual(values = test_linetypes) +
  scale_shape_manual(values = test_shapes) +
  scale_x_continuous(breaks = seq(0.8, 1.5, by = 0.1)) +
  labs(
    title = "",
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





#_______________________________________
#Figur 2, varierende antal frihedsgrader
#_______________________________________

df_values <- c(100,10,5,3)
x_vals <- seq(-6, 6, length.out = 500)
plot_data <- data.frame()
for (df in df_values) {
  density_vals <- dt(x_vals,df) 
  plot_data <- rbind(plot_data, data.frame(x = x_vals, density = density_vals, df = as.factor(df)))
}
ggplot(plot_data, aes(x = x, y = density, color = df)) +
  labs(title = "", 
       x = "", y = "", color = "Frihedsgrader") +
  geom_line(size = 1) +
  theme_minimal() +
  theme(legend.position = "top")


# Define gamma values to loop over
set.seed(1234)
df_seq <- c(3,4,5,6,7,8,9,10,12,15,20,30,50,100)
M0 <- 10000  # Under H0
M1 <- 2500   # Under H1
T <- 250

# Storage for results
power_list <- list()

# Simulate under H0 once 
Z2_H0 <- generate_Z2_df(T, 100, 100,1, 0.025, 2.378497, M0)
Z3_H0 <- generate_Z3_df(T, 100, 100,1, 0.025, M0)
Z4_H0 <- generate_Z4_df(T, 100, 100,1, 0.025, 2.378497, M0)
VaR1_H0 <- generate_VaR1(T, 100,1, 0.01, 2.36422, M0)
T2_H0 <- generate_T2_df(T, 100, 100,1, 0.025, 2.378497, M0)

# Critical values at 4.2% level
crit_Z2 <- quantile(Z2_H0, probs = 0.042, type = 1)
crit_Z3 <- quantile(Z3_H0, probs = 0.042, type = 1)
crit_Z4 <- quantile(Z4_H0, probs = 0.042, type = 1)
crit_VaR1 <- quantile(VaR1_H0, probs = 0.042, type = 1)
crit_T2 <- 0.042 #Obs brug  0.042, alternativt  quantile(T2_H0, probs = 0.042, type = 1)

# Loop over df values
for (d in df_seq) {
  
  Z2_sim <- generate_Z2_df(T,100,d, 1, 0.025, 2.378497, M1)
  Z3_sim <- generate_Z3_df(T,100,d, 1, 0.025, M1)
  Z4_sim <- generate_Z4_df(T,100,d, 1, 0.025, 2.378497, M1)
  VaR1_sim <- generate_VaR1(T, d,1, 0.01, 2.36422, M1)
  T2_sim <- generate_T2_df(T,100,d, 1, 0.025, 2.378497, M1)
  
  power_list[[length(power_list) + 1]] <- data.frame(
    df = d,  
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

# Combine into one data frame
power_df <- do.call(rbind, power_list)

power_df$df_label <- as.character(power_df$df)
power_df$df_label[power_df$df == 100000] <- "normal"

# Plot
ggplot(power_df, aes(x = factor(df_label, levels = c(as.character(3:10), "12", "15", "20", "30", "50","100")),
                     y = Power, color = Test, shape= Test, group = Test)) +
  geom_line(aes(group = Test), linewidth = 1) +
  scale_color_manual(values = test_colors) +
  scale_linetype_manual(values = test_linetypes) +
  scale_shape_manual(values = test_shapes) +
  geom_point(size = 2) +
  labs(
    title = "",
    x = expression(df),
    y = "Power",
    color = "Test"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  )






#Figur 3
dfs <- c(100, 10, 5, 3)
x_range <- seq(-6, 6, length.out = 1000)
q_2.5_ref <- qt(0.025, df = 100)
shifts <- qt(0.025, df = dfs) - q_2.5_ref


data <- do.call(rbind, lapply(1:length(dfs), function(i) 
  data.frame(x = x_range, density = dt(x_range + shifts[i], df = dfs[i]), df = factor(dfs[i]))))

print(data.frame(df = dfs, Original_2.5_Quantile = qt(0.025, df = dfs), Shifted_2.5_Quantile = q_2.5_ref))

ggplot(data, aes(x = x, y = density, color = df)) +
  geom_line(linewidth = 1) + geom_vline(xintercept = q_2.5_ref, linetype = "dashed", color = "black") +
  labs(title = "", x = "x", y = "", color = "Frihedsgrader") +
  theme_minimal()


# Degrees of freedom
set.seed(1234)
q_2.5_ref <- qt(0.025, df = 100)
gamma <- 1
dfs <- c(3,4,5,6,7,8,9,10,12,15,20,30,50,100)
M0 <- 10000  
M1 <- 2500 

# Storage for results
power_list <- list()

# Simulate under H0 (df = 100)
Z1_H0 <- generate_Z1(250, 100, 0.025, 2.37, M0)
Z2_H0 <- generate_Z2_c(250, 100, 0.025, 2.37, M0)
Z3_H0 <- generate_Z3_c(250, 100, 0.025, M0)
Z4_H0 <- generate_Z4_c(250, 100, 1, 0.025, 2.37, M0)
VaR1_H0 <- generate_VaR1_c(250, 100, 2.36, M0)
T2_H0 <- generate_T2_c(250, 100, 0.025, 2.37, M0)

# Critical values at 4.2% level
crit_Z1 <- quantile(Z1_H0, probs = 0.042, type = 1)
crit_Z2 <- quantile(Z2_H0, probs = 0.042, type = 1)
crit_Z3 <- quantile(Z3_H0, probs = 0.042, type = 1)
crit_Z4 <- quantile(Z4_H0, probs = 0.042, type = 1)
crit_VaR1 <- quantile(VaR1_H0, probs = 0.042, type = 1)
crit_T2 <- 0.042 

# Loop over dfs
for (d in dfs) {
  
  Z1_sim <- generate_Z1(250, d, 0.025, 2.37, M1)
  Z2_sim <- generate_Z2_c(250, d, 0.025, 2.37, M1)
  Z3_sim <- generate_Z3_c(250, d, 0.025, M1)
  Z4_sim <- generate_Z4_c(250, d, 1, 0.025, 2.37, M1)
  VaR1_sim <- generate_VaR1_c(250, d, 2.36, M1)
  T2_sim <- generate_T2_c(250, d, 0.025, 2.37, M1)
  
  power_list[[length(power_list) + 1]] <- data.frame(
    df = d,
    Test = c("Z1", "Z2", "Z3", "Z4", "VaR1","T2"),
    Power = c(
      mean(Z1_sim <= crit_Z1),
      mean(Z2_sim <= crit_Z2),
      mean(Z3_sim <= crit_Z3),
      mean(Z4_sim <= crit_Z4),
      mean(VaR1_sim <= crit_VaR1),
      mean(2*apply(1-pnorm(T2_sim),1,min) <= crit_T2)
    )
  )
}

# Combine into one data frame
power_df3 <- do.call(rbind, power_list)


# Plot
ggplot(power_df3, aes(x = factor(df, levels = c(3,4,5,6,7,8,9,10,12,15,20,30,50,100)), y = Power, color = Test,shape = Test, group = Test)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = test_colors) +
  scale_linetype_manual(values = test_linetypes) +
  scale_shape_manual(values = test_shapes) +
  labs(
    title = "",
    x = expression(df),
    y = "Power",
    color = "Test"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  )



