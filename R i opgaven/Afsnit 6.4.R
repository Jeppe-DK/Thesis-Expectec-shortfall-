library(ggplot2)
library(stats)
library(Rlab)
library(fGarch)
library(VaRES)
library(dplyr)
library(gridExtra)
library(knitr)



location_values <- c(-1, 0, 1)
df_values <- c(3, 5, 10, 100)
alpha_values <- c(0.025, 0.01)
options(digits=5)

compute_ES <- function(location, df, alpha) {
  integrand <- function(q) location + qt(q, df)
  -integrate(integrand, 0, alpha)$value / alpha
}

results <- data.frame()

for (location in location_values) {
  for (df in df_values) {
    for (alpha in alpha_values) {
      VaR <- -(location + qt(alpha, df))
      ES <- compute_ES(location, df, alpha)
      results <- rbind(results, data.frame(
        Location = location,
        DF = df,
        Alpha = alpha,
        VaR = VaR,
        ES = ES
      ))
    }
  }
}


print(results)
#__________________________________________________
#Undersøg for Z1
set.seed(12345)
generate_Z1 <- function(T, df,location,alpha,ES, M) {
  Z_values <- numeric(M)
  for (i in 1:M) {
    X_t <- rt(T,df)+location
    I_t <- ifelse(X_t -(qt(alpha,df)+location) < 0, 1, 0)
    N_t <- sum(I_t)
    Z_values[i] <- sum(X_t*I_t)/(ES*(ifelse(N_t<1,1,N_t)))+1
  }
  return(Z_values)
}

df3_1 <- generate_Z1(250,3,-1,0.025,6.0396,10000)
df3_0 <- generate_Z1(250,3,0,0.025,5.0396,10000)
df3__1 <- generate_Z1(250,3,1,0.025,4.0396,10000)

df5_1 <- generate_Z1(250,5,-1,0.025,4.52158,10000)
df5_0 <- generate_Z1(250,5,0,0.025,3.52158,10000)
df5__1 <- generate_Z1(250,5,1,0.025,2.521582,10000)

df10_1 <- generate_Z1(250,10,-1,0.025,3.8190,10000)
df10_0 <- generate_Z1(250,10,0,0.025,2.8190,10000)
df10__1 <- generate_Z1(250,10,1,0.025,1.8190,10000)

df100_1 <- generate_Z1(250,100,-1,0.025,3.37850,10000)
df100_0 <- generate_Z1(250,100,0,0.025,2.37850,10000)
df100__1 <- generate_Z1(250,100,1,0.025,1.37850,10000)



#Undersøg om den kritiske værdi er stabil for 5% 
options(digits=2)
results_matrix <- data.frame(
  DF = c(3, 5, 10, 100),
  Location_minus1 = c(
    quantile(df3_1, 0.05),
    quantile(df5_1, 0.05),
    quantile(df10_1, 0.05),
    quantile(df100_1, 0.05)
  ),
  Location_0 = c(
    quantile(df3_0, 0.05),
    quantile(df5_0, 0.05),
    quantile(df10_0, 0.05),
    quantile(df100_0, 0.05)
  ),
  Location_1 = c(
    quantile(df3__1, 0.05),
    quantile(df5__1, 0.05),
    quantile(df10__1, 0.05),
    quantile(df100__1, 0.05)
  )
)
print(results_matrix)

#__________________________________________________
#Undersøg for Z2
set.seed(12345)
generate_Z2 <- function(T, df, location, alpha,ES, M) {
  Z_values <- numeric(M)
  for (i in 1:M) {
    X_t <- location+rt(T, df)
    I_t <- ifelse(X_t+(-(qt(alpha,df)+location)) < 0, 1, 0)
    Z_values[i] <- sum(X_t*I_t)/(T*0.025*ES)+1
  }
  return(Z_values)
}

df3_1 <- generate_Z2(250,3,-1,0.025,6.0396,10000)
df3_0 <- generate_Z2(250,3,0,0.025,5.0396,10000)
df3__1 <- generate_Z2(250,3,1,0.025,4.0396,10000)

df5_1 <- generate_Z2(250,5,-1,0.025,4.52158,10000)
df5_0 <- generate_Z2(250,5,0,0.025,3.52158,10000)
df5__1 <- generate_Z2(250,5,1,0.025,2.521582,10000)

df10_1 <- generate_Z2(250,10,-1,0.025,3.8190,10000)
df10_0 <- generate_Z2(250,10,0,0.025,2.8190,10000)
df10__1 <- generate_Z2(250,10,1,0.025,1.8190,10000)

df100_1 <- generate_Z2(250,100,-1,0.025,3.37850,10000)
df100_0 <- generate_Z2(250,100,0,0.025,2.37850,10000)
df100__1 <- generate_Z2(250,100,1,0.025,1.37850,10000)

#Undersøg om den kritiske værdi er stabil for 5% 
options(digits=2)
results_matrix <- data.frame(
  DF = c(3, 5, 10, 100),
  Location_minus1 = c(
    quantile(df3_1, 0.05),
    quantile(df5_1, 0.05),
    quantile(df10_1, 0.05),
    quantile(df100_1, 0.05)
  ),
  Location_0 = c(
    quantile(df3_0, 0.05),
    quantile(df5_0, 0.05),
    quantile(df10_0, 0.05),
    quantile(df100_0, 0.05)
  ),
  Location_1 = c(
    quantile(df3__1, 0.05),
    quantile(df5__1, 0.05),
    quantile(df10__1, 0.05),
    quantile(df100__1, 0.05)
  )
)
print(results_matrix)

#__________________________________________________
#Undersøg for Z3
set.seed(12345)
generate_Z3 <- function(T, df, location, alpha, M) {
  Z_values <- numeric(M)
  for (i in 1:M){
    Xdata <- rt(T, df)+location
    EST <- -1/floor(alpha*T)*sum(sort(Xdata)[1:floor(alpha*T)])
    EV <- -T/floor(alpha*T)*integrate(function(p) pbeta(1-p,T-floor(T*alpha),floor(T*alpha))*
                                        (location+qt(p,df)),
                                      lower = 0,
                                      upper = 1)$value
    
    Z_values[i] <- -EST/EV+1
  }
  return(Z_values)
}

df3_1 <- generate_Z3(250,3,-1,0.025,10000)
df3_0 <- generate_Z3(250,3,0,0.025,10000)
df3__1 <- generate_Z3(250,3,1,0.025,10000)

df5_1 <- generate_Z3(250,5,-1,0.025,10000)
df5_0 <-  generate_Z3(250,5,0,0.025,10000)
df5__1 <-  generate_Z3(250,5,1,0.025,10000)

df10_1 <- generate_Z3(250,10,-1,0.025,10000)
df10_0 <- generate_Z3(250,10,0,0.025,10000)
df10__1 <- generate_Z3(250,10,1,0.025,10000)

df100_1 <- generate_Z3(250,100,-1,0.025,10000)
df100_0 <- generate_Z3(250,100,0,0.025,10000)
df100__1 <- generate_Z3(250,100,1,0.025,10000)

#Undersøg om den kritiske værdi er stabil for 5% 
options(digits=2)
results_matrix <- data.frame(
  DF = c(3, 5, 10, 100),
  Location_minus1 = c(
    quantile(df3_1, 0.05),
    quantile(df5_1, 0.05),
    quantile(df10_1, 0.05),
    quantile(df100_1, 0.05)
  ),
  Location_0 = c(
    quantile(df3_0, 0.05),
    quantile(df5_0, 0.05),
    quantile(df10_0, 0.05),
    quantile(df100_0, 0.05)
  ),
  Location_1 = c(
    quantile(df3__1, 0.05),
    quantile(df5__1, 0.05),
    quantile(df10__1, 0.05),
    quantile(df100__1, 0.05)
  )
)
print(results_matrix)



#__________________________________________________
#Undersøg for Z4
set.seed(12345)
generate_Z4 <- function(T, df, location, alpha, ES, M) {
  Z_values <- numeric(M)
  for (i in 1:M) {
    X_t <- location+rt(T, df)
    VaR_t <- -(qt(alpha,df)+location)
    I_t <- ifelse(X_t +VaR_t< 0, 1, 0)
    Z_values[i] <- sum(alpha*(ES-VaR_t)+(X_t+VaR_t)*I_t)/(T*0.025*ES)
  }
  return(Z_values)
}

df3_1 <- generate_Z4(250,3,-1,0.025,6.04,10000)
df3_0 <- generate_Z4(250,3,0,0.025,5.04,10000)
df3__1 <- generate_Z4(250,3,1,0.025,4.04,10000)

df5_1 <- generate_Z4(250,5,-1,0.025,4.52158,10000)
df5_0 <- generate_Z4(250,5,0,0.025,3.52158,10000)
df5__1 <- generate_Z4(250,5,1,0.025,2.521582,10000)

df10_1 <- generate_Z4(250,10,-1,0.025,3.8190,10000)
df10_0 <- generate_Z4(250,10,0,0.025,2.8190,10000)
df10__1 <- generate_Z4(250,10,1,0.025,1.8190,10000)

df100_1 <- generate_Z4(250,100,-1,0.025,3.37850,10000)
df100_0 <- generate_Z4(250,100,0,0.025,2.37850,10000)
df100__1 <- generate_Z4(250,100,1,0.025,1.37850,10000)

#Undersøg om den kritiske værdi er stabil for 5% 
options(digits=2)
results_matrix <- data.frame(
  DF = c(3, 5, 10, 100),
  Location_minus1 = c(
    quantile(df3_1, 0.05),
    quantile(df5_1, 0.05),
    quantile(df10_1, 0.05),
    quantile(df100_1, 0.05)
  ),
  Location_0 = c(
    quantile(df3_0, 0.05),
    quantile(df5_0, 0.05),
    quantile(df10_0, 0.05),
    quantile(df100_0, 0.05)
  ),
  Location_1 = c(
    quantile(df3__1, 0.05),
    quantile(df5__1, 0.05),
    quantile(df10__1, 0.05),
    quantile(df100__1, 0.05)
  )
)
print(results_matrix)









#_____________________________________________________________________________
#Undersøger for ændringer i DF 
set.seed(1234)

Z1 <- function(T, ES_matrix, Var_matrix, X_matrix) {
  I_t <- ifelse(X_matrix + Var_matrix < 0, 1, 0)
  N_t <- rowSums(I_t)
  Z_values <- rowSums((X_matrix * I_t) / (ES_matrix * ifelse(N_t < 1, 1, N_t))) + 1
  return(Z_values)
}
Z2 <- function(T, ES_matrix, Var_matrix, X_matrix) {
  I_t <- ifelse(X_matrix + Var_matrix < 0, 1, 0)  
  Z_values <- rowSums((X_matrix * I_t) / (T * 0.025 * ES_matrix)) + 1  
  return(Z_values)
}
Z3 <- function(T, X_matrix, alpha, df_val, var_sd_post) {
  eps <- 1e-6  # to avoid boundary issues
  
  if (df_val == "Normal") stop("This logic requires a t-distribution df_val.")
  
  df_num <- as.numeric(df_val)
  if (is.na(df_num)) stop(paste("Invalid df_val:", df_val))
  
  alpha_floor <- floor(alpha * T)
  
  # Compute integral value once
  
  # Compute EST and EV for each simulation
  results <- sapply(1:nrow(X_matrix), function(i) {
    x <- X_matrix[i, ]
    sd_row <- var_sd_post[i, ]
    
    # Step 1: Compute Y_i^(t)
    pt_scaled <- pstd(x,0, sd_row, nu = df_num) 
    # Step 2: For each time t, get EST[t] from sorted values
    EST_vec <- rep(NA, T)
    for (t in 1:T) {
      EST_vec[t] <- -mean(sort(qstd(pt_scaled,0,sd_row[t], nu = df_num))[1:alpha_floor])
    }
    
    # Step 3: EV is a vector
    integral_val <- integrate(function(p) {
      pbeta(1 - p, T - alpha_floor, alpha_floor) * qstd(p,0,sd_row[t], nu = df_num)
    }, lower = 0, upper = 1 )$value
    EV_vec <- (-T/alpha_floor)*integral_val
    
    # Step 4: Sum over EST[t] / EV[t]
    mean(EST_vec / EV_vec)
  })
  
  # Final Z3 value
  Z3_value <- 1 - results
  
  return(Z3_value)
}
Z4 <- function(T, ES_matrix, Var_matrix, X_matrix){
  I_t <- ifelse(X_matrix + Var_matrix < 0, 1, 0) 
  Z_values <- rowMeans(1-(Var_matrix/ES_matrix))+rowSums((X_matrix+Var_matrix)*I_t/(0.025*T*ES_matrix))
  return(Z_values)
}
VaR1 <- function(X_matrix,Var_matrix){
  I_t <- ifelse(X_matrix + Var_matrix < 0, 1, 0)*(-1) 
  VaR1_values <-rowSums(I_t)
  return(VaR1_values)
}




#Undersøger for ændringer i DF 
# Parameters
T <- 250
M <- 1000
burn_in <- 250
T_total <- T + burn_in
alpha <- 0.025

# GARCH(1,1) parameters
w <- 0.01
a <- 0.1
b <- 0.85

# Degrees of freedom to test under H1
df_values <- c(3, 5,10, 20, 50)
df_h0 <- 10

# Storage lists
X_list <- list()
Var_list <- list()
ES_list <- list()
var_sd_list <- list()  # NEW: store var_sd_post here

set.seed(123)

# Function to simulate GARCH paths
simulate_garch_paths <- function(M, T_total, epsilons) {
  X <- matrix(0, nrow = M, ncol = T_total)
  sigma2 <- matrix(0, nrow = M, ncol = T_total)
  
  # Initial conditions
  sigma2[,1] <- 0.05
  X[,1] <- 0
  
  for (t in 2:T_total) {
    sigma2[,t] <- w + a * X[,t-1]^2 + b * sigma2[,t-1]
    X[,t] <- sqrt(sigma2[,t]) * epsilons[,t]
  }
  
  list(X = X, sigma2 = sigma2)
}

### ----- H0 Simulation ----- ###
cat("Simulating H0 (t-distribution with df=10)\n")

epsilons <- matrix(rstd(M * T_total,0,1,nu = df_h0), nrow = M, ncol = T_total)
sim_result <- simulate_garch_paths(M, T_total, epsilons)

# Only keep post-burn-in
X_post <- sim_result$X[, (burn_in + 1):T_total]
sigma2_post <- sim_result$sigma2[, (burn_in + 1):T_total]

# Compute var_sd directly
var_sd_post <- sqrt(w + a * sim_result$X[, (burn_in):(T_total-1)]^2 + b * sim_result$sigma2[, (burn_in):(T_total-1)])

# Compute VaR and ES
Var_post <- -qstd(alpha,0,var_sd_post,nu = df_h0)
ES_post_mat <- matrix(0, nrow = M, ncol = T)
for (i in 1:M) {
  for (j in 1:T) {
    sd_val <- var_sd_post[i, j]
    ES_post_mat[i, j] <- -1 / alpha * integrate(
      function(q) qstd(q, 0, sd_val, nu = df_h0),
      lower = 0, upper = alpha
    )$value
  }
}



# Store
X_list[["H0"]] <- X_post
Var_list[["H0"]] <- Var_post
ES_list[["H0"]] <- ES_post_mat
var_sd_list[["H0"]] <- var_sd_post   # Store var_sd

### ----- H1 Simulations ----- ###
for (df_val in df_values) {
  cat("Simulating H1 with df =", df_val, "\n")
  
  if (df_val == "Normal") {
    epsilons <- matrix(rnorm(M * T_total), nrow = M, ncol = T_total)
  } else {
    epsilons <- matrix(rstd(M * T_total, 0,1, nu = as.numeric(df_val)), nrow = M, ncol = T_total)
  }
  
  sim_result <- simulate_garch_paths(M, T_total, epsilons)
  
  X_post <- sim_result$X[, (burn_in + 1):T_total]
  sigma2_post <- sim_result$sigma2[, (burn_in + 1):T_total]
  
  var_sd_post <- sqrt(w + a * sim_result$X[, (burn_in):(T_total-1)]^2 + b * sim_result$sigma2[, (burn_in):(T_total-1)])
  
  Var_post <- -qstd(alpha,0,var_sd_post, nu = as.numeric(df_val))
  ES_post_mat <- matrix(0, nrow = M, ncol = T)
  for (i in 1:M) {
    for (j in 1:T) {
      sd_val <- var_sd_post[i, j]
      ES_post_mat[i, j] <- -1 / alpha * integrate(
        function(q) qstd(q, 0, sd_val, nu = df_h0),
        lower = 0, upper = alpha
      )$value
    }
  }
  df_label <- ifelse(df_val == "Normal", "Normal", paste0("df", df_val))
  X_list[[paste0("H1_", df_label)]] <- X_post
  Var_list[[paste0("H1_", df_label)]] <- Var_post
  ES_list[[paste0("H1_", df_label)]] <- ES_post_mat
  var_sd_list[[paste0("H1_", df_label)]] <- var_sd_post   # Store var_sd
}

# Assume level = 0.05
level <- 0.05

# Initialize list to store quantile results
quantile_results <- list()

# List of hypothesis names (keys in your X_list, Var_list, ES_list)
hypotheses <- names(X_list)

# Loop through each hypothesisz
# Inside your loop:
for (H in hypotheses) {
  cat("Processing", H, "\n")
  
  # --- NEW: Extract correct df_val
  if (H == "H0") {
    df_val <- df_h0
  } else if (grepl("Normal", H)) {
    df_val <- "Normal"
  } else {
    df_val <- sub("H1_df", "", H)
  }
  
  # Then call Z2 
  
  Z2_val <- quantile(Z2(T, ES_list[[H]], Var_list[[H]], X_list[[H]]), probs = level, type = 1)
  quantile_results[[H]] <- list(
    Z2 = Z2_val
  )
}

# Convert list to data frame
quantile_df <- bind_rows(quantile_results, .id = "Hypothesis")

# Round columns 
quantile_df <- quantile_df %>%
  mutate(across(Z2, round, 4))

# View result
print(quantile_df)

# Plot PDF of Z3 for each hypothesis and print quantile value
par(mfrow = c(2, ceiling(length(hypotheses)/2)))  # layout for multiple plots

for (H in hypotheses) {
  cat("Plotting Z2 PDF for", H, "\n")
  
  # Extract correct df_val from hypothesis name
  if (H == "H0") {
    df_val <- df_h0
  } else if (grepl("Normal", H)) {
    df_val <- "Normal"
  } else {
    df_val <- sub("H1_df", "", H)
  }
  
  # Compute Z2 values
  Z2_vals <- Z2(T, ES_list[[H]], Var_list[[H]], X_list[[H]])
  
  # Compute the quantile
  q_val <- quantile(Z2_vals, probs = level, type = 1)
  
  # Plot the density
  plot(density(Z2_vals), 
       main = paste("Z2 PDF -", H, "\n", level*100, "% quantile =", round(q_val, 4)),
       xlab = "Z2 value",
       ylab = "Density",
       col = "steelblue",
       lwd = 2)
  
  # Add vertical line at quantile
  abline(v = q_val, col = "red", lty = 2)
}




#Forskel i power for t=250 og t=1000.

generate_T2 <- function(T, df, gamma, alpha, ES, M) {
  T_values <- matrix(0, nrow = M, ncol = 2)
  for (j in 1:M) {
    X_t <- gamma * rt(T, df)
    VaR <- -qt(alpha, df)  
    
    I_t <- ifelse(X_t + VaR < 0, 1, 0)
    v1 <- alpha - I_t
    v2 <- VaR- ES + (1 / alpha) * (-VaR - X_t) * I_t
    v <- cbind(v1, v2)
    
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

set.seed(66)
t250 <- generate_T2(250,100,1.1,0.025,2.378,1000)
power_250_gamma <- mean(2*apply(1-pnorm(t250),1,min) <= 0.05)
set.seed(66)
t1000 <- generate_T2(1000,100,1.1,0.025,2.378,1000)
power_1000_gamma <- mean(2*apply(1-pnorm(t1000),1,min) <= 0.05)


c(power_250_gamma,power_1000_gamma)


generate_Z2 <- function(T, df, gamma, alpha,ES, M) {
  Z_values <- numeric(M)
  for (i in 1:M) {
    X_t <- gamma*rt(T, df)
    I_t <- ifelse(X_t -qt(alpha,df) < 0, 1, 0)
    Z_values[i] <- sum(X_t*I_t)/(T*0.025*ES)+1
  }
  return(Z_values)
}

set.seed(66)
critz2_250 <- quantile(generate_Z2(250,100,1,0.025,2.378,1000),0.05,type=1)
set.seed(66)
critz2_1000 <-quantile(generate_Z2(1000,100,1,0.025,2.378,1000),0.05,type=1)
z250 <- generate_Z2(250,100,1.1,0.025,2.378,1000)
power_250_gamma <- mean(z250 <= critz2_250)
set.seed(66)
z1000 <- generate_Z2(1000,100,1.1,0.025,2.378,1000)
power_1000_gamma <- mean(z1000 <= critz2_1000)

c(power_250_gamma,power_1000_gamma)
