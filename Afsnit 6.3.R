library(quantreg)
library(dplyr)
library(stats)
library(ggplot2)
library(tidyr)
library(gridExtra)

Z2 <- function(T, ES_matrix, Var_matrix, X_matrix) {
  I_t <- ifelse(X_matrix + Var_matrix < 0, 1, 0) 
  Z_values <- rowSums((X_matrix * I_t) / (T * 0.025 * ES_matrix)) + 1  
  return(Z_values)
}
Z3 <- function(T, X_matrix, alpha, df_val, var_sd_post) {
  df_num <- as.numeric(df_val)
  if (is.na(df_num)) stop(paste("Invalid df_val:", df_val))
  
  alpha_floor <- floor(alpha * T)
  
  
 
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
    EV_vec <- rep(NA, T)
    for (t in 1:T) {
      EV_vec[t] <- (-T/alpha_floor)*integrate(function(p) {
        pbeta(1 - p, T - alpha_floor, alpha_floor) * qstd(p,0,sd_row[t], nu = 10)
      }, lower = 0, upper = 1 )$value}
    
    
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
df_values <- c(3, 4, 5, 6, 7, 8, 9, 10,100)
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


VaR1_H0 <- VaR1(X_list[["H0"]], Var_list[["H0"]])
beta <- seq(0.03,0.07,0.001)
quantile(VaR1_H0, beta, type = 1)



# --- First: Calculate power ---
level <- 0.058
power_df <- data.frame()

# Critical values under H0
crit_Z2 <- quantile(Z2(T, ES_list[["H0"]], Var_list[["H0"]], X_list[["H0"]]), probs = level, type = 1)
crit_Z3 <- quantile(Z3(T, X_list[["H0"]],0.025,df_h0,var_sd_list[["H0"]]), probs = level, type = 1)
crit_Z4 <- quantile(Z4(T, ES_list[["H0"]], Var_list[["H0"]], X_list[["H0"]]), probs = level, type = 1)
crit_VaR1 <- quantile(VaR1(X_list[["H0"]], Var_list[["H0"]]), probs = level, type = 1)
crit_T2 <-  level 

# Save order
df_labels <- c(3, 4, 5, 6, 7, 8, 9, 10,100)
df_indices <- 1:length(df_labels)

i <- 1  # make sure to initialize i

for (df_name in names(X_list)) {
  if (df_name == "H0") next
  
  # Extract numeric df from df_name like "H1_df5"
  df_val <- as.numeric(gsub("H1_df", "", df_name))
  
  # Calculate test stats
  Z2_vals <- Z2(T, ES_list[[df_name]], Var_list[[df_name]], X_list[[df_name]])
  Z3_vals <- Z3(T, X_list[[df_name]], alpha, df_val, var_sd_list[[df_name]])
  Z4_vals <- Z4(T, ES_list[[df_name]], Var_list[[df_name]], X_list[[df_name]])
  VaR1_vals <- VaR1(X_list[[df_name]], Var_list[[df_name]])
  T2_vals <-  T2(T, ES_list[[df_name]], Var_list[[df_name]], X_list[[df_name]])
  
  # Save
  power_df <- rbind(
    power_df,
    data.frame(df_index = i, df_label = df_labels[i], Test = "Z2", Power = mean(Z2_vals < crit_Z2)),
    data.frame(df_index = i, df_label = df_labels[i], Test = "Z3", Power = mean(Z3_vals < crit_Z3)),
    data.frame(df_index = i, df_label = df_labels[i], Test = "Z4", Power = mean(Z4_vals < crit_Z4)),
    data.frame(df_index = i, df_label = df_labels[i], Test = "VaR1", Power = mean(VaR1_vals < crit_VaR1)),
    data.frame(df_index = i, df_label = df_labels[i], Test = "T2", Power =  mean(2*apply(1-pnorm(T2_vals),1,min) <= crit_T2))
  )
  
  i <- i + 1
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
  "Z2" = "solid",
  "Z3" = "solid",
  "Z4" = "solid",
  "VaR1" = "solid",
  "T2" = "solid"
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

# Your plot
ggplot(power_df, aes(x = df_index, y = Power, color = Test, group = Test, linetype = Test, shape = Test)) +
  geom_line(position = position_dodge(width = 0.1), size = 0.5) +
  geom_point(position = position_dodge(width = 0.1), size = 2.5) +
  scale_x_continuous(
    breaks = df_indices,
    labels = df_labels
  ) +
  scale_color_manual(values = test_colors) +   # <-- Use your color vector
  scale_shape_manual(values = test_shapes) +   # <-- Use your shape vector
  scale_linetype_manual(values = test_linetypes) +  # assuming you already defined linetypes
  ylim(0, 1) +
  labs(
    title = "",
    x = "df",
    y = "Power",
    color = "Test",         # <-- fixed legend titles
    shape = "Test",
    linetype = "Test"
  ) +
  geom_hline(yintercept = 0.053, linetype = "longdash", color = "gray40") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
power_df

#Lav plot af dag til dag afkast for en "vej" 

# Extract your simulated 250 returns
returns <- X_list[["H0"]][1, 1:250]

# Create a time index for days
days <- 1:250

# Basic plot
plot(days, returns, type = "l",
     xlab = "Dag",
     ylab = "Afkast",
     main = "",
     col = "blue",
     lwd = 1)

# Optional: add a horizontal line at 0 to make positive/negative returns clearer
abline(h = 0, col = "black", lty = 3)






#_________________________________________________
#Undersøg for ændringer i variansen af epsilon 
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
df <- 10 
# Degrees of freedom to test under H1
sd0 <- 1 
sd_values <- c(0.8,1.0,1.1,1.2,1.3,1.4)


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
cat("Simulating H0 (t-distribution with sd=1)\n")

epsilons <- matrix(rstd(M * T_total,0,sd0,nu = df), nrow = M, ncol = T_total)
sim_result <- simulate_garch_paths(M, T_total, epsilons)

# Only keep post-burn-in
X_post <- sim_result$X[, (burn_in + 1):T_total]
sigma2_post <- sim_result$sigma2[, (burn_in + 1):T_total]

# Compute var_sd directly
var_sd_post <- sqrt(w + a * sim_result$X[, (burn_in):(T_total-1)]^2 + b * sim_result$sigma2[, (burn_in):(T_total-1)])

# Compute VaR and ES
Var_post <- -qstd(alpha,0,var_sd_post*sd0,nu = df)
ES_post_mat <- matrix(0, nrow = M, ncol = T)
for (i in 1:M) {
  for (j in 1:T) {
    ES_post_mat[i, j] <- -1/alpha * integrate(
      function(q) qstd(q, 0,    var_sd_post[i, j], nu = df),
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
for (sd_val in sd_values) {
  cat("Simulating H1 with sd =", sd_val, "\n")
  
  if (sd_val == "Normal") {
    epsilons <- matrix(rnorm(M * T_total), nrow = M, ncol = T_total)
  } else {
    epsilons <- matrix(rstd(M * T_total, 0,as.numeric(sd_val), nu = df), nrow = M, ncol = T_total)
  }
  
  sim_result <- simulate_garch_paths(M, T_total, epsilons)
  
  X_post <- sim_result$X[, (burn_in + 1):T_total]
  sigma2_post <- sim_result$sigma2[, (burn_in + 1):T_total]
  
  var_sd_post <- sqrt(w + a * sim_result$X[, (burn_in):(T_total-1)]^2 + b * sim_result$sigma2[, (burn_in):(T_total-1)])
  
    Var_post <- -qstd(alpha,0,var_sd_post, nu = df)
    ES_post_mat <- matrix(0, nrow = M, ncol = T)
    for (i in 1:M) {
      for (j in 1:T) {
        ES_post_mat[i, j] <- -1/alpha * integrate(
          function(q) qstd(q, 0, var_sd_post[i, j], nu = df),
          lower = 0, upper = alpha
        )$value
      }
    }
  
  
  sd_label <- gsub("\\.0$", "", paste0("sd", sprintf("%.1f", sd_val)))
  X_list[[paste0("H1_", sd_label)]] <- X_post
  Var_list[[paste0("H1_",  sd_label)]] <- Var_post
  ES_list[[paste0("H1_",  sd_label)]] <- ES_post_mat
  var_sd_list[[paste0("H1_",  sd_label)]] <- var_sd_post   # Store var_sd
}





VaR1_H0 <- VaR1(X_list[["H0"]], Var_list[["H0"]])
beta <- seq(0.03,0.07,0.001)
quantile(VaR1_H0, beta, type = 1)



# --- First: Calculate power ---
level <- 0.058
power_sd <- data.frame()

# Critical values under H0
crit_Z2 <- quantile(Z2(T, ES_list[["H0"]], Var_list[["H0"]], X_list[["H0"]]), probs = level, type = 1)
crit_Z3 <- quantile(Z3(T, X_list[["H0"]],alpha,df,var_sd_list[["H0"]]), probs = level, type = 1)
crit_Z4 <- quantile(Z4(T, ES_list[["H0"]], Var_list[["H0"]], X_list[["H0"]]), probs = level, type = 1)
crit_VaR1 <- quantile(VaR1(X_list[["H0"]], Var_list[["H0"]]), probs = level, type = 1)
crit_T2 <-  level #Obs brug qnorm(0.025), alternativt ---


# Save order
sd_labels <- c(0.8,1.0,1.1,1.2,1.3,1.4)
sd_indices <- 1:length(sd_labels)

i <- 1  # make sure to initialize i

for (sd_name in names(X_list)) {
  if (sd_name == "H0") next
  
  # Extract numeric sd from sd_name like "H1_sd5"
  sd_val <- as.numeric(gsub("H1_sd", "", sd_name))
  
  # Calculate test stats
  Z2_vals <- Z2(T, ES_list[[sd_name]], Var_list[[sd_name]], X_list[[sd_name]])
  Z3_vals <- Z3(T, X_list[[sd_name]], alpha, df, var_sd_list[[sd_name]])
  Z4_vals <- Z4(T, ES_list[[sd_name]], Var_list[[sd_name]], X_list[[sd_name]])
  VaR1_vals <- VaR1(X_list[[sd_name]], Var_list[[sd_name]])
  T2_vals <- T2(T, ES_list[[sd_name]], Var_list[[sd_name]], X_list[[sd_name]])
  
  # Save
  power_sd <- rbind(
    power_sd,
    data.frame(sd_index = i, sd_label = sd_labels[i], Test = "Z2", Power = mean(Z2_vals < crit_Z2)),
    data.frame(sd_index = i, sd_label = sd_labels[i], Test = "Z3", Power = mean(Z3_vals < crit_Z3)),
    data.frame(sd_index = i, sd_label = sd_labels[i], Test = "Z4", Power = mean(Z4_vals < crit_Z4)),
    data.frame(sd_index = i, sd_label = sd_labels[i], Test = "VaR1", Power = mean(VaR1_vals < crit_VaR1)),
    data.frame(sd_index = i, sd_label = sd_labels[i], Test = "T2", Power =  mean(2*apply(1-pnorm(T2_vals),1,min) <= crit_T2))
  )
  
  i <- i + 1
}


# --- Now Plot ---
ggplot(power_sd, aes(x = sd_index, y = Power, color = Test, group = Test, linetype = Test, shape = Test)) +
  geom_line(position = position_dodge(width = 0.1), size = 0.5) +
  geom_point(position = position_dodge(width = 0.1), size = 2.5) +
  scale_x_continuous(
    breaks = sd_indices,
    labels = sd_labels
  ) +
  scale_color_manual(values = test_colors) +   # <-- Use your color vector
  scale_shape_manual(values = test_shapes) +   # <-- Use your shape vector
  scale_linetype_manual(values = test_linetypes) +  # assuming you already defined linetypes
  ylim(0, 1) +
  labs(
    title = "",
    x = "sigma",
    y = "Power",
    color = "Test",         # <-- fixed legend titles
    shape = "Test",
    linetype = "Test"
  ) +
  geom_hline(yintercept = level, linetype = "longdash", color = "gray40") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )








#_________________________________________________
#Undersøg for ændringer i den langsigtede varians 
# Parameters
T <- 250
M <- 1000
burn_in <- 250
T_total <- T + burn_in
alpha <- 0.025

# GARCH(1,1) parameters
a <- 0.1
b <- 0.85
df <- 10 
sd <- 1
# Degrees of freedom to test under H1

wh0 <- 0.03 
w_values <- c(0.01,0.03,0.04,0.05,0.06,0.1)



# Storage lists
X_list <- list()
Var_list <- list()
ES_list <- list()
var_sd_list <- list()  

set.seed(123)

# Function to simulate GARCH paths
simulate_garch_paths <- function(M, T_total, epsilons,w) {
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
cat("Simulating H0 (t-distribution with sd=1)\n")

epsilons <- matrix(rstd(M * T_total,0,sd,nu = df), nrow = M, ncol = T_total)
sim_result <- simulate_garch_paths(M, T_total, epsilons,wh0)

# Only keep post-burn-in
X_post <- sim_result$X[, (burn_in + 1):T_total]
sigma2_post <- sim_result$sigma2[, (burn_in + 1):T_total]

# Compute var_sd directly
X_hist <- sim_result$X[, (burn_in):(T_total - 1)]
M_rows <- nrow(X_hist)
T_cols <- ncol(X_hist)

var_sd_post <- matrix(0, nrow = M_rows, ncol = T_cols)
var_sd_post[, 1] <- sqrt(wh0 + a * X_hist[, 1]^2 + b * wh0)  # Initial sigma^2



for (t in 2:T_cols) {
  var_sd_post[, t] <- sqrt(wh0 + a * X_hist[, t]^2 + b * var_sd_post[, t - 1]^2)
}

# Compute VaR and ES
Var_post <- -qstd(alpha,0,var_sd_post,nu = df)
ES_post_mat <- matrix(0, nrow = M, ncol = T)
for (i in 1:M) {
  for (j in 1:T) {
    ES_post_mat[i, j] <- -1/alpha * integrate(
      function(q) qstd(q, 0,    var_sd_post[i, j], nu = df),
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
set.seed(123)
for (w_val in w_values) {
  cat("Simulating H1 with w =", w_val, "\n")
  
  epsilons <- matrix(rstd(M * T_total, 0, sd, nu = df), nrow = M, ncol = T_total)
  
  
  # Use current w_val for GARCH simulation
  sim_result <- simulate_garch_paths(M, T_total, epsilons, w_val)
  
  X_post <- sim_result$X[, (burn_in + 1):T_total]
  sigma2_post <- sim_result$sigma2[, (burn_in + 1):T_total]
  
  # Compute var_sd_post using wh0, not w_val
  X_hist <- sim_result$X[, (burn_in):(T_total - 1)]
  M_rows <- nrow(X_hist)
  T_cols <- ncol(X_hist)
  
  var_sd_post <- matrix(0, nrow = M_rows, ncol = T_cols)
  var_sd_post[, 1] <- sqrt(wh0 + a * X_hist[, 1]^2 + b * wh0)  # initialize with sigma2_0 = wh0
  
  for (t in 2:T_cols) {
    var_sd_post[, t] <- sqrt(wh0 + a * X_hist[, t]^2 + b * var_sd_post[, t - 1]^2)
  }
  
  Var_post <- -qstd(alpha, 0, var_sd_post, nu = df)
  ES_post_mat <- matrix(0, nrow = M, ncol = T)
  for (i in 1:M) {
    for (j in 1:T) {
      ES_post_mat[i, j] <- -1 / alpha * integrate(
        function(q) qstd(q, 0, var_sd_post[i, j], nu = df),
        lower = 0, upper = alpha
      )$value
    }
  }
  
  w_label <- paste0("w_", gsub("\\.", "_", formatC(w_val, format = "f", digits = 2)))
  X_list[[paste0("H1_", w_label)]] <- X_post
  Var_list[[paste0("H1_", w_label)]] <- Var_post
  ES_list[[paste0("H1_", w_label)]] <- ES_post_mat
  var_sd_list[[paste0("H1_", w_label)]] <- var_sd_post  # Store var_sd using wh0
}





VaR1_H0 <- VaR1(X_list[["H0"]], Var_list[["H0"]])
beta <- seq(0.03,0.07,0.001)
quantile(VaR1_H0, beta, type = 1)



# --- First: Calculate power ---
level <- 0.051
power_w <- data.frame()

# Critical values under H0
crit_Z2 <- quantile(Z2(T, ES_list[["H0"]], Var_list[["H0"]], X_list[["H0"]]), probs = level, type = 1)
crit_Z3 <- quantile(Z3(T, X_list[["H0"]],alpha,df,var_sd_list[["H0"]]), probs = level, type = 1)
crit_Z4 <- quantile(Z4(T, ES_list[["H0"]], Var_list[["H0"]], X_list[["H0"]]), probs = level, type = 1)
crit_VaR1 <- quantile(VaR1(X_list[["H0"]], Var_list[["H0"]]), probs = level, type = 1)
crit_T2 <-  level #Obs brug qnorm(0.025), alternativt ---




# Save order
w_labels <- c(0.01,0.03,0.04,0.05,0.06,0.1)
w_indices <- 1:length(w_labels)

i <- 1  

for (w_name in names(X_list)) {
  if (w_name == "H0") next
  
  # Extract numeric sd from sd_name like "H1_sd5"
  w_val <- as.numeric(gsub("_", ".", gsub("H1_w", "", w_name)))  # NEW
  
  # Calculate test stats
  Z2_vals <- Z2(T, ES_list[[w_name]], Var_list[[w_name]], X_list[[w_name]])
  Z3_vals <- Z3(T, X_list[[w_name]], alpha, df, var_sd_list[[w_name]])
  Z4_vals <- Z4(T, ES_list[[w_name]], Var_list[[w_name]], X_list[[w_name]])
  VaR1_vals <- VaR1(X_list[[w_name]], Var_list[[w_name]])
  T2_vals <- T2(T, ES_list[[w_name]], Var_list[[w_name]], X_list[[w_name]])
  
  # Save
  power_w <- rbind(
    power_w,
    data.frame(w_index = i, w_label = w_labels[i], Test = "Z2", Power = mean(Z2_vals < crit_Z2)),
    data.frame(w_index = i, w_label = w_labels[i], Test = "Z3", Power = mean(Z3_vals < crit_Z3)),
    data.frame(w_index = i, w_label = w_labels[i], Test = "Z4", Power = mean(Z4_vals < crit_Z4)),
    data.frame(w_index = i, w_label = w_labels[i], Test = "VaR1", Power = mean(VaR1_vals < crit_VaR1)),
    data.frame(w_index = i, w_label = w_labels[i], Test = "T2", Power =  mean(2*apply(1-pnorm(T2_vals),1,min) <= crit_T2))
  )
  
  i <- i + 1
}



# --- Now Plot ---
ggplot(power_w, aes(x = w_index, y = Power, color = Test, group = Test, linetype = Test, shape = Test)) +
  geom_line(position = position_dodge(width = 0.1), size = 0.5) +
  geom_point(position = position_dodge(width = 0.1), size = 2.5) +
  scale_x_continuous(
    breaks = w_indices,
    labels = w_labels
  ) +
  scale_color_manual(values = test_colors) +   # <-- Use your color vector
  scale_shape_manual(values = test_shapes) +   # <-- Use your shape vector
  scale_linetype_manual(values = test_linetypes) +  # assuming you already defined linetypes
  ylim(0, 1) +
  labs(
    title = "",
    x = "w",
    y = "Power",
    color = "Test",         # <-- fixed legend titles
    shape = "Test",
    linetype = "Test"
  ) +
  geom_hline(yintercept = level, linetype = "longdash", color = "gray40") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )









#_________________________________________________
#Undersøger for ændringer i a og sætter b=0.95-a
# Parameters
T <- 250
M <- 1000
burn_in <- 250
T_total <- T + burn_in
alpha <- 0.025

# GARCH(1,1) parameters
b <- 0.95 
df <- 10 
sd <- 1
w <- 0.01 
# Degrees of freedom to test under H1

a0 <- 0.2
a_values <- c(0,0.05,0.1,0.15,0.2,0.3)




# Storage lists
X_list <- list()
Var_list <- list()
ES_list <- list()
var_sd_list <- list()  # NEW: store var_sd_post here

set.seed(123)

# Function to simulate GARCH paths
simulate_garch_paths <- function(M, T_total, epsilons,a) {
  X <- matrix(0, nrow = M, ncol = T_total)
  sigma2 <- matrix(0, nrow = M, ncol = T_total)
  
  # Initial conditions
  sigma2[,1] <- 0.05
  X[,1] <- 0
  
  for (t in 2:T_total) {
    sigma2[,t] <- w + a*X[,t-1]^2 + (b-a)*sigma2[,t-1]
    X[,t] <- sqrt(sigma2[,t]) * epsilons[,t]
  }
  
  list(X = X, sigma2 = sigma2)
}

### ----- H0 Simulation ----- ###
cat("Simulating H0 (t-distribution with sd=1)\n")

epsilons <- matrix(rstd(M * T_total,0,sd,nu = df), nrow = M, ncol = T_total)
sim_result <- simulate_garch_paths(M, T_total, epsilons,a0)

# Only keep post-burn-in
X_post <- sim_result$X[, (burn_in + 1):T_total]
sigma2_post <- sim_result$sigma2[, (burn_in + 1):T_total]

# Compute var_sd directly
X_hist <- sim_result$X[, (burn_in):(T_total - 1)]
M_rows <- nrow(X_hist)
T_cols <- ncol(X_hist)

var_sd_post <- matrix(0, nrow = M_rows, ncol = T_cols)
var_sd_post[, 1] <- sqrt(w + (a0)*X_hist[, 1]^2 + (b-a0)*w)  # Initial sigma^2

for (t in 2:T_cols) {
  var_sd_post[, t] <- sqrt(w + (a0)*X_hist[, t]^2 + (b-a0)*var_sd_post[, t - 1]^2)
}

# Compute VaR and ES
Var_post <- -qstd(alpha,0,var_sd_post,nu = df)
ES_post_mat <- matrix(0, nrow = M, ncol = T)
for (i in 1:M) {
  for (j in 1:T) {
    ES_post_mat[i, j] <- -1/alpha * integrate(
      function(q) qstd(q, 0,    var_sd_post[i, j], nu = df),
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
set.seed(123)
for (a_val in a_values) {
  cat("Simulating H1 with a =", a_val, "\n")
  
  epsilons <- matrix(rstd(M * T_total, 0, sd, nu = df), nrow = M, ncol = T_total)
  
  
  # Use current w_val for GARCH simulation
  sim_result <- simulate_garch_paths(M, T_total, epsilons, a_val)
  
  X_post <- sim_result$X[, (burn_in + 1):T_total]
  sigma2_post <- sim_result$sigma2[, (burn_in + 1):T_total]
  
  # Compute var_sd_post using wh0, not w_val
  X_hist <- sim_result$X[, (burn_in):(T_total - 1)]
  M_rows <- nrow(X_hist)
  T_cols <- ncol(X_hist)
  
  var_sd_post <- matrix(0, nrow = M_rows, ncol = T_cols)
  var_sd_post[, 1] <- sqrt(w + (a0) * X_hist[, 1]^2 + (b-a0)*w)  # initialize with sigma2_0 = wh0
  
  for (t in 2:T_cols) {
    var_sd_post[, t] <- sqrt(w + (a0) * X_hist[, t]^2 + (b-a0)* var_sd_post[, t - 1]^2)
  }
  
  Var_post <- -qstd(alpha, 0, var_sd_post, nu = df)
  ES_post_mat <- matrix(0, nrow = M, ncol = T)
  for (i in 1:M) {
    for (j in 1:T) {
      ES_post_mat[i, j] <- -1 / alpha * integrate(
        function(q) qstd(q, 0, var_sd_post[i, j], nu = df),
        lower = 0, upper = alpha
      )$value
    }
  }
  
  a_label <- paste0("a_", gsub("\\.", "_", formatC(a_val, format = "f", digits = 2)))
  X_list[[paste0("H1_", a_label)]] <- X_post
  Var_list[[paste0("H1_", a_label)]] <- Var_post
  ES_list[[paste0("H1_", a_label)]] <- ES_post_mat
  var_sd_list[[paste0("H1_", a_label)]] <- var_sd_post  # Store var_sd using wh0
}





VaR1_H0 <- VaR1(X_list[["H0"]], Var_list[["H0"]])
beta <- seq(0.03,0.07,0.001)
quantile(VaR1_H0, beta, type = 1)



# --- First: Calculate power ---
level <- 0.04
power_a <- data.frame()

# Critical values under H0
crit_Z2 <- quantile(Z2(T, ES_list[["H0"]], Var_list[["H0"]], X_list[["H0"]]), probs = level, type = 1)
crit_Z3 <- quantile(Z3(T, X_list[["H0"]],alpha,df,var_sd_list[["H0"]]), probs = level, type = 1)
crit_Z4 <- quantile(Z4(T, ES_list[["H0"]], Var_list[["H0"]], X_list[["H0"]]), probs = level, type = 1)
crit_VaR1 <- quantile(VaR1(X_list[["H0"]], Var_list[["H0"]]), probs = level, type = 1)
crit_T2 <-  level #Obs brug qnorm(0.025), alternativt ---




# Save order
a_labels <- c(0,0.05,0.1,0.15,0.2,0.3)
a_indices <- 1:length(a_labels)

i <- 1  # make sure to initialize i

for (a_name in names(X_list)) {
  if (a_name == "H0") next
  
  # Extract numeric sd from sd_name like "H1_sd5"
  a_val <- as.numeric(gsub("_", ".", sub("H1_a_", "", a_name)))
  
  # Calculate test stats
  Z2_vals <- Z2(T, ES_list[[a_name]], Var_list[[a_name]], X_list[[a_name]])
  Z3_vals <- Z3(T, X_list[[a_name]], alpha, df, var_sd_list[[a_name]])
  Z4_vals <- Z4(T, ES_list[[a_name]], Var_list[[a_name]], X_list[[a_name]])
  VaR1_vals <- VaR1(X_list[[a_name]], Var_list[[a_name]])
  T2_vals <- T2(T, ES_list[[a_name]], Var_list[[a_name]], X_list[[a_name]])
  
  # Save
  power_a <- rbind(
    power_a,
    data.frame(a_index = i, a_label = a_labels[i], Test = "Z2", Power = mean(Z2_vals < crit_Z2)),
    data.frame(a_index = i, a_label = a_labels[i], Test = "Z3", Power = mean(Z3_vals < crit_Z3)),
    data.frame(a_index = i, a_label = a_labels[i], Test = "Z4", Power = mean(Z4_vals < crit_Z4)),
    data.frame(a_index = i, a_label = a_labels[i], Test = "VaR1", Power = mean(VaR1_vals < crit_VaR1)),
    data.frame(a_index = i, a_label = a_labels[i], Test = "T2", Power = mean(2*apply(1-pnorm(T2_vals),1,min) <= crit_T2))
  )
  
  i <- i + 1
}


# --- Now Plot --- 
ggplot(power_a, aes(x = a_index, y = Power, color = Test, group = Test, linetype = Test, shape = Test)) +
  geom_line(position = position_dodge(width = 0.1), size = 0.5) +
  geom_point(position = position_dodge(width = 0.1), size = 2.5) +
  scale_x_continuous(
    breaks = a_indices,
    labels = a_labels
  ) +
  scale_color_manual(values = test_colors) +   # <-- Use your color vector
  scale_shape_manual(values = test_shapes) +   # <-- Use your shape vector
  scale_linetype_manual(values = test_linetypes) +  # assuming you already defined linetypes
  ylim(0, 1) +
  labs(
    title = "",
    x = "alpha",
    y = "Power",
    color = "Test",     
    shape = "Test",
    linetype = "Test"
  ) +
  geom_hline(yintercept = level, linetype = "longdash", color = "gray40") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )








