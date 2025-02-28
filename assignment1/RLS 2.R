#setwd("/Users/MortenHelsoe/Desktop/Skole & Arbejde/MSc. HCAI/1. Semester/Timeseries/Assignment1")
setwd("/Users/ceciliehvilsted/Desktop/DTU-MSc/Sem2/TimeSeries/Assignment1")

#install.packages("gridExtra")
library(ggplot2)
library(reshape2)
library(gridExtra)

# Load data
D <- read.csv("DST_BIL54.csv")

# Convert time to POSIXct format
D$time <- as.POSIXct(paste0(D$time, "-01"), format="%Y-%m-%d", tz="UTC")
D$year <- as.numeric(format(D$time, "%Y")) + as.numeric(format(D$time, "%m"))/12
D$total <- as.numeric(D$total) / 1E6  # Convert to floating point

# Train-test split
teststart <- as.POSIXct("2024-01-01", tz="UTC")
Dtrain <- D[D$time < teststart, ]
Dtest <- D[D$time >= teststart, ]

# Prepare design matrices
X <- cbind(1, Dtest$year)
Xtest <- cbind(1, Dtest$year)

y <- matrix(Dtest$total, ncol=1)
n <- nrow(X)  # Number of observations
p <- 2  # Number of parameters

# -----------------------------------------------------------------------------------

# sanity check:
# using RLS notation for computing OLS with all 26 observations:
R26 <- t(X)%*%X
print(R26)

h26 <- t(X)%*%y
print(h26)

RLS26 <- solve(R26)%*%h26
print(RLS26)
# this corresponds exactly to the OLS estimates (as it should)

# (back to slides)


###########################
# Now iterate through data:
# For each step:
# - calculate R and theta
# - calculate one-step prediction

lambda <- 0.9

# we will use the entire dataset (not only the training data from before):

n <- length(X[,1])

# initialise containers for parameter estimates (Theta) and one step predictions:
Theta <- matrix(NA, nrow=n, ncol=p)
OneStepPred <- matrix(NA, nrow=n)

# 1 # very first step:
x1 <- X[1,]
R_1 <- x1%*%t(x1) # R is a pxp matrix
h_1 <- x1*y[1]    # h is a px1 vector (but R prints it in a row..)

# to estimate theta we need to invert R:
#solve(R_1)
# in this very first step R cannot be inverted - too soon to estimate parameters!
# (we cannot estimate p parameters drom only one datapoint)


# 2 # second step - first time to estimate parameters and make prediction
x2 <- X[2,]
R_2 <- lambda*R_1 + x2 %*% t(x2)
h_2 <- lambda*h_1 + x2 * y[2]

print(kappa(R_2))

# Add a small regularization term to the diagonal of R_2
epsilon <- 1e-8
R_2_reg <- R_2 + epsilon * diag(nrow(R_2))

# Now try to invert the regularized matrix
solve(R_2_reg)
#solve(R_2)
# R is now invertible (we can estimate p parameters from p observations)

# we estimate theta (for the first time - so not yet using "update" formula):
R_2 <- R_2_reg
Theta[2,] <- solve(R_2) %*% h_2

# we predict one step ahead:
OneStepPred[2+1] <- X[2+1,]%*%Theta[2,]

# 3 # third step - first time to use update formula
x3 <- X[3,]
R_3 <- lambda*R_2 + x3 %*% t(x3)
Theta[3,] <- Theta[2,] + solve(R_3) %*% x3 %*% (y[3] - t(x3) %*% Theta[2,])

# we predict one step ahead:
OneStepPred[3+1] <- X[3+1,]%*%Theta[3,]

# next many steps # - update and predict

R <- R_3

for(i in 4:n){
  x <- X[i, ]
  # Update
  R <- lambda*R + x %*% t(x)
  Theta[i, ] <- Theta[i-1, ] + solve(R) %*% x %*% (y[i] - t(x) %*% Theta[i-1, ])
}

# predict
for(i in 4:n-1){
  OneStepPred[i+1] <- X[i+1, ]%*%Theta[i, ]
}

# Plot estimate of intercept:
plot(Theta[,1])

# Plot estimate of slope:
plot(Theta[,2])

# Printing last theta estimates
print(Theta[n,])

Theta_df <- data.frame(
  Time = 1:n,
  Intercept = Theta[,1],
  Slope = Theta[,2]
)

# Plot one step predictions:
ggplot(Dtrain, aes(x=year, y=total)) +
  geom_point(col="black") +
  geom_point(aes(y=OneStepPred), col="blue", size=1) + 
  geom_line(aes(y=OneStepPred), col="blue")

# Horizon

# Function to calculate k-step RMSE for a given lambda
calculate_k_step_rmse <- function(lambda, k) {
  # Initialize containers for parameter estimates (Theta) and k-step predictions
  Theta <- matrix(NA, nrow=n, ncol=p)
  kStepPred <- matrix(NA, nrow=n, ncol=k)
  
  # 1. Very first step
  x1 <- X[1,]
  R_1 <- x1 %*% t(x1)  # R is a pxp matrix
  h_1 <- x1 * y[1]     # h is a px1 vector (but R prints it in a row..)
  
  # 2. Second step - first time to estimate parameters and make prediction
  x2 <- X[2,]
  R_2 <- lambda * R_1 + x2 %*% t(x2)
  h_2 <- lambda * h_1 + x2 * y[2]
  
  # Add a small regularization term to the diagonal of R_2
  epsilon <- 1e-8
  R_2_reg <- R_2 + epsilon * diag(nrow(R_2))
  
  # Now try to invert the regularized matrix
  Theta[2,] <- solve(R_2_reg) %*% h_2
  
  # Predict k steps ahead
  for (j in 1:k) {
    if (2 + j <= n) {
      kStepPred[2 + j, j] <- X[2 + j,] %*% Theta[2,]
    }
  }
  
  # 3. Third step - first time to use update formula
  x3 <- X[3,]
  R_3 <- lambda * R_2 + x3 %*% t(x3)
  Theta[3,] <- Theta[2,] + solve(R_3 + epsilon * diag(nrow(R_3))) %*% x3 %*% (y[3] - t(x3) %*% Theta[2,])
  
  # Predict k steps ahead
  for (j in 1:k) {
    if (3 + j <= n) {
      kStepPred[3 + j, j] <- X[3 + j,] %*% Theta[3,]
    }
  }
  
  # Next many steps - update and predict
  R <- R_3
  
  for (i in 4:n) {
    x <- X[i, ]
    # Update
    R <- lambda * R + x %*% t(x)
    Theta[i, ] <- Theta[i-1, ] + solve(R + epsilon * diag(nrow(R))) %*% x %*% (y[i] - t(x) %*% Theta[i-1, ])
    
    # Predict k steps ahead
    for (j in 1:k) {
      if (i + j <= n) {
        kStepPred[i + j, j] <- X[i + j,] %*% Theta[i,]
      }
    }
  }
  
  # Calculate k-step residuals and RMSE
  kStepResiduals <- matrix(NA, nrow=n, ncol=k)
  for (j in 1:k) {
    kStepResiduals[, j] <- y - kStepPred[, j]
  }
  kStepRMSE <- sqrt(colMeans(kStepResiduals^2, na.rm=TRUE))

  if (length(kStepRMSE) < 12) {
    kStepRMSE <- c(kStepRMSE, rep(NA, 12 - length(kStepRMSE)))
  }
  
  return(kStepRMSE)
}

# Sequence of lambda values
lambda_seq <- seq(0.5, 0.99, by=0.01)

rmse_results <- sapply(lambda_seq, function(lambda) {
  calculate_k_step_rmse(lambda, 12)  # Ensure a fixed length return
})

print(str(rmse_results))
print(dim(rmse_results))

# Convert results to a data frame for plotting
rmse_df <- as.data.frame(t(rmse_results))
colnames(rmse_df) <- paste0("k=", 1:12)
rmse_df$lambda <- lambda_seq


# Function to plot RMSE for each k
plotting_rmse_k <- function(rmse_matrix, lambda_seq) {
  plots <- list()
  
  for (k in 1:12) {
    rmse_values <- rmse_matrix[k, ]  # Correct indexing (k-th row, all lambda values)
    
    # Filter out NA values
    valid_indices <- !is.na(rmse_values)
    lambda_numeric <- lambda_seq[valid_indices]
    rmse_values <- rmse_values[valid_indices]
    
    # Create a data frame for plotting
    plot_df <- data.frame(lambda = lambda_numeric, RMSE = rmse_values)
    
    # Create the plot
    p <- ggplot(plot_df, aes(x=lambda, y=RMSE)) +
      geom_line(color="blue") +  # Line plot to show trends
      geom_point(color="red") +  # Add points for clarity
      labs(title=sprintf("RMSE for k = %d", k), x="Lambda", y="RMSE") +
      theme_minimal()
    
    plots[[k]] <- p
  }
  
  # Arrange plots in a grid
  do.call("grid.arrange", c(plots, ncol=3))
}

# Run the function
plotting_rmse_k(rmse_results, lambda_seq)

# Print minimum RMSE values and lambda for each horizon k
print_min_rmse_lambda <- function(rmse_df) {
  for (k in 1:12) {
    rmse_values <- as.numeric(unlist(rmse_df[, k]))
    min_rmse <- min(rmse_values, na.rm = TRUE)
    min_lambda <- rmse_df$lambda[which.min(rmse_values)]
    cat(sprintf("Horizon k = %d: Minimum RMSE = %.4f at lambda = %.2f\n", k, min_rmse, min_lambda))
  }
}

print_min_rmse_lambda(rmse_df)