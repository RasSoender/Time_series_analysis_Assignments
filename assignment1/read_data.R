### Read training data
#! Perhaps you need to set the working directory!?
setwd("/Users/ceciliehvilsted/Desktop/DTU-MSc/Sem2/TimeSeries/Assignment1")
D <- read.csv("DST_BIL54.csv")
str(D)

# See the help
?strftime
D$time <- as.POSIXct(paste0(D$time,"-01"), "%Y-%m-%d", tz="UTC")
D$time
class(D$time)

## Year to month for each of them
D$year <- 1900 + as.POSIXlt(D$time)$year + as.POSIXlt(D$time)$mon / 12

## Make the output variable a floating point (i.e.\ decimal number)
D$total <- as.numeric(D$total) / 1E6

## Divide intro train and test set
teststart <- as.POSIXct("2024-01-01", tz="UTC")
Dtrain <- D[D$time < teststart, ]
Dtest <- D[D$time >= teststart, ]

### 1. Plot data ###

#install.packages("ggplot2")
library(ggplot2)

# 1.1

Dtrain$x <- seq(2018, by = 1/12, length.out = nrow(Dtrain))

# Plot the training data
ggplot(Dtrain, aes(x = year, y = total)) +
  geom_line(color = "blue") +
  labs(
    title = "Number of registered vehicles in Denmark",
    x = "Year",
    y = "Total vehicles (in millions)"
  ) +
  theme_minimal()

# 1.2

# The number of registered vehicles shows a general upward trend over time, meaning more vehicles are being registered each month. This can be due to population growth and economic expansion.
# The rate of increase is not uniform, as there are be periods of faster or slower growth. For instance it stabilizes around year 2020, but then greatly increases in the second half of 2020 and until 2022.

### 2. Linear trend model ###

Dtrain$total[1:3]


### 3. OLS - global linear trend model ###

# 3.1

# Formula: theta_hat = (X^T @ X)^(-1) @ X^T @ Y
# This minimizes the sum of squared errors
# Design Matrix X: first column is all ones (for intercept), second column is x_t (time variable)
model <- lm(total ~ x, data = Dtrain)
summary(model)

# 3.2

# Extract parameter estimates
theta_hat <- coef(model)

# Extract standard errors
se_theta_hat <- summary(model)$coefficients[, "Std. Error"]

# Print results
cat("Estimated parameters:\n")
cat("θ1̂ (Intercept) =", round(theta_hat[1], 3), " (SE =", round(se_theta_hat[1], 3), ")\n")
cat("θ2̂ (Slope) =", round(theta_hat[2], 3), " (SE =", round(se_theta_hat[2], 3), ")\n")

# Small SE means more confidence in estimates.
# Large SE means estimates are more uncertain.

# Create the plot: estimated mean as a line with the observations as points
ggplot(Dtrain, aes(x = x, y = total)) +
  geom_point(color = "blue") +  # Observed data as points
  geom_abline(intercept = theta_hat[1], slope = theta_hat[2], color = "red", linewidth = 1.2) +  # Fitted trend line
  labs(
    title = "Linear trend model for vehicle registrations",
    x = "Time (years)",
    y = "Total vehicles (in millions)"
  ) +
  theme_minimal()

# 3.3

# Create new data frame for 2024
future_x <- seq(max(Dtrain$x) + 1/12, by = 1/12, length.out = 12)
future_data <- data.frame(x = future_x)

# Predict values with 95% prediction intervals
predictions <- predict(model, newdata = future_data, interval = "prediction", level = 0.95)

# Combine predictions into a table
forecast_table <- data.frame(
  Month = format(seq(as.Date("2024-01-01"), by = "month", length.out = 12), "%Y-%b"),
  Predicted_Value = round(predictions[, "fit"], 3),
  Lower_Bound = round(predictions[, "lwr"], 3),
  Upper_Bound = round(predictions[, "upr"], 3)
)

print(forecast_table)


# 3.4

# Blue points: Actual training data.
# Red line: Fitted linear model for the training set.
# Black line: Forecasted values for 2024.
# Shaded black region: 95% prediction interval, showing the uncertainty in predictions.

# Create future time points for prediction (2024)
future_x <- seq(max(Dtrain$x) + 1/12, by = 1/12, length.out = 12)
future_data <- data.frame(x = future_x)

# Predict values with 95% prediction intervals
predictions <- predict(model, newdata = future_data, interval = "prediction", level = 0.95)

# Combine forecast data with predictions
future_data$Predicted <- predictions[, "fit"]
future_data$Lower_Bound <- predictions[, "lwr"]
future_data$Upper_Bound <- predictions[, "upr"]

# Plot the training data, fitted trend, and forecast
ggplot() +
  # Training data points
  geom_point(data = Dtrain, aes(x = x, y = total), color = "blue", size = 2) +
  
  # Fitted linear model line
  geom_abline(intercept = coef(model)[1], slope = coef(model)[2], color = "red", linewidth = 1.2) +
  
  # Forecasted trend line
  geom_line(data = future_data, aes(x = x, y = Predicted), color = "black", linewidth = 1.2) +
  
  # Prediction interval shading
  geom_ribbon(data = future_data, aes(x = x, ymin = Lower_Bound, ymax = Upper_Bound), fill = "black", alpha = 0.2) +
  
  # Labels and theme
  labs(
    title = "Vehicle registration trend and forecast for 2024",
    x = "Time (years)",
    y = "Total vehicles (in millions)"
  ) +
  theme_minimal()


# 3.5 


# 3.6

# Compute residuals
Dtrain$residuals <- residuals(model)

ggplot(Dtrain, aes(x = fitted(model), y = residuals)) +
  geom_point(color = "blue") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Residuals vs. fitted values",
    x = "Fitted values",
    y = "Residuals"
  ) +
  theme_minimal()

# If the points are randomly scattered, homoscedasticity is satisfied. They are not since there is a wave trend.
# If a pattern exists (e.g., funnel shape), the variance is not constant. The spread of residuals appears somewhat uneven, which violates the assumption.

# Histogram of residuals
ggplot(Dtrain, aes(x = residuals)) +
  geom_histogram(bins = 20, fill = "blue", alpha = 0.5, color = "black") +
  labs(title = "Histogram of residuals", x = "Residuals", y = "Frequency") +
  theme_minimal()

# The histogram should resemble a normal distribution. It doesnt really.

# Q-Q plot (Normality check)
qqnorm(Dtrain$residuals)
qqline(Dtrain$residuals, col = "red")

# The Q-Q plot should have points roughly along the red line, which is not really the case.

