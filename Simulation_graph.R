library(ggplot2)
library(tidyr)
library(dplyr)

# Sample sizes
n_vals <- c(50,100,200,500,1000)

# Updated MSE values from latest simulation
mse_alpha <- c(108.2436, 33.7425, 1.4852, 0.5758, 0.1857)
mse_b     <- c(218.0672, 75.2569, 18.9527, 2.4207, 0.5733)
mse_eta   <- c(20.0982, 12.6440, 2.2783, 0.6922, 0.1197)

# Create dataframe
df <- data.frame(
  SampleSize = rep(n_vals, 3),
  Parameter  = factor(rep(c("alpha","b","eta"), each = 5)),
  MSE = c(mse_alpha, mse_b, mse_eta)
)

# Plot
p <- ggplot(df, aes(x = SampleSize, y = MSE,
                    color = Parameter, linetype = Parameter)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  scale_y_log10() +
  labs(x = "Sample Size",
       y = "MSE") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = c(0.8, 0.80),
    legend.text = element_text(size = 10)
  )

p

# Mean estimates
alpha_hat <- c(3.1777, 2.0266, 1.0733, 0.9814, 0.8757)
b_hat     <- c(4.7862, 3.2855, 2.4207, 1.5652, 1.3567)
eta_hat   <- c(3.3920, 3.0005, 2.4765, 2.2779, 2.0772)

# Create dataframe
df_est <- data.frame(
  SampleSize = rep(n_vals, 3),
  Parameter  = factor(rep(c("alpha","b","eta"), each = 5)),
  Estimate   = c(alpha_hat, b_hat, eta_hat)
)

# Plot
p_est <- ggplot(df_est, aes(x = SampleSize, y = Estimate,
                            color = Parameter, linetype = Parameter)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  labs(x = "Sample Size",
       y = "Estimate") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = c(0.8, 0.80),
    legend.text = element_text(size = 10)
  )

p_est