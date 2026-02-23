library(ggplot2)
library(dplyr)
library(gridExtra)

# -----------------------------
# PG CDF (4 parameters)
# -----------------------------
pg_cdf <- function(x, alpha, b, eta, xm) {
  ifelse(
    x < 0, 0,
    1 - (xm / (xm + eta * (exp(b * x) - 1)))^alpha
  )
}

# -----------------------------
# PG PDF (4 parameters)
# -----------------------------
pg_pdf <- function(x, alpha, b, eta, xm) {
  ifelse(
    x < 0, 0,
    alpha * b * eta * xm^alpha * exp(b * x) *
      (xm + eta * (exp(b * x) - 1))^(-alpha - 1)
  )
}
xm <- 1   # fixed

'param_configs <- data.frame(
  alpha = c(1.0, 1.5, 0.5, 2.0, 1.0, 3.0, 1.5,1.2,1.6,1.0,1.4,0.4,2.0, 2.0, 2.0, 2.0, 2.0),
  b     = c(1.0, 1.5, 1.0, 0.5, 2.0, 1.0, 0.5,1.0,0.9,1.4,1.1,2,0.7, 0.8, 0.9, 1.0, 1.1),
  eta   = c(1.0, 1.5, 2.0, 1.0, 0.5, 1.5, 2.0,1.0,1.1,0.8,1.4,1,0.02, 0.03, 0.04, 0.05, .06)
)'

param_configs <- data.frame(
  x_m =c(1,1,1,1,1,1),
  alpha = c(1.0, 1, 1.5, 2 ,2.5, 2),
  b     = c(.9, 1, 1.1, .8, .5, 1),
  eta   = c(.6, .01, .001, .9, .1, 1)
)

param_configs$label <- paste0(
  "(x_m=", param_configs$x_m,
  ", α=", param_configs$alpha,
  ", b=", param_configs$b,
  ", η=", param_configs$eta,
  ")"
)

x_pdf <- seq(0, 10, length.out = 600)
x_cdf <- seq(0, 10, length.out = 600)



pdf_data <- expand.grid(x = x_pdf, id = 1:nrow(param_configs)) %>%
  mutate(
    alpha = param_configs$alpha[id],
    b     = param_configs$b[id],
    eta   = param_configs$eta[id],
    label = param_configs$label[id],
    value = pg_pdf(x, alpha, b, eta, xm)
  )

cdf_data <- expand.grid(x = x_cdf, id = 1:nrow(param_configs)) %>%
  mutate(
    alpha = param_configs$alpha[id],
    b     = param_configs$b[id],
    eta   = param_configs$eta[id],
    label = param_configs$label[id],
    value = pg_cdf(x, alpha, b, eta, xm)
  )

# 1. Define the specific color palette from the image
# Adding a 6th color (orange) to cover all 6 configs in your data
my_colors <- c("black", "red", "blue", "green3", "pink", "orange")

'p_pdf <- ggplot(pdf_data,
                aes(x = x, y = value,
                    color = label, linetype = label)) +
  geom_line(linewidth = 1, linetype = "solid") +
  labs(
    title = "Pareto–Gompertz PDF (xm = 1)",
    x = "x",
    y = "f(x)",
    color = "Parameters",
    linetype = "Parameters"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "inside"
  )


p_pdf'
# 2. Updated PDF Plot
p_pdf <- ggplot(pdf_data, 
                aes(x = x, y = value, color = label)) +
  geom_line(linewidth = 1) +
  # Apply the manual color scale here
  scale_color_manual(values = my_colors,name=NULL) + 
  labs(
    title = "PDF of PG models",
    x = "x",
    y = "f(x)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = c(0.72, 0.75),
    legend.text = element_text(size = 8)
  )
'p_cdf <- ggplot(cdf_data,
                aes(x = x, y = value,
                    color = label, linetype = label)) +
  geom_line(linewidth = 1, linetype = "solid") +
  labs(
    title = "Pareto–Gompertz CDF (xm = 1)",
    x = "x",
    y = "F(x)",
    color = "Parameters",
    linetype = "Parameters"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = c(0.72, 0.35),
    legend.justification = c(0, 1),
    )
    legend.text = element_text(size = 8)
p_cdf'
# 3. Updated CDF Plot
p_cdf <- ggplot(cdf_data, 
                aes(x = x, y = value, color = label)) +
  geom_line(linewidth = 1) +
  # Apply the manual color scale here
  scale_color_manual(values = my_colors,name=NULL) +
  labs(
    title = "CDF of PG models",
    x = "x",
    y = "F(x)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = c(.75,.20),
    legend.text = element_text(size = 8)
  )

grid.arrange(p_cdf, ncol = 1)
grid.arrange(p_pdf, ncol = 1)

library(cowplot)

plot_grid(
  p_pdf, p_cdf,
  label_size = 12,
  ncol = 2,
  align = "hv"
)