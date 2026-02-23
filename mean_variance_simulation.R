# =====================================================
# Load library
# =====================================================
library(plotly)

# =====================================================
# Pareto–Gompertz PDF
# =====================================================
pg_pdf <- function(x, alpha, eta, xm = 1, b = 1) {
  bx <- pmin(b * x, 50)  # numerical stabilization
  num <- alpha * eta * b * xm^alpha * exp(bx)
  den <- (xm + eta * (exp(bx) - 1))^(alpha + 1)
  num / den
}

# =====================================================
# Numerical Mean
# =====================================================
pg_mean_num <- function(alpha, eta, xm = 1, b = 1) {
  if (alpha <= 1) return(NA)
  integrate(
    function(x) x * pg_pdf(x, alpha, b, xm, eta),
    lower = 0,
    upper = 20 / eta,
    rel.tol = 1e-10
  )$value
}

# =====================================================
# Numerical Variance
# =====================================================
pg_var_num <- function(alpha, eta, xm = 1, b = 1) {
  if (alpha <= 2) return(NA)
  m1 <- pg_mean_num(alpha, b, xm, eta)
  integrate(
    function(x) (x - m1)^2 * pg_pdf(x, alpha, b, xm, eta),
    lower = 0,
    upper = 20 / eta,
    rel.tol = 1e-10
  )$value
}

# =====================================================
# Parameter Grid
# =====================================================
alpha_vals <- c(2.1, 2.2, 2.3, 2.4, 2.5)
eta_vals     <- c(0.5, 1, 1.5, 2, 2.5)

# =====================================================
# Compute Mean & Variance Table
# =====================================================
results <- expand.grid(alpha = alpha_vals, eta = eta_vals)
results$Mean     <- mapply(pg_mean_num, results$alpha, results$eta)
results$Variance <- mapply(pg_var_num, results$alpha, results$eta)

# Print table
print(results)

# =====================================================
# Surface Matrices for 3D Plot
# =====================================================
Z_mean <- outer(alpha_vals, eta_vals, Vectorize(pg_mean_num))
Z_var  <- outer(alpha_vals, eta_vals, Vectorize(pg_var_num))

# =====================================================
# Surface Plot Function
# =====================================================
plot_surface <- function(z_matrix, x_vals, y_vals, z_title) {
  plot_ly(
    x = x_vals,
    y = y_vals,
    z = z_matrix,
    type = "surface"
  ) %>%
    layout(
      scene = list(
        xaxis = list(title = "eta"),
        yaxis = list(title = "Alpha"),
        zaxis = list(title = z_title),
        camera = list(eye = list(x = 1.7, y = 1.7, z = 1.2))  # diagonal view
      )
    )
}

# =====================================================
# Plot 3D Surfaces
# =====================================================
fig_mean <- plot_surface(Z_mean, eta_vals, alpha_vals, "Mean")
fig_mean

fig_var <- plot_surface(Z_var, eta_vals, alpha_vals, "Variance")
fig_var
