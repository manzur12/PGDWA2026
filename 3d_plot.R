# =====================================================
# Pareto–Gompertz PDF
# =====================================================
pg_pdf <- function(x, alpha, eta, xm = 1, b = 1) {
  bx <- pmin(b * x, 50)  # avoid overflow
  num <- alpha * eta * b * xm^alpha * exp(bx)
  den <- (xm + eta * (exp(bx) - 1))^(alpha + 1)
  num / den
}

# =====================================================
# Numerical Mean
# =====================================================
pg_mean_num <- function(alpha, eta, xm = 1, b = 1) {
  integrate(
    function(x) x * pg_pdf(x, alpha, b, xm, eta),
    lower = 0,
    upper = 50,  # sufficiently large to approximate infinity
    rel.tol = 1e-10
  )$value
}

# =====================================================
# Grid of parameters
# =====================================================
alpha_vals <- seq(2.1, 2.5, length.out = 20)
eta_vals     <- seq(0.5, 2.5, length.out = 20)

# Initialize matrix for means
mean_matrix <- matrix(0, nrow = length(alpha_vals), ncol = length(eta_vals))

# Compute numerical mean over grid
for (i in seq_along(alpha_vals)) {
  for (j in seq_along(eta_vals)) {
    mean_matrix[i, j] <- pg_mean_num(alpha_vals[i], eta_vals[j])
  }
}

mean_matrix
# =====================================================
# 3D Surface Plot (Diagonal/angled view)
# =====================================================
library(plotly)

fig <- plot_ly(
  x = ~eta_vals, 
  y = ~alpha_vals, 
  z = ~mean_matrix, 
  type = 'surface',
  colorscale = "Viridis",
  showscale = FALSE
  
)

fig <- fig %>% layout(
  scene = list(
    xaxis = list(title = "eta"),
    yaxis = list(title = "alpha"),
    zaxis = list(title = "Mean"),
    camera = list(eye = list(x = 1.5, y = 1.5, z = 1.2)) # diagonal/angled view
  )
)

fig


# =====================================================
# Pareto–Gompertz PDF
# =====================================================
pg_pdf <- function(x, alpha, eta, xm = 1, b = 1) {
  bx <- pmin(b * x, 50)   # numerical stabilization
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
alpha_vals <- seq(2.1, 2.5, length.out = 200)
eta_vals     <- seq(0.5, 2.5, length.out = 200)

# =====================================================
# Surface Matrix (Variance)
# =====================================================
Z_var <- outer(
  alpha_vals,
  eta_vals,
  Vectorize(function(a, eta) pg_var_num(a, eta, xm = 1, b = 1))
)
Z_var
# =====================================================
# 3D Surface Plot (Variance)
# =====================================================
library(plotly)

fig_var <- plot_ly(
  x = eta_vals,
  y = alpha_vals,
  z = Z_var,
  type = "surface",
  showscale = FALSE
)

fig_var <- fig_var %>% layout(
  scene = list(
    xaxis = list(title = "eta"),
    yaxis = list(title = "Alpha"),
    zaxis = list(title = "Variance"),
    camera = list(
      eye = list(x = 1.7, y = 1.7, z = 1.2)  # diagonal angle
    )
  )
)

fig_var








