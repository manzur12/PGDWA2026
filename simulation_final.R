# -----------------------------
# PG Distribution Simulation & MLE
# -----------------------------

library(bbmle)

set.seed(123)

# True parameters
alpha_true <- 0.8
b_true     <- 1.2
eta_true   <- 2
xm         <- 1

# Monte Carlo settings
R  <- 1000                # number of replications
ns <- c(50, 100, 200, 500, 1000)  # sample sizes

# Function to simulate PG data
rpg <- function(n, alpha, b, eta, xm=1){
  u <- runif(n)
  (1/b) * log(1 + (xm/eta) * ((1 - u)^(-1/alpha) - 1))
}

# Negative log-likelihood for PG
nll_pg <- function(alpha, b, eta, x) {
  if (alpha <= 0 || b <= 0 || eta <= 0) return(Inf)
  
  bx <- b * x
  temp <- 1 - eta + eta * exp(bx)
  
  if (any(temp <= 0) || any(!is.finite(temp))) return(Inf)
  
  -sum(log(alpha) + log(b) + log(eta) + bx - (alpha + 1) * log(temp))
}

# Storage for results
out <- data.frame()

# Loop over different sample sizes
for(n in ns){
  
  alpha_hat <- b_hat <- eta_hat <- numeric(R)
  
  for(i in 1:R){
    
    # Simulate data
    x <- rpg(n, alpha_true, b_true, eta_true, xm)
    
    # Fit MLE
    fit <- try(
      mle2(
        nll_pg,
        start = list(alpha=0.6, b=1, eta=2),
        method = "L-BFGS-B",
        lower = c(1e-6, 1e-6, 1e-6),
        data = list(x = x)
      ),
      silent = TRUE
    )
    
    if(inherits(fit, "try-error")) next  # skip failed iteration
    
    est <- coef(fit)
    alpha_hat[i] <- est["alpha"]
    b_hat[i]     <- est["b"]
    eta_hat[i]   <- est["eta"]
  }
  
  # Remove zeros if any iterations failed
  alpha_hat <- alpha_hat[alpha_hat > 0]
  b_hat     <- b_hat[b_hat > 0]
  eta_hat   <- eta_hat[eta_hat > 0]
  
  # Store mean estimates and MSE
  out <- rbind(out,
               data.frame(
                 n = n,
                 alpha_hat = mean(alpha_hat),
                 b_hat     = mean(b_hat),
                 eta_hat   = mean(eta_hat),
                 MSE_alpha = mean((alpha_hat - alpha_true)^2),
                 MSE_b     = mean((b_hat - b_true)^2),
                 MSE_eta   = mean((eta_hat - eta_true)^2)
               )
  )
}

# Round and display
Table4 <- round(out, 4)
Table4
