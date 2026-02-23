##########################################
# Maximum Likelihood Estimation
# PG, Pareto, Gompertz, PE, PW, Exp, Weibull
##########################################

library(bbmle)
library(e1071)
library(goftest)
library(ggplot2)


# =====================================================
#52 weekly extreme wind speed data in m/s.
# =====================================================

'x <- c(24.52, 69.15, 33.21, 30.79, 18.12, 94.34, 55.40, 28.28,
       49.73, 15.92, 21.81, 21.80, 30.56, 18.56, 87.63, 40.53,
       43.31, 20.10, 27.56, 90.16, 57.14, 28.98, 43.03, 49.13,
       50.47, 31.09, 27.17, 20.59, 65.99, 18.05, 19.44, 54.70,
       31.34, 48.32, 30.11, 46.22, 40.60, 24.76, 46.64, 75.50,
       45.78, 26.14, 23.06, 33.01, 51.56, 55.22, 67.04, 34.59,
       24.32, 55.88, 31.74, 76.17)'
'Chiodo, E.; De Falco, P. The Inverse Burr Distribution for Extreme Wind Speed Prediction: Genesis,
Identification and Estimation. Electr. Pow. Syst. Res. 2016, 141, 549–561. [CrossRef]'

# =====================================================








# =====================================================
#Data set 1
#The data set given in Table A.1 is taken from [1] Murthy and his colleagues (2004), p180. It represents the failure times of 50 components(per 1000h).

# =====================================================
'x <- c(0.036,0.058,0.061,0.074,0.078,0.086,0.102,
       0.103,0.114,0.116,0.148,0.183,0.192,0.254,
       0.262,0.379,0.381,0.538,0.570,0.574,0.590,
       0.618,0.645,0.961,1.228,1.600,2.006,2.054,
       2.804,3.058,3.076,3.147,3.625,3.704,3.931,
       4.073,4.393,4.534,4.893,6.274,6.816,7.896,
       7.904,8.022,9.337,10.940,11.020,13.880,
       14.730,15.080)'


# =====================================================
# Rainfall Data 
# https://www.researchgate.net/publication/317825442_The_Extended_Erlang-Truncated_Exponential_distribution_Properties_and_application_to_rainfall_data
# =====================================================

x <- c(24.8, 140.9, 54.1 ,153.5, 47.9, 165.5, 68.5, 153.1, 254.7 ,175.3 ,87.6, 150.6,
       147.9,354.7 ,128.5, 150.4 ,119.2 ,69.7 ,185.1, 153.4 ,121.7 ,99.3, 126.9 ,150.1,
       149.1, 143, 125.2 ,97.2 ,179.3 ,125.8, 101 ,89.8 ,54.6, 283.9 ,94.3 ,165.4,
       48.3, 69.2 ,147.1 ,114.2 ,159.4 ,114.9, 58.5, 76.6, 20.7 ,107.1, 244.5, 126,
       122.2, 219.9 ,153.2 ,145.3, 101.9 ,135.3, 103.1 ,74.7 ,174 ,126, 144.9, 226.3,
       96.2 ,149.3 ,122.3, 164.8 ,188.6, 273.2 ,61.2 ,84.3, 130.5 ,96.2 ,155.8, 194.6,
       92 ,131 ,137, 106.8 ,131.6, 268.2 ,124.5, 147.8 ,294.6, 101.6 ,103.1 ,247.5,
       140.2 ,153.3, 91.8, 79.4 ,149.2, 168.6, 127.7 ,332.8 ,261.6 ,122.9 ,273.4 ,178,
       177 ,108.5, 115 ,241, 76 ,127.5 ,190 ,259.5 ,301.5)

#33Taxes revenue data
'x<- c(5.9,20.4,14.9,16.2,17.2,7.8,6.1,9.2,10.2,9.6,
      13.3,8.5,21.6,18.5,5.1,6.7,17,8.6,9.7,39.2,
      35.7,15.7,9.7,10,4.1,36,8.5,8,9.2,26.2,21.9,
      16.7,21.3,35.4,14.3,8.5,10.6,19.1,20.5,7.1,
      7.1,18.1,16.5,11.9,7,8.6,12.5,10.3,11.2,6.1,
      8.4,11,11.6,11.9,5.2,6.8,8.9,7.1,10.8)'


summary(x)
n  <- length(x)
xm <- min(x)
xm

# =====================================================
# 1. Pareto–Gompertz Distribution
# =====================================================
nll_pg <- function(alpha, b, eta) {
  
  if (alpha <= 0 || b <= 0 || eta <= 0) return(Inf)
  
  bx <- pmin(b * x, 50)
  temp <- 1 - eta + eta * exp(bx)
  if (any(temp <= 0) || any(!is.finite(temp))) return(Inf)
  
  -( sum(log(alpha) + log(b) + log(eta) +
           bx - (alpha + 1) * log(temp)) )
}

mle_pg <- mle2(
  nll_pg,
  start  = list(alpha=0.5, b=0.001, eta=1),
  method = "L-BFGS-B",
  lower  = c(1e-6, 1e-6, 1e-6)
)
'nll_pg <- function(lalpha, lb, leta) {
  alpha <- exp(lalpha)
  b     <- exp(lb)
  eta   <- exp(leta)
  
  bx <- pmin(b*x, 50)
  temp <- 1 - eta + eta*exp(bx)
  if (any(temp <= 0)) return(Inf)
  
  -sum(log(alpha) + log(b) + log(eta) +
         bx - (alpha+1)*log(temp))
}

mle_pg <- mle2(
  nll_pg,
  start=list(lalpha=log(0.5), lb=log(0.001), leta=log(0.1)),
  method="BFGS"
)'


# =====================================================
# 3. Gompertz Distribution
# =====================================================
nll_gompertz <- function(b, eta) {
  
  if (b <= 0 || eta <= 0) return(Inf)
  
  bx <- pmin(b * x, 20)   # overflow protection
  
  ll <- sum(
    log(b) + log(eta) +
      bx -
      eta * (exp(bx) - 1)
  )
  
  if (!is.finite(ll)) return(Inf)
  
  -ll
}




mle_g <- mle2(
  nll_gompertz,
  start = list(b = 0.01, eta = 0.5),
  method = "L-BFGS-B",
  lower = c(b = 1e-8, eta = 1e-8)
)

# =====================================================
# 4. Pareto–Exponential Distribution
# =====================================================
nll_pe <- function(lambda, alpha) {
  if (lambda <= 0 || alpha <= 0) return(Inf)
  -( n*log(alpha) + n*alpha*log(xm) -
       n*log(lambda) -
       (alpha+1)*sum(log(xm + x/lambda)) )
}

mle_pe <- mle2(
  nll_pe,
  start=list(lambda=1000, alpha=1),
  method="L-BFGS-B",
  lower=c(1e-6,1e-6)
)


# =====================================================
# 7. Weibull Distribution
# =====================================================
nll_weibull <- function(k, lambda) {
  if (k <= 0 || lambda <= 0) return(Inf)
  -( sum(log(k) - log(lambda) +
           (k-1)*log(x/lambda) - (x/lambda)^k) )
}

mle_w <- mle2(
  nll_weibull,
  start=list(k=1, lambda=mean(x)),
  method="L-BFGS-B",
  lower=c(1e-6,1e-6)
)


#Summary of MLE
summary(mle_pg)
summary(mle_g)
summary(mle_pe)
summary(mle_w)


profile(mle_pg)
confint(mle_pg)



# =====================================================
# AICc Function
# =====================================================
AICc_manual <- function(fit, n) {
  k <- length(coef(fit))
  AIC(fit) + (2*k*(k+1))/(n-k-1)
}

BIC_manual <- function(fit, n) {
  k <- length(coef(fit))
  ll <- as.numeric(logLik(fit))
  -2 * ll + k * log(n)
}

# Pareto–Gompertz CDF
cdf_pg <- function(x, alpha, b, eta) {
  1 - (1 - eta + eta * exp(b*x))^(-alpha)
}

# Gompertz CDF
cdf_gompertz <- function(x, b, eta) {
  1 - exp(-eta * (exp(b*x) - 1))
}

# Pareto–Exponential CDF
cdf_pe <- function(x, lambda, alpha) {
  1 - (xm / (xm + x/lambda))^alpha
}

# Weibull CDF
cdf_weibull <- function(x, k, lambda) {
  1 - exp(-(x/lambda)^k)
}
# =====================================================
# Parameter Estimates
# =====================================================
pg_hat <- coef(mle_pg)
g_hat  <- coef(mle_g)
pe_hat <- coef(mle_pe)
w_hat  <- coef(mle_w)

# =====================================================
# Model Comparison
# =====================================================
info_criteria <- rbind(
  PG = c(logLik(mle_pg), AIC(mle_pg), AICc_manual(mle_pg,n), BIC_manual(mle_pg,n)),
  Gompertz = c(logLik(mle_g), AIC(mle_g), AICc_manual(mle_g,n), BIC_manual(mle_g,n)),
  ParetoExponential = c(logLik(mle_pe), AIC(mle_pe), AICc_manual(mle_pe,n), BIC_manual(mle_pe,n)),
  Weibull = c(logLik(mle_w), AIC(mle_w), AICc_manual(mle_w,n), BIC_manual(mle_w,n))
)

colnames(info_criteria) <- c("LogLik","AIC","AICc","BIC")
round(info_criteria,3)

gof_tests <- rbind(
  PG = c(
    ks.test(jitter(x), function(z) cdf_pg(z, pg_hat[1], pg_hat[2], pg_hat[3]))$statistic,
    cvm.test(x, function(z) cdf_pg(z, pg_hat[1], pg_hat[2], pg_hat[3]))$statistic,
    ad.test(x, function(z) cdf_pg(z, pg_hat[1], pg_hat[2], pg_hat[3]))$statistic
  ),
  Gompertz = c(
    ks.test(jitter(x), function(z) cdf_gompertz(z, g_hat[1], g_hat[2]))$statistic,
    cvm.test(x, function(z) cdf_gompertz(z, g_hat[1], g_hat[2]))$statistic,
    ad.test(x, function(z) cdf_gompertz(z, g_hat[1], g_hat[2]))$statistic
  ),
  ParetoExponential = c(
    ks.test(jitter(x), function(z) cdf_pe(z, pe_hat[1], pe_hat[2]))$statistic,
    cvm.test(x, function(z) cdf_pe(z, pe_hat[1], pe_hat[2]))$statistic,
    ad.test(x, function(z) cdf_pe(z, pe_hat[1], pe_hat[2]))$statistic
  ),
  Weibull = c(
    ks.test(jitter(x), function(z) cdf_weibull(z, w_hat[1], w_hat[2]))$statistic,
    cvm.test(x, function(z) cdf_weibull(z, w_hat[1], w_hat[2]))$statistic,
    ad.test(x, function(z) cdf_weibull(z, w_hat[1], w_hat[2]))$statistic
  )
)

colnames(gof_tests) <- c("KS","CvM","AD")
round(gof_tests,3)

final_comparison <- cbind(
  info_criteria[,c("LogLik","AIC","AICc","BIC")],
  gof_tests
)

round(final_comparison,4)

# ECDF VD Fitted CDF
'Fn <- ecdf(x)
xx <- seq(min(x), max(x), length.out = 500)

plot(Fn,
     main = "ECDF vs Fitted CDFs",
     xlab = "x",
     ylab = "Distribution Function",
     verticals = TRUE,
     do.points = FALSE,
     lwd = 2,
     col = "black")

lines(xx, cdf_pg(xx, pg_hat[1], pg_hat[2], pg_hat[3]),
      col = "red", lwd = 2)

lines(xx, cdf_gompertz(xx, g_hat[1], g_hat[2]),
      col = "blue", lwd = 2, lty = 2)

lines(xx, cdf_pe(xx, pe_hat[1], pe_hat[2]),
      col = "darkgreen", lwd = 2, lty = 3)

lines(xx, cdf_weibull(xx, w_hat[1], w_hat[2]),
      col = "purple", lwd = 2, lty = 4)

legend("bottomright",
       legend = c("ECDF",
                  "Pareto-Gompertz",
                  "Gompertz",
                  "Pareto-Exponential",
                  "Weibull"),
       col = c("black", "red", "blue", "darkgreen", "purple"),
       lty = c(1,1,2,3,4),
       lwd = 2,
       bty = "n")'
# Order corresponds to:
# 1 = ECDF
# 2 = Pareto–Gompertz
# 3 = Gompertz
# 4 = Pareto–Exponential
# 5 = Weibull
# 6 = (extra if needed later)

my_colors <- c("black", "green3", "blue1", "red", "purple", "orange")

xx <- seq(min(x), max(x), length.out = 500)

cdf_data <- rbind(
  data.frame(x = xx,
             value = cdf_pg(xx, pg_hat[1], pg_hat[2], pg_hat[3]),
             label = "Pareto–Gompertz"),
  data.frame(x = xx,
             value = cdf_gompertz(xx, g_hat[1], g_hat[2]),
             label = "Gompertz"),
  data.frame(x = xx,
             value = cdf_pe(xx, pe_hat[1], pe_hat[2]),
             label = "Pareto–Exponential"),
  data.frame(x = xx,
             value = cdf_weibull(xx, w_hat[1], w_hat[2]),
             label = "Weibull")
)

ecdf_data <- data.frame(
  x = sort(x),
  value = ecdf(x)(sort(x)),
  label = "Empirical CDF"
)

plot_data <- rbind(ecdf_data, cdf_data)

# Force legend & color order
plot_data$label <- factor(
  plot_data$label,
  levels = c("Empirical CDF",
             "Pareto–Gompertz",
             "Gompertz",
             "Pareto–Exponential",
             "Weibull")
)

p_ecdf <- ggplot(plot_data,
                 aes(x = x, y = value, color = label)) +
  geom_step(data = subset(plot_data, label == "ECDF"),
            linewidth = .8) +
  geom_line(data = subset(plot_data, label != "ECDF"),
            linewidth = .9) +
  scale_color_manual(values = my_colors, name = NULL) +
  labs(
    x = "x",
    y = "F(x)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = c(.75,.20),
    legend.text = element_text(size = 8)
  )

p_ecdf

curve_col <- "tomato"     # red/orange curve
ref_col   <- "purple"     # diagonal reference line


#p-p plot
pp_df <- data.frame(
  Empirical = (1:length(x) - 0.5) / length(x),
  Theoretical = cdf_pg(sort(x),
                       pg_hat[1], pg_hat[2], pg_hat[3])
)

ggplot(pp_df, aes(x = Empirical, y = Theoretical)) +
  geom_line(color = curve_col, linewidth = .8) +
  geom_abline(intercept = 0, slope = 1,
              color = ref_col, linewidth = 0.9) +
  labs(
    x = "Empirical CDF",
    y = "Theoretical CDF"
  ) +
  theme_bw()
#Q-Q

p <- (1:length(x) - 0.5) / length(x)
# Quantile function of Pareto–Gompertz distribution
q_pg <- function(p, alpha, b, eta) {
  (1 / b) * log( ((1 - p)^(-1 / alpha) - 1 + eta) / eta )
}


qq_df <- data.frame(
  Data = sort(x),
  Theoretical = q_pg(p, pg_hat[1], pg_hat[2], pg_hat[3])
)

ggplot(qq_df, aes(x = Data, y = Theoretical)) +
  geom_point(color = curve_col, size = 2) +
  geom_abline(intercept = 0, slope = 1,
              color = ref_col, linewidth = 0.9) +
  labs(
    x = "Data quantiles",
    y = "Theoretical quantiles"
  ) +
  theme_bw()

#ttt plot

ttt_df <- function(x) {
  x <- sort(x)
  n <- length(x)
  data.frame(
    r_n = (1:n) / n,
    TTT = cumsum(x) / sum(x)
  )
}

df_ttt <- ttt_df(x)

ggplot(df_ttt, aes(x = r_n, y = TTT)) +
  geom_line(color = curve_col, linewidth = .8) +
  geom_abline(intercept = 0, slope = 1,
              color = ref_col, linewidth = 0.9) +
  labs(
    x = "r/n",
    y = "TTT"
  ) +
  ggtitle("TTT Plot") +
  theme_bw()





# Hazard function using MLE estimates
h_pg <- function(x, alpha, b, eta){
  (alpha * eta * b * exp(b*x)) /
    (1 - eta + eta * exp(b*x))
}

xx <- seq(min(x), max(x), length.out = 500)

df_hazard <- data.frame(
  x = xx,
  hazard = h_pg(xx, pg_hat[1], pg_hat[2], pg_hat[3])
)

ggplot(df_hazard, aes(x = x, y = hazard)) +
  geom_line(color = curve_col, linewidth = .8) +
  labs(
    x = "x",
    y = "h(x)"
  ) +
  ggtitle("Hazard plot") +
  theme_bw()





