source("scripts/00_setup.R")

# Variance analytic solution (with fixed sigma^2 = 1)

sigma2 <- 1
s <- 20
n <- 10
var_y_hat <- NA

for (i in 1:n_sims) {
  x_i <- runif(10, 10, 20)
  var_y_hat[i] <- (sigma2 / (n*sum(x_i^2)  - sum(x_i)^2)) * (sum(x_i^2) - 2*s*sum(x_i) + s^2*n)
}



# dgp parameters (from Rainey 2017)
n <- 10
b_cons <- 2.5
b_edu <- 0.1
sigma2 <- 1

# number of simuations
n_sims <- 10000  # mc simulations
n_ktw <- 10000  # number of simulations of qi

# do simulations
df <- data.frame(edu)
prog <- progress_estimated(n_sims)
qi_mle <- qi_avg <- qi_mle_c <- qi_avg_c <- qi_log <- qi_median <- var_qi_log <- var_y_hat_analytic <- numeric(n_sims)

sigma_h <- NA
sigma_2_hat <- NA

for (sim in 1:n_sims) {
  edu <- runif(n, 10, 20)
  Xbeta <- b_cons + b_edu*edu
  error <- rnorm(n, 0, sqrt(sigma2))      # draw error
  log_income <- Xbeta + error             # income on log scale
  df$income <- exp(log_income)            # income
  fit <- lm(log(income) ~ edu, data = df) # regression estimation
  beta_hat <- coef(fit)                   # coefs
  Sigma_hat <- vcov(fit)                  # vcov
  beta_tilde <- MASS::mvrnorm(n_ktw, mu = beta_hat, Sigma = Sigma_hat) # simulation
  
  var_qi_log[sim] <- var((beta_tilde[, 1] + beta_tilde[, 2]*20))
  
  sigma_2_hat <- (summary(fit)$sigma^2)
  var_y_hat_analytic[sim] <- (sigma_2_hat / (n*sum(df$edu^2)  - sum(df$edu)^2)) * (sum(df$edu^2) - 2*s*sum(df$edu) + s^2*n)
  
  prog$tick()$print()
  
}


#pdf("figures/var_distr.pdf")
plot(density(var_y_hat),
     xlab = "",
     ylab = "",
     #yaxt = "n",
     bty = "n",
     main = "Variance y_hat",
     xlim = c(0, 2.5),
     col = "black")
lines(density(var_qi_log), col = "red")
lines(density(var_y_hat_analytic), col = "blue", lty = "dashed")
abline(v = 0.308)
legend("topright",
       legend = c("Analytic solution over samples \nwith fixed sigma^2", 
                  "Calculate variance after step 3",
                  "Analytic after 3 steps \n(estimated sigma^2)"),
       lty = 1,
       col = c("black", "red", "blue"),
       bty = "n")
#dev.off()
