library(viridis)
library(MASS)

n_draws <- 100
X <- cbind(1, rnorm(n_draws), rbinom(n_draws, size = 1, prob = 0.5))
betas <- c(0, 1, 2)
orig_mean <- X %*% betas
orig_sd <- 0.5
ln_y <- rnorm(n_draws, mean = orig_mean, sd = orig_sd)

plot(X[,2], ln_y, col = X[,3]+1)

y <- exp(ln_y)

plot(X[,2], y, col = X[,3] + 1)

dens_ln_y <- density(ln_y, from = min(ln_y), to = max(ln_y))
dens_y <- density(y, from = min(y), to = max(y))

col_vec <- viridis(4)

plot(
  dens_ln_y,
  xlim = c(min(c(y, ln_y)), max(c(y, ln_y))),
  ylim = c(-0.05, max(c(
    dens_y$y + 0.7, dens_ln_y$y
  ))),
  bty = "n",
  ylab = "",
  yaxt = "n",
  type = "n",
  xlab = "",
  main = paste0(
    n_draws,
    " draws from ln(y) and y for ln(y) = N(",
    mean(orig_mean),
    ", ",
    orig_sd^2,
    ")"
  )
)

points(y,
       rep(0.45, length(y)),
       pch = 16,
       cex = 0.5,
       col = col_vec[2])
points(
  ln_y,
  rep(-0.05, length(ln_y)),
  pch = 16,
  cex = 0.5,
  col = col_vec[1]
)
segments(
  y0 = rep(-0.04, length(ln_y)),
  y1 = rep(0.44, length(y)),
  x0 = ln_y,
  x1 = y,
  col = adjustcolor("grey", alpha = 0.2)
)
lines(dens_y$x, dens_y$y + 0.5, lwd = 2, col = col_vec[2])
lines(dens_ln_y$x, dens_ln_y$y, lwd = 2, col = col_vec[1])


points(median(y),
       0.45,
       pch = 16,
       cex = 0.5,
       col = col_vec[3])
points(median(ln_y),
       -0.05 ,
       pch = 16,
       cex = 0.5,
       col = col_vec[3])
segments(
  y0 = -0.04,
  y1 = 0.44,
  x0 = median(ln_y),
  x1 = median(y),
  col = col_vec[3]
)
segments(
  y0 = -0.04,
  y1 = 0.44,
  x0 = mean(ln_y),
  x1 = mean(y),
  col = col_vec[4],
  lty = "dashed"
)

legend(
  "topright",
  legend = c("Log Scale", "Original Scale", "Median", "Mean"),
  col = c(col_vec[1], col_vec[2], col_vec[3], col_vec[4]),
  lty = c("solid", "solid", "solid", "dashed"),
  lwd = 2,
  bty = "n"
)


betas <- solve(t(X) %*% X) %*% t(X) %*% ln_y
e <- (diag(1, n_draws) - X %*%  solve(t(X) %*% X) %*% t(X)) %*% ln_y
sigma_sq <- t(e) %*% e / (n_draws - ncol(X))

cov <- as.numeric(sigma_sq) * solve(t(X) %*% X)

nsim <- 1000

S <- mvrnorm(nsim, betas, cov)
Xb <- S %*% t(X)
pred_ln_y <- Xb + rnorm(nsim, 0, sd = sqrt(sigma_sq))
pred_y <- exp(pred_ln_y)

points(rowMeans(S %*% t(X)), rep(-0.02, nsim), pch = 16, cex = 0.3)
points(pred_ln_y, rep(0, nsim*nrow(X)), pch = 16, cex = 0.3)

points(pred_y, rep(0.47, nsim*nrow(X)), pch = 16, cex = 0.3)
points(rowMeans(pred_y), rep(0.49, nsim), pch = 16, cex = 0.3)
points(apply(pred_y, 1, median), rep(0.51, nsim), pch = 16, cex = 0.3)
points(rowMeans(exp(Xb)), rep(0.53, nsim), pch = 16, cex = 0.3)

simulate_exp_ci <- function(nsim, distribution, N = 1000) {
  pred_y <- list()
  for (n in 1:N) {
    sim_ln_y <-
      rnorm(
        nsim,
        mean = mean(distribution),
        sd = sd(distribution) / sqrt(length(distribution))
      )
    pred_ln_y <- sim_ln_y + rnorm(nsim, 0, sd = sd(ln_y))
    pred_y[[paste0(n)]] <- exp(pred_ln_y)
    cat(n, "\n")
  }
  
  return(pred_y)
}


pred_y <- simulate_exp_ci(1000, ln_y)

points(sim_ln_y, rep(0, nsim), pch = 16, cex = 0.3)

pred_ln_y <- sim_ln_y + rnorm(nsim, 0, sd = sd(ln_y))
points(pred_ln_y, rep(-0.02, nsim), pch = 16, cex = 0.3)


points(pred_y$'1', rep(0.47, nsim), pch = 16, cex = 0.3)

points(sapply(pred_y, mean),
       rep(0.5, 1000),
       pch = 16,
       cex = 0.3)
points(sapply(pred_y, median),
       rep(0.52, 1000),
       pch = 16,
       cex = 0.3)
points(exp(sim_ln_y), rep(0.54, 1000), pch = 16, cex = 0.4)

mean(exp(sim_ln_y))
mean(sapply(pred_y, mean))



reg_int <- lm(ln_y ~ 1)

sqrt(vcov(reg_int))
sd(reg_int$residuals)
sd(ln_y) / sqrt(length(ln_y))
sd(ln_y)
mean(ln_y)

X <- matrix(rep(1, 1000), nrow = 1000, ncol = 1)

beta <- solve(t(X) %*% X) %*% t(X) %*% ln_y

e <- ln_y - X %*% beta

sigma_sq <- t(e) %*% e / (nrow(e) - ncol(X))
var(ln_y)

beta_se <- as.numeric(sigma_sq) * solve(t(X) %*% X)
summary(reg_int)
sqrt(beta_se)
sqrt(var(ln_y) / 1000)
