

pop_size <- 10000

# means
mus <- c(0, 1)

# SD
sds <- c(1, 1)
sds_mat <- diag(sds)

# Correlation matrix (we assume some multicollinearity)
cor_mat <- matrix(c(1, 0.2, 0.2, 1), nrow = 2, ncol = 2)

# Convert to variance covariance matrix
varcov <- sds_mat %*% cor_mat %*% sds_mat

# true coefficients
b <- c(0, 2, 1)

# Choose sigma.est
sigma <- 0.5

gen_pop <- function(m,          # population size
                    mus,        # means
                    varcov,     # variance-covariance matrix
                    b,          # true coefficients
                    sigma       # sigma est
) {
  X <- mvrnorm(n = m, mu = mus, Sigma = varcov)
  mu_y <-
    cbind(1, X) %*% b + rnorm(m, 0, sigma)
  Y <- exp(mu_y)
  population <- data.frame("Y" = Y,
                           "X1" = X[, 1],
                           "X2" = X[, 2])
  return(population)
}

m <- 1000

pop <- gen_pop(
  m = m,
  mus = mus,
  varcov = varcov,
  b = b,
  sigma_est = sigma_est
)

dta <- pop[sample(1:m, size = 100),]


reg1 <- lm(log(Y) ~ X1 + X2, data = dta)


# Simulation

X1_range <- seq(min(dta$X1), max(dta$X1), length.out = 100)

scenario <- cbind(1, X1_range, mean(dta$X2))

S <- mvrnorm(1000, coef(reg1), vcov(reg1))

EV_logged <- S %*% t(scenario)

# unlog (this is incorrect!)

EV_wrong <- exp(EV_logged)

EV_wrong_mean <- apply(EV_wrong, 2, mean)
EV_wrong_ci <- apply(EV_wrong, 2, quantile, probs = c(0.025, 0.975))

# Step 4 & 5

# Get sigma

sigma_est <- summary(reg1)$sigma

# empty dataframes to store unlogged expected values

EV_unlogged <- matrix(nrow = dim(EV_logged)[1],
                      ncol = dim(EV_logged)[2])


pb = txtProgressBar(min = 0, max = length(dim(EV1_unlogged)[1]), initial = 0) 

for (i in 1:dim(EV_unlogged)[1]) {
  for (k in 1:dim(EV_unlogged)[2]){
    EV_unlogged[i, k] <- mean(exp(EV_logged[i, k] + rnorm(1000, 0, sigma_est)))
  }
  progress(i, dim(EV_unlogged)[1])
}


## Summarize

EV_mean <- apply(EV_unlogged, 2, mean)
EV_ci <- apply(EV_unlogged, 2, quantile, probs = c(0.025, 0.975))


# Wooldridge Fix


EV_W <- exp(EV_logged) * exp(sigma_est^2 / 2)
EV_W_mean <- apply(EV_W, 2, mean)
EV_W_ci <- apply(EV_W, 2, quantile, probs = c(0.025, 0.975))

plot(X1_range,
     EV_ci[1,],
     type = "n",
     ylim = c(0, 35),
     xlab = "X1",
     ylab = "Y")
polygon(c(rev(X1_range), X1_range), 
        c(rev(EV_ci[2,]), EV_ci[1,]),
        col = "gray80",
        border = NA)
lines(X1_range,
      EV_mean)

# wrong 
polygon(c(rev(X1_range), X1_range), 
        c(rev(EV_wrong_ci[2,]), EV_wrong_ci[1,]),
        col = alpha("red", 0.6),
        border = NA)
lines(X1_range,
      EV_wrong_mean)

# Wooldridge
lines(X1_range,
      EV_W_mean,
      col = "blue")
lines(X1_range,
      EV_W_ci[1,],
      col = "blue",
      lty = "dashed")
lines(X1_range,
      EV_W_ci[2,],
      col = "blue",
      lty = "dashed")



