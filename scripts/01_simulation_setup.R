source("scripts/00_setup.R")

# population size
m <- 10000

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
b <- c(2, 0.3, -0.5)

# Choose sigma.est
sigma_est <- 0.5

gen_pop <- function(m,          # population size
                    mus,        # means
                    varcov,     # variance-covariance matrix
                    b,          # true coefficients
                    sigma_est   # sigma est
                    ) {
                    X <- mvrnorm(n = m, mu = mus, Sigma = varcov)
                    mu_y <-
                      cbind(1, X) %*% b + rnorm(m, 0, sigma_est)
                    Y <- exp(mu_y)
                    population <- data.frame("Y" = Y,
                                             "X1" = X[, 1],
                                             "X2" = X[, 2])
                    return(population)
                    }

pop <- gen_pop(
  m = m,
  mus = mus,
  varcov = varcov,
  b = b,
  sigma_est = sigma_est
)

dta <- pop[sample(1:m, size = 100),]

# Plot the data

plot(dta$X1, log(dta$Y))
plot(dta$X1, dta$Y)

plot(dta$X2, log(dta$Y))
plot(dta$X2, dta$Y)

reg1 <- lm(log(Y) ~ X1 + X2, data = dta)

mean(dta$Y)
median(dta$Y)
mean(log(dta$Y))

mu_Y_hat <- cbind(1, dta$X1, dta$X2) %*% coef(reg1)
mean(mu_Y_hat)

mu_Y_hat <-
  cbind(1, mean(dta$X1), mean(dta$X2)) %*% coef(reg1) + var(reg1$residuals) /
  2
mu_Y_hat

median(log(dta$Y))

median(dta$Y)
exp(mu_Y_hat)
mean(dta$Y)

mm <- cbind(1, dta$X1, dta$X2)

M <- 1000
nsim <- 1000

sim_list <- list()

for (sim in 1:nsim) {
  S <- mvrnorm(M, coef(reg1), vcov(reg1))
  
  mu_Y_ova <- S %*% t(mm)
  
  # Predicted Values
  
  fundamental_uncertainty <-
    matrix(rnorm(prod(dim(mu_Y_ova)), 0, sd = sd(reg1$residuals)),
           nrow = dim(mu_Y_ova)[1],
           ncol = dim(mu_Y_ova)[2])
  
  pred_mu_Y_ova <- mu_Y_ova + fundamental_uncertainty
  
  sim_list[[sim]] <- pred_mu_Y_ova
  cat(sim, "\n")
}

mean_y <- sapply(sim_list, function(x)
  mean((exp(x))))

median_y <- sapply(sim_list, function(x)
  median((exp(x))))

robust_mean_y <-  sapply(sim_list, function(x)
  mean((exp(x[x<quantile(x, 0.999)]))))

mean_col_median_y <- sapply(sim_list, function(x) (apply(exp(x), 2, mean)))

mean(dta$Y)
hist(mean_col_median_y)

robust_mean_y

pred_Y <- sapply(sim_list, function(x) exp(x))

plot(density(dta$Y, from = 0), xlim = c(0, 200), ylim = c(0, 0.3))
lines(density(mean_col_median_y, from = 0), col = "red")

mean(dta$Y)
hist(mean_y)
hist(robust_mean_y)

median(dta$Y)
mean(log(dta$Y))
hist(sapply(sim_list, mean))

hist(median_y)

pred_Y_ova <- exp(pred_mu_Y_ova)

# pred_Y <- rowMeans(pred_Y_ova)

hist(pred_Y_ova)
hist(dta$Y)
mean(pred_Y_ova)
median(pred_Y_ova)
median(dta$Y)
mean(dta$Y)



# Standard Simulation step 1-3

sim_function <- function(lm_obj, nsim = 1000, scenario){
  # Step 1: Get the regression coefficients
  beta_hat <- coef(lm_obj)
  # Step 2: Generate sampling distribution
  # Step 2.1: Get the variance-covariance matrix.
  V_hat <- vcov(lm_obj)
  # Step 2.2: Draw from the multivariate normal distribution.
  library(MASS)
  S <- mvrnorm(nsim, beta_hat, V_hat)
  # Step 3: Choose interesting covariate values.
  # Make sure the matrix multiplication also works for single scenarios
  if(is.null(nrow(scenario))){
    scenario <- matrix(scenario, nrow = 1)
  }
  # Print a message if the scenario does not fit the regression.
  if(ncol(scenario) != length(lm_obj$coefficients)){
    return(cat("The scenario has the wrong number of variables."))
  }
  # Step 4: Calculate Quantities of Interest -
  # Expected Values
  EV <- S %*% t(scenario)
  
  return(EV)
}


# average case approach

ac_EV <- sim_function(lm_obj = reg1,
                      nsim = 1000,
                      scenario = cbind(1, mean(dta$X1), mean(dta$X2)))


mean(exp(ac_EV))
median(exp(ac_EV))
exp(reg1$coefficients %*% c(1, mean(dta$X1), mean(dta$X2)))


# observed value approach

ovs_EV <- sim_function(lm_obj = reg1,
                       nsim = 1000,
                       scenario = cbind(1, dta$X1, dta$X2))

mean(dta$Y)
mean(exp(ovs_EV))
median(exp(ovs_EV))
median(dta$Y)




### Investigate Error (disregarding estimation uncertainty) ###

# variation in theta, hold sigma^2 constant

theta <- seq(0, 100, 1)

bias_data <- data.frame("theta" = theta,
                        "mean_Y" = NA,
                        "median_Y" = NA)

mean_median_Y <- function(theta, sigma_est, m){
  mean_Y <- mean(exp(rnorm(m, theta, sigma_est)))
  median_Y <- median(exp(rnorm(m, theta, sigma_est)))
  
  return(c(mean_Y, median_Y))
}

for (i in 1:length(theta)) {
  mean_median <- mean_median_Y(theta[i], sigma_est = 0.5, m = 10000)
  bias_data$mean_Y[i] <- mean_median[1]
  bias_data$median_Y[i] <- mean_median[2]
}


pdf("figures/error_exp_theta.pdf")
plot(x = exp(theta), y = exp(theta),
     type = "n",
     col = "red",
     xlab = expression(e^theta),
     ylab = "",
     xaxt = "n", yaxt = "n",
     bty = "n",
     ylim = c(min(exp(theta)), max(bias_data$mean_Y)),
     asp = 1,
     main = expression("Difference between " ~ e^theta ~ "and E(Y) conditional on" ~ theta ~ "," ~ sigma^2 == 0.5))

# Grid
segments(x0 = seq(min(exp(theta)), max(exp(theta)), length.out = 10),
         x1 = seq(min(exp(theta)), max(exp(theta)), length.out = 10),
         y0 = min(exp(theta)), y1 = max(exp(theta)),
         col = "grey", lwd = 1.2,
         lty = "dashed")
segments(y0 = seq(min(exp(theta)), max(exp(theta)), length.out = 10),
         y1 = seq(min(exp(theta)), max(exp(theta)), length.out = 10),
         x0 = min(exp(theta)), x1 = max(exp(theta)),
         col = "grey", lwd = 1.2,
         lty = "dashed")

axis(1,
     at = c(exp(min(theta)), exp(max(theta))),
     labels = c(expression(e^0), expression(e^100)),
     pos = 0, lwd = 1.2)
axis(2,
     las = 2,
     at = c(exp(min(theta)), exp(max(theta))),
     labels = c(expression(e^0), expression(e^100)),
     pos = 0, lwd = 1.2)



# exp(theta) == exp(theta)
lines(x = exp(theta), y = exp(theta),
      col = "red")


lines(x = exp(theta),
      y = bias_data$mean_Y,
      col = "black")
#lines(x = exp(theta),
#      y = bias_data$median_Y,
#      col = "blue",
#      lty = "dashed")

lines(x = exp(theta),
      y = exp(theta) * exp((0.5^2) / 2),
      col = "green",
      lty = "dashed")


midscale <- (exp(min(theta)) + exp(max(theta))) / 2

text(x = midscale + 1/5*midscale,
     y = midscale,
     labels = expression(e^theta == e^theta),
     col = "red")
text(x = midscale,
     y = midscale + 1/3*midscale,
     labels = expression(E(Y)),
     col = "black")
arrows(x0 = exp(bias_data$theta[101]),
       y0 = exp(bias_data$theta[101]),
       x1 = exp(bias_data$theta[101]),
       y1 = bias_data$mean_Y[101],
       code = 3,
       length = 0.1,
       lwd = 1.5)
text(x = exp(bias_data$theta[101] + 0.09),
     y = exp(bias_data$theta[101]) + (bias_data$mean_Y[101] - exp(bias_data$theta[101])) / 2,
     labels = expression(E(Y) - e^theta))
dev.off()


# variation in sigma^2, holding constant theta

sigma2 <- seq(0, 5, length.out = 100)
theta_fix <- 1

bias_data2 <- data.frame("sigma2" = sigma2,
                         "mean_Y" = NA,
                         "median_Y" = NA)



for (i in 1:nrow(bias_data2)) {
  mean_median <- mean_median_Y(theta = theta_fix, sigma_est = bias_data2$sigma2[i], m = 10000000)
  bias_data2$mean_Y[i] <- mean_median[1]
  bias_data2$median_Y[i] <- mean_median[2]
}

bias_data2$dif <- bias_data2$mean_Y - rep(exp(2), 100)

pdf("figures/error_sigma2.pdf")
plot(x = bias_data2$sigma2, 
     y = bias_data2$dif,
     type = "l",
     col = "black",
     xlab = expression(sigma^2),
     ylab = expression(E(Y) - e^theta),
     ylim = c(-1/5*abs(min(bias_data2$dif)),
             max(bias_data2$dif)),
     xlim = c(0, 5),
     xaxt = "n", yaxt = "n",
     bty = "n",
     main = expression("Difference between" ~ e^theta ~ "and E(Y) conditional on" ~ sigma^2 ~ "," ~ theta == 1))  
abline(h = 0,
       lty = "dashed",
       col = "black")
axis(1,
     at = min(bias_data2$sigma2):max(bias_data2$sigma2))
lines(x = bias_data2$sigma2,
      y = exp(theta_fix) * exp((bias_data2$sigma2^2)/2),
      lty = "dashed",
      col = "green")
dev.off()