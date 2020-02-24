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

mean_median_Y <- function(theta, sigma_est){
  mean_Y <- mean(exp(rnorm(10000, theta, sigma_est)))
  median_Y <- median(exp(rnorm(10000, theta, sigma_est)))
  
  return(c(mean_Y, median_Y))
}

mean_median_Y(theta[1], sigma_est = 0.5)

for (i in 1:length(theta)) {
  mean_median <- mean_median_Y(theta[i], sigma_est = 0.5)
  bias_data$mean_Y[i] <- mean_median[1]
  bias_data$median_Y[i] <- mean_median[2]
}


pdf("error_exp_theta.pdf")
plot(x = exp(theta), y = exp(theta),
     type = "l",
     col = "red",
     xlab = expression(e^theta),
     ylab = "",
     xaxt = "n", yaxt = "n",
     bty = "n",
     ylim = c(min(exp(theta)), max(bias_data$mean_Y)),
     asp = 1,
     main = expression("Difference between " ~ e^theta ~ "and E(Y)"))  
axis(1,
     at = c(exp(min(theta)), exp(max(theta))),
     labels = c(expression(e^0), expression(e^100)))
axis(2,
     las = 2,
     at = c(exp(min(theta)), exp(max(theta))),
     labels = c(expression(e^0), expression(e^100)))

lines(x = exp(theta),
      y = bias_data$mean_Y,
      col = "black")
#lines(x = exp(theta),
#      y = bias_data$median_Y,
#      col = "blue",
#      lty = "dashed")

midscale <- (exp(min(theta)) + exp(max(theta))) / 2

text(x = midscale,
     y = midscale - 1/6*midscale,
     labels = expression(e^theta),
     col = "red")
text(x = midscale,
     y = midscale + 1/3*midscale,
     labels = expression(E(Y)),
     col = "black")
dev.off()

