source("scripts/00_setup.R")
source("scripts/03_sim_function.R")
set.seed(200226)
# population size
pop_size <- 100000

n_obs <- 1000
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
b <- c(2, 0.3,-0.5)

# Choose sigma.est
sigma_est <- 1

gen_pop <- function(pop_size,
                    mus,
                    varcov,
                    b,
                    sigma_est) {
  X <- mvrnorm(n = pop_size, mu = mus, Sigma = varcov)
  mu_y <-
    cbind(1, X) %*% b + rnorm(pop_size, 0, sigma_est)
  Y <- exp(mu_y)
  population <- data.frame("Y" = Y,
                           "X1" = X[, 1],
                           "X2" = X[, 2])
  return(population)
}

pop <- gen_pop(
  pop_size = pop_size,
  mus = mus,
  varcov = varcov,
  b = b,
  sigma_est = sigma_est
)

mc_nsim <- 1000
aca_mean_coverage <-
  aca_mean_error <-
  aca_median_coverage <-
  aca_median_error <-
  ova_mean_coverage <-
  ova_mean_error <-
  ova_median_coverage <- ova_median_error <- rep(NA, mc_nsim)

ova <- F

for (mc_run in 1:mc_nsim) {
  cat("\n Monte Carlo Simulation. Run: ", mc_run, "\n")
  dta <- pop[sample(1:pop_size, size = n_obs), ]
  
  ff <- log(Y) ~ X1 + X2
  tmp_reg <- lm(ff, data = dta)
  
  mf <- model.frame(ff, dta)
  mm <- model.matrix(mf, dta)
  
  pv_aca <- pred_values(coef(tmp_reg),
                        vcov(tmp_reg),
                        apply(mm, 2, mean),
                        model_sigma = sd(tmp_reg$residuals)) %>%
    transform_predicted_values()
  
  
  qi_mean_aca <- pv_aca %>% summarize_sim_values(qi = mean)
  
  true_mean_aca <-
    exp(apply(mm, 2, mean) %*% b +  0.5 * sigma_est ^ 2)
  
  ci_mean_aca <- quantile(qi_mean_aca, c(0.025, 0.975))
  
  
  aca_mean_coverage[mc_run] <-
    true_mean_aca >= ci_mean_aca[1] &
    true_mean_aca <= ci_mean_aca[2]
  aca_mean_error[mc_run] <- true_mean_aca - mean(qi_mean_aca)
  
  qi_median_aca <- pv_aca %>% summarize_sim_values(qi = median)
  
  true_median_aca <- exp(apply(mm, 2, mean) %*% b)
  
  ci_median_aca <- quantile(qi_median_aca, c(0.025, 0.975))
  
  aca_median_coverage[mc_run] <-
    true_median_aca >= ci_median_aca[1] &
    true_median_aca <= ci_median_aca[2]
  aca_median_error[mc_run] <-
    true_median_aca - median(qi_median_aca)
  
  
  if(ova){
  pv_ova <- pred_values(coef(tmp_reg),
                        vcov(tmp_reg),
                        mm,
                        model_sigma = sd(tmp_reg$residuals)) %>%
    transform_predicted_values()
  
  
  qi_mean_ova <- pv_ova %>% summarize_sim_values(qi = mean)
  
  true_mean_ova <- mean(exp(mm %*% b +  0.5 * sigma_est ^ 2))
  
  ci_mean_ova <-
    quantile(apply(qi_mean_ova, 1, mean), c(0.025, 0.975))
  
  
  ova_mean_coverage[mc_run] <-
    true_mean_ova >= ci_mean_ova[1] &
    true_mean_ova <= ci_mean_ova[2]
  ova_mean_error[mc_run] <- true_mean_ova - mean(qi_mean_ova)
  
  qi_median_ova <- pv_ova %>% summarize_sim_values(qi = median)
  
  true_median_ova <- median(exp(mm %*% b))
  
  ci_median_ova <-
    quantile(apply(qi_median_ova, 1, median), c(0.025, 0.975))
  
  ova_median_coverage[mc_run] <-
    true_median_ova >= ci_median_ova[1] &
    true_median_ova <= ci_median_ova[2]
  ova_median_error[mc_run] <-
    true_median_ova - median(qi_median_ova)
  }
}


res_error <- list(aca_mean_error,
                  aca_median_error,
                  ova_mean_error,
                  ova_median_error)

par(mfrow = c(2, 2))
lapply(res_error, function(x){
  hist(x)
  abline(v = mean(x))
})


res_coverage <- list(aca_mean_coverage,
                     aca_median_coverage,
                     ova_mean_coverage,
                     ova_median_coverage)

sapply(res_coverage, mean)



res_error <- list(aca_mean_error,
     aca_median_error,
     ova_mean_error,
     ova_median_error)

par(mfrow = c(2, 2))
lapply(res_error, function(x){
  hist(x)
  abline(v = mean(x))
  abline(v = median(x), lty = "dashed")
})


res_coverage <- list(aca_mean_coverage,
                  aca_median_coverage,
                  ova_mean_coverage,
                  ova_median_coverage)

sapply(res_coverage, mean)
