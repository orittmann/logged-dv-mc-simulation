setup_pop <-
  function(n,
           ff,
           type = "num",
           coefs = 1,
           zero_centered = T,
           noise_variance = 1) {
    variables <- all.vars(as.formula(ff))
    variable_list <- list()
    
    if (length(type) == length(variables) - 1) {
      type_vec <- type
    } else {
      type_vec <- rep(type[1], length(variables) - 1)
    }
    
    if (zero_centered == T) {
      for (variable in variables[2:length(variables)]) {
        variable_list[[paste0(variable)]] <-
          if (type_vec[variable == variables[2:length(variables)]] == "bin") {
            rbinom(n, size = 1, prob = 0.5)
          } else {
            rnorm(n, 0, 1)
          }
      }
    } else {
      for (variable in variables[2:length(variables)]) {
        variable_list[[paste0(variable)]] <-
          if (type_vec[variable == variables[2:length(variables)]] == "bin") {
            rbinom(n, size = 1, prob = runif(1))
          } else {
            rnorm(n, rnorm(1, 0, sd = 10), runif(1, min = 1, max = 10))
          }
      }
    }
    tmp_df <- do.call(data.frame, variable_list)
    tmp_ff <- as.formula(sub(".*~", "~", ff))
    mf <- model.frame(tmp_ff, tmp_df)
    mm <- model.matrix(tmp_ff, data = mf)
    # If the number of provided coefficients is not the same number as the number
    # of actual coefficients we take the first value of the vector and set all
    # coefficients to that value.
    if (length(coefs) == ncol(mm)) {
      coef_vec <- coefs
    } else {
      coef_vec <- rep(coefs[1], ncol(mm))
    }
    # Create the dependent variable based on the model matrix and coefficients.
    # Add some random noise.
    variable_list[[paste0(variables[1])]] <-
      mm %*% coef_vec + rnorm(n, 0, noise_variance)
    df <- do.call(data.frame, variable_list)
    return(df)
  }

ff <- "y~x1*x2"
n <- 1000


pop_df <-
  setup_pop(
    n,
    ff,
    type = c("num", "bin"),
    coefs = c(0, 1,-2, 0.3),
    noise_variance = 0.4
  )
dev.off()
plot(pop_df$x1, exp(pop_df$y), col = pop_df$x2 + 1)
plot(pop_df$x1, (pop_df$y), col = pop_df$x2 + 1)

reg_ia <- lm(ff, data = pop_df)

tmp_ff <- as.formula(sub(".*~", "~", ff))
mf <- model.frame(tmp_ff, pop_df)
mm <- model.matrix(tmp_ff, data = mf)

mf <- model.frame(tmp_ff, pop_df)
mm <- model.matrix(tmp_ff, data = mf)

source("scripts/03_sim_function.R")

scenario_1 <-
  cbind(1, seq(min(pop_df$x1), max(pop_df$x1), by = 0.2), 1, seq(min(pop_df$x1), max(pop_df$x1), by = 0.2))

pred_1 <-
  pred_values(
    model_theta = coef(reg_ia),
    model_vcov = vcov(reg_ia),
    scenario = scenario_1,
    model_sigma = sd(reg_ia$residuals)
  )

exp_pred_1 <- transform_predicted_values(pred_1)

sum_pred_1 <- summarize_sim_values(pred_1)
sum_exp_pred_1 <- summarize_sim_values(exp_pred_1)

scenario_0 <-
  cbind(1, seq(min(pop_df$x1), max(pop_df$x1), by = 0.2), 0, 0)

pred_0 <-
  pred_values(
    model_theta = coef(reg_ia),
    model_vcov = vcov(reg_ia),
    scenario = scenario_0,
    model_sigma = sd(reg_ia$residuals)
  )

exp_pred_0 <- transform_predicted_values(pred_0)

sum_pred_0 <- summarize_sim_values(pred_0)
sum_exp_pred_0 <- summarize_sim_values(exp_pred_0)





par(mfrow = c(2, 2))

plot(
  pop_df$x1,
  (pop_df$y),
  col = adjustcolor(pop_df$x2 + 1, alpha = 0.5),
  bty = "n",
  las = 1,
  main = "Log-Scale",
  ylab = "log(Y)",
  xlab = "X1",
  pch = 16
)
lines(seq(min(pop_df$x1), max(pop_df$x1), by = 0.2),
      apply(sum_pred_1, 2, mean),
      col = "red",
      lwd = 2)
lines(
  seq(min(pop_df$x1), max(pop_df$x1), by = 0.2),
  apply(sum_pred_1, 2, quantile, 0.025),
  col = "red",
  lwd = 1,
  lty = "dashed"
)
lines(
  seq(min(pop_df$x1), max(pop_df$x1), by = 0.2),
  apply(sum_pred_1, 2, quantile, 0.975),
  col = "red",
  lwd = 1,
  lty = "dashed"
)
lines(seq(min(pop_df$x1), max(pop_df$x1), by = 0.2),
      apply(sum_pred_0, 2, mean),
      col = "black",
      lwd = 2)
lines(
  seq(min(pop_df$x1), max(pop_df$x1), by = 0.2),
  apply(sum_pred_0, 2, quantile, 0.025),
  col = "black",
  lwd = 1,
  lty = "dashed"
)
lines(
  seq(min(pop_df$x1), max(pop_df$x1), by = 0.2),
  apply(sum_pred_0, 2, quantile, 0.975),
  col = "black",
  lwd = 1,
  lty = "dashed"
)




plot(
  pop_df$x1,
  exp(pop_df$y),
  col = adjustcolor(pop_df$x2 + 1, alpha = 0.5),
  bty = "n",
  las = 1,
  main = "Original-Scale",
  ylab = "log(Y)",
  xlab = "X1",
  pch = 16
)
lines(
  seq(min(pop_df$x1), max(pop_df$x1), by = 0.2),
  apply(sum_exp_pred_1, 2, mean),
  col = "red",
  lwd = 2
)
lines(
  seq(min(pop_df$x1), max(pop_df$x1), by = 0.2),
  apply(sum_exp_pred_1, 2, quantile, 0.025),
  col = "red",
  lwd = 1,
  lty = "dashed"
)
lines(
  seq(min(pop_df$x1), max(pop_df$x1), by = 0.2),
  apply(sum_exp_pred_1, 2, quantile, 0.975),
  col = "red",
  lwd = 1,
  lty = "dashed"
)
lines(
  seq(min(pop_df$x1), max(pop_df$x1), by = 0.2),
  apply(sum_exp_pred_0, 2, mean),
  col = "black",
  lwd = 2
)
lines(
  seq(min(pop_df$x1), max(pop_df$x1), by = 0.2),
  apply(sum_exp_pred_0, 2, quantile, 0.025),
  col = "black",
  lwd = 1,
  lty = "dashed"
)
lines(
  seq(min(pop_df$x1), max(pop_df$x1), by = 0.2),
  apply(sum_exp_pred_0, 2, quantile, 0.975),
  col = "black",
  lwd = 1,
  lty = "dashed"
)



fd <- sum_pred_0 - sum_pred_1

plot(
  seq(min(pop_df$x1), max(pop_df$x1), by = 0.2),
  apply(fd, 2, mean),
  type = "l",
  bty = "n",
  las = 1,
  main = "First Difference Log-Scale",
  ylab = "FD E(log(Y)|X2 = 0) - E(log(Y)|X2 = 1)",
  xlab = "X1"
)
lines(
  seq(min(pop_df$x1), max(pop_df$x1), by = 0.2),
  apply(fd, 2, quantile, 0.025),
  col = "black",
  lwd = 1,
  lty = "dashed"
)
lines(
  seq(min(pop_df$x1), max(pop_df$x1), by = 0.2),
  apply(fd, 2, quantile, 0.975),
  col = "black",
  lwd = 1,
  lty = "dashed"
)


exp_fd <- sum_exp_pred_0 - sum_exp_pred_1

plot(
  seq(min(pop_df$x1), max(pop_df$x1), by = 0.2),
  apply(exp_fd, 2, mean),
  type = "l",
  bty = "n",
  las = 1,
  main = "First Difference Original-Scale",
  ylab = "FD E(Y|X2 = 0) - E(Y|X2 = 1)",
  xlab = "X1"
)
lines(
  seq(min(pop_df$x1), max(pop_df$x1), by = 0.2),
  apply(exp_fd, 2, quantile, 0.025),
  col = "black",
  lwd = 1,
  lty = "dashed"
)
lines(
  seq(min(pop_df$x1), max(pop_df$x1), by = 0.2),
  apply(exp_fd, 2, quantile, 0.975),
  col = "black",
  lwd = 1,
  lty = "dashed"
)
