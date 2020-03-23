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
            rbinom(n, size = 1, prob = 0.5)
          } else {
            rnorm(n, 7, 1.5)
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

ff <- "y~x1+x2"
n <- 1000


pop_df <-
  setup_pop(
    n,
    ff,
    zero_centered = F,
    type = c("num", "bin"),
    coefs = c(3, 0.7,-2),
    noise_variance = 1.5
  )
dev.off()
plot(pop_df[,1], exp(pop_df[,3]), col = pop_df[,2] + 1)
plot(pop_df[,1], (pop_df[,3]), col = pop_df[,2] + 1)


hist(exp(pop_df[,3]) / 1000, breaks = 50)
reg_ia <- lm(ff, data = pop_df)

tmp_ff <- as.formula(sub(".*~", "~", ff))
mf <- model.frame(tmp_ff, pop_df)
mm <- model.matrix(tmp_ff, data = mf)

mf <- model.frame(tmp_ff, pop_df)
mm <- model.matrix(tmp_ff, data = mf)

source("scripts/03_sim_function.R")

scenario_1 <-
  cbind(1, seq(min(pop_df$x1), max(pop_df$x1), by = 0.2), 1)

pred_1 <-
  pred_values(M = 1000,
    model_theta = coef(reg_ia),
    model_vcov = vcov(reg_ia),
    scenario = scenario_1,
    model_sigma = sd(reg_ia$residuals)
  )

exp_pred_1 <- transform_predicted_values(pred_1)

sum_pred_1 <- summarize_sim_values(pred_1)
sum_exp_pred_1 <- summarize_sim_values(exp_pred_1)

scenario_0 <-
  cbind(1, seq(min(pop_df$x1), max(pop_df$x1), by = 0.2), 0)

pred_0 <-
  pred_values(M = 1000,
    model_theta = coef(reg_ia),
    model_vcov = vcov(reg_ia),
    scenario = scenario_0,
    model_sigma = sd(reg_ia$residuals)
  )

exp_pred_0 <- transform_predicted_values(pred_0)

sum_pred_0 <- summarize_sim_values(pred_0)
sum_exp_pred_0 <- summarize_sim_values(exp_pred_0)








#pdf("figures/dummy_plot.pdf", width = 16)
m <-matrix(c(1, 2, 3, 4), ncol = 2, byrow = T)
layout(m,
       widths = c(1,1))

# Plot 1: Regressions on log scale
par(mar = c(5,     # bottom (5)
            6.5,   # left   (4)
            2,     # top    (4)
            0.5    # right  (2)
),
xpd = T)


col_vec <- viridis::viridis(2)

plot(
  pop_df$x1,
  (pop_df$y),
  col = adjustcolor(ifelse(pop_df$x2 == 0, col_vec[1], col_vec[2]), alpha = 0.5),
  las = 1,
  main = "",
  ylab = "",
  xlab = "",
  pch = 16,
  yaxt = "n",
  xaxt = "n",
  bty = "n"
)
box(lty = "solid", col = "grey50")
axis(2, las = 2, 
     col.axis = "grey50", col = "grey50")
axis(1, 
     col.axis = "grey50", col = "grey50")
mtext(2, 
      text = "log(INCOME)", 
      line = 4.1)
mtext(3, 
      text = expression(bold("Log-Scale")), 
      line = 0.4)
mtext(1, 
      text = "EDUCATION",
      line = 3.5)
lines(seq(min(pop_df$x1), max(pop_df$x1), by = 0.2),
      apply(sum_pred_1, 2, mean),
      col = col_vec[2],
      lwd = 2)
lines(
  seq(min(pop_df$x1), max(pop_df$x1), by = 0.2),
  apply(sum_pred_1, 2, quantile, 0.025),
  col = col_vec[2],
  lwd = 1,
  lty = "dashed"
)
lines(
  seq(min(pop_df$x1), max(pop_df$x1), by = 0.2),
  apply(sum_pred_1, 2, quantile, 0.975),
  col = col_vec[2],
  lwd = 1,
  lty = "dashed"
)
lines(seq(min(pop_df$x1), max(pop_df$x1), by = 0.2),
      apply(sum_pred_0, 2, mean),
      col = col_vec[1],
      lwd = 2)
lines(
  seq(min(pop_df$x1), max(pop_df$x1), by = 0.2),
  apply(sum_pred_0, 2, quantile, 0.025),
  col = col_vec[1],
  lwd = 1,
  lty = "dashed"
)
lines(
  seq(min(pop_df$x1), max(pop_df$x1), by = 0.2),
  apply(sum_pred_0, 2, quantile, 0.975),
  col = col_vec[1],
  lwd = 1,
  lty = "dashed"
)
text(10, 3, labels = "FEMALE = 1", col = col_vec[2], font = 2)
text(3, 10, labels = "FEMALE = 0", col = col_vec[1], font = 2)

par(mar = c(5,     # bottom (5)
            6.5,   # left   (4)
            2,     # top    (4)
            0.5    # right  (2)
),
xpd = T)

ylimits <- c(0, max(apply(cbind(sum_exp_pred_1, sum_exp_pred_0), 2, quantile, c(0.025, 0.5, 0.975))))

plot(
  pop_df$x1,
  exp(pop_df$y),
  col = adjustcolor(ifelse(pop_df$x2 == 0, col_vec[1], col_vec[2]), alpha = 0.5),
  las = 1,
  main = "",
  ylab = "",
  xlab = "",
  pch = 16,
  yaxt = "n",
  xaxt = "n",
  bty = "n",
  ylim = ylimits
)
box(lty = "solid", col = "grey50")

marks <- formatC(round(seq(0, max(apply(cbind(sum_exp_pred_1, sum_exp_pred_0), 2, quantile, c(0.025, 0.5, 0.975))), length.out = 5), 0), "d")
axis(2, las = 2, at = marks, labels = marks,
     col.axis = "grey50", col = "grey50")

axis(1,
     col.axis = "grey50", col = "grey50")

mtext(2,
      text = "INCOME", 
      line = 4.1)
mtext(3, 
      text = expression(bold("Original-Scale")), 
      line = 0.4)
mtext(1, 
      text = "EDUCATION",
      line = 3.5)
lines(
  seq(min(pop_df$x1), max(pop_df$x1), by = 0.2),
  apply(sum_exp_pred_1, 2, mean),
  col = col_vec[2],
  lwd = 2
)
lines(
  seq(min(pop_df$x1), max(pop_df$x1), by = 0.2),
  apply(sum_exp_pred_1, 2, quantile, 0.025),
  col = col_vec[2],
  lwd = 1,
  lty = "dashed"
)
lines(
  seq(min(pop_df$x1), max(pop_df$x1), by = 0.2),
  apply(sum_exp_pred_1, 2, quantile, 0.975),
  col = col_vec[2],
  lwd = 1,
  lty = "dashed"
)
lines(
  seq(min(pop_df$x1), max(pop_df$x1), by = 0.2),
  apply(sum_exp_pred_0, 2, mean),
  col = col_vec[1],
  lwd = 2
)
lines(
  seq(min(pop_df$x1), max(pop_df$x1), by = 0.2),
  apply(sum_exp_pred_0, 2, quantile, 0.025),
  col = col_vec[1],
  lwd = 1,
  lty = "dashed"
)
lines(
  seq(min(pop_df$x1), max(pop_df$x1), by = 0.2),
  apply(sum_exp_pred_0, 2, quantile, 0.975),
  col = col_vec[1],
  lwd = 1,
  lty = "dashed"
)


fd <- sum_pred_0 - sum_pred_1

plot(
  seq(min(pop_df$x1), max(pop_df$x1), by = 0.2),
  apply(fd, 2, mean),
  col = col_vec[1],
  las = 1,
  ylim = c(0, 4),
  lwd = 2, 
  main = "",
  ylab = "",
  xlab = "",
  pch = 16,
  yaxt = "n",
  xaxt = "n",
  bty = "n",
  type = "l"
)
box(lty = "solid", col = "grey50")

axis(2, las = 2, 
     col.axis = "grey50", col = "grey50")

axis(1,
     col.axis = "grey50", col = "grey50")

mtext(2,
      text = "E(log(INCOME)|FEMALE = 0) - \n  E(log(INCOME)|FEMALE = 1)", 
      line = 4.1)
mtext(3, 
      text = expression(bold("First Difference on Log-Scale")), 
      line = 0.4)
mtext(1, 
      text = "EDUCATION",
      line = 3.5)

lines(
  seq(min(pop_df$x1), max(pop_df$x1), by = 0.2),
  apply(fd, 2, quantile, 0.025),
  col = col_vec[1],
  lwd = 1,
  lty = "dashed"
)
lines(
  seq(min(pop_df$x1), max(pop_df$x1), by = 0.2),
  apply(fd, 2, quantile, 0.975),
  col = col_vec[1],
  lwd = 1,
  lty = "dashed"
)


exp_fd <- sum_exp_pred_0 - sum_exp_pred_1
ylimits <- c(0, max(apply(exp_fd, 2, quantile, c(0.025, 0.5, 0.975))))
plot(
  seq(min(pop_df$x1), max(pop_df$x1), by = 0.2),
  apply(exp_fd, 2, mean),
  col = col_vec[1],
  las = 1,
  lwd = 2, 
  main = "",
  ylab = "",
  xlab = "",
  pch = 16,
  yaxt = "n",
  xaxt = "n",
  bty = "n",
  type = "l",
  ylim = ylimits
)
box(lty = "solid", col = "grey50")

axis(2, las = 2, 
     col.axis = "grey50", col = "grey50")

axis(1,
     col.axis = "grey50", col = "grey50")

mtext(2,
      text = "E(INCOME|FEMALE = 0) - \n  E(INCOME|FEMALE = 1)", 
      line = 4.1)
mtext(3, 
      text = expression(bold("First Difference on Original-Scale")), 
      line = 0.4)
mtext(1, 
      text = "EDUCATION",
      line = 3.5)

lines(
  seq(min(pop_df$x1), max(pop_df$x1), by = 0.2),
  apply(exp_fd, 2, quantile, 0.025),
  col = col_vec[1],
  lwd = 1,
  lty = "dashed"
)
lines(
  seq(min(pop_df$x1), max(pop_df$x1), by = 0.2),
  apply(exp_fd, 2, quantile, 0.975),
  col = col_vec[1],
  lwd = 1,
  lty = "dashed"
)
dev.copy(pdf, "figures/dummy_example.pdf", width = 16, height = 9)
dev.off()

