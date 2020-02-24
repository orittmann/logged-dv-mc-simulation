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


#pdf("figures/error_exp_theta.pdf")
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
#dev.off()


# variation in sigma^2, holding constant theta

sigma2 <- seq(0, 5, length.out = 100)
theta_fix <- 1

bias_data2 <- data.frame("sigma2" = sigma2,
                         "mean_Y" = NA,
                         "median_Y" = NA)



for (i in 1:nrow(bias_data2)) {
  mean_median <- mean_median_Y(theta = theta_fix, sigma_est = bias_data2$sigma2[i], m = 1000000)
  bias_data2$mean_Y[i] <- mean_median[1]
  bias_data2$median_Y[i] <- mean_median[2]
}

bias_data2$dif <- bias_data2$mean_Y - rep(exp(2), 100)

# pdf("figures/error_sigma2.pdf")
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
# dev.off()