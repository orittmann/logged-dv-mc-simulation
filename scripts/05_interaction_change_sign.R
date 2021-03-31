
library(viridis)

b0 <- 0
b1 <- 2
b2 <- 3
b3 <- -1
beta <- c(b0, b1, b2, b3)

n <- 1000

x <- rnorm(n, 0, 1)
#x <- runif(n, -3, 3)
z <- rbinom(n, 1, 0.5)

X <- cbind(1, x, z, x*z)

sigma2 <- 1 
lny <- X %*% beta + rnorm(n, 0, sqrt(sigma2))

y <- exp(lny)


plot(x, lny,
     col = ifelse(z == 1,
                  viridis(2, 0.5)[1],
                  viridis(2, 0.5)[2]),
     pch = 19)
xseq <- seq(-3, 3, 0.1)
y_pred_0 <- (b0 + b1*xseq) 
y_pred_1 <- (b0 + b1*xseq + b2*1 + b3*1*xseq)

lines(xseq, y_pred_0, 
      col = viridis(2, 1)[2],
      lwd = 2)
lines(xseq, y_pred_1, 
      col = viridis(2, 1)[1],
      lwd = 2)


plot(x, y,
     col = ifelse(z == 1,
                  viridis(2, 0.5)[1],
                  viridis(2, 0.5)[2]),
     pch = 19)

# andy's approach
t <- (1 / (b1+b3)) - (b2/b3)
t <- (log(b1) - log(b1+b3) - b2) / b3



# E(y)
var_exp <- (sigma2 / (n * sum(x^2) - sum(x)^2)) * 
  (sum(x^2) - 2*xseq * sum(x) + xseq^2*n)


xseq <- seq(-3, 3, 0.1)
y_pred_0 <- exp(b0 + b1*xseq) * exp(sigma2/2) * exp(var_exp/2)
y_pred_1 <- exp(b0 + b1*xseq + b2*1 + b3*1*xseq) * exp(sigma2/2) * exp(var_exp/2)

y_pred_0 <- exp(b0 + b1*xseq + sigma2/2 + var_exp/2)
y_pred_1 <- exp(b0 + b1*xseq + b2*1 + b3*1*xseq + sigma2/2 + var_exp/2)


lines(xseq, y_pred_0, 
      col = viridis(2, 1)[2],
      lwd = 2)
lines(xseq, y_pred_1, 
      col = viridis(2, 1)[1],
      lwd = 2)
abline(v = t)



# derivative with respect to z (equivalent to first difference) (log-scale)
xseq <- seq(-3, 3, 0.1)
y_pred_0 <- b0 + b1*xseq + b2*0 + b3*0*xseq
y_pred_1 <- b0 + b1*xseq + b2*1 + b3*1*xseq

derivative_z <- b2 + b3*xseq

plot(xseq, derivative_z,
     type = "l")
lines(x = xseq, 
      y = y_pred_1 - y_pred_0, 
      type = "l", col = "red")

# derivative with respect to z (equivalent to first difference)
xseq <- seq(-3, 3, 0.1)

a <- sigma2 / (n * sum(x^2) - sum(x)^2)
sigma2_sampl <- a * (sum(x^2) - 2*xseq * sum(x) + xseq^2*n)

y_pred_0 <- exp(b0 + b1*xseq + b2*0 + b3*0*xseq + sigma2/2 + sigma2_sampl/2)
y_pred_1 <- exp(b0 + b1*xseq + b2*1 + b3*1*xseq + sigma2/2 + sigma2_sampl/2)

derivative_z <- exp(b0 + b1*xseq + b2*0 + b3*0*xseq + sigma2/2 + sigma2_sampl/2) *
  (b2 + b3*xseq + 
     (a * (2*n*xseq - 2*sum(x)) / 2))

derivative_z0 <- exp(b0 + b1*xseq + b2*0 + b3*0*xseq + sigma2/2 + sigma2_sampl/2) *
  (b2 + b3*xseq + 
     (a * (n*xseq - sum(x))))
derivative_z1 <- exp(b0 + b1*xseq + b2*1 + b3*1*xseq + sigma2/2 + sigma2_sampl/2) *
  (b2 + b3*xseq + 
     (a * (n*xseq - sum(x))))
derivative_z05 <- exp(b0 + b1*xseq + b2*0.5 + b3*0.5*xseq + sigma2/2 + sigma2_sampl/2) *
  (b2 + b3*xseq + 
     (a * (n*xseq - sum(x))))


plot(xseq, derivative_z0,
     type = "l",
     ylim = c(0, 300))
lines(xseq, derivative_z1,
     type = "l")
lines(xseq, derivative_z05,
      type = "l")
lines(x = xseq, 
      y = y_pred_1 - y_pred_0, 
      type = "l", col = "red")
t <- (1 / (b1+b3)) - (b2/b3)
t <- (log(b1) - log(b1+b3) - b2) / b3
abline(v = t)








# derivative with respect to z and x (equivalent to b3)
derivative_zx0 <- 
  exp(b0 + b1*xseq + b2*0 + b3*0*xseq + sigma2/2 + sigma2_sampl/2) * 
  (b3*0 + b1) * (b2 + b3*xseq + a*(n*xseq - sum(x))) *
  exp(b0 + b1*xseq + b2*0 + b3*0*xseq + sigma2/2 + sigma2_sampl/2) *
  (b3 + a*n)

derivative_zx1 <- 
  exp(b0 + b1*xseq + b2*1 + b3*1*xseq + sigma2/2 + sigma2_sampl/2) * 
  (b3*1 + b1) * (b2 + b3*xseq + a*(n*xseq - sum(x))) *
  exp(b0 + b1*xseq + b2*1 + b3*1*xseq + sigma2/2 + sigma2_sampl/2) *
  (b3 + a*n)

derivative_zx05 <- 
  exp(b0 + b1*xseq + b2*0.5 + b3*0.5*xseq + sigma2/2 + sigma2_sampl/2) * 
  (b3*0.5 + b1) * (b2 + b3*xseq + a*(n*xseq - sum(x))) *
  exp(b0 + b1*xseq + b2*0.5 + b3*0.5*xseq + sigma2/2 + sigma2_sampl/2) *
  (b3 + a*n)


plot(xseq, derivative_zx05,
     type = "p",
     pch = 19,
     cex = 0.7,
     ylim = c(-100000,500000),
     xlim = c(-3,3),
     col = ifelse(derivative_zx0 < 0, "red", "blue"))
lines(xseq, derivative_zx05)

t <- (1 / (b1+b3)) - (b2/b3)
t <- (log(b1) - log(b1+b3) - b2) / b3
abline(v = t)




b0 <- 0
b1 <- 2
b2 <- 3
b3 <- -1
beta <- c(b0, b1, b2, b3)
n <- 100000

#x <- rnorm(n, 0, 1)
x <- runif(n, -3, 3)
z <- rbinom(n, 1, 0.5)

X <- cbind(1, x, z, x*z)

sigma2 <- 1 
lny <- X %*% beta + rnorm(n, 0, sqrt(sigma2))

y <- exp(lny)

plot(x = x,
     y = y,
     type = "n",
     ylim = c(0,1000))
lines(ksmooth(x = x[z == 0],
              y = y[z == 0],
              "normal",
              bandwidth = 0.5),
      col = "red", lwd = 2)
lines(ksmooth(x = x[z == 1],
              y = y[z == 1],
              "normal",
              bandwidth = 0.5),
      col = "red", lwd = 2)
lines(xseq, 
      y_pred_0, type = "l", col = "blue")

lines(xseq, 
      y_pred_1, type = "l", col = "blue")


y_pred_01 <- exp(b0 + b1*xseq + b2*0 + b3*0*xseq + sigma2/2)
y_pred_11 <- exp(b0 + b1*xseq + b2*1 + b3*1*xseq + sigma2/2)
lines(xseq, 
      y_pred_01, type = "l", col = "orange")

lines(xseq, 
      y_pred_11, type = "l", col = "orange")




der_logy <- b3
xseq <- seq(1, 2.2, 0.1)
der_y <- exp(b0 + b1*xseq + b2*0 + b3*0*xseq + sigma2/2)*(b3*0 + b1) +
  exp(b0 + b1*xseq + b2*0 + b3*0*xseq + sigma2/2)*b3

plot(xseq, der_y, type = "l")
abline(h = der_logy)






