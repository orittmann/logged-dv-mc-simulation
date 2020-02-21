source("scripts/00_setup.R")

# population size
m <- 10000

# means
mus <- c(5, 10)

# SD
sds <- c(4, 3)
sds_mat <- diag(sds)

# Correlation matrix (we assume some multicollinearity)
cor_mat <- matrix(c(1, 0.2, 0.2, 1), nrow = 2, ncol = 2)

# Convert to variance covariance matrix
varcov <- sds_mat %*% cor_mat %*% sds_mat 

# true coefficients
b <- c(2, 0.3, -0.5)

# Choose sigma.est
sigma_est <- 1

gen_pop <- function(m,          # population size
                    mus,        # means
                    varcov,     # variance-covariance matrix
                    b,          # true coefficients
                    sigma_est   # sigma est
                    ){
  X <- mvrnorm(n = m, mu = mus, Sigma = varcov) 
  mu_y <- cbind(1, X) %*% b
  Y <- exp(rnorm(m, mu_y, sigma_est))
  population <- data.frame("Y" = Y,
                           "X1" = X[,1],
                           "X2" = X[,2])
  return(population)
}

pop <- gen_pop(m = m,
               mus = mus, 
               varcov = varcov, 
               b = b, 
               sigma_est = sigma_est)

dta <- pop[sample(1:m, size = 1000),]

# Plot the data

plot(dta$X1, log(dta$Y))
plot(dta$X1, dta$Y)

plot(dta$X2, log(dta$Y))
plot(dta$X2, dta$Y)

