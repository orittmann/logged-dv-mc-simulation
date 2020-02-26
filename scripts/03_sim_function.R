# betas
# vcov
# scenario
# (stochastic component)
# (parameters for sc)
# methods (ova, aca, median?, vector, matrix)

pred_values <-
  function(model_theta,
           model_vcov,
           scenario,
           M = 1000,
           nsim = 1000,
           model_sigma) {
    require(MASS)
    cat("\n Simulating Predicted Values: \n")
    if (is.null(nrow(scenario))) {
      sim_array <- array(NA, dim = c(nsim, 1, M))
    } else {
      sim_array <- array(NA, dim = c(nsim, nrow(scenario), M))
    }
    S <-  mvrnorm(nsim, model_theta, model_vcov)
    if (is.null(nrow(scenario))) {
      scenario <- t(scenario)
    }
    Xb <- S %*% t(scenario)
    
    pb <- txtProgressBar(min = 0, max = M, style = 3)
    for (m in 1:M) {
      # Predicted Values
      
      fundamental_uncertainty <-
        matrix(rnorm(prod(dim(Xb)), 0, sd = model_sigma),
               nrow = dim(Xb)[1],
               ncol = dim(Xb)[2])
      
      pred_Xb <- Xb + fundamental_uncertainty
      
      sim_array[, , m] <- pred_Xb
      
      setTxtProgressBar(pb, m)
    }
    close(pb)
    return(sim_array)
  }


transform_predicted_values <- function(obj, invlink = exp) {
  invlink(obj)
}

summarize_sim_values <- function(obj, qi = mean) {
  apply(obj, c(1, 2), qi)
}