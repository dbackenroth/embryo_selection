library(data.table)
library(dplyr)
library(mvtnorm)

Simulation <- function(x_sigma, c_sigma, e_sigma, n, thresholds) {
  n_diseases <- nrow(x_sigma)
  x <- rmvnorm(n = n, mean = rep(0, n_diseases), sigma = x_sigma)
  e <- rmvnorm(n = n, mean = rep(0, n_diseases), sigma = e_sigma)
  sum <- x + e
  s <- matrix(NA, nrow = n, ncol = n_diseases)
  for (i in 1:n_diseases) {
    s[, i] <- sum[, i] > thresholds[i]
  }
  p <- numeric(length = n)
  for (j in 1:nrow(x)) {
    p[j] <- 1 - pmvnorm(lower = rep(-Inf, n_diseases), upper = thresholds,
                                    mean = x[j, ], 
                                    sigma = e_sigma)
    
  }
  # x are PRSs
  # e are noise effects
  # we assume c (shared genetic term) = 0
  # s is disease status
  # p is probability of one or more diseases
  return(list(x = x, e = e, s = s, p = p))
}

EvaluateBySim <- function(x_sigma, e_sigma, n, prev1, prev2) {
  sims <- Simulation(x_sigma = x_sigma, c_sigma = diag(0, 2), e_sigma = e_sigma, n = n, 
                     thresholds = c(qnorm(1 - prev1), qnorm(1 - prev2)))
  
  p <- sims$p
  s <- sims$s
  has_any_disease <- s[, 1] | s[, 2]
  
  res <- map_dfr(1:15, function(f_size) {
    families <- rep(1:(n %/% f_size), each = f_size)
    
    df <- data.frame(p = p, d = has_any_disease) %>%
      slice(1:length(families))
    df$f <- families
    summ <- df %>%
      dplyr::group_by(f) %>%
      dplyr::summarise(random_has_disease = d[1], 
                       best_has_disease = d[which.min(p)]) %>%
      mutate(f_size = f_size)
  }) 
  res2 <- res %>%
    dplyr::group_by(f_size) %>%
    dplyr::summarise(r = mean(random_has_disease), 
                     b = mean(best_has_disease))
  plot(res2$f_size, 1 - res2$b / mean(res2$r), type = 'l', xlab = "# embryos", ylab = "% risk reduction (of 1 or 2 diseases)", 
       ylim = c(0, 1))
  return(res2)
}

# assume x, c and e are independent
r1var <- 0.2
r2var <- 0.2
prev1 <- 0.01
prev2 <- 0.01
rho <- 0.05
n <- 250000
x_sigma <- matrix(c(r1var / 2, sqrt(r1var* r2var / 4) * rho, 
                    sqrt(r1var * r2var / 4) * rho, r2var / 2), nrow = 2)
c_sigma <- x_sigma
e_sigma <- diag(c(1 - r1var, 1 - r2var))

sim_res <- EvaluateBySim(x_sigma = x_sigma, e_sigma = e_sigma, n = n, prev1 = prev1, prev2 = prev2)

EvaluateByIntegration <- function(x_sigma = x_sigma, e_sigma = e_sigma, prev1 = prev1, prev2 = prev2) {
  
  prob_finder <- function(mean_1, mean_2, t_1, t_2, sigma){
      1 - pmvnorm(lower = c(-Inf, -Inf), upper = c(t_1, t_2), 
              mean = c(mean_1, mean_2), 
              sigma = sigma)
  }
  grid_size <- 0.05
  s <- seq(-5, 5, by = grid_size)
  grid <- expand.grid(mean_1 = s, mean_2 = s)
  
  r <- pmap(grid, prob_finder, t_1 = qnorm(1 - prev1), t_2 = qnorm(1 - prev2), sigma = e_sigma) %>% unlist()
  print("finished setting up grid")
  grid$prob <- r
  grid <- as.data.table(grid)
  res <- map_dfr(c(1:15), function(f_size) {
    integrand <- function(p, grid, grid_size, f_size) {
      map(p, function(x) {
        this_grid <- grid[prob >= x, ]
        densities <- dmvnorm(as.matrix(this_grid[, .(mean_1, mean_2)]), mean = c(0, 0), sigma = x_sigma) 
        ans <- (sum(densities) * grid_size ^ 2) ^ f_size
      }) %>% 
        unlist()
    }
    r <- integrate(integrand, lower = 0, upper = 1, grid = grid, grid_size = grid_size, f_size = f_size)
    data.frame(f_size = f_size, int = r$value)
  })
  return(res)
}

int_res <- EvaluateByIntegration(x_sigma = x_sigma, e_sigma = e_sigma, prev1 = prev1, prev2 = prev1)
lines(int_res$f_size, 1 - int_res$int / max(int_res$int), col = "red")
legend("bottomright", c("simulation", "integration"), col = c("black", "red"), lwd = c(2,2), bty = 'n')
# 
# 
# rho <- 0.05
# r1var <- 0.15
# r2var <- 0.15
# e1var <- 0.85
# e2var <- 0.85
# 
# m <- diag(c(r1var / 2, r2var / 2, r1var / 2, r2var / 2, e1var, e2var))
# m[1, 2] <- m[2, 1] <- m[3, 4] <- m[4, 3] <- sqrt(r1var / 2) * sqrt(r2var / 2) * rho
# y_mat <- matrix(c(1, 0, 1, 0, 1, 0, 
#                   0, 1, 0, 1, 0, 1), nrow = 2, byrow = T)
# y_mean <- rep(0, 6)
# y_var <- y_mat %*% m %*% t(y_mat)