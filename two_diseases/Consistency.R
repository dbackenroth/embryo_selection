source("two_disease_sims.R")

n_family <- 100
r1 <- 0.05
r2 <- 0.05
rho <- 0.03
n_embryo <- 5
Do <- function(seed) {
  set.seed(seed)
  ll <- SimEmbryosAndParents(n_family = n_family, r1 = r1, r2 = r2, rho = rho, n_embryo = n_embryo)
  children <- ll$embryo_ind
  parents <- ll$parents_raw
  res <- matrix(NA, nrow = 100, ncol = 100)
  for (i in 1:100) {
    ch <- children + parents[i, ]
    trait1_p <- ProbTraitGivenOnePS(ch[, 1], K = 0.01, R = r1)
    trait2_p <- ProbTraitGivenOnePS(ch[, 2], K = 0.01, R = r2)
    ss <- trait1_p + trait2_p
    dd <- data.frame(ss = ss, f = rep(1:100, each = 5))
    df <- group_by(dd, f) %>%
      summarise(m = which.min(ss))
    res[, i] <- df$m
  }
  browser()
}