library(mvtnorm)
library(dplyr)
library(purrr)

Test1 <- function() {
  k1 <- 0.01
  k2 <- 0.02
  r1 <- sqrt(0.05)
  r2 <- sqrt(0.05)
  rho <- 0.1
  s <- SimEmbryosAndParents(n_family = 1000000, r1 = r1, r2 = r2, rho = rho, n_embryo = 1)
  p1 <- ProbTraitGivenOnePS(s$children$ps2, k = k2, r = r2)
  # we have R^2 = h^4 / (h^2 + v^2)
  # if the error variance is m percent of h^2
  # then R^2 = h^4 / ((1 + m)h^2)
  # so h^2 = R^2 (1 + m)
  m1 <- 0.1
  m2 <- 0.1
  hsnp1 <- R1 * sqrt(1 + m1)
  hsnp2 <- R2 * sqrt(1 + m2)
  p2 <- ProbTrait2GivenPS1and2(ps1 = s$children$ps1, ps2 = s$children$ps2, hsnp1 = hsnp1, 
                               hsnp2 = hsnp2, r1 = r1, r2 = r2, k2 = k2)
  browser()
}

ProbTraitGivenOnePS <- function(ps, k, r) {
  thresh <- qnorm(1 - k)
  1 - pnorm((thresh - ps) / sqrt(1 - r^2))
}

ProbTrait2GivenPS1and2 <- function(ps1, ps2, hsnp1, hsnp2, r1, r2, k2) {
  thresh <- qnorm(1 - k2)
  e <- ExpectedValueTrait2(hsnp1 = hsnp1, hsnp2 = hsnp2, r1 = r1, r2 = r2, ps1 = ps1, ps2 = ps2)
  v <- VarianceTrait2(r1 = r1, r2 = r2, rho = rho, hsnp1 = hsnp1, hsnp2 = hsnp2)
  1 - pnorm((thresh - e) / sqrt(v))
}

ExpectedValueTrait2 <- function(hsnp1, hsnp2, r1, r2, ps1, ps2) {
  # hsnp1, hsnp2: square roots of SNP-based heritabilities
  # r1, r2: square root of proportion of variance in liability for disease
  #         explained by score for that disease
  # ps1, ps2: scaled propensity score (with variance equal to r1, r2)
  num <- hsnp1 * hsnp2 * (hsnp2 ^ 2 - r2 ^ 2) * ps1 + hsnp2 ^ 2 * (hsnp1 ^ 2 - r1 ^ 2) * ps2
  den <- hsnp1 ^ 2 * hsnp2 ^ 2 - r1 ^ 2 * r2 ^ 2
  num / den
}

VarianceTrait2 <- function(r1, r2, rho, hsnp1, hsnp2) {
  num <- rho ^ 2 * r1 ^ 2 * (hsnp2 ^ 2 - r2 ^ 2) + r2 ^ 2 * (hsnp1 ^ 2 - rho ^ 2 * r1 ^ 2)
  den <- hsnp1 ^ 2 * hsnp2 ^ 2 - rho ^ 2 * r1 ^ 2 * r2 ^ 2
  hsnp2 ^ 2 * (1 - num / den) + 1 - hsnp2 ^ 2
}

SimEmbryosAndParents <- function(n_family, r1, r2, rho, n_embryo) {
  # rho: correlation between scores of two diseases
  # r1, r2: square root of proportion of variance in liability for disease
  #         explained by score for that disease
  mean <- c(0, 0)
  off_diag <- rho * r1 * r2
  sigma <- 1/2 * matrix(c(r1^2,     off_diag, 
                          off_diag, r2^2), nrow = 2)
  # for each rep we need one set of parents and n_embryo embryos
  n <- n_family * (1 + n_embryo)
  variates <- rmvnorm(n, mean = mean, sigma = sigma)
  parent_rows <- seq(1, n, by = n_embryo + 1)
  other_rows <- setdiff(1:n, parent_rows)
  final <- variates[other_rows, ] + variates[rep(parent_rows, each = n_embryo), ]
  family_num <- rep(1:n_family, each = n_embryo)
  embryo_num <- rep(1:n_embryo, times = n_family)
  children <- as.data.frame(final) %>%
    set_names(c("ps1", "ps2")) %>%
    mutate(family = family_num, 
           embryo = embryo_num)
  parents <- variates[parent_rows, ] %>%
    as.data.frame() %>%
    set_names(c("ps1", "p2")) %>%
    mutate(family = 1:n_family)
  list(children = children, parents = parents, 
       embryo_ind = variates[other_rows, ], 
       parents_raw = variates[parent_rows, ])
}