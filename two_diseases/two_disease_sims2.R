library(mvtnorm)

Simulate <- function(hsnp1sq = 0.05, v1sq = 0.01, h1sq = 0.3, 
                     hsnp2sq = 0.05, v2sq = 0.01, h2sq = 0.3, 
                     rho = 0.05, n = 1000, plot = T) {
  
  Amu <- matrix(0, nrow = 8, ncol = 1)
  dimnames <- c("ps1", "ps2", "s1", "s2", "g1p", "g2p", "e1", "e2")
  Asig <- matrix(0, nrow = 8, ncol = 8, dimnames = list(dimnames, dimnames))
  Asig[1, 1] <- hsnp1sq + v1sq
  Asig[2, 2] <- hsnp2sq + v2sq
  Asig[1, 2] <- rho * sqrt(hsnp1sq * hsnp2sq)
  Asig[3, 3] <- Asig[1, 3] <- v1sq
  Asig[4, 4] <- Asig[2, 4] <- v2sq
  Asig[5, 5] <- h1sq - hsnp1sq
  Asig[6, 6] <- h2sq - hsnp2sq
  Asig[5, 6] <- rho * sqrt((h1sq - hsnp1sq) * (h2sq - hsnp2sq))
  Asig[7, 7] <- 1 - h1sq
  Asig[8, 8] <- 1 - h2sq
  Asig[lower.tri(Asig)] = t(Asig)[lower.tri(Asig)]
  
  # get joint distribution of y1, ps1, ps2
  B <- matrix(c(1, 0, -1, 0, 1, 0, 1, 0,    # y1
                1, 0, 0, 0, 0, 0, 0, 0,     # ps1
                0, 1, 0, 0, 0, 0, 0, 0),    # ps2
              nrow = 3, byrow = T)
  Bsig <- B %*% Asig %*% t(B)
  Bmu <- B %*% Amu
  
  #Ymu_cond_both <- Bmu[1, ] + 
  Y1sig_cond <- Bsig[1, 1] - Bsig[1, 2:3, drop = F] %*% solve(Bsig[2:3, 2:3]) %*% Bsig[2:3, 1, drop = F]
  Y1mu_premult <- Bsig[1, 2:3, drop = F] %*% solve(Bsig[2:3, 2:3])
  k1 <- 0.01
  k2 <- 0.02
  ps <- rmvnorm(n, mean = Amu[1:2, 1], sigma = Asig[1:2, 1:2])
  
  thresh <- qnorm(1 - k1)
  p1 <- 1 - pnorm((thresh - ps[, 1]) / sqrt(1 - Asig[1, 1]))
  p2 <- 1 - pnorm((thresh - as.numeric(ps %*% t(Y1mu_premult))) / sqrt(as.numeric(Y1sig_cond)))
  
  summary(abs(p2 - p1) / p2) %>%
    print()
  res <- data.frame(p1 = p1, p2 = p2)
  if (plot) {
    print(ggplot(res, aes(x = p1, y = p2)) + 
      geom_point() + 
      geom_abline() + 
      theme_bw() + 
      xlab("prob. trait 1 based on PS1") + 
      ylab("prob. trait 1 based on PS1 and PS2"))
  }
  print(mean(p1))
  print(mean(p2))
}
