library(purrr)
library(Hmisc)

source("explore_score_variance.R")
source("helpers.R")
source("score_analysis.R")

GetRiskReduction <- function(dat, strategy = "lowest_risk", n = NULL, 
                             q_exclude = NULL, threshold = NULL, prevalence, ...) {
  first_child <- dat[cnum == 1, .(risk, p1, p2, p1_pheno, p2_pheno)]
  first_child_by_parent_status <- first_child[, .(mean_risk = mean(risk)),
                                              by = .(p1_pheno, p2_pheno)]
  if (strategy == "lowest_risk") {
    selected_child <- dat[cnum <= n, .(risk = min(risk)), 
                          by = .(p1, p2, p1_pheno, p2_pheno)]
  } else if (strategy == "hre") {
    dat[, `:=`(above_threshold = score >= threshold)]
    selected_child <- dat[cnum <= 5, .(all_above = all(above_threshold),
                                       risk = if_else(all(above_threshold), risk[1], 
                                                      risk[!above_threshold][1])), 
                          by = .(p1, p2, p1_pheno, p2_pheno)]
  }
  strategy_risk <- selected_child[, .(mean_risk = mean(risk)), 
                                  by = .(p1_pheno, p2_pheno)]
  
  probs <- c((1 - prevalence)^2, 2 * (1 - prevalence) * prevalence, prevalence^2)
  #child1_by_parent_status$prob_couple <- probs
  strategy_risk_ave <- sum(strategy_risk$mean_risk * probs)
  first_child_risk_ave <- sum(first_child_by_parent_status$mean_risk * probs)
  
  #RRR <- 1 - sum(strategy_risk$mean_risk * probs / first_child_by_parent_status$mean_risk)
  RRR <- 1 - strategy_risk_ave / first_child_risk_ave
  ret <- data.frame(RRR = RRR) %>%
    mutate(n = n, strategy = strategy, q_exclude = q_exclude, threshold = threshold)
  return(ret)
}

if (F) {
  schiz_r2 <- SchizR2()
  ScoreAnalysis(disease = SCHIZ, prevalence = SCHIZ_PREVALENCE, r2_liability = schiz_r2, norm_approximation = T)
  crohns_r2 <- CrohnsR2()
  ScoreAnalysis(disease = CROHNS, prevalence = CROHNS_PREVALENCE, r2_liability = crohns_r2, norm_approximation = T)
}

ScoreAnalysis <- function(disease, prevalence, r2_liability, norm_approximation) {
  dat <- get_data(disease = disease, prevalence = prevalence)

  sim <- dat$sim_scores
  parents <- dat$p_info %>%
    as.data.frame() %>%
    mutate(prev = if_else(pheno == 0, 1 - prevalence, prevalence)) %>%
    group_by(pheno) %>%
    mutate(weight = prev / n())
  
  wmean <- wtd.mean(x = parents$score, weights = parents$weight, normwt = T)
  wvar <- wtd.var(x = parents$score, weights = parents$weight, normwt = T)
  
  quantiles <- data.frame(prob = seq(0.6, 0.99, by = 0.01))
  quantiles$q_exclude <- 1 - quantiles$prob
  if (norm_approximation) {
    quantiles$q_val <- qnorm(quantiles$prob, mean = wmean, sd = sqrt(wvar))
  } else {
    quantiles$q_val <- wtd.quantile(x = parents$score, prob = quantiles$prob, 
                                    weights = parents$weight, normwt = T)
  }
  #type=c('quantile','(i-1)/(n-1)','i/(n+1)','i/n')
  
  grid_hre <- data.frame(strategy = rep("hre", nrow(quantiles)), 
                         q_exclude = quantiles$q_exclude, 
                         threshold = quantiles$q_val)
  res_hre <- pmap_dfr(grid_hre, GetRiskReduction, dat = sim, prevalence = prevalence)
  res_hre$theoretical <- NA
  for (i in 1:nrow(res_hre)) {
    res_hre$theoretical[i] <- risk_reduction_exclude(r = sqrt(r2_liability), K = prevalence, n = 5, 
                                                     q = res_hre$q_exclude[i])
  }
  
  grid_lr <- data.frame(strategy = rep("lowest_risk", 20), 
                        n = 1:20)
  res_lr <- pmap_dfr(grid_lr, GetRiskReduction, dat = sim, prevalence = prevalence)
  res_lr$theoretical <- NA
  for (i in 1:nrow(res_lr)) {
    res_lr$theoretical[i] <- risk_reduction_lowest(r = sqrt(r2_liability), K = prevalence, n = res_lr$n[i])
  }
  plot(res_lr$n, res_lr$RRR, type = 'l', ylim = c(0, 1))
  lines(res_lr$n, res_lr$theoretical, col = "blue")
  
  p_hre <- ggplot(res_hre %>% 
                    mutate(RRR = RRR * 100, 
                           q_exclude = q_exclude * 100,
                           theoretical = theoretical * 100),
                  aes(x = q_exclude, y = RRR)) + 
    geom_point() + 
    geom_line(aes(x = q_exclude, y = theoretical)) + 
    theme_bw() + 
    ylab("Risk reduction (%)") + 
    xlab("Percentile PRS to exclude")
  ggsave(glue("Results/{disease}_hre_r2_{r2_liability}_norm_approx{norm_approximation}.pdf"), height = 3, width = 3)
  p_lr <- ggplot(res_lr %>% 
                   mutate(RRR = RRR * 100, 
                          theoretical = theoretical * 100),
                 aes(x = n, y = RRR)) + 
    geom_point() + 
    geom_line(aes(x = n, y = theoretical)) + 
    theme_bw() + 
    ylab("Risk reduction (%)") + 
    xlab("Number of embryos")
  ggsave(glue("Results/{disease}_lr_r2_{r2_liability}_norm_approx{norm_approximation}.pdf"), height = 3, width = 3)
}
#best_of_5 <- GetRiskReduction(dat, strategy = "lowest_risk", n = 5)
#best_of_5 <- GetRisk(dat, )


if (F) {
  #rr_a <- 1 / (1 + exp(-rr))
  #rr2 <- predict(mod, newdata = sim_scores, type = "response")
  
  child1 <- sim_scores[cnum == 1, .(risk, p1, p2, p1_pheno, p2_pheno)]
  child1_by_parent_status <- child1[, .(mean_risk = mean(risk)), 
                                    by = .(p1_pheno, p2_pheno)]
  best_child_of_5 <- sim_scores[cnum <= 5, .(risk = min(risk)), 
                                by = .(p1, p2, p1_pheno, p2_pheno)]
  best_child_of_5_by_parent_status <- best_child_of_5[, .(mean_risk = mean(risk)), 
                                                      by = .(p1_pheno, p2_pheno)]
  probs <- c(0.99^2, NA, 0.01^2)
  probs[2] <- 1 - sum(probs, na.rm = T)
  #child1_by_parent_status$prob_couple <- probs
  RR <- 1 - sum(best_child_of_5_by_parent_status$mean_risk * probs / child1_by_parent_status$mean_risk)
  
  risk_reduction_lowest(r = 0.0583, K = 0.01, n = 5)
}
# adjust logistic regression model
#   https://core.ac.uk/download/pdf/61320341.pdf
# add q / (1 - q) to intercept


#5) I think we would need to simulate between 1k-10k families. 
#   I think we should simulate two scenarios. 
#    (1) Random parents. Here we need to sample controls with 
#        probability 99%, and cases with probability 1% (the prevalence). 
#        The disease status is in the *fam file. 
#    (2) Conditional on parental disease status. Here we can simply 
#        choose pairs where either both are healthy (will probably 
#        be covered by part (1)), or one healthy and one sick, 
#        or both sick.
#6) In each family, we compute the PRS for each child. 
#   Then use the logistic regression model learned in (2) 
#   to compute the risk of the child. (I guess using "average" sex 
#   if sex was a covariate).
#7) Then report: (A) the average risk across families when using the 
#   first offspring in each family. (B) The average risk across 
#   families when selecting the offspring based on the score 
#   (according to either of the two strategies). I think we 
#   can use the same simulated families for all experiments. 
#   Basically once you have the families, all information you 
#   need to save is the disease status of the parents, and the 
#   scores of the parents and all children.
#8) Then we can compare to the theory. There's a bit of thinking to 
#   do on how to compute the variance explained on the liability 
#   scale, though we know the rough number (around 10%).