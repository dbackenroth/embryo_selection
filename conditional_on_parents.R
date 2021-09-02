source("explore_score_variance.R")
source("helpers.R")
source("score_analysis.R")

ConditionalAnalysis <- function() {
  cc <- map_dfr(as.list(0:2), function(x){
    res <- ScoreAnalysis(disease = SCHIZ, prevalence = SCHIZ_PREVALENCE, r2_liability = SchizR2(SCHIZ_PREVALENCE), 
                norm_approximation = T, num_couples = NULL, n_parents_diseased = x)
    both <- bind_rows(res$res_lr, res$res_hre) %>%
      mutate(num_parents_diseased = x)
  }) %>%
    mutate(pschiz = paste0("Parents with Schiz.=", num_parents_diseased))
  
  ggplot(cc %>% filter(strategy == "lowest_risk"), 
         aes(x = n, y = RRR, 
             col = pschiz)) + 
    geom_line() + 
    theme_bw()
  ggsave("~/Desktop/lre_schiz.pdf")
  ggplot(cc %>% filter(strategy == "hre"), 
         aes(x = q_exclude, y = RRR,
             col = pschiz)) + 
    geom_line() + 
    theme_bw()
  ggsave("~/Desktop/hre_schiz.pdf")
  browser()
}
