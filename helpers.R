library(MASS)
library(glue)
library(dplyr)
library(data.table)
library(ggplot2)
library(latex2exp)

get_data_helper <- function() {
  DIR <- "~/Documents/Docs/embryos/"
  sim_scores <- fread(glue("gunzip -c {DIR}/sim_scores.csv_old.gz"))
  sim_scores[, c("p1", "p2", "cnum") := tstrsplit(id, "_", fixed=TRUE, type.convert = T)]
  
  parents_info <- fread(glue("{DIR}/LIJMC37_.1.sample"))
  p_info <- parents_info[, .(ID_1, plink_pheno)]
  
  parents_scores <- fread(glue("{DIR}/LIJMC_scores.sscore"))
  return(list(sim_scores = sim_scores, p_info = p_info, parents_scores = parents_scores))
}

compare_variance <- function() {
  l <- get_data_helper()
  children <- l$sim_scores
  var_children <- children %>%
    group_by(p1, p2) %>%
    summarise(var = var(score), 
              mean = mean(score))
  p_info <- l$p_info
  setnames(p_info, c("#IID", "pheno"))
  parents <- l$parents_scores
  parents <- merge(parents, p_info, by = "#IID")
  
  parent_scores0 <- unique(children[, .(p1, p2)])
  parents1 <- parents %>%
    transmute(p1 = `#IID`, p1_score = SCORE1_AVG, p1_pheno = pheno)
  parents2 <- parents %>%
    transmute(p2 = `#IID`, p2_score = SCORE1_AVG, p2_pheno = pheno)
  parent_scores <- reduce(list(parent_scores0, parents1, parents2), left_join) %>%
    mutate(average_score = (p1_score + p2_score) / 2) 
  bins_data <- left_join(var_children, parent_scores, by = c("p1", "p2")) %>%
    mutate(num_parents_schiz = (p1_pheno == 2) + (p2_pheno == 2), 
           parent_score_diff = abs(p1_score - p2_score)) %>%
    ungroup()
  # variance depends on difference between parent scores as well as average score
  # difference depends on average score
  mod <- lm(var ~ poly(parent_score_diff, 3), bins_data)
  bd <- bins_data %>%
    ungroup() %>%
    mutate(resid = residuals(mod))
  ggplot(bd, aes(x = average_score, y = resid)) + geom_smooth()
  browser()
  plot(bins_data$average_score, residuals(mod))
  
  p <- ggplot(bins_data, aes(x = average_score, y = var)) + 
    geom_smooth() #+ 
    #facet_wrap(~num_parents_schiz)
  
  browser()
  
  
  
  parents <- parents %>%
    group_by(pheno) %>%
    mutate(weight = if_else(pheno == 1, 0.99 / n(), 0.01 / n()))
  wtd_var <- wtd.var(parents$SCORE1_AVG, weights = parents$weight, normwt = T)
  cat("Weighted variance parents:", wtd_var, "\n")
  cat("Unweighted variance parents:", var(parents$SCORE1_AVG), "\n")
  cat("Mean variance children conditional on parents:", mean(var_children$var), "\n")
  cat("Variance parents by phenotype\n")
  parents %>%
    ungroup() %>%
    mutate(pheno = if_else(pheno==2, "Schiz.", "Not schiz.")) %>%
    group_by(pheno) %>%
    summarise(v = var(SCORE1_AVG))
  
  var <- wtd_var
  n <- 20
  sims <- rchisq(600000, df = n - 1) / (n-1) * (wtd_var / 2)
  var_children$v2 <-  var_children$var
  aa <- bind_rows(var_children %>% transmute(type = "observed", var = v2), 
                  data.frame(var = sims) %>%
                    mutate(type = "theoretical"))
  ggplot(aa, aes(x = var, col = type)) + 
    geom_density() + 
    xlab("Sampling variance") +
    geom_vline(xintercept = wtd_var / 2)
  
  browser()
}

get_data <- function() {
  l <- get_data_helper()
  sim_scores <- l$sim_scores
  parents_scores <- l$parents_scores
  p_info <- l$p_info
  
  setnames(p_info, c("p1", "p1_pheno"))
  sim_scores <- merge(sim_scores, p_info)
  setnames(p_info, c("p2", "p2_pheno"))
  sim_scores <- merge(sim_scores, p_info, by = "p2")
  sim_scores[, `:=`(id = NULL)]
  browser()
  
  
  couples <- unique(sim_scores[, .(p1, p2, p1_pheno, p2_pheno)])
  
  #parents_scores <- fread(glue("{DIR}/LIJMC_scores.sscore"))
  setnames(p_info, c("#IID", "pheno"))
  parents_scores <- merge(parents_scores, p_info, by = "#IID")
  parents_scores[, `:=`(pheno = as.numeric(pheno) - 1)]
  setnames(parents_scores, "SCORE1_AVG", "score")
  mod <- glm(pheno ~ score, family = binomial(), data = parents_scores)
  
  #sim_scores[, `:=`(risk = predict(mod, newdata = sim_scores, type = "response"))]
  SAMPLE_PREV <- mean(parents_scores$pheno == 1)
  adjusted_link <- predict(mod, newdata = sim_scores, type = "link") + 
    log(PREV / (1 - PREV)) - 
    log(SAMPLE_PREV / (1 - SAMPLE_PREV))
  adjusted_risk <- 1 / (1 + exp(-adjusted_link))
  sim_scores[, `:=`(risk = adjusted_risk, unadj_risk = predict(mod, newdata = sim_scores, type = "response"))]
  parents_scores[, `:=`(risk = predict(mod, newdata = parents_scores, type = "response"))]
  return(list(sim = sim_scores, parents = parents_scores))
}

risk_reduction_lowest = function(r,K,n)
{
  zk = qnorm(K, lower.tail=F)
  integrand_lowest = function(t)
    return(dnorm(t)*pnorm((zk-t*sqrt(1-r^2/2)) / (r/sqrt(2)), lower.tail=F)^n)
  risk = integrate(integrand_lowest,-Inf,Inf)$value
  return((K-risk)/K)
}

risk_reduction_exclude = function(r,K,q,n)
{
  zk = qnorm(K, lower.tail=F)
  zq = qnorm(q, lower.tail=F)
  integrand_t = function(t,u)
    return(dnorm(t)*pnorm((zk-r/sqrt(2)*(u+t))/sqrt(1-r^2),lower.tail=F))
  integrand_u = function(us)
  {
    y = numeric(length(us))
    for (i in seq_along(us))
    {
      u = us[i]
      beta = zq*sqrt(2)-u
      internal_int1 = integrate(integrand_t,-Inf,beta,u)$value
      denom1 = pnorm(beta)
      if (denom1==0) {denom1=1e-300} # Avoid dividing by zero
      numer1 = 1-pnorm(beta,lower.tail=F)^n
      internal_int2 = integrate(integrand_t,beta,Inf,u)$value
      denom2 = pnorm(beta,lower.tail=F)
      if (denom2==0) {denom2=1e-300} # Avoid dividing by zero
      numer2 = pnorm(beta,lower.tail=F)^n
      y[i] = dnorm(u) * (numer1/denom1*internal_int1 + numer2/denom2*internal_int2)
    }
    return(y)
  }
  risk = integrate(integrand_u,-Inf,Inf)$value
  return((K-risk)/K)
}