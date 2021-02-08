library(MASS)
library(glue)
library(dplyr)
library(data.table)

get_data <- function() {
  
  DIR <- "~/Documents/Docs/embryos/"
  sim_scores <- fread(glue("gunzip -c {DIR}/sim_scores.csv.gz"))
  sim_scores[, c("p1", "p2", "cnum") := tstrsplit(id, "_", fixed=TRUE, type.convert = T)]
  
  parents_info <- fread(glue("{DIR}/LIJMC37_.1.sample"))
  p_info <- parents_info[, .(ID_1, plink_pheno)]
  setnames(p_info, c("p1", "p1_pheno"))
  sim_scores <- merge(sim_scores, p_info)
  setnames(p_info, c("p2", "p2_pheno"))
  sim_scores <- merge(sim_scores, p_info, by = "p2")
  sim_scores[, `:=`(id = NULL)]
  
  couples <- unique(sim_scores[, .(p1, p2, p1_pheno, p2_pheno)])
  
  parents_scores <- fread(glue("{DIR}/LIJMC_scores.sscore"))
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