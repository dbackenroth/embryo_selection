library(MASS)
library(glue)
library(dplyr)
library(data.table)
library(ggplot2)
library(latex2exp)
library(purrr)
library(readxl)
library(lubridate)

# ID, AGE, plink_pheno, sex (1, 2)
GetAges <- function(disease = "schiz") {
  if (disease == "crohns") {
    DIR <- "~/Documents/Docs/embryos/crohns/"
    Process <- function(x) {
      x %>%
        mutate(ID = paste(FIDplink, IIDplink, sep = "_"), 
               plink_pheno = if_else(diagnosis=="Control", 1, 2), 
               sex = if_else(PhenoGender == "Male", 1, 2), 
               AGE = runif(n(), min = 25, max = 100))
    }
    cases <- read_excel(glue("{DIR}/SampleInfo.xlsx"), sheet = "Cases") %>%
      Process()
    controls <- read_excel(glue("{DIR}/SampleInfo.xlsx"), sheet = "Controls") %>%
      Process()
    both <- bind_rows(cases, controls) %>%
      dplyr::select(ID, AGE, plink_pheno, sex)
    return(both)
  }
  if (disease == "schiz") {
    DIR <- "~/Documents/Docs/embryos/fwdexternalreagesoftheajszsamples/"
    ids <- fread("~/Documents/Docs/embryos/LIJMC37_.1.sample")
    control_ids <- ids[plink_pheno == 1, ID_1]
    controls <- read_excel(glue("{DIR}/control samples data.xls"))
    control_data <- controls %>%
      filter(INL_ID %in% control_ids) %>%
      dplyr::select(ID = INL_ID, AGE = AGE)
    
    #Attached are 3 files relating to my schizophrenia case-control cohort:
    # 1) age data for the controls (straightforward);
    
    # 2 & 3) age data for the cases, divided by prefix (SZP and SZA) - 
    # ages need to be calculated by subtracting date of blood draw (tab 1, column B) 
    # from month/year of birth (tab 2, columns C&D).
    
    CaseAge <- function(f, num1, num2) {
      s1 <- read_excel(glue("{DIR}/{f}.xls"), sheet = num1) %>%
        dplyr::select(ID = `Book Number`, blood_draw_date = `1-2`) %>%
        mutate(blood_draw_date = as.Date(blood_draw_date, "%d/%m/%y"))
      s2 <- read_excel(glue("{DIR}/{f}.xls"), sheet = num2) %>%
        dplyr::select(ID = `Book Number`, month_birth = `2-9-1`, year_birth = `2-9-2`) %>%
        mutate(month_birth = if_else(is.na(month_birth), "06", month_birth)) %>%
        mutate(date_birth = as.Date(glue("{year_birth}-{month_birth}-01")))
      joined <- left_join(s1, s2, by = "ID") %>%
        mutate(AGE = round(as.numeric(blood_draw_date - date_birth) / 365)) %>%
        dplyr::select(ID, AGE)
      return(joined)
    }
    cases1 <- CaseAge("SZA", num1 = 12, num2 = 11)
    cases2 <- CaseAge("SZP", num1 = 13, num2 = 12)
    
    all <- bind_rows(control_data, cases1, cases2) %>%
      distinct()
    
    ids_joined <- inner_join(ids %>% mutate(ID = ID_1), all, by = "ID")
    return(ids_joined[, .(ID, AGE, plink_pheno, sex)])
  }
}

get_data_helper <- function(disease) {
  DIR <- "~/Documents/Docs/embryos/"
  if (disease == "schiz") {
    sim_scores <- fread(glue("gunzip -c {DIR}/sim_scores.csv.gz"))
    sim_scores[, c("p1", "p2", "cnum") := tstrsplit(id, "_", fixed=TRUE, type.convert = T)]
    p_info <- GetAges(disease = disease)
    parents_scores <- fread(glue("{DIR}/LIJMC_scores.sscore"))
    gwas_info <- fread(glue("{DIR}/daner_flipped.tsv"))
  } 
  if (disease == "crohns") {
    dir <- glue("{DIR}/{disease}/")
    p_info <- GetAges(disease = disease)
    sim_scores <- fread(glue("gunzip -c {dir}/sim_scores.csv.gz"))
    sim_scores[, c("p1", "p2", "cnum") := tstrsplit(id, "%", fixed=TRUE, type.convert = T)]
    gwas_info <- fread(glue("{dir}/daner_flipped.tsv"))
    parents_scores <- fread(glue("{dir}/parents_scores.sscore"))
  }
  sim_scores[, `:=`(id = NULL)]
  p_info <- as.data.table(p_info)
  p_info[, `:=`(pheno = as.numeric(plink_pheno) - 1)]
  p_info[, `:=`(plink_pheno = NULL)]
  
  setnames(p_info, "ID", "IID")
  setnames(parents_scores, "#IID", "IID")
  p_info <- merge(p_info, parents_scores, by = "IID")
  setnames(p_info, "SCORE1_AVG", "score")
  p1 <- p_info[, .(IID, pheno, score)]
  
  p2 <- copy(p1)
  cols <- c("IID", "pheno", "score")
  setnames(p1, cols, c("p1", "p1_pheno", "p1_score"))
  setnames(p2, cols, c("p2", "p2_pheno", "p2_score"))
  s1 <- merge(sim_scores, p1, by = "p1")
  s2 <- merge(s1, p2, by = "p2")
  setnames(gwas_info, "new_BETA", "BETA")
  return(list(sim_scores = s2, p_info = p_info, gwas_info = gwas_info))
}

get_data <- function(disease = "schiz", prevalence = 0.01, num_couples = 5000) {
  l <- get_data_helper(disease = disease)
  sim_scores <- l$sim_scores %>%
    mutate(num_parents_cases = p1_pheno + p2_pheno, 
           couple = paste(p1, p2))
  if (!is.null(num_couples)) {
    num_aval_couples <- sim_scores %>%
      group_by(num_parents_cases) %>%
      count() %>%
      ungroup()
    stopifnot(n_distinct(num_aval_couples$n) == 1)
    num_to_sample <- num_aval_couples %>%
      mutate(num_to_sample = case_when(num_parents_cases == 0 ~ round(num_couples * (1 - prevalence) ^ 2), 
                                       num_parents_cases == 1 ~ round(num_couples * 2 * prevalence * (1 - prevalence)), 
                                       num_parents_cases == 2 ~ round(num_couples * prevalence ^ 2))) %>%
      dplyr::select(-n)
    num_to_sample$num_to_sample[num_to_sample$num_parents_cases == 2] <- num_couples - sum(num_to_sample$num_to_sample[num_to_sample$num_parents_cases < 2])
    set.seed(1)
    sel_couples <- sim_scores %>%
      left_join(num_to_sample, by = "num_parents_cases") %>%
      dplyr::select(couple, num_parents_cases, num_to_sample) %>%
      distinct() %>%
      split(.$num_parents_cases) %>%
      map_dfr(function(x){
        rows_to_sample <- sample(1:nrow(x), x$num_to_sample[1])
        x[rows_to_sample, ]
      })
    sim_scores_sel <- sim_scores %>%
      filter(couple %in% sel_couples$couple)
    sim_scores <- as.data.table(sim_scores_sel)
  }
  p_info <- l$p_info
  mod <- glm(pheno ~ score, family = binomial(), data = p_info)
  
  #sim_scores[, `:=`(risk = predict(mod, newdata = sim_scores, type = "response"))]
  SAMPLE_PREV <- mean(p_info$pheno == 1)
  
  adjusted_link <- predict(mod, newdata = sim_scores, type = "link") + 
    log(prevalence / (1 - prevalence)) - 
    log(SAMPLE_PREV / (1 - SAMPLE_PREV))
  adjusted_risk <- 1 / (1 + exp(-adjusted_link))
  sim_scores[, `:=`(risk = adjusted_risk, unadj_risk = predict(mod, newdata = sim_scores, type = "response"))]
  p_info[, `:=`(risk = predict(mod, newdata = p_info, type = "response"))]
  return(list(sim_scores = sim_scores, p_info = p_info, gwas_info = l$gwas_info))
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