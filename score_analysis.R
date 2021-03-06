library(glue)
library(data.table)
library(pROC)
source("helpers.R")
#source("get_ages.R")

#dir <- "/vol/sci/bio/data/shai.carmi/db2175/embryo_selection/"

CrohnsR2 <- function(prevalence = 0.013) {
  dir <- "~/Documents/Docs/embryos/crohns/"
  scores_file <- glue("{dir}/parents_scores.sscore")
  demog <- GetAges("crohns")
  GetR2(scores_file, demog, K = prevalence, label = "crohns") %>%
    return()
}

#scores_file <- glue("{dir}/LIJMC_scores.sscore_OLD.sscore")
SchizR2 <- function(prevalence = 0.01) {
  dir <- "~/Documents/Docs/embryos/"
  scores_file <- glue("{dir}/LIJMC_scores.sscore")
  demog <- GetAges("schiz")
  GetR2(scores_file, demog, K = prevalence, label = "schizophrenia") %>%
    return()
}

GetR2 <- function(scores_file, demog, K, label) {
  
  scores <- fread(scores_file)
  setnames(scores, "#IID", "ID")
  merged <- merge(scores, demog)
  regr_data <- merged[, .(ID, SCORE1_AVG, sex, AGE, plink_pheno)]
  regr_data$pheno <- as.numeric(regr_data$plink_pheno) - 1
  
  no_nas <- regr_data %>%
    filter(!is.na(AGE))
  
  mod1 <- glm(pheno ~ SCORE1_AVG, data = no_nas, family = binomial())
  #mod2 <- glm(pheno ~ SCORE1_AVG + sex, data = regr_data)
  #mod3 <- glm(pheno ~ SCORE1_AVG + sex + AGE, data = regr_data)
  #mod4 <- glm(pheno ~ SCORE1_AVG + AGE, data = regr_data)
  mod2 <- glm(pheno ~ SCORE1_AVG + sex + poly(AGE, 2), data = no_nas, family = binomial())
  mod3 <- glm(pheno ~ SCORE1_AVG + sex, data = no_nas, family = binomial())
  no_nas$pred1 <- predict(mod1, newdata = no_nas, type = "response")
  no_nas$pred2 <- predict(mod2, newdata = no_nas, type = "response")
  no_nas$pred3 <- predict(mod3, newdata = no_nas, type = "response")
  roc1 <- roc(response = no_nas$pheno, predictor = no_nas$pred1)
  roc2 <- roc(response = no_nas$pheno, predictor = no_nas$pred2)
  roc3 <- roc(response = no_nas$pheno, predictor = no_nas$pred3)
  auc1 <- auc(roc1)
  auc2 <- auc(roc2)
  auc3 <- auc(roc3)
  cat("AUC of model with just score:", auc1, "\n")
  cat("AUC of model with age and sex:", auc2, "\n")
  cat("AUC of model with sex:", auc3, "\n")
  
  ggplot(no_nas, aes(x = pred1, y = as.numeric(plink_pheno) - 1)) + 
    geom_smooth() + 
    geom_abline() + 
    theme_bw() +
    xlab("Model prediction") + 
    xlim(0, 1) + 
    ylab("% cases (smoothed)")
  ggsave(glue("Results/{label}_model_calibration.pdf"))
  ggplot(no_nas, aes(x = pred1)) + 
    geom_density() + 
    xlab("Model prediction") + 
    theme_bw() + 
    xlim(0, 1)
  ggsave(glue("Results/{label}_predictions_in_training_set.pdf"))
  
  
  corr <- cor(regr_data$SCORE1_AVG, regr_data$pheno)
  
  nc <- sum(regr_data$pheno == 1)
  nt <- sum(regr_data$pheno == 0)
  r2 <- corr^2
  #print(cor(preds, merged$pheno))
  #print(cor(merged$SCORE1_AVG, merged$pheno))
  #print(table(merged$pheno))
  
  #K = 0.01
  #nc = 897
  #nt = 1629
  #r2 = 0.09
  
  P = nc/(nc+nt)
  t = qnorm(1-K)
  z = dnorm(t)
  m = z/K
  
  C = K^2*(1-K)^2 / (z^2*P*(1-P))
  theta = m * (P-K)/(1-K) * (m*(P-K)/(1-K)-t)
  r2_liab = r2*C / (1+r2*theta*C)
  
  cat(sprintf("R^2 on the liability scale: %g\n", r2_liab))
  return(r2_liab)
}