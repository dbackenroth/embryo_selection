library(glue)
library(data.table)

dir <- "/vol/sci/bio/data/shai.carmi/db2175/embryo_selection/"

scores_file <- glue("{dir}/LIJMC_scores.sscore")
sample_file <- glue("{dir}/LIJMC/LIJMC37_.1.sample")

scores <- fread(scores_file)
samples <- fread(sample_file)

setnames(scores, "#IID", "ID_1")
merged <- merge(scores, samples)

merged$pheno <- as.numeric(merged$plink_pheno) - 1

mod <- glm(pheno ~ SCORE1_AVG, data = merged)
preds <- predict(mod, newdata = merged, type = "response")
print(cor(preds, merged$pheno))
print(cor(merged$SCORE1_AVG, merged$pheno))
print(table(merged$pheno))
