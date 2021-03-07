library(data.table)
library(dplyr)
set.seed(5)
sample <- fread("/vol/sci/bio/data/shai.carmi/db2175/embryo_selection/LIJMC/LIJMC37_.1.sample")
samples <- sample[, .(ID_1, plink_pheno)]
samples1 <- copy(samples)
samples2 <- copy(samples)
setnames(samples1, c("ID.1", "pheno.1"))
setnames(samples2, c("ID.2", "pheno.2"))
merged <- merge(as.data.frame(samples1), as.data.frame(samples2))

subset <- merged %>% filter(!ID.1 == 0, !ID.2 == 0, ID.1 < ID.2)
selected <- subset %>% 
  group_by(pheno.1, pheno.2) %>% 
	slice(sample(1:n(), 10000)) %>%
	ungroup()

write.csv(selected, "selected_couples.csv", row.names = F)
