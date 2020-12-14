library(glue)
library(purrr)
library(data.table)
library(dplyr)

dir <- "/vol/sci/bio/data/shai.carmi/db2175/embryo_selection/"
peds_dir <- glue("{dir}/Peds/")
files <- Sys.glob(glue("{peds_dir}/*sscore"))
all <- map_dfr(files, function(x){fread(x)})
all_sel <- all %>%
  transmute(id = `#IID`, score = SCORE1_AVG)
write.csv(all_sel, glue("{dir}/sim_scores.csv"), row.names = F)
