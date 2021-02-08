library(readxl)
library(glue)
library(data.table)
library(dplyr)
library(lubridate)
DIR <- "~/Documents/Docs/embryos/fwdexternalreagesoftheajszsamples/"

GetAges <- function() {
  
  ids <- fread("~/Documents/Docs/embryos/LIJMC37_.1.sample")
  control_ids <- ids[plink_pheno == 1, ID_1]
  controls <- read_excel(glue("{DIR}/control samples data.xls"))
  control_data <- controls %>%
    filter(INL_ID %in% control_ids) %>%
    select(ID = INL_ID, AGE = AGE)
  
  #Attached are 3 files relating to my schizophrenia case-control cohort:
  # 1) age data for the controls (straightforward);
  
  # 2 & 3) age data for the cases, divided by prefix (SZP and SZA) - 
  # ages need to be calculated by subtracting date of blood draw (tab 1, column B) 
  # from month/year of birth (tab 2, columns C&D).
  
  CaseAge <- function(f, num1, num2) {
    s1 <- read_excel(glue("{DIR}/{f}.xls"), sheet = num1) %>%
      select(ID = `Book Number`, blood_draw_date = `1-2`) %>%
      mutate(blood_draw_date = as.Date(blood_draw_date, "%d/%m/%y"))
    s2 <- read_excel(glue("{DIR}/{f}.xls"), sheet = num2) %>%
      select(ID = `Book Number`, month_birth = `2-9-1`, year_birth = `2-9-2`) %>%
      mutate(month_birth = if_else(is.na(month_birth), "06", month_birth)) %>%
      mutate(date_birth = as.Date(glue("{year_birth}-{month_birth}-01")))
    joined <- left_join(s1, s2, by = "ID") %>%
      mutate(AGE = round(as.numeric(blood_draw_date - date_birth) / 365)) %>%
      select(ID, AGE)
    return(joined)
  }
  cases1 <- CaseAge("SZA", num1 = 12, num2 = 11)
  cases2 <- CaseAge("SZP", num1 = 13, num2 = 12)
  
  all <- bind_rows(control_data, cases1, cases2) %>%
    distinct()
  
  ids_joined <- inner_join(ids %>% mutate(ID = ID_1), all, by = "ID")
  return(ids_joined)
}