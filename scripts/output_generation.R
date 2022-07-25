rm(list = ls())
library(tidyverse)
library(readxl)
library(dplyr)

# Import simulated output -------------------------------------------------
path <- "output/m_based_simulation_tables/"
files <- list.files(path)
cov_files <- files[str_detect(files, "cov")]
mse_files <- files[str_detect(files, "mse")]

d_cov <-
  map(cov_files, ~ 
      read_excel(paste0(path, .x))) %>%
    setNames(cov_files) %>%
    imap(~ .x %>%
          mutate(case = .y))  %>%
  bind_rows() %>%
  separate(case, c("case", NA, "train_data", NA), sep = "[_]") # %>%
  # separate(case, c("case", NA, "train_data"), sep = "[_]")  %>%
  # separate("train_data", c("train_data", NA), sep = "[.]") %>%
  # mutate(outcome = ifelse(is.na(`Train with full data`), `Train with minority data`, `Train with full data`)) %>%
  #select(case, outcome, train_data, `CMMP & Lmer*`, "Regrouped CMMP & Regrouped Lmer*")

writexl::write_xlsx(d_cov, "output/final_table/m_based_all_cov.xlsx")

d_mse <-
  map(mse_files, ~ 
        read_excel(paste0(path, .x))) %>%
  setNames(mse_files) %>%
  imap(~ .x %>%
         mutate(case = .y))  %>%
  bind_rows() %>%
  separate(case, c("case", NA, "train_data", NA), sep = "[_]") # %>%
  # separate(case, c("case", NA, "train_data"), sep = "[_]")  %>%
  # separate("train_data", c("train_data", NA), sep = "[.]") %>%
  # select("case", "train_data","rowname", "RP", "Random Forest", "Boosted Reg Tree", "CMMP", 
  #        "Lmer*", "Regrouped Lmer*", "Regrouped CMMP", "PURE", "ULMM CMMP left", 
  #        "ULMM CMMP right", "Average" )
writexl::write_xlsx(d_mse, "output/final_table/mbased_all_mse.xlsx")
