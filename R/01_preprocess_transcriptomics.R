library(tidyverse)

# directory containing the files
data_dir <- "data/transcriptomics/FBA_RNA_windows_TPMs"

# list files
files <- c(
  "fba_ready_iCHO3K_TPM.csv",
  "fba_ready_iCHO1776_TPM.csv",
  "fba_ready_iCHO2441_TPM.csv"
)

process_file <- function(file) {
  
  df <- read_csv(file.path(data_dir, file))
  
  # extract model name from filename
  model_raw <- stringr::str_extract(file, "iCHO[0-9A-Za-z]+")
  
  # map to final model names used in Python
  model_name <- case_when(
    model_raw == "iCHO1776" ~ "iCHOv1",
    model_raw == "iCHO2441" ~ "iCHO2441",
    model_raw == "iCHO3K" ~ "iCHO3K_cho_prod_generic_unblocked"
  )
  
  df_processed <- df %>%
    
    pivot_longer(
      cols = -EntrezID,
      names_to = "exp_window",
      values_to = "expression"
    ) %>%
    
    separate(exp_window, into = c("experiment", "window"), sep = "_") %>%
    
    mutate(
      condition = case_when(
        experiment %in% c("E13","E15","E17","E19") ~ "Constant",
        experiment %in% c("E14","E16","E18","E20") ~ "Temp. shifted"
      )
    ) %>%
    
    group_by(EntrezID, condition, window) %>%
    summarise(expression = mean(expression), .groups = "drop") %>%
    
    mutate(context = paste(condition, window, sep = "_")) %>%
    
    select(EntrezID, context, expression) %>%
    
    pivot_wider(
      names_from = context,
      values_from = expression
    ) %>%
    
    mutate(model = model_name) %>%
    
    select(model, everything())
  
  return(df_processed)
}

expr_combined <- purrr::map_dfr(files, process_file)

write_csv(expr_combined, "data/transcriptomics/FBA_RNA_windows_TPMs/charrun_FBA_RNA3_all_V01_20260717_VS.csv")

#quick qc
n_distinct(expr_combined$EntrezID)
sum(is.na(expr_combined))
expr_combined %>%
  pivot_longer(-c(EntrezID,model)) %>%
  ggplot(aes(value)) +
  geom_histogram(bins = 50)
