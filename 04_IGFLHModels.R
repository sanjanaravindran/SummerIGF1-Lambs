set.seed(2025)
# Load libraries
if (!require("pacman")) install.packages("pacman")
pacman::p_load(see, tidybayes, rstanarm, rethinking, cmdstanr, parameters, tidyverse, lmerTest, ggeffects, plyr, reshape2, rptR, viridis, cowplot, bayesplot, patchwork, ggpubr, mgcv)

# Set wd
setwd("~/Desktop/IGF1_MS/Analysis/R_Scripts")
source("00_functions.R")

# Read Data
igf_lh_data <- read.csv('2025Mar_LambIGF1_Anon.csv',  header = T, stringsAsFactors = F, fileEncoding="UTF-8-BOM")

# # Temp dataset for stan models 
temp <- igf_lh_data

# Conver Sex to female or not
temp <- temp %>%
  mutate(SexF = case_when(Sex == 1 ~ 1, 
                          Sex == 2 ~ 0))

# Mean imputing mum age
temp <- temp %>%
  mutate(MumAge = case_when(is.na(MumAge) ~ mean(MumAge, na.rm=T),
                            TRUE ~ MumAge))

# Subset relevant columns
temp <- temp %>%
  select(ID, Weight, ForeLeg, HornLen, Horn, Survival, BredAsAYearling, IGF1, SexF, Twin, PopSize, BirthYear, MumID, MumAge, ELISARunDate, PlateNumber, DaysSinceBirth, BirthWt) 

# Split data into different subsets (weight, foreleg, survival, repro)
temp_s <- temp %>% select(-ForeLeg, -Weight, -HornLen, -Horn, -BredAsAYearling, -BirthWt, -DaysSinceBirth) %>% drop_na()
temp_sw <- temp %>% select(-ForeLeg, -HornLen, -Horn, -BredAsAYearling, -BirthWt, -DaysSinceBirth) %>% drop_na()
temp_r <- temp %>% select(-ForeLeg, -Weight, -HornLen, -Horn, -Survival, -BirthWt, -DaysSinceBirth) %>% drop_na()
temp_rw <- temp %>% select(-ForeLeg, -HornLen, -Horn, -Survival, -BirthWt, -DaysSinceBirth) %>% drop_na()

# # Convert cols to numeric
cols_n <- c("ELISARunDate", "PlateNumber", "BirthYear", "MumID")
temp_s <- convert_to_num_fac(temp_s, cols_n)
temp_sw <- convert_to_num_fac(temp_sw, cols_n)
temp_r <- convert_to_num_fac(temp_r, cols_n)
temp_rw <- convert_to_num_fac(temp_rw, cols_n)

# Rescale variables
temp_s <- standardize_columns(temp_s, c("IGF1", "PopSize", "MumAge"))
temp_sw <- standardize_columns(temp_sw, c("IGF1", "Weight", "PopSize", "MumAge"))
temp_r <- standardize_columns(temp_r, c("IGF1", "PopSize", "MumAge"))
temp_rw <- standardize_columns(temp_rw, c("IGF1", "Weight", "PopSize", "MumAge"))

# Prepare data list to pass to stan model
data_survival <- prepare_data_list(temp_s)
data_survwt <- prepare_data_list(temp_sw)
data_repro <- prepare_data_list(temp_r)
data_reprowt <- prepare_data_list(temp_rw)

# Run the model using cmdstanr
file_s <- c("./StanModels/IGF_Survival.stan")
file_sw <- c("./StanModels/IGF_SurvivalCtrlWt.stan")
file_r <- c("./StanModels/IGF_Reproduction.stan")
file_rw <- c("./StanModels/IGF_ReproductionCtrlWt.stan")

fit_mod_surv <- run_stan_model(file_s, data_survival)
fit_mod_survwt <- run_stan_model(file_sw, data_survwt)
fit_mod_repro <- run_stan_model(file_r, data_repro)
fit_mod_reprowt <- run_stan_model(file_rw, data_reprowt)

# Check the model
fit_mod_list <- list(
  fit_mod_surv,
  fit_mod_survwt,
  fit_mod_repro,
  fit_mod_reprowt
)
run_cmdstan_diagnose(fit_mod_list)

# Check posterior dist
color_scheme_set("mix-teal-pink")
(p8 <- generate_and_plot_ppc(fit_mod_surv, "Survival_rep", temp_s, temp_s$Survival, "\nFirst-Year Survival"))
(p5 <- generate_and_plot_ppc(fit_mod_survwt, "Survival_rep", temp_sw, temp_sw$Survival, "Survival (Model controlling for weight)"))
(p9 <- generate_and_plot_ppc(fit_mod_repro, "BredAsAYearling_rep", temp_r, temp_r$BredAsAYearling, "\nFirst-Year Breeding Probability"))
(p7 <- generate_and_plot_ppc(fit_mod_reprowt, "BredAsAYearling_rep", temp_rw, temp_rw$BredAsAYearling, "Reproduction (Model controlling for weight)"))


# Exrtact samples 
# Add "true" IGF values estimated in model to original dataset with observed values to compare against
temp_s <- process_fit_mod(fit_mod_surv, temp_s)
temp_sw <- process_fit_mod(fit_mod_survwt, temp_sw)
temp_r <- process_fit_mod(fit_mod_repro, temp_r)
temp_rw <- process_fit_mod(fit_mod_reprowt, temp_rw)

# Plotting
# Scatter-plot of IGF_true and observed IGF values 
(violplot_s <- plot_violin(temp_s, "First-Year Survival"))
(violplot_sw <- plot_violin(temp_sw, "Survival (Model controlling for weight)"))
(violplot_r <- plot_violin(temp_r, "First-Year Breeding Probability"))
(violplot_rw <- plot_violin(temp_rw, "First-Year Breeding Probability (Model controlling for weight)"))

# Compare means and variances of IGF_obs vs IGF_true
var(temp_s$IGF1_sc)
var(temp_s$IGF_true)
cor.test(temp_s$IGF1_sc, temp_s$IGF_true)

var(temp_sw$IGF1_sc)
var(temp_sw$IGF_true)
cor.test(temp_sw$IGF1_sc, temp_sw$IGF_true)

var(temp_r$IGF1_sc)
var(temp_r$IGF_true)
cor.test(temp_r$IGF1_sc, temp_r$IGF_true)

var(temp_rw$IGF1_sc)
var(temp_rw$IGF_true)
cor.test(temp_rw$IGF1_sc, temp_rw$IGF_true)


# Model summary
# List the predictors for each model
predictors_surv <- c("alpha", "beta_IGF", "beta_SexF", "sigma_e1", "sigma_t", "sigma_u", "sigma_v", "sigma_w")
predictors_survwt <- c("alpha", "beta_IGF", "beta_SexF", "beta_Weight", "sigma_e1", "sigma_t", "sigma_u", "sigma_v", "sigma_w")
predictors_repro <- c("alpha", "beta_IGF", "beta_SexF",  "sigma_e1", "sigma_t", "sigma_u", "sigma_v", "sigma_w")
predictors_reprowt <- c("alpha", "beta_IGF", "beta_SexF", "beta_Weight","sigma_e1", "sigma_t", "sigma_u", "sigma_v", "sigma_w")

# Call the function for each model
full_post_surv <- summarize_model(fit_mod_surv, predictors_surv)
full_post_surv$Trait <- c("Surv")
full_post_survwt <- summarize_model(fit_mod_survwt, predictors_survwt)
full_post_survwt$Trait <- c("SurvWt")
full_post_repro <- summarize_model(fit_mod_repro, predictors_repro)
full_post_repro$Trait <- c("Repro")
full_post_reprowt <- summarize_model(fit_mod_reprowt, predictors_reprowt)
full_post_reprowt$Trait <- c("ReproWt")


full_post_df <- bind_rows(full_post_surv, full_post_survwt, full_post_repro, full_post_reprowt)
full_post_df
# write.table(full_post_df, file = "./Aug2024_IGFDataAnalysis/TableLHNoPreds.csv", sep = ",",row.names = FALSE)
