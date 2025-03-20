set.seed(2025)
# Load libraries
if (!require("pacman")) install.packages("pacman")
pacman::p_load(see, tidybayes, rstanarm, rethinking, cmdstanr, parameters, tidyverse, lmerTest, ggeffects, plyr, reshape2, rptR, viridis, cowplot, bayesplot, patchwork, ggpubr, mgcv)

# Set wd
setwd("~/Desktop/IGF1_MS/Analysis/R_Scripts")
source("00_functions.R")

# Read Data
igf_lh_data <- read.csv('2025Mar_LambIGF1.csv',  header = T, stringsAsFactors = F, fileEncoding="UTF-8-BOM")

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

# Split data into different subsets (weight, foreleg, growth)
temp_w <- temp %>% select(-ForeLeg, -HornLen, -Horn, -Survival, -BredAsAYearling, -BirthWt, -DaysSinceBirth) %>% drop_na()
temp_f <- temp %>% select(-Weight, -HornLen, -Horn, -Survival, -BredAsAYearling, -BirthWt, -DaysSinceBirth) %>% drop_na()
temp_fw <- temp %>% select(-HornLen, -Survival, -Horn, -BredAsAYearling, -BirthWt, -DaysSinceBirth) %>% drop_na()
temp_h <- temp %>% filter(Horn == 3) %>% select(-Weight, -ForeLeg, -Survival, -BredAsAYearling, -BirthWt, -DaysSinceBirth) %>% drop_na()
temp_g <- temp %>% select(-ForeLeg, -Horn,  -HornLen, -Survival, -BredAsAYearling) %>% drop_na()
temp_hw <- temp %>% filter(Horn == 3) %>% select(-ForeLeg, -Survival, -BredAsAYearling, -BirthWt, -DaysSinceBirth) %>% drop_na()

# # Convert cols to numeric
cols_n <- c("ELISARunDate", "PlateNumber", "BirthYear", "MumID")
temp_w <- convert_to_num_fac(temp_w, cols_n)
temp_f <- convert_to_num_fac(temp_f, cols_n)
temp_fw <- convert_to_num_fac(temp_fw, cols_n)
temp_h <- convert_to_num_fac(temp_h, cols_n)
temp_g <- convert_to_num_fac(temp_g, cols_n)
temp_hw <- convert_to_num_fac(temp_hw, cols_n)

# Rescale variables
temp_w <- standardize_columns(temp_w, c("Weight", "IGF1", "PopSize", "MumAge"))
temp_f <- standardize_columns(temp_f, c("ForeLeg", "IGF1", "PopSize", "MumAge"))
temp_fw <- standardize_columns(temp_fw, c("ForeLeg", "IGF1", "PopSize", "Weight", "MumAge"))
temp_h <- standardize_columns(temp_h, c("HornLen", "IGF1", "PopSize", "MumAge"))
temp_g <- standardize_columns(temp_g, c("Weight", "IGF1", "PopSize", "DaysSinceBirth", "BirthWt", "MumAge"))
temp_hw <- standardize_columns(temp_hw, c("HornLen", "IGF1", "PopSize", "Weight", "MumAge"))


# Prepare data list to pass to stan model
data_weight <- prepare_data_list(temp_w)
data_foreleg <- prepare_data_list(temp_f)
data_forelegw<- prepare_data_list(temp_fw)
data_hornlen <- prepare_data_list(temp_h)
data_growth <- prepare_data_list(temp_g)
data_hornlenw <- prepare_data_list(temp_hw)

# Run the model using cmdstanr
# Foreleg Model
file_f <- c("./StanModels/IGF_Foreleg.stan")
fit_mod_foreleg <- run_stan_model(file_f, data_foreleg)

# Weight Model
file_w <- c("./StanModels/IGF_Weight.stan")
fit_mod_weight <- run_stan_model(file_w, data_weight)

# Foreleg ctrling for Weight Model
file_fw <- c("./StanModels/IGF_ForelegCtrlWt.stan")
fit_mod_forelegw <- run_stan_model(file_fw, data_forelegw)

# Growth model
file_g <- c("./StanModels/IGF_Growth.stan")
fit_mod_growth <- run_stan_model(file_g, data_growth)

# Horn Model
file_h <- c("./StanModels/IGF_HornLen.stan")
fit_mod_hornlen <- run_stan_model(file_h, data_hornlen)

# Horn Model ctrling for Weight 
file_hw <- c("./StanModels/IGF_HornLenCtrlWt.stan")
fit_mod_hornlenw <- run_stan_model(file_hw, data_hornlenw)

# Check the model diagnostics
fit_mod_list <- list(
  fit_mod_weight, 
  fit_mod_foreleg, 
  fit_mod_growth,
  fit_mod_forelegw,
  fit_mod_hornlen,
  fit_mod_hornlenw
)
run_cmdstan_diagnose(fit_mod_list)

# Check posterior dist
color_scheme_set("mix-teal-pink")
(p2 <- generate_and_plot_ppc(fit_mod_weight, "Weight_rep", temp_w, temp_w$Weight, "\nAugust Weight"))
(p1 <- generate_and_plot_ppc(fit_mod_foreleg, "ForeLeg_rep", temp_f, temp_f$ForeLeg, "\nAugust Foreleg Length"))
(p5 <- generate_and_plot_ppc(fit_mod_hornlen, "HornLen_rep", temp_h, temp_h$HornLen, "\nAugust Horn Length"))
(p4 <- generate_and_plot_ppc(fit_mod_growth, "Weight_rep", temp_g, temp_g$Weight, "\nAugust Weight\n(Model accounting for lamb age and birth weight)"))
(p3 <- generate_and_plot_ppc(fit_mod_forelegw, "ForeLeg_rep", temp_fw, temp_fw$ForeLeg, "August Foreleg Length\n(Model accounting for August body weight)"))
(p6 <- generate_and_plot_ppc(fit_mod_hornlenw, "HornLen_rep", temp_hw, temp_hw$HornLen, "August Horn Length\n(Model accounting for August body weight)"))

# Exrtact samples 
# Add "true" IGF values estimated in model to original dataset with observed values to compare against
temp_w <- process_fit_mod(fit_mod_weight, temp_w)
temp_f <- process_fit_mod(fit_mod_foreleg, temp_f)
temp_h <- process_fit_mod(fit_mod_hornlen, temp_h)
temp_g <- process_fit_mod(fit_mod_growth, temp_g)
temp_fw <- process_fit_mod(fit_mod_forelegw, temp_fw)
temp_hw <- process_fit_mod(fit_mod_hornlenw, temp_hw)

# Plotting
# Scatter-plot of IGF_true and observed IGF values 
(violplot_f <- plot_violin(temp_f, "August Foreleg Length"))
(violplot_h <- plot_violin(temp_h, "August Horn Length"))
(violplot_w <- plot_violin(temp_w, "August Weight"))
(violplot_g <- plot_violin(temp_g, "August Weight\n(Model accounting for lamb age and birth weight)"))

(violplot_fw <- plot_violin(temp_fw, "ForeLeg Ctrl Wt"))
(violplot_hw <- plot_violin(temp_hw, "HornLen Ctrl Wt"))

# Compare means and variances of IGF_obs vs IGF_true
var(temp_w$IGF1)
var(temp_w$IGF_true)
cor.test(temp_w$IGF_true, temp_w$IGF1)

var(temp_f$IGF1_sc)
var(temp_f$IGF_true)
cor.test(temp_f$IGF_true, temp_f$IGF1)

var(temp_h$IGF1_sc)
var(temp_h$IGF_true)
cor.test(temp_h$IGF_true, temp_h$IGF1)

var(temp_g$IGF1_sc)
var(temp_g$IGF_true)
cor.test(temp_g$IGF_true, temp_g$IGF1)

var(temp_fw$IGF1_sc)
var(temp_fw$IGF_true)
cor.test(temp_fw$IGF_true, temp_fw$IGF1)

var(temp_hw$IGF1_sc)
var(temp_hw$IGF_true)
cor.test(temp_hw$IGF_true, temp_hw$IGF1)


# Model summary
# List the predictors for each model
predictors_lin <- c("alpha", "beta_IGF","beta_SexF", "sigma_e", "sigma_e1", "sigma_t", "sigma_u", "sigma_v", "sigma_w")
predictors_linw <- c("alpha", "beta_IGF", "beta_SexF","beta_Weight", "sigma_e", "sigma_e1", "sigma_t", "sigma_u", "sigma_v", "sigma_w")
predictors_growth <- c("alpha", "beta_IGF", "beta_SexF", "beta_DSB", "beta_BirthWt", "sigma_e", "sigma_e1", "sigma_t", "sigma_u", "sigma_v", "sigma_w")

# Call the function for each model
full_post_weight <- summarize_model(fit_mod_weight, predictors_lin)
full_post_weight$Trait <- c("Weight")
full_post_foreleg <- summarize_model(fit_mod_foreleg, predictors_lin)
full_post_foreleg$Trait <- c("Foreleg")
full_post_hornlen <- summarize_model(fit_mod_hornlen, predictors_lin)
full_post_hornlen$Trait <- c("HornLen")
full_post_growth <- summarize_model(fit_mod_growth, predictors_growth)
full_post_growth$Trait <- c("Growth")
full_post_forelegw <- summarize_model(fit_mod_forelegw, predictors_linw)
full_post_forelegw$Trait <- c("ForelegWt")
full_post_hornlenw <- summarize_model(fit_mod_hornlenw, predictors_linw)
full_post_hornlenw$Trait <- c("HornLenWt")

full_post_df <- bind_rows(full_post_weight, full_post_foreleg, full_post_hornlen, full_post_growth, full_post_forelegw, full_post_hornlenw)
# write.table(full_post_df, file = "./Aug2024_IGFDataAnalysis/Table2A.csv", sep = ",",row.names = FALSE)
