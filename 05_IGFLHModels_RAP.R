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
temp_sw <- standardize_columns(temp_sw, c("IGF1", "PopSize", "Weight", "MumAge"))
temp_r <- standardize_columns(temp_r, c("IGF1", "PopSize", "MumAge"))
temp_rw <- standardize_columns(temp_rw, c("IGF1", "PopSize", "Weight", "MumAge"))


# Prepare data list to pass to stan model
data_survival <- prepare_data_list(temp_s)
data_survwt <- prepare_data_list(temp_sw)
data_repro <- prepare_data_list(temp_r)
data_reprowt <- prepare_data_list(temp_rw)

# Run the model using cmdstanr
file_s <- c("./StanModels/IGF_Survival_RAP.stan")
file_sw <- c("./StanModels/IGF_SurvivalCtrlWt_RAP.stan")
file_r <- c("./StanModels/IGF_Reproduction_RAP.stan")
file_rw <- c("./StanModels/IGF_ReproductionCtrlWt_RAP.stan")

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
(p4 <- generate_and_plot_ppc(fit_mod_surv, "Survival_rep", temp_s, temp_s$Survival, "Survival"))
(p5 <- generate_and_plot_ppc(fit_mod_survwt, "Survival_rep", temp_sw, temp_sw$Survival, "Survival (Model controlling for weight)"))
(p6 <- generate_and_plot_ppc(fit_mod_repro, "BredAsAYearling_rep", temp_r, temp_r$BredAsAYearling, "Reproduction"))
(p7 <- generate_and_plot_ppc(fit_mod_reprowt, "BredAsAYearling_rep", temp_rw, temp_rw$BredAsAYearling, "Reproduction (Model controlling for weight)"))


# Exrtact samples 
# Add "true" IGF values estimated in model to original dataset with observed values to compare against
temp_s <- process_fit_mod(fit_mod_surv, temp_s)
temp_sw <- process_fit_mod(fit_mod_survwt, temp_sw)
temp_r <- process_fit_mod(fit_mod_repro, temp_r)
temp_rw <- process_fit_mod(fit_mod_reprowt, temp_rw)

# Plotting
# Scatter-plot of IGF_true and observed IGF values 
(violplot_s <- plot_violin(temp_s, "Survival"))
(violplot_sw <- plot_violin(temp_sw, "Survival (Model controlling for weight)"))
(violplot_r <- plot_violin(temp_r, "Reproduction"))
(violplot_rw <- plot_violin(temp_rw, "Reproduction (Model controlling for weight)"))

# Compare means and variances of IGF_obs vs IGF_true
var(temp_s$IGF1_sc)
var(temp_s$IGF_true)
cor.test(temp_s$IGF1_sc, temp_s$IGF_true)

var(temp_sw$IGF1_sc)
var(temp_sw$IGF_true)
cor.test(temp_sw$IGF1_sc, temp_sw$IGF_true)
# 
var(temp_r$IGF1_sc)
var(temp_r$IGF_true)
cor.test(temp_r$IGF1_sc, temp_r$IGF_true)

var(temp_rw$IGF1_sc)
var(temp_rw$IGF_true)
cor.test(temp_rw$IGF1_sc, temp_rw$IGF_true)

#
# Get posterior draws of weight and compare with observed
# bayesplot::mcmc_intervals(fit_mod_weight$draws(), pars = c("beta_IGF", "beta_SexF"))
post_survival <- process_fit_bernoulli_model_ma(fit_mod_surv, temp_s, 
                                                c("alpha", "beta_IGF", "beta_SexF", "beta_Twin", "beta_PopSize", "beta_MumAge1", "beta_MumAge2"), 
                                                "Survival")

post_survwt <- process_fit_bernwt_model_ma(fit_mod_survwt, temp_sw, 
                                           c("alpha", "beta_IGF", "beta_SexF", "beta_Twin", "beta_PopSize", "beta_Weight", "beta_MumAge1", "beta_MumAge2"), 
                                           "Survival")

post_repro <- process_fit_bernoulli_model_ma(fit_mod_repro, temp_r, 
                                             c("alpha", "beta_IGF", "beta_SexF", "beta_Twin", "beta_PopSize", "beta_MumAge1", "beta_MumAge2"), 
                                             "BredAsAYearling")

post_reprowt <- process_fit_bernwt_model_ma(fit_mod_reprowt, temp_rw, 
                                            c("alpha", "beta_IGF", "beta_SexF", "beta_Twin", "beta_PopSize", "beta_Weight", "beta_MumAge1", "beta_MumAge2"), 
                                            "BredAsAYearling")


#Plotting
(plot_surv <- ggplot(data= temp_s, aes(x=IGF_true, y=Survival)) +
    geom_point(colour="#5e8d83", alpha=0.6) +
    geom_line(data=post_survival$summary_500draws_mu,  aes(y=value, group=draw), colour="#339E66FF", alpha=1/15) +
    geom_line(data=post_survival$mu_median,  aes(y=value), size=1.2, colour="#078282FF", alpha=0.8) +
    theme_cowplot() +
    ylab("First-Year Overwinter Survival") + 
    xlab("") +
    panel_bg(fill = "gray95", color = NA) +
    grid_lines(color = "white") +
    xlab("Normalized & Corrected IGF-1 Concentration") +
    annotate("text", x=1.75, y=0.95, label= expression(beta[IGF-1]~"= 0.270 (95% CI: -0.034"~-~"0.595)")) +
    theme(legend.text = element_text(size=14),
          legend.title = element_text(size=14))) 


(plot_survwt <- ggplot(data= temp_sw, aes(x=IGF_true, y=Survival)) +
    geom_point(colour="#ccb0be", alpha=0.6) +
    geom_line(data=post_survwt$summary_500draws_mu,  aes(y=value, group=draw), colour="#72668a", alpha=1/15) +
    geom_line(data=post_survwt$mu_median,  aes(y=value), size=1.2, colour="#374971", alpha=0.8) +
    theme_cowplot() +
    ylab("First-Year Overwinter Survival") + 
    xlab("") +
    panel_bg(fill = "gray95", color = NA) +
    grid_lines(color = "white") +
    xlab("Normalized & Corrected IGF-1 Concentration") +
    annotate("text", x=1.75, y=0.95, label= expression(beta[IGF-1]~"= 0.107 (95% CI: -0.251"~-~"0.448)")) +
    ggtitle("Model accounting for August body mass") +
    theme(legend.text = element_text(size=14),
          legend.title = element_text(size=14),
          plot.title = element_text(size = 16, face="plain")))

(plot_repro <- ggplot(data= temp_r, aes(x=IGF_true, y=BredAsAYearling)) +
    geom_point(colour="#a67d65", alpha=0.6) +
    geom_line(data=post_repro$summary_500draws_mu,  aes(y=value, group=draw), colour="#a65e58", alpha=1/15) +
    geom_line(data=post_repro$mu_median,  aes(y=value), size=1.2, colour="#400101", alpha=0.8) +
    theme_cowplot() +
    ylab("First-Year Reproduction") + 
    xlab("") +
    panel_bg(fill = "gray95", color = NA) +
    grid_lines(color = "white") +
    annotate("text", x=1.75, y=0.95, label= expression(beta[IGF-1]~"= 0.384 (95% CI: 0.082"~-~"0.694)")) +
    theme(legend.text = element_text(size=14),
          legend.title = element_text(size=14)))

(plot_reprowt <- ggplot(data= temp_rw, aes(x=IGF_true, y=BredAsAYearling)) +
    geom_point(colour="#ff9b54", alpha=0.6) +
    geom_line(data=post_reprowt$summary_500draws_mu,  aes(y=value, group=draw), colour="#fa9500", alpha=1/15) +
    geom_line(data=post_reprowt$mu_median,  aes(y=value), size=1.2, colour="#ff7f51", alpha=0.8) +
    theme_cowplot() +
    ylab("First-Year Reproduction") + 
    xlab("") +
    panel_bg(fill = "gray95", color = NA) +
    grid_lines(color = "white") +
    annotate("text", x=1.75, y=0.95, label= expression(beta[IGF-1]~"= 0.270 (95% CI: -0.110"~-~"0.630)")) +
    ggtitle("Model accounting for August body mass") +
    theme(legend.text = element_text(size=14),
          legend.title = element_text(size=14),
          plot.title = element_text(size = 16, face="plain")))

# Save both survival and repro plots
plot_fin <- (plot_repro + plot_reprowt) /(plot_surv + plot_survwt) + plot_annotation(tag_levels = 'A')
plot_fin
ggsave("May25_IGFReproAndSurv_Fig4.tiff", plot_fin, dpi=600, height = 10, bg="white" )


# Model summary
# List the predictors for each model
predictors_surv <- c("alpha", "beta_IGF", "beta_SexF", "beta_Twin", "beta_PopSize", "beta_MumAge1", "beta_MumAge2", "sigma_e1", "sigma_t", "sigma_u", "sigma_v", "sigma_w")
predictors_survwt <- c("alpha", "beta_IGF", "beta_SexF", "beta_Twin", "beta_PopSize", "beta_Weight", "beta_MumAge1", "beta_MumAge2", "sigma_e1", "sigma_t", "sigma_u", "sigma_v", "sigma_w")
predictors_repro <- c("alpha", "beta_IGF", "beta_SexF", "beta_Twin", "beta_PopSize", "beta_MumAge1", "beta_MumAge2", "sigma_e1", "sigma_t", "sigma_u", "sigma_v", "sigma_w")
predictors_reprowt <- c("alpha", "beta_IGF", "beta_SexF", "beta_Twin", "beta_PopSize", "beta_Weight", "beta_MumAge1", "beta_MumAge2", "sigma_e1", "sigma_t", "sigma_u", "sigma_v", "sigma_w")

# Call the function for each model
full_post_surv <- summarize_model(fit_mod_surv, predictors_surv)
full_post_survwt <- summarize_model(fit_mod_survwt, predictors_survwt)
full_post_repro <- summarize_model(fit_mod_repro, predictors_repro)
full_post_reprowt <- summarize_model(fit_mod_reprowt, predictors_reprowt)

full_post_df <- bind_rows(full_post_surv, full_post_survwt,
                          full_post_repro, full_post_reprowt
                          )
full_post_df
