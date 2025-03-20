set.seed(2025)
# Load libraries
if (!require("pacman")) install.packages("pacman")
pacman::p_load(see, tidybayes, rstanarm, rethinking, cmdstanr, parameters, tidyverse, lmerTest, ggeffects, plyr, reshape2, rptR, viridis, cowplot, bayesplot, patchwork, ggpubr, mgcv)

# Set wd
setwd("~/Desktop/IGF1_MS/Analysis/R_Scripts")
source("00_functions.R")

# Read Data
igf_lh_data <- read.csv('2025Mar_LambIGF1.csv',  header = T, stringsAsFactors = F, fileEncoding="UTF-8-BOM")

# Recoding Sex as female or not
igf_lh_data <- igf_lh_data %>%
  mutate(SexF = case_when(Sex == 1 ~ "1", 
                          Sex == 2 ~ "0"))

# Specifying horntype as normal, polled or scurred
igf_lh_data <- igf_lh_data %>%
  mutate(HornType = case_when(Horn == 1 ~ "Scurred", 
                              Horn == 2 ~ "Polled",
                              Horn == 3 ~ "Normal"))

# Convert relevant columns to factor 
cols_f <- c("SexF", "Twin", "HornType")
igf_lh_data <- convert_to_factor(igf_lh_data, cols_f)

# Mean imputing mum age
temp <- igf_lh_data %>%
  mutate(MumAge = case_when(is.na(MumAge) ~ mean(MumAge, na.rm=T),
                            TRUE ~ MumAge))

# Rescale variables 
cols_sc <- c("PopSize", "MumAge")
temp <- standardize_columns(temp, cols_sc)

# Plot IGF-1 histogram
igf_hist <- ggplot(temp, aes(x=IGF1)) +
  geom_histogram() +
  labs(subtitle = "Plasma IGF-1 Concentration (ng/ml)") +
  xlab("")
igf_hist

############
# Analyse  #
############

# IGF-1 Model
mod_igf <- rstanarm::stan_lmer(IGF1 ~  SexF + Twin + PopSize_sc + poly(MumAge_sc, degree = 2, raw = TRUE) + (1|ELISARunDate) + (1|PlateNumber) + (1|BirthYear) + (1|MumID),
                               cores=4, 
                               seed=12345,
                               data=temp)

# IGF-1 Model (continuous predictors are unscaled)
mod_igf_unscaled <- rstanarm::stan_lmer(IGF1 ~  SexF + Twin + PopSize + poly(MumAge, degree = 2, raw = TRUE) + (1|ELISARunDate) + (1|PlateNumber) + (1|BirthYear) + (1|MumID),
                               cores=4, 
                               seed=12345,
                               data=temp)

# Model diagnostics
color_scheme_set("mix-teal-purple")
pp_check(mod_igf) + labs(subtitle="IGF-1")

rhat(mod_igf) %>%
  as.data.frame() %>%
  filter(.[[1]] > 1.1 & .[[1]] < 0.99)

mcmc_neff(neff_ratio(mod_igf), size = 2)

# Extract posterior draws for later use
posterior_igf <- as.array(mod_igf)

# Extract variance components (median; note: estimates differ a fair bit when you take mean instead)
posterior_igf_df <- as.data.frame(mod_igf)
med_igf_df <- median_hdci(posterior_igf_df)
p_value(mod_igf, method="hdi")

# Plot predictors
predictors <- as.array(mod_igf, pars = c("SexF1", "Twin1", "PopSize_sc", "poly(MumAge_sc, degree = 2, raw = TRUE)1", "poly(MumAge_sc, degree = 2, raw = TRUE)2"
                                         # "residBW"
                                         ))
p1 <- bayesplot::mcmc_areas(predictors, point_est = c("median"), prob = 0.5, prob_outer = 0.95) +
  theme_cowplot() +
  vline_0(color = "darkgray", linetype = 2) +
  scale_x_continuous(limits=c(-120,50)) +
  scale_y_discrete(labels=c("Maternal Age (2)", "Maternal Age(1)","Population Size", "Twin:1", "Sex:Female"), limits = rev) 
p1  

# Getting model predictions (conditional effects)
# Sex effect
sex_eff <- temp %>%
  modelr::data_grid(PopSize_sc = mean(PopSize_sc), 
                    MumAge_sc = mean(MumAge_sc),
                    HornType = c("Normal"),
                    SexF = unique(SexF), 
                    Twin=0) 

sex_eff <- convert_to_factor(sex_eff, cols_f)

sex_preds <- add_epred_draws(mod_igf, newdata = sex_eff, re_formula = NA) %>% 
  group_by(SexF)  %>% 
  median_hdci(.epred, .width = .95) 

p2 <- ggplot(sex_preds, aes(x = SexF)) +
  geom_jitter(data = temp, aes(x=SexF, y = IGF1, colour=SexF), size = 1.5,  alpha=0.3, height=0) +
  geom_pointrange(aes(ymin = .lower, ymax=.upper, y=.epred)) +
  theme_cowplot() +
  scale_colour_brewer(name = "Sex=Female", palette="Set2") +
  ylab("Plasma IGF-1 Concentration\n(ng/ml)") + 
  xlab("") +
  scale_x_discrete(labels=c("Male", "Female")) +
  panel_bg(fill = "gray95", color = NA) +
  grid_lines(color = "white") +
  theme(legend.text = element_text(size=18),
        legend.title = element_text(size=18)) +
  theme(legend.position="none")+
  theme(axis.text.x = element_text(size = 18)) +
  theme(axis.title.y = element_text(size = 18)) 

p2 


# Twin effect
twin_eff <- temp %>%
  modelr::data_grid(PopSize_sc = mean(PopSize_sc), 
                    MumAge_sc = mean(MumAge_sc),
                    HornType = c("Normal"),
                    Twin = unique(Twin), 
                    SexF=1) 

twin_eff <- convert_to_factor(twin_eff, cols_f)

twin_preds <- add_epred_draws(mod_igf, newdata = twin_eff, re_formula = NA) %>% 
  group_by(Twin)  %>% 
  median_hdci(.epred, .width = .95) 

p3 <- ggplot(twin_preds, aes(x = Twin)) +
  geom_jitter(data = temp, aes(x=Twin, y = IGF1, colour=Twin), size = 1.5,  alpha=1, height=0) +
  geom_pointrange(aes(ymin = .lower, ymax=.upper, y=.epred)) +
  theme_cowplot() +
  scale_colour_brewer(name = "Twin", palette="Pastel1") +
  ylab("") + 
  xlab("") +
  scale_x_discrete(labels=c("Singleton", "Twin")) +
  panel_bg(fill = "gray95", color = NA) +
  grid_lines(color = "white") +
  theme(legend.text = element_text(size=10),
        legend.title = element_text(size=10))+
  theme(legend.position="none") +
  theme(axis.text.x = element_text(size = 18)) +
  theme(axis.title.y = element_text(size = 18)) 


p3


# Pop effect
pop_eff <- temp %>%
  modelr::data_grid(PopSize = modelr::seq_range(PopSize, n = 200), 
                    MumAge = mean(MumAge),
                    HornType=c("Normal"),
                    SexF = 1, 
                    Twin=0) 

pop_eff <- convert_to_factor(pop_eff, cols_f)

pop_preds <- add_epred_draws(mod_igf_unscaled, newdata = pop_eff, re_formula = NA, ndraws=500) %>% 
  group_by(PopSize) 

pop_preds_med <- add_epred_draws(mod_igf_unscaled, newdata = pop_eff, re_formula = NA) %>% 
  group_by(PopSize) %>%
  median_hdci(.epred, .width = .95) 

p5 <- ggplot(pop_preds, aes(x = PopSize, y = IGF1)) +
  geom_jitter(data = temp, size = 1.5,  alpha=0.2, height=0, width=0.1, color = "#350") +
  geom_line(aes(y = .epred, group = .draw), alpha=1/15, color = "#298")  +
  geom_line(data=pop_preds_med,  aes(y=.epred), linewidth=1.2, colour="#355", alpha=0.8) +
  theme_cowplot() +
  scale_colour_brewer(name = "Population Size") +
  # ylab("") + 
  ylab("Plasma IGF-1 Concentration\n(ng/ml)") + 
  xlab("Population Size") +
  panel_bg(fill = "gray95", color = NA) +
  grid_lines(color = "white") +
  theme(legend.text = element_text(size=10),
        legend.title = element_text(size=10)) +
  theme(axis.text.x = element_text(size = 18)) +
  theme(axis.title.y = element_text(size = 18)) +
  theme(axis.title.x = element_text(size = 18)) 

p5

# # MumAge effect 
matage_eff <- temp %>%
  modelr::data_grid(MumAge = modelr::seq_range(MumAge, n = 200), 
                    PopSize = mean(PopSize),
                    HornType=c("Normal"),
                    SexF = 1, 
                    Twin=0) 

matage_eff <- convert_to_factor(matage_eff, cols_f)

matage_preds <- add_epred_draws(mod_igf_unscaled, newdata = matage_eff, re_formula = NA, ndraws=500) %>% 
  group_by(MumAge) 

matage_preds_med <- add_epred_draws(mod_igf_unscaled, newdata = matage_eff, re_formula = NA) %>% 
  group_by(MumAge) %>%
  median_hdci(.epred, .width = .95) 


p6 <- ggplot(matage_preds, aes(x = MumAge, y = IGF1)) +
  geom_jitter(data = temp, size = 1.5,  alpha=0.2, height=0, width=0.1, color = "#843") +
  geom_line(aes(y = .epred, group = .draw), alpha=1/15, color = "#950")  +
  geom_line(data=matage_preds_med,  aes(y=.epred), linewidth=1.2, colour="#900", alpha=0.8) +
  theme_cowplot() +
  scale_colour_brewer(name = "Maternal Age") +
  ylab("") + 
  xlab("Maternal Age (in years)") +
  panel_bg(fill = "gray95", color = NA) +
  grid_lines(color = "white") +
  theme(legend.text = element_text(size=10),
        legend.title = element_text(size=10)) +
  scale_x_continuous(breaks=c(1,3,5,7,9,11,13))+
  theme(axis.text.x = element_text(size = 18)) +
  theme(axis.title.y = element_text(size = 18)) +
  theme(axis.title.x = element_text(size = 18)) 

p6


# Arrange the plots
plot1 <- cowplot::plot_grid(p1 , 
                            nrow = 1,
                            labels = c("A"),
                            align = "h")

plot1
# ggsave("./Aug2024_IGFDataAnalysis/Fig1.tiff", plot1, dpi=600, width=6, height=4, bg="white" )


plot2 <- cowplot::plot_grid(p2, 
                            p3 , 
                            p5 ,
                            p6,
                            nrow = 2, 
                            ncol = 2,
                            labels = c("A", "B", "C", "D"),
                            align = "h")

plot2
# ggsave("./Fig2.tiff", plot2, dpi=600, width=12, height=9, bg="white" )


color_scheme_set("viridisA")
p1_pp <- pp_check(mod_igf, nreps=200) 

plot3 <- cowplot::plot_grid(igf_hist ,
                            p1_pp,
                            nrow = 2,
                            labels = c("A", "B"),
                            align = "vh")

plot3

# ggsave("./Aug2024_IGFDataAnalysis/FigS1_IGFHistPostDist.tiff", plot3, dpi=600, width=5, height=6, bg="white" )


# Extract variance components (median; note: estimates differ a fair bit when you take mean instead)
posterior_igf_df <- as.data.frame(mod_igf)
med_igf_df <- t(median_hdci(posterior_igf_df))
p_value(mod_igf, method="hdi")
vec_m <- as.data.frame(matrix(med_igf_df[c(1:18, 1030:1044),], ncol=3, nrow=11, byrow = TRUE))

vec_m$Parameter <- c("Intercept", "SexF", "Twin", "PopSize", "MumAge1", "MumAge2", "Resid_SD", "MumID_var", "Plate_var", "BirthYear_var", "ELISADate_var")
# write.table(vec_m, file = "./Aug2024_IGFDataAnalysis/Table1.csv", sep = ",",row.names = FALSE)

# Extract variance components (median; note: estimates differ a fair bit when you take mean instead)
posterior_igf_df <- as.data.frame(mod_igf_unscaled)
med_igf_df <- t(median_hdci(posterior_igf_df))
p_value(mod_igf_unscaled, method="hdi")
vec_m <- as.data.frame(matrix(med_igf_df[c(1:18, 1030:1044),], ncol=3, nrow=11, byrow = TRUE))

vec_m$Parameter <- c("Intercept", "SexF", "Twin", "PopSize", "MumAge1", "MumAge2", "Resid_SD", "MumID_var", "Plate_var", "BirthYear_var", "ELISADate_var")
# write.table(vec_m, file = "./Aug2024_IGFDataAnalysis/Table1_unscaled.csv", sep = ",",row.names = FALSE)

# Extract variance components (mean)
variance_table <- as.data.frame(VarCorr(mod_igf))
variance_table <- add_prop_var(variance_table)

# Model summary table 
modelsummary(mod_igf, statistic = "conf.int")
print(mod_igf, digits=2)
summary(mod_igf, probs=c(.025, .975), digits=2)

#######################################################
# Model looking at whether horn type predicts igf1
# Horn effect
# HornType model
temp_horns <- temp %>%
  filter(!is.na(HornType))

temp_horns <- standardize_columns(temp_horns, cols_sc)

mod_igf_horns <- rstanarm::stan_lmer(IGF1 ~  SexF + Twin + HornType + PopSize_sc + poly(MumAge_sc, degree = 2, raw = TRUE) + (1|ELISARunDate) + (1|PlateNumber) + (1|BirthYear) + (1|MumID),
                                     cores=4, 
                                     seed=12345,
                                     data=temp_horns)

horn_eff <- temp_horns %>%
  modelr::data_grid(PopSize_sc = mean(PopSize_sc), 
                    MumAge_sc = mean(MumAge_sc),
                    HornType = unique(HornType),
                    Twin = 0, 
                    SexF=1) 

horn_eff <- convert_to_factor(horn_eff, cols_f)

horn_preds <- add_epred_draws(mod_igf_horns, newdata = horn_eff, re_formula = NA) %>% 
  group_by(HornType)  %>% 
  median_hdci(.epred, .width = .95) 

p4 <- ggplot(horn_preds, aes(x = HornType)) +
  geom_jitter(data = temp_horns, aes(x=HornType, y = IGF1, colour=HornType), size = 1.5,  alpha=0.2, height=0) +
  geom_pointrange(aes(ymin = .lower, ymax=.upper, y=.epred)) +
  theme_cowplot() +
  scale_color_manual(values=c("#089", "#450", "#620")) +
  # scale_colour_brewer(name = "HornType") +
  ylab("Plasma IGF-1 Concentration\n(ng/ml)") + 
  xlab("Horn Type") +
  # scale_x_discrete(labels=c("Normal", "Polled", "Scurred")) +
  panel_bg(fill = "gray95", color = NA) +
  grid_lines(color = "white") +
  theme(legend.text = element_text(size=10),
        legend.title = element_text(size=10))+
  theme(legend.position="none") +
  theme(axis.text.x = element_text(size = 18)) +
  theme(axis.title.y = element_text(size = 18)) +
  theme(axis.title.x = element_text(size = 18)) 


p4
# ggsave("./Aug2024_IGFDataAnalysis/FigS1_Horns.tiff", p4, dpi=600, width=6, height=4, bg="white" )

# Extract variance components (median; note: estimates differ a fair bit when you take mean instead)
posterior_igf_df <- as.data.frame(mod_igf_horns)
med_igf_df <- t(median_hdci(posterior_igf_df))
p_value(mod_igf_horns, method="hdi")
vec_m <- as.data.frame(matrix(med_igf_df[c(1:24, 1036:1050),], ncol=3, nrow=13, byrow = TRUE))

vec_m$Parameter <- c("Intercept", "SexF", "Twin", "HTPolled", "HTScurred", "PopSize", "MumAge1", "MumAge2", "Resid_SD", "MumID_var", "Plate_var", "BirthYear_var", "ELISADate_var")
# write.table(vec_m, file = "./Aug2024_IGFDataAnalysis/Table1_horntype.csv", sep = ",",row.names = FALSE)

########################################################
# Add models explroing storage time effects 
# IGF-1 Model
# IGF-1 Model (continuous predictors are unscaled)
# Rescale variables 
cols_st <- c("StorageTime")
temp <- standardize_columns(temp, cols_st)

mod_igf_unscaled_st <- rstanarm::stan_lmer(IGF1 ~  SexF + Twin + PopSize + StorageTime + poly(MumAge, 2) + (1|ELISARunDate) + (1|PlateNumber) + (1|BirthYear) + (1|MumID),
                                           cores=4, 
                                           seed=12345,
                                           data=temp)

mod_igf_st <- rstanarm::stan_lmer(IGF1 ~  SexF + Twin + PopSize_sc + StorageTime_sc + poly(MumAge_sc, 2) + (1|ELISARunDate) + (1|PlateNumber) + (1|BirthYear) + (1|MumID),
                                        cores=4, 
                                        seed=12345,
                                        data=temp)

# Pop effect
storagetime_eff <- temp %>%
  modelr::data_grid(StorageTime = modelr::seq_range(StorageTime, n = 200), 
                    PopSize = mean(PopSize),
                    MumAge = mean(MumAge),
                    HornType=c("Normal"),
                    SexF = 1, 
                    Twin=0) 

storagetime_eff <- convert_to_factor(storagetime_eff, cols_f)

st_preds <- add_epred_draws(mod_igf_unscaled_st, newdata = storagetime_eff, re_formula = NA, ndraws=500) %>% 
  group_by(StorageTime) 

st_preds_med <- add_epred_draws(mod_igf_unscaled_st, newdata = storagetime_eff, re_formula = NA) %>% 
  group_by(StorageTime) %>%
  median_hdci(.epred, .width = .95) 

p7 <- ggplot(st_preds, aes(x = StorageTime, y = IGF1)) +
  geom_jitter(data = temp, size = 1.5,  alpha=0.2, height=0, width=0.1, color = "#EB984E") +
  geom_line(aes(y = .epred, group = .draw), alpha=1/15, color = "#EC7063")  +
  geom_line(data=st_preds_med,  aes(y=.epred), linewidth=1.2, colour="#CB4335", alpha=0.8) +
  theme_cowplot() +
  scale_colour_brewer(name = "Storage Time (in years)") +
  # ylab("") + 
  ylab("Plasma IGF-1 Concentration\n(ng/ml)") + 
  xlab("Storage Time (in years)") +
  # annotate("text", x=7, y=1200, label= expression(beta[StorageTime]~"= 7.559 (95% CI: 0.788"~-~"14.025)")) +
  panel_bg(fill = "gray95", color = NA) +
  grid_lines(color = "white") +
  scale_x_continuous(n.breaks=10) +
  theme(legend.text = element_text(size=10),
        legend.title = element_text(size=10))

p7

# ggsave("Dec24_IGFStorageTime.tiff", p7, dpi=600, height = 6, width=8, bg="white" )

# Extract posterior draws for later use
posterior_igf_st <- as.array(mod_igf_st)

# Extract variance components (median; note: estimates differ a fair bit when you take mean instead)
posterior_igf_st_df <- as.data.frame(mod_igf_st)
med_igf_st_df <- as.data.frame(t(median_hdci(posterior_igf_st_df)))
p_value(posterior_igf_st_df, method="hdi")
med_igf_st_df$V1 <- as.numeric(med_igf_st_df$V1)
med_igf_st_df$squared <- med_igf_st_df$V1*med_igf_st_df$V1

# Extract posterior draws for later use
posterior_igf_st <- as.array(mod_igf_unscaled_st)

# Extract variance components (median; note: estimates differ a fair bit when you take mean instead)
posterior_igf_st_df <- as.data.frame(mod_igf_unscaled_st)
med_igf_st_df1 <- as.data.frame(t(median_hdci(posterior_igf_st_df)))

# Model summary table 
modelsummary(mod_igf_unscaled_st, statistic = "conf.int")
print(mod_igf_unscaled_st, digits=2)
summary(mod_igf_unscaled_st, probs=c(.025, .975), digits=2)
summary(mod_igf_st, probs=c(.025, .975), digits=3)
# summary(mod_igf, probs=c(.005, .995), digits=3)
p_value( mod_igf_st, method="hdi")


########################################################
# Add models explroing sample collection time effects 
# IGF-1 Model
# IGF-1 Model (continuous predictors are unscaled)
temp$SampleCollectionTimeofDay <- temp$CapHour_Tertiles
temp$SampleCollectionTimeofDay <- as.factor(temp$SampleCollectionTimeofDay)
levels(temp$SampleCollectionTimeofDay)
temp_toc <- temp %>%
  filter(!is.na(SampleCollectionTimeofDay))
temp_toc <- standardize_columns(temp_toc, cols_sc)

mod_igf_unscaled_toc <- rstanarm::stan_lmer(IGF1 ~  SexF + Twin + PopSize + SampleCollectionTimeofDay + poly(MumAge, 2) + (1|ELISARunDate) + (1|PlateNumber) + (1|BirthYear) + (1|MumID),
                                           cores=4, 
                                           seed=12345,
                                           data=temp_toc)

mod_igf_toc <- rstanarm::stan_lmer(IGF1 ~  SexF + Twin + PopSize_sc + SampleCollectionTimeofDay + poly(MumAge_sc, 2) + (1|ELISARunDate) + (1|PlateNumber) + (1|BirthYear) + (1|MumID),
                                            cores=4, 
                                            seed=12345,
                                            data=temp_toc)

# Extract variance components (median; note: estimates differ a fair bit when you take mean instead)
posterior_igf_toc_df <- as.data.frame(mod_igf_toc)
med_igf_toc_df <- as.data.frame(t(median_hdci(posterior_igf_toc_df)))
p_value(posterior_igf_toc_df, method="hdi")

# Horn effect
toc_eff <- temp_toc %>%
  modelr::data_grid(PopSize = mean(PopSize), 
                    MumAge = mean(MumAge),
                    SampleCollectionTimeofDay = unique(SampleCollectionTimeofDay),
                    HornType = c("Normal"),
                    Twin = 0, 
                    SexF=1) 

toc_eff <- convert_to_factor(toc_eff, cols_f)

toc_preds <- add_epred_draws(mod_igf_unscaled_toc, newdata = toc_eff, re_formula = NA) %>% 
  group_by(SampleCollectionTimeofDay)  %>% 
  median_hdci(.epred, .width = .95) 

p8 <- ggplot(toc_preds, aes(x = SampleCollectionTimeofDay)) +
  geom_jitter(data = temp_toc, aes(x=SampleCollectionTimeofDay, y = IGF1, colour=SampleCollectionTimeofDay), size = 1.5,  alpha=0.2, height=0) +
  geom_pointrange(aes(ymin = .lower, ymax=.upper, y=.epred)) +
  theme_cowplot() +
  scale_color_manual(values=c("#089", "#450", "#620")) +
  # scale_colour_brewer(name = "HornType") +
  ylab("Plasma IGF-1 Concentration\n(ng/ml)") + 
  xlab("Sample Collection Time of Day") +
  # scale_x_discrete(labels=c("Normal", "Polled", "Scurred")) +
  panel_bg(fill = "gray95", color = NA) +
  grid_lines(color = "white") +
  theme(legend.text = element_text(size=10),
        legend.title = element_text(size=10))+
  theme(legend.position="none") +
  theme(axis.text.x = element_text(size = 14))

p8
# ggsave("Dec24_IGFSampleCollectionTOC.tiff", p8, dpi=600, height = 6, width=8, bg="white" )

# Extract variance components (median; note: estimates differ a fair bit when you take mean instead)
posterior_igf_toc_df <- as.data.frame(mod_igf_unscaled_toc)
med_igf_toc_df <- median_hdci(posterior_igf_toc_df)
p_value(posterior_igf_st_df, method="hdi")

# Model summary table 
modelsummary(mod_igf_unscaled_st, statistic = "conf.int")
print(mod_igf_unscaled_st, digits=2)
summary(mod_igf_unscaled_st, probs=c(.025, .975), digits=2)