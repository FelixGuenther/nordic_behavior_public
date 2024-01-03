# _targets.R file
library(targets)
library(stantargets)
library(tarchetypes)
library(tidyverse)
library(lubridate)
library(DescTools)
library(readxl)

source("R/functions.R")
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("tidyverse", "lubridate", "bayesplot", 
                            "DescTools", "splines", "zoo", "readxl"),
               error = "abridge",
               memory = "transient",
               garbage_collection = TRUE)

library(future)
library(future.callr)
plan(callr)

list(
  # Fit model
  # Germany
  tar_stan_mcmc(
    name = ger_comb_redu_imp,
    stan_files = c("./stan/mod_add_ger_redu_imp.stan"),
    data = prepare_data_ger(regions_ger = 1:2,
                            gamma=1/6,
                            num_week = 44,
                            categories = c("grocery_and_pharmacy",  
                                           "retail_and_recreation",
                                           "transit_stations",
                                           "workplaces"),
                            mult = c("Temperature")),
    parallel_chains = 4,
    seed = 12213,
    iter_warmup = 200, iter_sampling = 800,
    adapt_delta = 0.99,
    max_treedepth = 15,
    refresh=250,
    init = list(list(u_rt=array(0.02, dim = c(2,75)), 
                     sd_u_rt=rep(0.09, 2),
                     phi_star=array(0.01, dim=c(2,4))),
                list(u_rt=array(0.015, dim = c(2,75)), 
                     sd_u_rt=rep(0.1, 2),
                     phi_star=array(0.005, dim=c(2,4))),
                list(u_rt=array(0.03, dim = c(2,75)), 
                     sd_u_rt=rep(0.07,2),
                     phi_star=array(0.0075, dim=c(2,4))),
                list(u_rt=array(0.01, dim = c(2,75)), 
                     sd_u_rt=rep(0.08,2),
                     phi_star=array(0.01, dim=c(2,4))))
  ),
  # Summary and plots
  # Model fit
  tar_target(plot_fit_bavaria_berlin,
             fit_plot_ger(mcmc_summary_ger = ger_comb_redu_imp_summary_mod_add_ger_redu_imp, 
                          data_ger = ger_comb_redu_imp_data)),
  # Estimated multiplicative effect
  tar_target(mult_eff,
             ger_comb_redu_imp_summary_mod_add_ger_redu_imp %>% 
               filter(variable == "phi_mult[1]") %>% 
               select(median, q5, q95, rhat) %>%
               mutate(median = exp(20*median),
                      q5 = exp(20*q5),
                      q95=exp(20*q95))),
  
  # Explained variability in contacts
  tar_target(post_rel_ger_redu_imp,
             post_rel_cont(mcmc = ger_comb_redu_imp_mcmc_mod_add_ger_redu_imp,
                           data = ger_comb_redu_imp_data, 
                           n_weeks = ger_comb_redu_imp_data$n_weeks,
                           n_days = ger_comb_redu_imp_data$n_days)),
  tar_target(post_rel_ger_redu_imp_tot,
             post_rel_cont_ctry(mcmc = ger_comb_redu_imp_mcmc_mod_add_ger_redu_imp,
                                data = ger_comb_redu_imp_data,
                                n_weeks = ger_comb_redu_imp_data$n_weeks,
                                n_days = ger_comb_redu_imp_data$n_days)),
  
  # Plot multiplicative effects
  tar_target(phi_kr_tab_ger_redu_imp,
             phi_kr_tab(mcmc = ger_comb_redu_imp_mcmc_mod_add_ger_redu_imp,
                        data = ger_comb_redu_imp_data, 
                        region_names = tibble(county_code = c("DE-BY", "DE-BE"),
                                              county_name = c("Bavaria", "Berlin")))),
  
  # Relative share of categories at BL Germany
  tar_target(rel_share_ger,
             share_bl_ger(mcmc = ger_comb_redu_imp_mcmc_mod_add_ger_redu_imp,
                          data = ger_comb_redu_imp_data)
             ),
  tar_target(fig_5,
             rel_share_ger %>% group_by(reg, cov) %>%
               summarise(share_med = median(share),
                         q05 = quantile(share, 0.05),
                         q95 = quantile(share, 0.95)) %>%
               mutate(cov = factor(cov, levels = c("workplaces", "transit_stations", 
                                                   "grocery_and_pharmacy", "retail_and_recreation"),
                                   labels = c("Work", 
                                              "Transit", 
                                              "Groc. & phar.",
                                              "Retail & recr.")),
                      reg = factor(reg, levels = c("DE-BY", "DE-BE"),
                                   labels = c("Bavaria", "Berlin"))) %>%
               ggplot() + 
               geom_col(aes(cov, share_med, fill=reg, group=reg),  position = position_dodge2(),
                        col = "black") +
               geom_errorbar(aes(cov, ymin=q05, ymax=q95, group=reg), position = position_dodge2()) +
               ylab("Relative contribution") +
               xlab("") +
               theme(legend.title = element_blank(), legend.position = "bottom"))
)  