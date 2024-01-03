#' Function that reads in prepared data files and returns data-list as 
#' expected by STAN model SIR_age_group_weekly
#' 
#' @param gamma inverse of infectious period
#' @param regions which regions to consider (integer index), defaults to NULL (equals all regions)
#' @param num_week number of weeks to consider in model fitting, defaults to NULL (all weeks)
#' @return list of data objects as expected by STAN model SIR_age_group_weekly_r



#' Function that reads in prepared data files and returns data-list as 
#' expected by STAN model SIR_age_group_weekly
#' 
#' @param gamma inverse of infectious period
#' @param regions which regions to consider (integer index), defaults to NULL (equals all regions)
#' @param num_week number of weeks to consider in model fitting, defaults to NULL (all weeks)
#' @return list of data objects as expected by STAN model SIR_age_group_weekly_r
prepare_data_ger = function(gamma=1/6,
                            num_week = 10,
                            regions_ger = 1:2,
                            categories = c("retail_and_recreation"),
                            mult = c("Temperature"),
                            vacc_eff = 0.9) {
  
  days = as.character(seq(lubridate::ymd("2020-02-24"),
                          lubridate::ymd("2020-12-27"),
                          "1 day"))
  # Variants
  variants_ger = readRDS("./data/VariantsGermany.rds")
  variants_ger = variants_ger[,days,]
  # Covariates
  covariates_ger = readRDS("./data/CovariatesGermany_wParks.rds")

  covariates_ger = covariates_ger[,days,]
  
  # Multiplicative effects design matrix
  # Temperature
  temp_ger = covariates_ger[,,"Temperature", drop=FALSE]
  for (i in 1:dim(temp_ger)[1]) {
    temp_ger[i,, "Temperature"] = (temp_ger[i,, "Temperature"] - temp_ger[i,1, "Temperature"])
  }

  mult_cov_ger = abind::abind(temp_ger, variants_ger)[,,mult, drop=FALSE]
  
  
  beh_cat_ger = covariates_ger[,,categories, drop=FALSE]/100
  
  for (i in 1:dim(beh_cat_ger)[1]) {
    for(j in 1:dim(beh_cat_ger)[3]) {
      c1 = beh_cat_ger[i,,j][1]
      cT = rev(beh_cat_ger[i,,j])[1]
      beh_cat_ger[i,,j] = zoo::rollmean(beh_cat_ger[i,,j], k = 7, align = "right", 
                                        fill = c(c1, NA, cT))
    }
  }
  
  
  pop_r_age_ger = readRDS("./data/PopulationGermany.rds")

  n_vaccinated_ger = readRDS("./data/VaccineGermany.rds")
  n_vaccinated_ger = n_vaccinated_ger[,days,]

  p_hosp_age_ger = readRDS("./data/HospRisk_Wildtype_Germany.rds")

  tau_hosp = readRDS("./data/TimeToHospitalGermany.rds")[1,]
  n_agegroup_ger=apply(pop_r_age_ger, 1, function(x) sum(!is.na(x)))

  
  hosp_dat_raw_ger = read_csv("./data/2023-04-29_Deutschland_COVID-19-Hospitalisierungen.csv")
  
  hosp_dat_raw_ger = hosp_dat_raw_ger %>% filter(Altersgruppe=="00+") %>% 
    group_by(Bundesland) %>%
    filter(Datum <= "2021-01-01") %>%
    filter(Datum %in% seq(lubridate::ymd("2020-03-01"), 
                          lubridate::ymd("2021-01-01"), 
                          "7 days")) %>%
    arrange(Bundesland, Datum) %>%
    ungroup() %>%
    filter(Bundesland %in% c("Berlin", "Bayern")) %>%
    group_by(Bundesland) %>%
    arrange(Datum) %>%
    mutate(week = 1:n()) %>%
    ungroup() %>%
    filter(week<=num_week) %>%
    select(Datum, Bundesland, 
           cases = `7T_Hospitalisierung_Faelle`) %>%
    mutate(Bundesland = factor(Bundesland, levels = c("Berlin", "Bayern"))) %>%
    arrange(Datum, Bundesland)
  
  
  hosp_data_ger = array(hosp_dat_raw_ger$cases,
                        dim = c(length(unique(regions_ger)),
                                num_week),
                        dimnames = list(unique(hosp_dat_raw_ger$Bundesland),
                                            1:num_week))
  
  # Restrict time
  mult_cov_ger = mult_cov_ger[,1:(num_week*7),,drop=FALSE]

  beh_cat_ger = beh_cat_ger[,1:(num_week*7),,drop=FALSE]
  # Remove negative vaccination counts, restrict to model period, apply vaccine efficacy
  n_vaccinated_ger = n_vaccinated_ger[,1:(num_week*7),,drop=FALSE] * vacc_eff
  
  # Hospital data
  hosp_data_ger = hosp_data_ger[,1:num_week,drop=FALSE]

  # Padded arrays
  n_vaccinated = n_vaccinated_ger
  
  beh_cat = beh_cat_ger
  
  mult_cov = mult_cov_ger
  
  n_agegroup = n_agegroup_ger
  pop_r_age = pop_r_age_ger
  
  p_hosp_age = p_hosp_age_ger
  hosp_data = hosp_data_ger
  
  # Select regions
  n_vaccinated = n_vaccinated[c(regions_ger),,,drop=FALSE]
  beh_cat = beh_cat[c(regions_ger),,,drop=FALSE]
  mult_cov = mult_cov[c(regions_ger),,,drop=FALSE]
  n_agegroup = n_agegroup[c(regions_ger),drop=FALSE]
  pop_r_age = pop_r_age[c(regions_ger),,drop=FALSE]
  p_hosp_age = p_hosp_age[c(regions_ger),,drop=FALSE]
  hosp_data = hosp_data[c(regions_ger),,drop=FALSE]
  
  list(n_days=dim(n_vaccinated)[2],
       n_weeks=dim(n_vaccinated)[2]/7,
       n_regions=dim(n_vaccinated)[1],
       n_cat = dim(beh_cat)[3],
       beh_cat=beh_cat,
       n_mult=dim(mult_cov)[3],
       mult_cov=mult_cov,
       n_agegroup=n_agegroup, 
       n_agegroup_max=max(n_agegroup),
       pop_r_age = pop_r_age[,1:max(n_agegroup), drop=FALSE], # Subset if padding is unneccessary
       n_vaccinated = n_vaccinated[,,1:max(n_agegroup), drop=FALSE],
       p_hosp_age = p_hosp_age[,1:max(n_agegroup), drop=FALSE],
       gamma = gamma,
       hosp_cutoff = length(tau_hosp)-1,
       tau_hosp=tau_hosp,
       hosp_data=hosp_data)
}

fit_plot_ger = function(mcmc_summary_ger, 
                        data_ger, 
                        reg_names = tibble(county_code = c("DE-BY", "DE-BE"), county_name=c("Bavaria", "Berlin")),
                        reg_sel = c("Bavaria", "Berlin")) {
  # Panel A Fit to data
  cust_theme = theme(legend.text = element_text(size=14), axis.text = element_text(size=14), axis.title = element_text(size=14), plot.title = element_text(size = 16))
  mcmc_summary = mcmc_summary_ger %>% filter(grepl("lambda_hosp_week", variable)) %>%
    mutate(indices =  gsub("\\]", "", gsub(pattern = "lambda_hosp_week\\[", "", variable)),
           region = sapply(indices, function(x) {str_split(x, ",")[[1]][1]}),
           week = as.numeric(sapply(indices, function(x) {str_split(x, ",")[[1]][2]})),
           mean = pmax(mean,0.5),
           q5=pmax(q5, 0.5),
           q95=pmax(q95,0.5),
           region=factor(region, levels=1:data_ger$n_regions, 
                         labels = dimnames(data_ger$hosp_data)[[1]]),
           region = factor(region, levels = c("Bayern", "Berlin"),
                           labels = c("DE-BY", "DE-BE"))) %>%
    left_join(reg_names %>% rename(region=county_code))
  
  mcmc_summary = mcmc_summary %>% filter(county_name %in% reg_sel) %>%
    mutate(date = ymd("2020-03-01")+(week-1)*7)
  
  hosp_dat = data_ger$hosp_data %>% as.data.frame() %>% 
    mutate(region=dimnames(data_ger$hosp_data)[[1]],
           region = factor(region, levels = c("Bayern", "Berlin"),
                           labels = c("DE-BY", "DE-BE"))) %>% 
    left_join(reg_names %>% rename(region=county_code)) %>%
    select(-region) %>%
    pivot_longer(-county_name, names_to = "week", values_to = "hosp_obs")  %>% 
    mutate(hosp_obs=pmax(hosp_obs, 0.5),
           week=as.numeric(week)) %>%
    filter(county_name %in% reg_sel) %>%
    mutate(date = ymd("2020-03-01")+(week-1)*7)
  
  fit_plot = mcmc_summary %>%
    ggplot() +
    geom_point(aes(date, mean), alpha=.75) +
    geom_errorbar(aes(date, ymin=q5, ymax=q95)) +
    geom_point(aes(date, hosp_obs), data = hosp_dat, alpha=.75, col = "darkgreen") + 
    xlab("Week") +
    ylab("New hosp. adm. (log2)") +
    scale_y_continuous(trans="log2", breaks=c(1,4,16,64, 256, 2^10)) +
    scale_x_date(date_breaks = "3 months", date_labels = "%b %y") +
    facet_wrap(~county_name) +
    theme_bw() +
    cust_theme
  # Panel B, R(t)
  r_t_sum = mcmc_summary_ger %>%
    filter(grepl("r_t", variable)) %>%
    mutate(region = as.numeric(str_replace(str_extract(variable, "\\[(\\d)+"), "\\[", "")),
           day = as.numeric(str_replace(str_extract(variable, "(\\d)+\\]"), "\\]", ""))) %>%
    select(mean, median, q5, q95, region, day) %>%
    mutate(date = ymd("2020-03-01") + day - 1 - 7,
           region=factor(region, levels=1:data_ger$n_regions, 
                         labels = dimnames(data_ger$hosp_data)[[1]]),
           region = factor(region, levels = c("Bayern", "Berlin"),
                           labels = c("DE-BY", "DE-BE"))) %>%
    left_join(reg_names %>% rename(region=county_code)) %>%
    filter(county_name %in% reg_sel)
  r_t_plot = r_t_sum %>%
    ggplot() +
    geom_line(aes(date, median, col = county_name)) +
    geom_ribbon(aes(date, ymin=q5, ymax=q95, fill = county_name), alpha=.3) +
    scale_y_continuous(trans="log2") +
    geom_hline(aes(yintercept=1), lty=2) +
    theme_bw() +
    theme(legend.position = "top", legend.title = element_blank()) +
    ylab("R(t)") +
    xlab("") +
    scale_x_date(breaks="3 month", date_labels = "%b %y") +
    cust_theme +
    scale_color_brewer(type = "qual", palette = 2) +
    scale_fill_brewer(type = "qual", palette = 2)
  
  # Panel C: R(t) contribution
  r_t_cont_time = mcmc_summary_ger %>%
    filter(grepl("share_beta_[1-4]_ex", variable)) %>%
    filter(!grepl("mean", variable)) %>%
    mutate(region = as.numeric(str_replace(str_extract(variable, "\\[(\\d)+"), "\\[", "")),
           day = as.numeric(str_replace(str_extract(variable, "(\\d)+\\]"), "\\]", "")),
           var = factor(substr(variable, 1, 15), 
                        levels = c("share_beta_1_ex", "share_beta_2_ex", "share_beta_3_ex", "share_beta_4_ex"),
                        labels = c("grocery_and_pharmacy",  "retail_and_recreation", "transit_stations",     
                                   "workplaces"))) %>%
    select(region, median, day, var) %>%
    mutate(date = ymd("2020-03-01") + day - 1 - 7,
           region=factor(region, levels=1:data_ger$n_regions, 
                         labels = dimnames(data_ger$hosp_data)[[1]]),
           region = factor(region, levels = c("Bayern", "Berlin"),
                           labels = c("DE-BY", "DE-BE"))) %>%
    pivot_wider(id_cols = c("region", "date"), values_from = "median", names_from = "var") %>%
    left_join(reg_names %>% rename(region=county_code)) %>%
    right_join(r_t_sum %>%
                 select(region, date, county_name, r_t = median)) %>%
    mutate(other = (1-grocery_and_pharmacy - retail_and_recreation - transit_stations - workplaces)*r_t,
           grocery_and_pharmacy = grocery_and_pharmacy*r_t,
           retail_and_recreation = retail_and_recreation*r_t,
           transit_stations = transit_stations*r_t,
           workplaces = workplaces*r_t) %>%
    pivot_longer(cols = c("grocery_and_pharmacy",  "retail_and_recreation", "transit_stations",     
                          "workplaces", "other"), names_to = "category", values_to = "r_t_cont") %>%
    mutate(category = factor(category, 
                             levels = rev(c("grocery_and_pharmacy",  "retail_and_recreation", "transit_stations",     
                                            "workplaces", "other")),
                             labels = c("Other", "Work", "Transit", "Retail & recr.", "Groc. & pharm.")))
  
  r_t_cont_plot = r_t_cont_time %>%
    ggplot() + 
    geom_col(aes(date, r_t_cont, fill = category)) +
    facet_wrap(~county_name, ncol=2) +
    geom_hline(aes(yintercept=1), lty=2) +
    ylab("R(t)") +
    xlab("") +
    scale_x_date(breaks="3 month", date_labels = "%b %y") +
    scale_fill_brewer(palette = 'Set1') +
    theme_bw() + 
    cust_theme +
    theme(legend.position = "top", legend.title = element_blank())
  
  plot = ggpubr::ggarrange(fit_plot, r_t_plot, r_t_cont_plot, nrow = 3, labels = "AUTO")
  plot
}

post_rel_cont = function(mcmc, 
                         data,
                         n_beh_cat = 4, 
                         n_weeks = 75,
                         n_days = 525,
                         cred_level = 0.9) {
  post_per_reg = function(region=1, 
                          mcmc_ob, 
                          data_ob, 
                          n_beh_cat = 4, 
                          n_weeks = 75,
                          n_days = 525,
                          cred_level = 0.9) {
    phi_rk_df = mcmc_ob$draws(
      variables = paste0("phi_star[", region, ",", 1:n_beh_cat,"]"), 
      format="draws_df") %>% tibble() %>% 
      select(-c(".chain", ".iteration", ".draw"))
    n_draws = dim(phi_rk_df)[1]
    
    u_t_df = mcmc_ob$draws(
      variables = paste0("u_rt[", region, ",", 1:n_weeks,"]"), 
      format="draws_df") %>% tibble() %>% 
      select(-c(".chain", ".iteration", ".draw"))
    
    eff_r_t = mcmc_ob$draws(
      variables = paste0("r_t[", region, ",", 1:n_days,"]"), 
      format="draws_df") %>% tibble() %>% 
      select(-c(".chain", ".iteration", ".draw"))
    
    i_t_df = mcmc_ob$draws(
      variables = paste0(rep(paste0("I[", region, ",", 1:n_days), each=7), ",",1:7,"]"), 
      format="draws_df") %>% tibble() %>% 
      select(-c(".chain", ".iteration", ".draw")) %>% 
      mutate(draw=1:n_draws) %>% 
      pivot_longer(-draw) %>%
      separate(name, sep=",", into = c("ind", "t", "a")) %>% 
      mutate(t=as.numeric(t)) %>%
      group_by(draw, t) %>% 
      summarise(value = sum(value)) %>%
      pivot_wider(id_cols = draw, names_from = t, names_prefix = "I_") %>%
      arrange(draw) %>% ungroup() %>%
      select(-draw)
    
    
    mean_pi_t = function(phi_rk, x_tk, u_t) {
      eta_t = (1+x_tk) %*% phi_rk %>% as.vector() 
      u_t = rep(u_t, each = 7)
      pi_t = eta_t / (eta_t + u_t)
      mean(pi_t)
    }
    
    share_sec_inf = function(phi_rk, x_tk, u_t, eff_r_t, i_t, gamma) {
      eta_t = (1+x_tk) %*% phi_rk %>% as.vector() 
      u_t = rep(u_t, each = 7)
      pi_t = eta_t / (eta_t + u_t)
      sum(pi_t * eff_r_t * i_t * gamma)/sum(eff_r_t * i_t * gamma) 
    }
    
    share_mcmc_ob_draws = sapply(1:n_draws,
                                 function(x) mean_pi_t(phi_rk = unlist(phi_rk_df[x,]),
                                                       x_tk = data_ob$beh_cat[region,,, drop = TRUE],
                                                       u_t = unlist(u_t_df[x,])))
    
    share_sec_inf_draws = sapply(1:n_draws,
                                 function(x) share_sec_inf(phi_rk = unlist(phi_rk_df[x,]),
                                                           x_tk = data_ob$beh_cat[region,,, drop = TRUE],
                                                           u_t = unlist(u_t_df[x,]),
                                                           eff_r_t = unlist(eff_r_t[x,]),
                                                           i_t = unlist(i_t_df[x,]),
                                                           gamma = data_ob$gamma))
    tibble(region = region,
           stat = c("mean_cont_day", "share_sec_inf"),
           mean = c(mean(share_mcmc_ob_draws), mean(share_sec_inf_draws)),
           median = c(median(share_mcmc_ob_draws), median(share_sec_inf_draws)),
           ci_lwr = c(quantile(share_mcmc_ob_draws, p = (1-cred_level)/2),
                      quantile(share_sec_inf_draws, p = (1-cred_level)/2)),
           ci_upr = c(quantile(share_mcmc_ob_draws, p = (1-(1-cred_level)/2)),
                      quantile(share_sec_inf_draws, p = (1-(1-cred_level)/2))))
  }
  
  res = purrr::map_df(1:data$n_regions, function(x) {
    post_per_reg(region=x, mcmc_ob = mcmc, 
                 data_ob = data,
                 n_beh_cat = n_beh_cat, 
                 n_weeks = n_weeks,
                 n_days = n_days,
                 cred_level = cred_level)
  })
  res %>% mutate(region = factor(region, levels = 1:data$n_regions,
                                 labels = dimnames(data$beh_cat)[[1]]))
}


post_rel_cont_ctry = function(mcmc, 
                              data,
                              n_beh_cat = 4, 
                              n_weeks = 75,
                              n_days = 525,
                              cred_level = 0.9) {
  rel_cont_reg_draw = function(region=1, 
                               mcmc_ob, 
                               data_ob, 
                               n_beh_cat = 4, 
                               n_weeks = 75,
                               n_days = 525,
                               cred_level = 0.9) {
    phi_rk_df = mcmc_ob$draws(
      variables = paste0("phi_star[", region, ",", 1:n_beh_cat,"]"), 
      format="draws_df") %>% tibble() %>% 
      select(-c(".chain", ".iteration", ".draw"))
    n_draws = dim(phi_rk_df)[1]
    
    u_t_df = mcmc_ob$draws(
      variables = paste0("u_rt[", region, ",", 1:n_weeks,"]"), 
      format="draws_df") %>% tibble() %>% 
      select(-c(".chain", ".iteration", ".draw"))
    
    eff_r_t = mcmc_ob$draws(
      variables = paste0("r_t[", region, ",", 1:n_days,"]"), 
      format="draws_df") %>% tibble() %>% 
      select(-c(".chain", ".iteration", ".draw"))
    
    i_t_df = mcmc_ob$draws(
      variables = paste0(rep(paste0("I[", region, ",", 1:n_days), each=7), ",",1:7,"]"), 
      format="draws_df") %>% tibble() %>% 
      select(-c(".chain", ".iteration", ".draw")) %>% 
      mutate(draw=1:n_draws) %>% 
      pivot_longer(-draw) %>%
      separate(name, sep=",", into = c("ind", "t", "a")) %>% 
      mutate(t=as.numeric(t)) %>%
      group_by(draw, t) %>% 
      summarise(value = sum(value)) %>%
      pivot_wider(id_cols = draw, names_from = t, names_prefix = "I_") %>%
      arrange(draw) %>% ungroup() %>%
      select(-draw)
    
    sec_inf = function(phi_rk, x_tk, u_t, eff_r_t, i_t, gamma) {
      eta_t = (1+x_tk) %*% phi_rk %>% as.vector() 
      u_t = rep(u_t, each = 7)
      pi_t = eta_t / (eta_t + u_t)
      tibble(sec_inf_cat = sum(pi_t * eff_r_t * i_t * gamma),
             sec_inf_tot = sum(eff_r_t * i_t * gamma)) 
    }
    
    
    sec_inf_draws = do.call(rbind, 
                            lapply(1:n_draws,
                                   function(x) 
                                     sec_inf(phi_rk = unlist(phi_rk_df[x,]),
                                             x_tk = data_ob$beh_cat[region,,, drop = TRUE],
                                             u_t = unlist(u_t_df[x,]),
                                             eff_r_t = unlist(eff_r_t[x,]),
                                             i_t = unlist(i_t_df[x,]),
                                             gamma = data_ob$gamma)))
    sec_inf_draws %>% mutate(region, draw = 1:n_draws)
  }
  
  res = purrr::map_df(1:data$n_regions, function(x) {
    rel_cont_reg_draw(region=x, 
                       mcmc_ob = mcmc, 
                       data_ob = data,
                       n_beh_cat = n_beh_cat, 
                       n_weeks = n_weeks,
                       n_days = n_days,
                       cred_level = cred_level)
  })
  res %>% group_by(draw) %>% 
    summarise(sec_inf_cat = sum(sec_inf_cat), 
              sec_inf_tot = sum(sec_inf_tot), 
              share_cat = sec_inf_cat/sec_inf_tot) %>% 
    summarise(share_cat_mean = mean(share_cat), 
              share_cat_med = median(share_cat),
              ci_lwr = quantile(share_cat, (1-cred_level)/2),
              ci_upr = quantile(share_cat, 1-((1-cred_level)/2)))
}


# phi_kr table
phi_kr_tab = function(mcmc, data, region_names) {
  phi_star_draws = mcmc$draws("phi_star", format="draws_df")
  
  draws_long = phi_star_draws %>% select(-c(.chain, .iteration)) %>% 
    pivot_longer(cols = starts_with("phi"))
  
  draws_long %>% 
    mutate(indices = gsub("\\]", "", gsub("phi_star\\[", "", name)),
           reg = purrr::map_chr(indices, function(x) strsplit(x, ",")[[1]][1]),
           cov = purrr::map_chr(indices, function(x) strsplit(x, ",")[[1]][2]),
           reg = factor(reg, levels = unique(reg), labels = dimnames(data$beh_cat)[[1]]),
           cov = factor(cov, levels = unique(cov), labels = dimnames(data$beh_cat)[[3]])) %>% 
    select(.draw, value, reg, cov) %>% 
    group_by(.draw, reg) %>% 
    mutate(share = value/sum(value)) %>% 
    ungroup() %>%
    group_by(reg, cov) %>%
    summarise(est_med = median(value),
              est_ci_lwr = quantile(value, 0.05),
              est_ci_upr = quantile(value, 0.95),
              share_med = median(share),
              share_ci_lwr = quantile(share, 0.05),
              share_ci_upr = quantile(share, 0.95)) %>%
    left_join(tibble(reg=names(rowSums(data$pop_r_age)), pop=rowSums(data$pop_r_age))) %>% 
    mutate(est_ci = paste0(format(est_med, digits = 2), " (", 
                           format(est_ci_lwr, digits = 2), "-",
                           format(est_ci_upr, digits = 2), ")"),
           share_ci = paste0(format(share_med, digits = 2), " (", 
                             format(share_ci_lwr, digits = 2), "-",
                             format(share_ci_upr, digits = 2), ")")) %>% 
    select(reg, cov, est_ci, share_ci) %>% 
    pivot_longer(cols = c("est_ci", "share_ci")) %>% 
    pivot_wider (id_cols=c("reg", "name"), names_from="cov", values_from="value") %>% 
    left_join(region_names, by = c("reg"="county_code"))
}

share_bl_ger = function(mcmc, data) {
  phi_star_draws = mcmc$draws("phi_star", format="draws_df")
  
  draws_long = phi_star_draws %>% select(-c(.chain, .iteration)) %>% 
    pivot_longer(cols = starts_with("phi"))
  
  share_per_draw = draws_long %>% 
    mutate(indices = gsub("\\]", "", gsub("phi_star\\[", "", name)),
           reg = purrr::map_chr(indices, function(x) strsplit(x, ",")[[1]][1]),
           cov = purrr::map_chr(indices, function(x) strsplit(x, ",")[[1]][2]),
           reg = factor(reg, levels = unique(reg), labels = dimnames(data$beh_cat)[[1]]),
           cov = factor(cov, levels = unique(cov), labels = dimnames(data$beh_cat)[[3]])) %>% 
    select(.draw, value, reg, cov) %>% 
    group_by(.draw, reg) %>% 
    mutate(share = value/sum(value)) %>% 
    ungroup()
  
  share_per_draw
}