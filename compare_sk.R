require(magrittr)
library(patchwork)
# Read in South Korea linelist data
gsheets_url <- paste0("https://docs.google.com/spreadsheets/d/",
                      "1nKRkOwnGV7RgsMnsYE6l96u4xxl3ZaNiTluPKEPaWm8")
url <- paste0(gsheets_url, "/export?format=csv&gid=0")

linelist <- suppressWarnings(
 suppressMessages(
  readr::read_csv(url) %>%
   dplyr::mutate(import_status =
                  dplyr::if_else(!is.na(source) & source == "import", "imported", "local"),
                 report_delay =
                  as.integer(as.Date(date_confirm) - as.Date(date_onset))) %>%
   ## remove annoying column
   dplyr::select(-`KCDC_no (https://www.cdc.go.kr/board/board.es?mid=a20501000000&bid=0015)`) %>%
   dplyr::mutate(case = 1:dplyr::n())
 )) %>%
 tidyr::drop_na(date_confirm)

formatted_linelist <- linelist %>%
 dplyr::rename(date_onset_symptoms = date_onset,
               date_confirmation = date_confirm,
               delay_confirmation = report_delay)

imported_cases <- formatted_linelist %>%
 dplyr::select(date_onset_symptoms, import_status) %>%
 dplyr::filter(!is.na(date_onset_symptoms), import_status == "imported") %>%
 dplyr::group_by(date_onset_symptoms) %>% 
 dplyr::summarise(confirm = dplyr::n()) %>%
 tidyr::complete(date_onset_symptoms = seq.Date(from = min(date_onset_symptoms), 
                          to = max(date_onset_symptoms), by = "day"),fill = list(confirm = 0))

local_cases <- formatted_linelist %>%
 dplyr::select(date_onset_symptoms, import_status) %>%
 dplyr::filter(!is.na(date_onset_symptoms), import_status == "local") %>%
 dplyr::group_by(date_onset_symptoms) %>% 
 dplyr::summarise(confirm = n()) %>%
 tidyr::complete(date_onset_symptoms = seq.Date(from = min(date_onset_symptoms), 
                                                to = max(date_onset_symptoms), by = "day"),fill = list(confirm = 0))

final_cases <- dplyr::full_join(imported_cases,local_cases,by = "date_onset_symptoms") %>% 
 dplyr::select(dates = date_onset_symptoms, imported = confirm.x, local = confirm.y) %>%
 dplyr::mutate(imported = tidyr::replace_na(imported,0),
               local = tidyr::replace_na(local,0))


ee_res <- EpiEstim::estimate_R(incid = final_cases, method = "parametric_si",
                     config = EpiEstim::make_config(list(mean_si = 4.7, 
                                                         std_si = 2.9,
                                                         t_start = 2:54,
                                                         t_end = 9:61)))

ee_df <- data.frame(med = ee_res$R$`Median(R)`,
                    date = ee_res$dates[ee_res$R$t_end],
                    LQ = ee_res$R$`Quantile.0.025(R)`,
                    UQ = ee_res$R$`Quantile.0.975(R)`,
                    model = "EpiEstim")

dat <- list(obs_local = final_cases$local,
            obs_imported = final_cases$imported,
            t = length(final_cases$local),
            tau = 7,
            si_mean = 4.7,
            si_sd = 2.9,
            w = EpiEstim::discr_si(0:59, mu = 4.7, sigma = 2.9))

options(mc.cores = parallel::detectCores())
mod <- rstan::stan_model(file = "r_fit_imports.stan")
fit <- rstan::sampling(mod, data = dat,chains = 4,iter = 1000, cores = 4)
mod_p <- rstan::stan_model(file = "r_fit_poisson.stan")
fit_p <- rstan::sampling(mod_p, data = dat,chains = 4,iter = 1000, cores = 4)
res_p <- rstan::extract(fit_p)
res <- rstan::extract(fit)

# log_lik <- loo::extract_log_lik(fit, merge_chains = FALSE)
# r_eff <- loo::relative_eff(exp(log_lik), cores = 4)
# loo::loo(log_lik, r_eff = r_eff, cores = 4)
# loo::waic(log_lik, r_eff = r_eff, cores = 6)

st_df <- data.frame(med = apply(res$R, 2, median),
                    LQ = apply(res$R, 2, FUN = function(x){quantile(x, prob=0.025)}),
                    UQ = apply(res$R, 2, FUN = function(x){quantile(x, prob=0.975)}),
                    date = final_cases$dates,
                    model = "stan_nb")

stp_df <- data.frame(med = apply(res_p$R, 2, median),
                    LQ = apply(res_p$R, 2, FUN = function(x){quantile(x, prob=0.025)}),
                    UQ = apply(res_p$R, 2, FUN = function(x){quantile(x, prob=0.975)}),
                    date = final_cases$dates,
                    model = "stan_pois")

res_sk <- rbind(rbind(rbind(ee_df,st_df),stp_df),res_full)
final_cases_sk <- final_cases
si_sk <- data.frame(EpiEstim = ee_res$si_distr[1:61], 
                    stan_pois = c(0,res_p$w2[1,]), 
                    stan_nb = c(0,res$w1[1,]),
                    t = 0:60)

saveRDS(final_cases_sk, "final_cases_sk.rds")
saveRDS(res_sk, "res_sk.rds")
saveRDS(si_sk, "si_sk.rds")

