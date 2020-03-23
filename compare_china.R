require(magrittr)
library(patchwork)

linelist <- NCoVUtils::get_international_linelist("China") %>%
 dplyr::mutate(import_status = "local")

final_cases <- linelist %>%
 dplyr::filter(!is.na(date_onset)) %>%
 dplyr::select(dates = date_onset) %>%
 dplyr::group_by(dates) %>%
 dplyr::summarise(local = dplyr::n()) %>%
 dplyr::arrange(dates) %>%
 tidyr::complete(dates = seq.Date(from = min(dates), 
                                                to = max(dates), by = "day"),fill = list(local = 0)) %>%
 dplyr::mutate(imported = 0)


final_cases[1,2:3] <- c(0,1)
final_cases

ee_res <- EpiEstim::estimate_R(incid = final_cases, method = "parametric_si",
                               config = EpiEstim::make_config(list(mean_si = 4.7, 
                                                                   std_si = 2.9)))
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

res_ch <- rbind(rbind(ee_df,st_df),stp_df)
final_cases_ch <- final_cases
si_ch <- data.frame(EpiEstim = ee_res$si_distr[1:59], 
                    stan_pois = c(0,res_p$w2[1,]), 
                    stan_nb = c(0,res$w1[1,]),
                    t = 0:58)

saveRDS(final_cases_ch, "final_cases_ch.rds")
saveRDS(res_ch, "res_ch.rds")
saveRDS(si_ch, "si_ch.rds")