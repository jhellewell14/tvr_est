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

p1 <- rbind(rbind(ee_df,st_df),stp_df) %>%
   dplyr::filter(date >= as.Date("2020-01-18")) %>%
 ggplot2::ggplot(ggplot2::aes(x=date,y=med)) +
 ggplot2::geom_line(ggplot2::aes(col = model)) +
 ggplot2::geom_ribbon(ggplot2::aes(fill = model, ymin = LQ, ymax = UQ),alpha=0.2) + 
 cowplot::theme_cowplot() +
 ggplot2::geom_hline(yintercept = 1, lty = 2) +
 ggplot2::ylab("Effective reproduction number") +
 ggplot2::xlab("Date") +
 ggplot2::scale_x_date(date_breaks = "1 week", date_labels = "%b %d")

p2 <- ggplot2::ggplot(data = final_cases) +
 ggplot2::geom_bar(ggplot2::aes(x = dates, y = imported+local), stat = "identity") + 
 cowplot::theme_cowplot() +
 ggplot2::ylab("Cases") + 
 ggplot2::xlab("") +
 ggplot2::scale_x_date(date_breaks = "1 week", date_labels = "%b %d")


p2 / p1

p3 <- data.frame(EpiEstim = ee_res$si_distr[1:61], 
                 stan_pois = c(0,res_p$w2[1,]), 
                 stan_nb = c(0,res$w1[1,]),
           t = 0:60) %>%
   tidyr::gather(key = "model", value = "freq", -t) %>%
   dplyr::filter(t < 20) %>%
   ggplot2::ggplot(ggplot2::aes(x = t, y = freq, col = model)) +
   ggplot2::geom_line() +
   cowplot::theme_cowplot() +
   ggplot2::xlab("Time") +
   ggplot2::ylab("Frequency") +
   ggplot2::scale_color_discrete(guide = "none")
 
p2 / p1 / p3 + patchwork::plot_layout(guides = "collect")

plot(res_p$w1[1,])
points(res_p$w2[1,],col = "green")
points(ee_res$si_distr,col="red")
points(EpiEstim::discr_si(1:60, mu = 4.7, sigma = 2.9),col = "blue")

