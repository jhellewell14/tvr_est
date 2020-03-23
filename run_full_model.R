library(magrittr)
# Compile stan model
mod_full <- rstan::stan_model(file = "r_fit_all_full.stan")

# Read in serial interval MCMC trace from Nishiura paper
# si_samps <- read.csv("nishiura-lognormal-truncated.csv") 

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

# Re-shape into onset to confirmation delay samples
df <- linelist %>% dplyr::filter(!is.na(date_onset)) %>%
 dplyr::mutate(days = as.numeric(date_confirm - date_onset)) %>%
 dplyr::select(days)

df$upper <- ifelse(df$days==0,1,df$days+1)
df$lower <- ifelse(df$days==0,0,df$days-1)
df$lower <- ifelse(df$lower == 0, 1e-06, df$lower)

# Read in imported_cases / local_cases for South Korea
formatted_linelist <- linelist %>%
 dplyr::rename(date_onset_symptoms = date_onset,
               date_confirmation = date_confirm,
               delay_confirmation = report_delay)


# Get WHO sit rep case counts ---------------------------------------------

total_cases <- TimeVaryingNCovR0::get_who_cases(country = "RepublicofKorea", daily = TRUE)


# Join imported and local cases -------------------------------------------

cases <- TimeVaryingNCovR0::get_local_import_case_counts(total_cases, linelist)

reported_cases <- cases %>%
 dplyr::rename(confirm = cases)

## Split cases into local and imported
local_cases <- reported_cases %>%
 dplyr::filter(import_status == "local")

imported_cases <- reported_cases %>%
 dplyr::filter(import_status == "imported")

dat_full <- list(obs_local = local_cases$confirm,
            obs_imported = imported_cases$confirm,
            t = length(local_cases$confirm),
            tau = 7,
            si_loc = 1.3780732,
            si_scale = 0.6184616,
            N=nrow(df),
            low=df$lower,
            up=df$upper,
            lam_mean=mean(df$days))

# mod_full <- rstan::stan_model(file = "r_fit_all_full.stan")
fit_full <- rstan::sampling(mod_full, data = dat_full,chains = 1,iter = 1000, cores = 4,
                          control = list(adapt_delta = 0.8))
res_full <- rstan::extract(fit_full)

res_full <- data.frame(med = apply(res_full$R,MARGIN = 2, median),
           UQ = apply(res_full$R, MARGIN = 2, 
                      FUN = function(x){quantile(x, prob=0.975)}),
           LQ = apply(res_full$R, MARGIN = 2, 
                      FUN = function(x){quantile(x, prob=0.025)}),
           date = local_cases$date,
           model = "stan_full")


p2 <- res_full %>%
 ggplot2::ggplot(ggplot2::aes(x=date, y = med, ymin = LQ, ymax = UQ)) + 
 ggplot2::geom_line() + 
 ggplot2::geom_ribbon(alpha = 0.2) + 
 cowplot::theme_cowplot() + 
 ggplot2::geom_hline(yintercept = 1, lty = 2) +
 ggplot2::ylab("Reproduction number") + 
 ggplot2::xlab("") +
 ggplot2::scale_x_date(breaks=seq(as.Date("2020/01/20"), as.Date("2020/03/16"),by = "week"),
                date_labels = "%b %d")


p1 <- data.frame(date = local_cases$date, confirm = local_cases$confirm + imported_cases$confirm) %>%
   ggplot2::ggplot(ggplot2::aes(x = date, y = confirm)) + 
   ggplot2::geom_bar(stat="identity") +
   ggplot2::scale_x_date(breaks=seq(as.Date("2020/01/20"), as.Date("2020/03/16"),by = "week"),
                date_labels = "%b %d") +
   cowplot::theme_cowplot()

library(patchwork)
p1 / p2


# plot(apply(res_sk$out_w_vec,2,median)[1:20])
# points(apply(res_sk$out_wv,2,median)[1:20],col="red")
# points(apply(res_sk$neg_wv,2, median), col = "blue")
# plot(apply(res_sk$out_i,2,median))
# points(apply(res_sk$out_i1,2,median), col = "red")
# points(apply(res_sk$out_i2,2,median), col = "blue")
# lines(dat_sk$obs_local+dat_sk$obs_imported,col = "red")
