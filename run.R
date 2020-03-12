imported_cases <- readRDS("imported_cases.rds")
local_cases <- readRDS("local_cases.rds")

dat <- list(obs_local = local_cases$confirm,
            obs_imported = imported_cases$confirm,
            t = length(local_cases$confirm),
            tau = 7,
            si_mean = 4.7,
            si_sd = 2.9)


mod <- rstan::stan_model(file = "r_fit_imports.stan")
fit <- rstan::sampling(mod, data = dat,chains = 4,iter = 1000, cores = 4)

rstan::summary(fit)
res <- rstan::extract(fit)
p1 <- data.frame(med = apply(res$R,MARGIN = 2, median),
           UQ = apply(res$R, MARGIN = 2, 
                      FUN = function(x){quantile(x, prob=0.975)}),
           LQ = apply(res$R, MARGIN = 2, 
                      FUN = function(x){quantile(x, prob=0.025)}),
           time = 1:length(dat$obs_local),
           date = local_cases$date) %>%
 ggplot2::ggplot(aes(x=date, y = med, ymin = LQ, ymax = UQ)) + 
 ggplot2::geom_line() + 
 ggplot2::geom_ribbon(alpha = 0.2) + 
 cowplot::theme_cowplot() + 
 ggplot2::geom_hline(yintercept = 1, lty = 2) +
 ggplot2::scale_x_date(date_breaks = "7 days") +
 ggplot2::ylab("Reproduction number") + 
 ggplot2::xlab("")


## PRIORS for serial interval method

mod1 <- rstan::stan_model("~/Downloads/r_fit_imports_sifit.stan")
fit2 <- rstan::sampling(mod1, data = dat,chains = 4,iter = 1000, cores = 4)
rstan::summary(fit2)
res2 <- rstan::extract(fit2)
p2 <- data.frame(med = apply(res2$R,MARGIN = 2, median),
                 UQ = apply(res2$R, MARGIN = 2, 
                            FUN = function(x){quantile(x, prob=0.975)}),
                 LQ = apply(res2$R, MARGIN = 2, 
                            FUN = function(x){quantile(x, prob=0.025)}),
                 time = 1:length(dat$obs_local),
                 date = local_cases$date) %>%
        ggplot2::ggplot(aes(x=date, y = med, ymin = LQ, ymax = UQ)) + 
        ggplot2::geom_line() + 
        ggplot2::geom_ribbon(alpha = 0.2) + 
        cowplot::theme_cowplot() + 
        ggplot2::geom_hline(yintercept = 1, lty = 2) +
        ggplot2::scale_x_date(date_breaks = "7 days") +
        ggplot2::ylab("Reproduction number") + 
        ggplot2::xlab("Hyperparameters")


### EXTREMELY SLOW ###

# dat2 <- list(obs_local = vapply(1:1000,
#                                FUN = function(x){local_cases$confirm + sample(0:1,
#                                                                               size = length(local_cases$confirm),
#                                                                               replace = TRUE,prob = c(0.9,0.1))},
#                                FUN.VALUE = rep(1,length(local_cases$confirm))),
#             obs_imported = vapply(1:1000,
#                                   FUN = function(x){imported_cases$confirm + sample(0:1,
#                                                                                  size = length(local_cases$confirm),
#                                                                                  replace = TRUE,prob = c(0.9,0.1))},
#                                   FUN.VALUE = rep(1,length(local_cases$confirm))),
#             t = length(local_cases$confirm),
#             tau = 7,
#             si_mean = rnorm(1000,4.7,0.01),
#             si_sd = rnorm(1000,2.9,0.01),
#             n_samp = 1000)

# mod2 <- rstan::stan_model("r_fit_sifit.stan")
# fit2 <- rstan::sampling(mod2,dat = dat2, chains = 1, iter = 1000)
