
library(foreach)
library(tidyverse)
library(MCMCvis)
library(rstan)
library(truncnorm)

empe_sites <- read.csv("Codes_EP/Codes_SAT/empe_sites_new.csv")
sites_update <- read.csv("Codes_EP/Codes_SAT/colony_attributes_update.csv")

data_pop <- readRDS("Codes_EP/Codes_SAT/data_pop_empe.rds")

theme_set(theme_bw())


# Load environmental data -------------------------------------------------

sites_all <- unique(readRDS("Codes_EP/Codes_SAT/data_env_empe.rds")$site_id)

data_esm <- readRDS("Codes_EP/Codes_SAT/data_env_empe.rds") %>%
  filter(site_id %in% unique(data_pop$sat$site_id)) %>%
  filter(year >= 2009 & year <= 2018) %>%
  select(site_id, year, contains(c("aice"))) %>%
  arrange(site_id)

data_env_avg <- data_esm %>%
  group_by(site_id) %>%
  summarise(across(-year, mean)) %>%
  ungroup()

env_mat <- data_env_avg[,-1] %>%
  as.matrix()

env_mat_std <- apply(env_mat, 2, function(x) (x - mean(x))/sd(x))


# Fit linear regression ---------------------------------------------------

#N_chains <- readRDS("Codes_EP/Codes_SAT/N_chains.rds")

## this is used only for R2 calculation
#N_mean_chains <- foreach(i = 1:50, .combine = "cbind") %do% {
  #z <- t(apply(N_chains[,seq(i, 500, 50)], 1, function(x) {
    #as.vector(mean(log(x[which(x != 0)])))
  #}))
  
  #as.vector(z)
#}

#y_mean <- apply(N_mean_chains, 2, mean)

#y <- apply(N_chains, 2, function(x) mean(log(x[x!=0])))
#y_sd <- apply(N_chains, 2, function(x) sd(log(x[x!=0])))
#N_data <- list(y = y, y_sd = y_sd, y_mean = y_mean)

N_data <- readRDS("Codes_EP/Codes_SAT/N_data.rds")

N_annual_log <- data.frame(y = N_data$y,
                           y_sd = N_data$y_sd,
                           site_id = rep(sites_update$site_id, 10))

N_annual_log <- N_annual_log[order(N_annual_log$site_id),]
y_mean <- N_data$y_mean[order(sites_update$site_id)]

x1 <- log(env_mat[,1])
x2 <- env_mat_std[,1]

idx <- which(!is.nan(N_annual_log$y))

dat_lm <- list(y = N_annual_log$y[idx], 
               y_sd = N_annual_log$y_sd[idx],
               y_mean = y_mean,
               N = 50,
               K = nrow(N_annual_log[idx,]),
               X = x1,
               site_no = as.numeric(as.factor(N_annual_log$site_id[idx])))

res_lm <- stan(file = "Codes_EP/Codes_SAT/lm5.stan", 
                 data = dat_lm,
                 iter = 5000,
                 thin = 2,
                 cores = 4,
                 control = list(adapt_delta = 0.99, 
                                max_treedepth = 15))

#if (!"results"  %in% list.files()) dir.create("results")
#saveRDS(res_lm, "Codes_EP/Codes_SAT/results_sat.rds")

# Model Checking
MCMCsummary(res_lm, params = c("mean_gt", "sd_gt"))
MCMCtrace(res_lm, params = "alpha", priors = rnorm(5000, 0, 20), pdf = F)
MCMCtrace(res_lm, params = "beta", priors = rnorm(5000, 0, 10), pdf = F)
MCMCtrace(res_lm, params = "mu_sigma", 
          priors = rtruncnorm(5000, 0, 5), pdf = F)
MCMCtrace(res_lm, params = "tau", 
          priors = rtruncnorm(5000, 0, 5), pdf = F)
MCMCtrace(res_lm, params = "sigma_site", 
          priors = rtruncnorm(5000, 0, 10), pdf = F)

# Create posterior csv files
param_chains <- MCMCchains(res_lm, params = c("alpha", "beta", "mu_sigma", 
                                              "sigma_time"))
write.csv(param_chains, "Codes_EP/Codes_SAT/param_chains.csv")

## Site effects
eps <- MCMCchains(res_lm, params = c("eps"), exact = F)
colnames(eps) <- unique(data_esm$site_id)

e_s <- matrix(0, ncol = 16, nrow = nrow(eps))
idx_site <- which(!sites_all %in% unique(data_esm$site_id))
colnames(e_s) <- sites_all[idx_site]

## Combine site effects
eps2 <- cbind(eps, e_s)
idx_site2 <- order(colnames(eps2)) 
eps2 <- eps2[,idx_site2]
write.csv(eps2, "Codes_EP/Codes_SAT/eps_chains.csv")

# Univariate plots
y1 <- MCMCsummary(res_lm, params = "mu")
y1_mean <- y1[,1]
y1_min <- y1[,3]
y1_max <- y1[,5]

ggplot() +
  geom_ribbon(aes(x = x1, ymin = y1_min, ymax = y1_max), alpha = 0.8, 
              fill = "grey") +
  #geom_segment(aes(x = env_mat[,1], xend = env_mat[,1], 
                   #y = N_annual_log$min, yend = N_annual_log$max)) +
  geom_line(aes(x = x1, y = y1_mean), col = "darkorange", size = 1.5, 
            linetype = 2) +
  #geom_point(aes(x = env_mat[,1], y = N_annual_log$mean), alpha= 0.8, 
             #color = "darkorange", size = 3) +
  #geom_point(aes(x = env_mat[,1], y = N_annual_log$mean), shape = 1, size = 3, 
             #stroke = 1.1) +
  labs(y = "Colony Abundance (log)", x = "SIC (Laying)") +
  theme(#axis.title.y = element_blank(),
    panel.border = element_blank(), 
    #panel.grid.major = element_blank()) +
    panel.grid.minor = element_blank())
