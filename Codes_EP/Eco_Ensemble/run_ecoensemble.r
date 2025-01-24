library(tidyverse)

## Calculate the percentage chance between 2009 and 2018 on the global pop size

## 2009 to 2018
glob_pop <- read.csv("N_glb.csv")
glob_pop$year <- c(2009:2018)

mobs <- lm(N_mean ~ year,
           data=glob_pop)

summary(mobs)
beta_obs <- mobs$coefficients[["year"]]

## calculate yearly % change
((glob_pop$N_mean[2]-glob_pop$N_mean[1])/glob_pop$N_mean[1])*100
yearly_perc_change_obs <- ((lead(glob_pop$N_mean)-glob_pop$N_mean)/glob_pop$N_mean)*100

################################################################################
## 1. Read in time series and turn into percent change


################################################################################
## import IPM preds
IPM_preds_dir <- "Ntotal/IPM_Fra/"
IPM_preds_ens <- dir(IPM_preds_dir, patter=".txt")
import_name <- c(paste(IPM_preds_dir, IPM_preds_ens, sep=""))

## make sure the order is correct
ens_numb <- sub('.*ens', '', import_name)
ens_numb <- sub("\\..*", "", ens_numb)
ens_numb <- as.numeric(ens_numb)
import_name <- import_name[order(ens_numb)]

length(1921:2100)

change_posterior <- c() ## based on lm slope
change_posterior_rmse <- c() ## based on yearly % change

IPM_preds_list <- list()

for(i in 1:length(import_name))
{
  tempi <- read.table(import_name[i], sep=",")

  tempi <- tempi %>%
    mutate_if(is.factor, as.character)

  tempi <- tempi %>%
    mutate_if(is.character, as.numeric)

  ## each column is a sim
  ## each row is a year

  tempi <- as.matrix(tempi)

  seq_yrs <- c(1921:2100)

  nsim=dim(tempi)[2]

  change_full <- c()
  change_full_v2 <- c()

  #tempi <- apply(tempi, 2, function(x) ((lead(x) - x)/(x) * 100) )
  for(t in 1:nsim)
  {
    ind_2009_18 <- which(seq_yrs %in% c(2009:2018))

    ## the mean abundance during the period 2009-2018
    vals_2009_18 <- mean(as.numeric(tempi[ind_2009_18,t]))

    ## the tseries 2009-2018
    vals_2009_18_tseries <- as.numeric(tempi[ind_2009_18,t])

    ## the yearly % change
    vals_2009_18_tseries_v2 <- ((lead(vals_2009_18_tseries)-vals_2009_18_tseries)/vals_2009_18_tseries)*100

    ## residuals in % change
    diffpreds <- vals_2009_18_tseries_v2 - yearly_perc_change_obs
    diffpreds <- diffpreds[-length(diffpreds)]

    if(length(which(is.na(vals_2009_18_tseries))) > 0)
    {
      beta_pred_t <- NA
    } else {
      mpreds <- lm(vals_2009_18_tseries ~ c(1:10))
      beta_pred_t <- as.numeric(mpreds$coefficients[2])
    }

    ## the slope value
    beta_pred_t

    ## the tseries = (Nt - mean2009-18)/mean2009-18
    val_t <- (tempi[,t]-vals_2009_18)/vals_2009_18 * 100
    tempi[,t] <- val_t

    ## 1 criterion. weigth by slope
    change_full <- c(change_full, beta_pred_t)

    ## 2 criterion. residuals in % change
    change_full_v2 <- cbind(change_full_v2, diffpreds)
  }

  ## 1 criterion
  change_posterior <- c(change_posterior, change_full)

  ## 2 criterion - need to calculate RMSE
  change_full_v2 <- as.matrix(change_full_v2)
  change_full_rmse_v2 <- as.numeric(apply(change_full_v2, 2, function(x) sqrt(sum(x^2)/ncol(change_full_v2))))
  ## 2 criterion - weight by RMSE
  change_posterior_rmse <- c(change_posterior_rmse, change_full_rmse_v2)

  rownames(tempi) <- c(1921:2100)
  IPM_preds_list[[i]] <- tempi
}

hist(change_posterior)
mean(change_posterior, na.rm=T) ## -5107.733
IPM_discrep <- mean(abs(change_posterior - beta_obs), na.rm=T) ##2475.213

## rmse % change
hist(change_posterior_rmse)
mean(change_posterior_rmse, na.rm=T) ## 1.687513
IPM_discrep_rmse <- mean(change_posterior_rmse, na.rm=T) ## 1.687513


##############################################################################################################
## import SAT preds
SAT_preds_dir <- "Ntotal/LE_CESM2/"
SAT_preds_ens <- dir(SAT_preds_dir, patter=".txt")
import_name <- c(paste(SAT_preds_dir, SAT_preds_ens, sep=""))

## make sure the order is correct
ens_numb <- sub('.*ens', '', import_name)
ens_numb <- sub("\\..*", "", ens_numb)
ens_numb <- as.numeric(ens_numb)
import_name <- import_name[order(ens_numb)]

length(1909:2100)

change_posterior <- c()
change_posterior_rmse <- c() ## based on yearly % change

SAT_preds_list <- list()

for(i in 1:length(import_name))
{
  tempi <- read.table(import_name[i], sep=",")

  tempi <- tempi %>%
    mutate_if(is.factor, as.character)

  ## each column is a sim
  ## each row is a year

  tempi <- as.matrix(tempi)

  #tempi <- apply(tempi, 2, function(x) ((lead(x) - x)/(x) * 100) )

  rownames(tempi) <- c(1909:2100)

  seq_yrs <- c(1909:2100)

  nsim=dim(tempi)[2]

  change_full <- c()
  change_full_v2 <- c()

  #tempi <- apply(tempi, 2, function(x) ((lead(x) - x)/(x) * 100) )
  for(t in 1:nsim)
  {
    ind_2009_18 <- which(seq_yrs %in% c(2009:2018))
    vals_2009_18 <- mean(as.numeric(tempi[ind_2009_18,t]))
    vals_2009_18_tseries <- as.numeric(tempi[ind_2009_18,t])

    vals_2009_18_tseries_v2 <- ((lead(vals_2009_18_tseries)-vals_2009_18_tseries)/vals_2009_18_tseries)*100
    #((vals_2009_18_tseries[4]-vals_2009_18_tseries[3])/vals_2009_18_tseries[3])*100
    ## residuals
    diffpreds <- vals_2009_18_tseries_v2 - yearly_perc_change_obs
    diffpreds <- diffpreds[-length(diffpreds)]

    if(length(which(is.na(vals_2009_18_tseries))) > 0)
    {
      beta_pred_t <- NA
    } else {
      mpreds <- lm(vals_2009_18_tseries ~ c(1:10))
      beta_pred_t <- as.numeric(mpreds$coefficients[2])
    }

    val_t <- (tempi[,t]-vals_2009_18)/vals_2009_18 * 100
    tempi[,t] <- val_t

    change_full <- c(change_full, beta_pred_t)
    change_full_v2 <- cbind(change_full_v2, diffpreds)

  }

  change_posterior <- c(change_posterior, change_full)

  change_full_v2 <- as.matrix(change_full_v2)
  change_full_rmse_v2 <- as.numeric(apply(change_full_v2, 2, function(x) sqrt(sum(x^2)/ncol(change_full_v2))))

  change_posterior_rmse <- c(change_posterior_rmse, change_full_rmse_v2)

  tempi <- tempi[which(rownames(tempi) %in% c(1921:2100)),]

  SAT_preds_list[[i]] <- tempi
}

length(SAT_preds_list)
hist(change_posterior)
mean(change_posterior)
which(is.na(change_posterior))
SAT_discrep <- mean(abs(change_posterior - beta_obs), na.rm=T) ## 2113.237

## rmse % change
hist(change_posterior_rmse)
mean(change_posterior_rmse, na.rm=T) ## 1.493522
SAT_discrep_rmse <- mean(change_posterior_rmse, na.rm=T) ## 1.493522


################################################################################
## import CMR preds
CMR_preds_dir <- "Ntotal/CMR_Stef/"
CMR_preds_ens <- dir(CMR_preds_dir, patter=".txt")

CMR_preds_ens_scen1 <- CMR_preds_ens[str_detect(CMR_preds_ens, c('_scen1'))]
CMR_preds_ens_scen2 <- CMR_preds_ens[str_detect(CMR_preds_ens, c('_scen2'))]
CMR_preds_ens_scen4 <- CMR_preds_ens[str_detect(CMR_preds_ens, c('_scen4'))]

import_name1 <- c(paste(CMR_preds_dir, CMR_preds_ens_scen1, sep=""))
## make sure the order is correct
ens_numb <- gsub(".*[_ens]([^.]+)[_].*", "\\1", import_name1)
ens_numb <- as.numeric(ens_numb)
import_name1 <- import_name1[order(ens_numb)]

import_name2 <- c(paste(CMR_preds_dir, CMR_preds_ens_scen2, sep=""))
ens_numb <- gsub(".*[_ens]([^.]+)[_].*", "\\1", import_name2)
ens_numb <- as.numeric(ens_numb)
import_name2 <- import_name2[order(ens_numb)]

import_name4 <- c(paste(CMR_preds_dir, CMR_preds_ens_scen4, sep=""))
ens_numb <- gsub(".*[_ens]([^.]+)[_].*", "\\1", import_name4)
ens_numb <- as.numeric(ens_numb)
import_name4 <- import_name4[order(ens_numb)]

cmr1_change_posterior <- c()
cmr2_change_posterior <- c()
cmr4_change_posterior <- c()

cmr1_change_posterior_rmse <- c()
cmr2_change_posterior_rmse <- c()
cmr4_change_posterior_rmse <- c()

CMR_scen1_preds_list <- list()
CMR_scen2_preds_list <- list()
CMR_scen4_preds_list <- list()

length(1900:2100)

for(i in 1:length(import_name1))
{
  ## each column is a sim
  ## each row is a year

  ## scenario1
  tempi <- read.table(import_name1[i], sep=",")

  tempi <- tempi %>%
    mutate_if(is.factor, as.character)

  tempi <- as.matrix(tempi)
  #tempi <- apply(tempi, 2, function(x) ((lead(x) - x)/(x) * 100) )
  rownames(tempi) <- c(1900:2100)

  seq_yrs <- c(1900:2100)

  nsim=dim(tempi)[2]
  #tempi <- apply(tempi, 2, function(x) ((lead(x) - x)/(x) * 100) )

  change_full <- c()
  change_full_v2 <- c()

  for(t in 1:nsim)
  {
    ind_2009_18 <- which(seq_yrs %in% c(2009:2018))
    vals_2009_18 <- mean(as.numeric(tempi[ind_2009_18,t]))
    vals_2009_18_tseries <- as.numeric(tempi[ind_2009_18,t])

    vals_2009_18_tseries_v2 <- ((lead(vals_2009_18_tseries)-vals_2009_18_tseries)/vals_2009_18_tseries)*100
    diffpreds <- vals_2009_18_tseries_v2 - yearly_perc_change_obs
    diffpreds <- diffpreds[-length(diffpreds)]

    if(length(which(is.na(vals_2009_18_tseries))) > 0)
    {
      beta_pred_t <- NA
    } else {
      mpreds <- lm(vals_2009_18_tseries ~ c(1:10))
      beta_pred_t <- as.numeric(mpreds$coefficients[2])
    }

    val_t <- (tempi[,t]-vals_2009_18)/vals_2009_18 * 100
    tempi[,t] <- val_t

    change_full <- c(change_full, beta_pred_t)
    change_full_v2 <- cbind(change_full_v2, diffpreds)
  }

  cmr1_change_posterior <- c(cmr1_change_posterior, change_full)

  change_full_v2 <- as.matrix(change_full_v2)
  change_full_rmse_v2 <- as.numeric(apply(change_full_v2, 2, function(x) sqrt(sum(x^2)/ncol(change_full_v2))))
  cmr1_change_posterior_rmse <- c(cmr1_change_posterior_rmse, change_full_rmse_v2)

  tempi <- tempi[which(rownames(tempi) %in% c(1921:2100)),]

  CMR_scen1_preds_list[[i]] <- tempi


  ## scenario2
  tempi <- read.table(import_name2[i], sep=",")
  tempi <- tempi %>%
    mutate_if(is.factor, as.character)

  tempi <- as.matrix(tempi)
  #tempi <- apply(tempi, 2, function(x) ((lead(x) - x)/(x) * 100) )
  rownames(tempi) <- c(1900:2100)

  seq_yrs <- c(1900:2100)

  nsim=dim(tempi)[2]
  #tempi <- apply(tempi, 2, function(x) ((lead(x) - x)/(x) * 100) )

  change_full <- c()
  change_full_v2 <- c()

  for(t in 1:nsim)
  {
    ind_2009_18 <- which(seq_yrs %in% c(2009:2018))
    vals_2009_18 <- mean(as.numeric(tempi[ind_2009_18,t]))
    vals_2009_18_tseries <- as.numeric(tempi[ind_2009_18,t])

    vals_2009_18_tseries_v2 <- ((lead(vals_2009_18_tseries)-vals_2009_18_tseries)/vals_2009_18_tseries)*100
    diffpreds <- vals_2009_18_tseries_v2 - yearly_perc_change_obs
    diffpreds <- diffpreds[-length(diffpreds)]

    if(length(which(is.na(vals_2009_18_tseries))) > 0)
    {
      beta_pred_t <- NA
    } else {
      mpreds <- lm(vals_2009_18_tseries ~ c(1:10))
      beta_pred_t <- as.numeric(mpreds$coefficients[2])
    }

    val_t <- (tempi[,t]-vals_2009_18)/vals_2009_18 * 100
    tempi[,t] <- val_t

    change_full <- c(change_full, beta_pred_t)
    change_full_v2 <- cbind(change_full_v2, diffpreds)
  }

  cmr2_change_posterior <- c(cmr2_change_posterior, change_full)

  change_full_v2 <- as.matrix(change_full_v2)
  change_full_rmse_v2 <- as.numeric(apply(change_full_v2, 2, function(x) sqrt(sum(x^2)/ncol(change_full_v2))))
  cmr2_change_posterior_rmse <- c(cmr2_change_posterior_rmse, change_full_rmse_v2)

  tempi <- tempi[which(rownames(tempi) %in% c(1921:2100)),]
  CMR_scen2_preds_list[[i]] <- tempi


  ## scenario4
  tempi <- read.table(import_name4[i], sep=",")
  tempi <- tempi %>%
    mutate_if(is.factor, as.character)

  tempi <- as.matrix(tempi)
  #tempi <- apply(tempi, 2, function(x) ((lead(x) - x)/(x) * 100) )
  rownames(tempi) <- c(1900:2100)

  seq_yrs <- c(1900:2100)

  nsim=dim(tempi)[2]
  #tempi <- apply(tempi, 2, function(x) ((lead(x) - x)/(x) * 100) )

  change_full <- c()
  change_full_v2 <- c()

  for(t in 1:nsim)
  {
    ind_2009_18 <- which(seq_yrs %in% c(2009:2018))
    vals_2009_18 <- mean(as.numeric(tempi[ind_2009_18,t]))
    vals_2009_18_tseries <- as.numeric(tempi[ind_2009_18,t])

    vals_2009_18_tseries_v2 <- ((lead(vals_2009_18_tseries)-vals_2009_18_tseries)/vals_2009_18_tseries)*100
    diffpreds <- vals_2009_18_tseries_v2 - yearly_perc_change_obs
    diffpreds <- diffpreds[-length(diffpreds)]

    if(length(which(is.na(vals_2009_18_tseries))) > 0)
    {
      beta_pred_t <- NA
    } else {
      mpreds <- lm(vals_2009_18_tseries ~ c(1:10))
      beta_pred_t <- as.numeric(mpreds$coefficients[2])
    }

    val_t <- (tempi[,t]-vals_2009_18)/vals_2009_18 * 100
    tempi[,t] <- val_t

    change_full <- c(change_full, beta_pred_t)
    change_full_v2 <- cbind(change_full_v2, diffpreds)
  }

  cmr4_change_posterior <- c(cmr4_change_posterior, change_full)

  change_full_v2 <- as.matrix(change_full_v2)
  change_full_rmse_v2 <- as.numeric(apply(change_full_v2, 2, function(x) sqrt(sum(x^2)/ncol(change_full_v2))))
  cmr4_change_posterior_rmse <- c(cmr4_change_posterior_rmse, change_full_rmse_v2)

  tempi <- tempi[which(rownames(tempi) %in% c(1921:2100)),]
  CMR_scen4_preds_list[[i]] <- tempi
}


length(CMR_scen1_preds_list)
length(CMR_scen2_preds_list)
length(CMR_scen4_preds_list)

hist(cmr1_change_posterior)
which(is.na(cmr1_change_posterior))

CMR_scen1_discrep <- mean(abs(cmr1_change_posterior - beta_obs), na.rm=T)
CMR_scen2_discrep <- mean(abs(cmr2_change_posterior - beta_obs), na.rm=T)
CMR_scen4_discrep <- mean(abs(cmr4_change_posterior - beta_obs), na.rm=T)

CMR_scen1_discrep_rmse <- mean(cmr1_change_posterior_rmse, na.rm=T) ## 1.504281
CMR_scen2_discrep_rmse <- mean(cmr2_change_posterior_rmse, na.rm=T) ## 1.626914
CMR_scen4_discrep_rmse <- mean(cmr4_change_posterior_rmse, na.rm=T) ## 1.480603

length(IPM_preds_list)
length(SAT_preds_list)
length(CMR_scen1_preds_list)


################################################################################
## calculate weights
################################################################################

## using RMSE
ipmw <- abs(1/IPM_discrep_rmse)
satw <- abs(1/SAT_discrep_rmse)
cmr1w <- abs(1/CMR_scen1_discrep_rmse)
cmr2w <- abs(1/CMR_scen2_discrep_rmse)
cmr4w <- abs(1/CMR_scen4_discrep_rmse)

wts <- c(ipmw, satw, cmr1w)
wts <- wts/sum(wts)
sum(wts)
round(wts, 3)
round(wts, 2)
round(wts, 1)
sum(round(wts, 2))
## 0.34 0.29 0.37


## using slope
ipmw <- abs(1/IPM_discrep)
satw <- abs(1/SAT_discrep)
cmr1w <- abs(1/CMR_scen1_discrep)
cmr2w <- abs(1/CMR_scen2_discrep)
cmr4w <- abs(1/CMR_scen4_discrep)

wts <- c(ipmw, satw, cmr1w)
wts <- wts/sum(wts)
sum(wts)
round(wts, 2)
sum(round(wts, 3))
## 0.29 0.32 0.40


################################################################################
## 2. Loop through ensamble

all_ens <- length(IPM_preds_list)

ens_preds_list <- list()

## Add a progress bar
bar <- txtProgressBar(min = 0, max = all_ens, initial = 0, style = 3)

str(IPM_preds_list[[1]])
str(SAT_preds_list[[1]])
str(CMR_scen1_preds_list[[1]])

for(e in 1:all_ens)
{
  IPM_e <- IPM_preds_list[[e]]
  SAT_e <- SAT_preds_list[[e]]

  ##  CMR output just considering scen1 here
  CMR_e <- CMR_scen1_preds_list[[e]]

  ##  Use different weights (IPM, SAT, CMR)
  wt <- wts

  nyrs <- nrow(IPM_e)
  nsims <- ncol(IPM_e)

  ens_avg <- matrix(nrow= nyrs,
                    ncol= nsims)

  ## loop through the years
  ## Each year: sample % change from posteriors of each models, based on respective weight
  ## here, we use same weight

  for(t in 1:nyrs)
  {
    IPM_e_yrt <- data.frame(ntot=as.vector(IPM_e[t,]),
                            prop=wt[1])
    SAT_e_yrt <- data.frame(ntot=as.vector(SAT_e[t,]),
                            prop=wt[2])
    CMR_e_yrt <- data.frame(ntot=as.vector(SAT_e[t,]),
                            prop=wt[3])

    all_e_yrt <- rbind.data.frame(IPM_e_yrt, SAT_e_yrt, CMR_e_yrt)

    ## sample without replacement
    out_avg_yrt <- sample(all_e_yrt[,1], size = 100, replace = F, prob=all_e_yrt[,2])

    ens_avg[t,] <- out_avg_yrt
  }

  ens_preds_list[[e]] <- ens_avg

  ## update the progress bar
  setTxtProgressBar(bar,e)
}

################################################################################
################################################################################

## extract ensamble mean and IQI (lci=0.25, uci=0.75 percentiles) from list

all_ens <- length(ens_preds_list)
nsims <- length(ens_preds_list[[1]][1,])
nyrs <- length(ens_preds_list[[1]][,1])

full_ens_mat <- matrix(nrow= nyrs,
                       ncol= nsims*all_ens)

for(t in 1:nyrs)
{
  full_vals <- c()

  for(e in 1:all_ens)
  {
    ens_e_yrt <- ens_preds_list[[e]][t,]
    full_vals <- c(full_vals, ens_e_yrt)
  }

  #length(full_vals)
  full_ens_mat[t,] <- full_vals
}

################################################################################

length(1921:2100)
str(full_ens_mat) ## rows are the years, cols are sims
i <- which(c(1921:2100) == 2073)
mean(full_ens_mat[i,])
length(which(as.numeric(full_ens_mat[i,]) <= c(-50)))/length(full_ens_mat[i,])
## 0.462 i.e. 2310/5000

## summarise
summ_df <- data.frame(Year = 1921:2100,
                      Mean_perc_change = as.numeric(apply(full_ens_mat, 1, mean,
                                                          na.rm=T)),
                      Median_perc_change = as.numeric(apply(full_ens_mat, 1, median,
                                                            na.rm=T)),
                      LCI_perc_change = as.numeric(apply(full_ens_mat, 1, quantile,
                                                         probs=0.025, na.rm=T)),
                      UCI_perc_change = as.numeric(apply(full_ens_mat, 1, quantile,
                                                         probs=0.975, na.rm=T))
)


ggplot() +
  theme_classic() +
  theme(panel.background = element_blank(),
        legend.position = c(1, 1),
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = "white"),
        legend.key = element_blank(),
        legend.title = element_text(colour="black", size=12),
        legend.text = element_text(colour="black", size = 12),
        plot.title = element_text(hjust = 0, vjust=0, face = "italic", size=12),
        text = element_text(family = "sans"),
        axis.title=element_text(size=12),
        axis.text.x=element_text(size=10)) +
  scale_x_continuous(breaks=seq(1921, 2100, by=10),
                     labels=seq(1921, 2100, by=10)) +
  geom_line(data=summ_df, aes(x=Year, y=Median_perc_change), col="black", linetype="solid") +
  geom_ribbon(data=summ_df, aes(x=Year, ymin = LCI_perc_change, ymax = UCI_perc_change),
              fill="red", alpha = .3) +
  ylab("Percent change")






