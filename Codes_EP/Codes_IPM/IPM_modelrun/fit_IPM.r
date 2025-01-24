################################################################################
##                                EPeng IPM                                   ##
################################################################################

library(lattice)
library(coda)
library(R2jags)
library(tidyverse)
library(ggplot2)
library(MCMCvis)
library(mcmcplots)
library(truncnorm)
library(dplyr)
library(zoo)

## Function to create m-array from encounter histories
## From Abadi, F., Barbraud, C., Gimenez, O. (2017). Integrated population modeling
##    reveals the impact of climate on the survival of juvenile emperor penguins.
##    Global Change Biology, 23, 1353-1359.
marray <- function(CH){
  nind <- dim(CH)[1]
  n.occasions <- dim(CH)[2]
  m.array <- matrix(data=0, ncol=n.occasions+1, nrow= n.occasions)
  # Calculate the number of released individuals at each time period
  for (t in 1:n.occasions){
    m.array[t,1] <- sum(CH[,t])
  }
  for (i in 1:nind){
    pos <- which(CH[i,] != 0)
    g <- length(pos)
    for (z in 1:(g-1)){
      m.array[pos[z], pos[z+1]] <- m.array[pos[z], pos[z+1]] +1
    } # z
  } # i
  # Calculate the number of individuals that is never recaptured
  for (t in 1:n.occasions){
    m.array[t,n.occasions +1] <- m.array[t,1] -
      sum(m.array[t,2:n.occasions])
  }
  out <- m.array[1:(n.occasions-1),2:(n.occasions+1)]
  return(out)
}

################################################################################
## Read in data
################################################################################

## count and productivity data
CountData <- read.csv("count_fecundity_updated.csv",sep=",",header=T)
colnames(CountData)[which(colnames(CountData) %in% c("Neggs"))] <- "Number.of.breeding.pairs"
colnames(CountData)[which(colnames(CountData) %in% c("Nfledged"))] <- "Number.of.chicks.fledged"

CountData$Number.of.chicks.fledged[which(CountData$Year %in% c(1972))]
CountData$Number.of.breeding.pairs[which(CountData$Year %in% c(1979))]

sel_yr <- 1979
end_yr <- 2020

CountData <- dplyr::filter(CountData, Year %in% c(1979:end_yr))
ind_yr <- which(CountData$Year==sel_yr)

BRAD <- CountData$Number.of.breeding.pairs[which(CountData$Year %in% c(sel_yr:end_yr))]
CHICK <- CountData$Number.of.chicks.fledged[which(CountData$Year %in% c(sel_yr:end_yr))]
BRAD <- as.vector(BRAD) ## N breeding pairs
CHICK <- as.vector(CHICK) ## N chicks
Time_counts <- length(BRAD)[1]

# Mark-recapture data (1971-1998), summarised as m-array
m.HC <- readRDS("m_array_EPENG.rds")

## two years with complete data loss (the vector will be used as a switch in the model)
data_loss_cmr <- rep(1, length(1979:1998))
index_loss <- which(seq(1979,1998, by=1) %in% c(1990, 1992))
data_loss_cmr[index_loss] <- 0
length(1979:1998)
length(data_loss_cmr)

################################################################################
## read environmental variables
################################################################################

env_data <- readRDS(file = "clim_data_byphase_1979_2021.rds")
env_data <- dplyr::filter(env_data, site_id%in% c("PGEO"))

ice_data <- readRDS(file = "SIC_byphase_1979_2020.rds")
ice_data <- dplyr::filter(ice_data, site_id%in% c("PGEO"))

## FAST ICE
fice <- read.csv("fast_ice_1979_2017_cut_for_climwin.csv")
fice$X <- NULL
fice$X.1 <- NULL
fice <- fice %>%
  mutate_if(is.factor, as.character)
fice$date_climwin <- as.Date(fice$date_climwin,
                             format = "%m/%d/%Y")
fice$year <- lubridate::year(fice$date_climwin)

## UPDATED FAST ICE FROM SATELLITE MEASUREMENTS
fice_1721 <- read.csv("NOW_2017_2021.csv")
fice_1721$date <- as.Date(fice_1721$date,
                          format = "%d.%m.%Y")
fice_1721 <- dplyr::filter(fice_1721, lubridate::month(date) %in% c(8,9,10,11))
fice_1721 <- fice_1721 %>%
  dplyr::group_by(lubridate::year(date),
                  lubridate::month(date)) %>%
  dplyr::summarise(NOW_dist_km = mean(NOW_distance_km))
colnames(fice_1721)[1] <- "year"
colnames(fice_1721)[2] <- "month"
fice_1721 <- fice_1721 %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(NOW_dist_km = mean(NOW_dist_km))
fice_1821 <- dplyr::filter(fice_1721, year %in% c(2018:2021))
fice_1821 <- data.frame(year=fice_1821$year,
                        NOW_rear=fice_1821$NOW_dist_km,
                        FIA_rear=NA)

## append datasets
NOW_FIA_rearing <- fice %>%
  dplyr::filter(lubridate::month(date_climwin) %in% c(08:12))
NOW_FIA_rearing <- NOW_FIA_rearing %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(NOW_rear=mean(NOW, na.rm=T),
                   FIA_rear=mean(FIA, na.rm=T))
NOW_FIA_rearing <- rbind.data.frame(NOW_FIA_rearing, fice_1821)
NOW_FIA_rearing <- dplyr::filter(NOW_FIA_rearing, year <= 2020)

## NOW & FIA during rearing (aug-dec)
NOW_rear_std <- (NOW_FIA_rearing$NOW_rear - mean(NOW_FIA_rearing$NOW_rear))/mean(NOW_FIA_rearing$NOW_rear)

## fecundity
## NOW rear
NOW_rear_std

## breeding probability:
## sic during rearing
SIC_rear <- ice_data$aice_rearing[which(env_data$year %in% c(1979:(end_yr)))]
SIC_rear_std <- (SIC_rear - mean(SIC_rear, na.rm=T))/mean(SIC_rear, na.rm=T)
## vwind during rearing
vwind_rear <- env_data$vatm_rearing[which(env_data$year %in% c(1979:(end_yr)))]
vwind_rear_std <- (vwind_rear - mean(vwind_rear))/mean(vwind_rear)
## sst during nonb
sst_nonbreed <- env_data$temp_nonbreed[which(env_data$year %in% c(1980:(end_yr+1)))]
sst_nonbreed_std <- (sst_nonbreed - mean(sst_nonbreed))/mean(sst_nonbreed)

## survival
## now during rear
NOW_rear_std
## npp during breeding
npp_breed <- env_data$phtc_breeding[which(env_data$year %in% c(1979:(end_yr)))]
npp_breed_std <- (npp_breed - mean(npp_breed))/mean(npp_breed)
## sic during nonb
SIC_nonbreed <- ice_data$aice_nonbreed[which(env_data$year %in% c(1980:(end_yr+1)))]
SIC_nonbreed_std <- (SIC_nonbreed - mean(SIC_nonbreed, na.rm=T))/mean(SIC_nonbreed, na.rm=T)
## sst during nonb
sst_nonbreed_std
## vwind rear
vwind_rear_std
## wspeed rear
wind_rear <- env_data$wind_rearing[which(env_data$year %in% c(1979:(end_yr)))]
wind_rear_std <- (wind_rear - mean(wind_rear))/mean(wind_rear)


################################################################################
## IPM
################################################################################

# sink("model/bestmodel_final.txt")
# cat("
# model {
#     #****************************************
#     # Define priors
#     #****************************************
#
#     # Initial population size
#
#     n1 ~ dnorm(1500, 0.0005)T(0,)
#     N1[1] <- round(n1)
#
#     n2 ~ dnorm(770, 0.0005)T(0,)
#     N2[1] <- round(n2)
#
#     n3 ~ dnorm(1600, 0.0001)T(0,)
#     N3[1] <- round(n3)
#
#     n4 ~ dnorm(1700, 0.00005)T(0,)
#     N4[1] <- round(n4)
#
#     n5 ~ dnorm(1800, 0.00005)T(0,)
#     N5[1] <- round(n5)
#
#     nbreed ~ dnorm(3431, 0.005)T(0,)
#     Nbreed[1] <- round(nbreed)
#
#     nnonbreed ~ dnorm(1700,0.001)T(0,)
#     Nonbreed[1] <- round(nnonbreed)
#
#     ## survival juvenile
#     mu.phi1 ~ dunif(0.5, 0.8)
#
#     ## adult survival, with env effect
#     beta0.phia ~ dunif(2, 3.5) # prior for intercept
#     betasst.phia ~ dunif(-2, 2) # prior for slope
#     betawind.phia ~ dunif(-2, 2) # prior for slope
#
#     ## probability of returning to breed, with env effect
#     beta0.breed ~ dnorm(1.25, 6.25)T(0.75, 1.75) # prior for intercept
#     betavwind.breed ~ dunif(-2, 2) # prior for slope
#
#     ## fecundity, with env effect
#     beta0.fec ~ dunif(-0.5, 1) # prior for intercept
#     betanow.fec ~ dunif(-2, 2) # prior for slope
#
#     ## recruiting
#     recr5 <- 0.22
#     recr6 <- 0.32
#
#     ## resight adults
#     mu.p ~ dunif(0,1)
#
#     ## process error
#     sigma2.phia ~  dnorm(0, 1)T(0,0.2)
#     tau.phia <- pow(sigma2.phia, -1)
#     sigma2.fec ~  dnorm(0, 1)T(0,)
#     tau.fec <- pow(sigma2.fec, -1)
#     sigma2.breed ~  dnorm(0, 1)T(0, 0.2)##T(0,)
#     tau.breed <- pow(sigma2.breed, -1)
#
#     for (t in 1:((Time_counts-1)))
#     {
#     # adult survival
#     logit(phia[t]) <- lphia.lim[t]
#     lphia.lim[t] <- min(999, max(-999,lphia[t]))
#     lphia[t] <- beta0.phia +
#                 betasst.phia*sst_nonbreed[t] +
#                 betawind.phia*wind_rear[t] +
#                 eps.phia[t]
#     eps.phia[t] ~ dnorm(0, tau.phia)
#
#     # juvenile survival
#     phi1[t] <- mu.phi1
#
#     ## breeding probability
#     logit(pbreed_B[t]) <- lpbreed_B.lim[t]
#     lpbreed_B.lim[t] <- min(999, max(-999,lpbreed_B[t]))
#     lpbreed_B[t] <- beta0.breed +
#                     betavwind.breed*vwind_rear[t] +
#                     eps.breed[t]
#     eps.breed[t] ~ dnorm(0, tau.breed)
#     }
#
#     for (t in 1:(Time_counts))
#     {
#     # productivity
#     logit(fec[t]) <- lfec.lim[t]
#     lfec.lim[t] <- min(999, max(-999,lfec[t]))
#     lfec[t] <- beta0.fec +
#                betanow.fec*now_rear[t] +
#                eps.fec[t]
#     eps.fec[t] ~ dnorm(0, tau.fec)
#     }
#
#     for (t in 1:(n.occasions-1))
#     {
#     # resighting, accounting for data loss
#     p[t] <- data_loss[t+1]*mu.p ## check indexing
#     } #t
#
#     #****************************************
#     # Likelihood for each data set
#     #****************************************
#
#     #****************************************
#     # Likelihood for population population count data
#     #****************************************
#
#     #****************************************
#     # System process
#     #****************************************
#
#     for (t in 2:(Time_counts))
#     {
#     # age0, t1 (fledglings)
#     # 0.5*Nbreed[t]*fec[t]
#
#     #age1
#     ##N1[t] <- round(phi1[t-1]*(0.5*Nbreed[t-1]*fec[t-1]))
#     mean1[t] <- round(0.5*Nbreed[t-1]*fec[t-1]*phi1[t-1])
#     N1[t] ~ dpois(mean1[t])
#
#     #age2
#     ##N2[t] <- round(phia[t-1]*(N1[t-1]))
#     N2[t] ~ dbin(phia[t-1],N1[t-1])
#
#     #age3
#     ##N3[t] <- round(phia[t-1]*(N2[t-1]))
#     N3[t] ~ dbin(phia[t-1],N2[t-1])
#
#     #age4 (prebreeders)
#     ##N4[t] <- round(phia[t-1]*N3[t-1])
#     N4[t] ~ dbin(phia[t-1],N3[t-1])
#
#     #age5 (pre breeders)
#     ##N5[t] <- round(phia[t-1]*((1-recr5)*N4[t-1] + (1-recr6)*N5[t-1]))
#     mean5[t] <- round(((1-recr5)*N4[t-1] + (1-recr6)*N5[t-1]))
#     N5[t] ~ dbin(phia[t-1], mean5[t])
#
#     #breeding adults
#     meanbreed[t] <-  round(recr5*N4[t-1] + recr6*N5[t-1] + pbreed_B[t-1]*Nbreed[t-1] +
#                            pbreed_B[t-1]*Nonbreed[t-1])
#     Nbreed[t] ~ dbin(phia[t-1], meanbreed[t])
#
#     #nonbreeding adults
#     meannonbreed[t] <- round(Nbreed[t-1]*(1-pbreed_B[t-1]) + Nonbreed[t-1]*(1-pbreed_B[t-1]))
#     Nonbreed[t] ~ dbin(phia[t-1], meannonbreed[t])
#     }
#
#     for (t in 1:(Time_counts))
#     {
#     Ntot[t] <- Nbreed[t] + Nonbreed[t] + N5[t] + N4[t] + N3[t] + N2[t] + N1[t]
#     } #t
#
#     #******************************
#     # Observation process
#     #******************************
#
#     for (t in 1:Time_counts){
#     Ka[t] ~ dpois(Nbreed[t])
#     }#t
#
#     #****************************************
#     # Likelihood for capture-recapture data:  m-array
#     #****************************************
#
#     # Define the multinomial likelihood
#     for (t in 1:(n.occasions-1)){
#     marr[t,1:n.occasions] ~ dmulti(pr[t, ],r[t])
#     }
#
#     # Main diagonal
#     for (t in 1:(n.occasions-1)){
#     q[t] <- 1-p[t] # Probability of non-recapture
#     pr[t,t] <- phia[t]*p[t] ## SHOULD ALSO ACCOUNT FOR PBREED?!
#     # Above main diagonal
#     for (j in (t+1):(n.occasions-1)){
#     pr[t,j] <- prod(phia[t:j])*prod(q[t:(j-1)])*p[j]
#     } # j
#     # Below main diagonal
#     for (j in 1:(t-1)){
#     pr[t,j] <- 0
#     } # j
#     } # t
#     # Last column: probability of non-recapture
#     for (t in 1:(n.occasions-1)){
#     pr[t,n.occasions] <- 1-sum(pr[t,1:(n.occasions-1)])
#     } # t
#
#     #**************************************************
#     # Likelihood for reproductive success data
#     #**************************************************
#     for (t in 1:(Time_counts))
#     {
#      ## fecundity, binomial
#      Kp[t] ~ dbin(fec[t],Ka[t])
#     }#t
# }#end model
# ",fill = TRUE)
# sink()


# Bundle data
ipm.data <- list(n.occasions=dim(m.HC)[2],
                 Ka=BRAD, Kp=CHICK,
                 Time_counts=length(BRAD),
                 marr=m.HC, r=rowSums(m.HC),
                 data_loss = data_loss_cmr,
                 now_rear = NOW_rear_std,
                 vwind_rear = vwind_rear_std,
                 sst_nonbreed = sst_nonbreed_std,
                 wind_rear = wind_rear_std)


# Parameters monitored
params <- c(
  "mu.phi1",
  "mu.p",

  "beta0.fec",
  "betanow.fec",
  "fec",
  "sigma2.fec",
  "eps.fec",

  "beta0.breed",
  "betavwind.breed",

  "pbreed_B",
  "sigma2.breed",
  "eps.breed",

  "beta0.phia",
  "betasst.phia",
  "betawind.phia",

  "phia",
  "sigma2.phia",
  "eps.phia",

  "N1", "N2", "N3", "N4", "N5", "Nbreed", "Nonbreed", "Ntot")

## MCMC settings, trial run
ni <- 1000 ## increase to 300k
nb <- 100 ## increase nb to 100k
nc <- 3 ## increase to 6 (see effect on neff)
nt <- 10 ## up to 1000
n.adapt <- 10

## MCMC settings, proper run
# ni <- 1000000
# nb <- 400000
# nc <- 6
# nt <- 1000
# n.adapt <- 10000

ipm.result <- NULL

# Initial values
random.inits <- function(){list(mu.phi1=runif(1, 0.5, 0.8),
                                mu.p=runif(1,0.5,1),

                                beta0.fec= runif(1,-0.5, 1), ##rnorm(1, 0, sqrt(2.78)),
                                betanow.fec=runif(1,-2,2), ##rnorm(1, 0, sqrt(2.78)),

                                sigma2.fec=rtruncnorm(1, a=0, mean=0, sd=1),

                                beta0.breed= rtruncnorm(1, a=0.75, b=1.75, mean=1.25, sd=0.4), ##rnorm(1, 0, sqrt(2.78)),
                                betavwind.breed=runif(1,-2,2), ##rnorm(1, 0, sqrt(2.78)),

                                sigma2.breed=truncnorm::rtruncnorm(1, a=0, b=0.2, mean=0, sd=1),

                                beta0.phia= runif(1, 2, 3.5),
                                betasst.phia=runif(1,-2,2), ##rnorm(1, 0, sqrt(2.78)),
                                betawind.phia=runif(1,-2,2), ##rnorm(1, 0, sqrt(2.78)),

                                sigma2.phia=rtruncnorm(1, a=0, b=0.2, mean=0, sd=1))}


#Running the model using JagsUI
start <- as.POSIXlt(Sys.time())
print(start)

## run once
ipm.result <- jags(data = ipm.data,
                   inits = random.inits,
                   params,
                   "model/bestmodel_final.txt",
                   n.chains=nc,
                   n.iter=ni,
                   n.burnin=nb,
                   n.thin=nt)

# ## export it
# saveRDS(ipm.result, file = "bestmodel_final_1979_2020.rds")


#################################################################################
## INSPECT RESULTS
#################################################################################

ipm.result <- readRDS(file = "bestmodel_final_1979_2020.rds")

library(MCMCvis)

MCMCvis::MCMCsummary(ipm.result) ## summary

MCMCvis::MCMCsummary(object = ipm.result, round = 3,
                               #Rhat = F, ## if failed to compute
                               #n.eff = F, ## if failed to compute
                               params = c(
                                 "mu.phi1"))

## plot, e.g. N breeders and fecundity
MCMCvis::MCMCplot(ipm.result, params="Nbreed", horiz = FALSE)
MCMCvis::MCMCplot(ipm.result, params="fec", horiz = FALSE)

# MCMCvis::MCMCtrace(object = ipm.result,
#                    #pdf = FALSE,
#                    ind = TRUE)

mean(as.vector(ipm.result$BUGSoutput$sims.list$fec))
quantile(as.vector(ipm.result$BUGSoutput$sims.list$fec), c(0.025, 0.975))
