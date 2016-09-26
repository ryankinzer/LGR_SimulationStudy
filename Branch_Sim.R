#------------------------------------------------------------------------------
# Simulate 9 branches and 1 black box data set for testing branching model.
#------------------------------------------------------------------------------
# Ryan N. Kinzer
# Created: 9/23/2016
# Modified: 9/23/2016
#------------------------------------------------------------------------------

library(MCMCpack)
library(FSA)
library(Rcapture)
library(boot)
library(msm)
library(lubridate)
library(magrittr)
library(plyr)
library(dplyr)
library(tidyr)
library(jagsUI)
library(ggplot2)

#source("./Script/lgd_run_timing.R") # generates the run timing parameters
# shown below.

# Branch run timing parameters
# A tibble: 3 × 6
#Group     n       mu    mu_sd    sd_mu    sd_sd
#<fctr> <int>    <dbl>    <dbl>    <dbl>    <dbl>
#1  Early    86 150.5457 19.91180 14.25784 6.954262
#2    Mid    42 170.9983 12.80495 22.21365 7.361752
#3   Late    21 182.9837 19.21272 22.16461 6.590507

#Black box run timing parameteres
# A tibble: 1 × 6 
#Species     n       mu    mu_sd    sd_mu    sd_sd
#<chr> <int>    <dbl>    <dbl>    <dbl>    <dbl>
#1 Chinook   149 160.8827 21.98749 17.52383 7.999656

source("SimFnc.R") # source LGR simulation function

my_trap_rate = data.frame(Week = 1:52,
                          trap.rate = 0.08,
                          trap.open = T)

dat <- SimulateLGRdata(N.lgr = 2500,
                       n.pops = c(1,1),
                       run.mu.mu = c(160, 1),
                       run.mu.sd = c(22, 1),
                       run.sd.mu = c(17,1),
                       run.sd.sd = c(8,1),
                       hatch.pop.prob = 0.1,
                       hnc.pop.prob = 0.1,
                       fallback.rate = 0.06,
                       reascension.rate = 1,
                       night.passage.rate = .06,
                       window.rate = 50/60,
                       marked.rate = 0.0,
                       ladder.det = 0.0,
                       start.date = ymd('20150101'),
                       trap.rate.df = my_trap_rate)

tmp <- SimulateLGRdata(N.lgr = 1000,
                trap.rate.df = my_trap_rate)

tmp <- SimulateLGRdata(N.lgr = 250, trap.rate.df = my_trap_rate)


