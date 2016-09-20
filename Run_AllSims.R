# Author: Kevin See
# Purpose: Run all simulation scenarios
# Created: 9/20/2016
# Last Modified: 9/20/2016
# Notes: 

#-----------------------------------------------------------------
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
library(ggplot2)
library(jagsUI)

theme_set(theme_bw())
setwd('IDFG_Comparison/SimulationStudy')

# this contains the function to simulate data
source('SimFnc.R')

#-----------------------------------------------------------------
# get file names of scenario scripts
#-----------------------------------------------------------------
script_nms = list.files('SimScenarios')
script_nms = script_nms[grepl('^RunSims', script_nms)]
length(script_nms)

for(i in 1:length(script_nms)) {
  source(paste0('SimScenarios/', script_nms[i]))
}

