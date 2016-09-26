# Author: Kevin See
# Purpose: Run all simulation scenarios
# Created: 9/20/2016
# Last Modified: 9/20/2016
# Notes: 

#-----------------------------------------------------------------
# set up AMI
#-----------------------------------------------------------------
# excludeSyncDropbox("*")
# includeSyncDropbox("LGR_SimulationStudy")
# 
# # set up RPushbullet
# install.packages(c('jsonlite', 'RPushbullet'))
# 
# library(jsonlite)
# library(RPushbullet)
# 
# cat(toJSON(list(key="dCj80OM9XSYx1xRUyvCp5KlN4aT456Kg", devices=c('ujx5LBhEzSusjAiVsKnSTs', 'ujx5LBhEzSusjzWIEVDzOK'), names=c('phone', 'Chrome'))), file='~/.rpushbullet.json')
# 
# detach('package:RPushbullet', unload=T)
# 
# library(RPushbullet)
# print(pbGetDevices())
# 
# # Test
# # if recipients = 1, this should go to phone (2 should go to Chrome)
# pbPost("note", "Test", "This came from R!", recipients=c(1))
# 
# install.packages(c('MCMCpack', 'FSA', 'Rcapture', 'boot', 'msm', 'lubridate', 'magrittr', 'dplyr', 'tidyr', 'ggplot2', 'jagsUI'))

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
library(RPushbullet)

theme_set(theme_bw())

# setwd('Dropbox/LGR_SimulationStudy')

# this contains the function to simulate data
source('SimFnc.R')

# this will send a pushbullet notification about any errors
options(error = function() {
  library(RPushbullet)
  pbPost("note", "Error", geterrmessage(), recipients=1)
})



#-----------------------------------------------------------------
# get file names of scenario scripts
#-----------------------------------------------------------------
script_nms = list.files('SimScenarios')
script_nms = script_nms[grepl('^RunSims', script_nms)]
length(script_nms)

#-----------------------------------------------------------------
# run scenario scripts and save results
#-----------------------------------------------------------------
for(i in 1:length(script_nms)) {
  set.seed(4)
  pbPost('note', paste('Starting', gsub('^RunSims_', '', script_nms[i])), recipients=1)
  ptm <- proc.time()
  source(paste0('SimScenarios/', script_nms[i]))
  pbPost('note', paste(gsub('^RunSims_', '', script_nms[i]), 'finished'), paste('It took', round(c(proc.time() - ptm)[3] / 3600, 1), 'hours to run.'), recipients=NA)
}

