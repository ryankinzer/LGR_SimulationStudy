
#---------------------------------------------
# source simulation functions
#--------------------------------------------
source("./Scripts/SimFnc.R")
source("./Scripts/SCOBIv2.R")
source("./Scripts/SCSrank.R")

#--------------------------------------------
# load library
#--------------------------------------------
library(SCOBI)
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
library(Hmisc) # needed for SCOBI model
#library(stats)
#library(stringr)
#library(car)

#library(jagsUI)

#----------------------------------------
# Raw LGTrappingDB data
#----------------------------------------
#View(exRawChnkAdultData)

#----------------------------------------
# Data formated for SCOBI function
#----------------------------------------
#View(chnkScobiInput)
#View(chnkWindowCounts)

#----------------------------------------
# Example SCOBI run
#----------------------------------------
#SCOBI(adultData = chnkScobiInput, windowData = chnkWindowCounts, Run = "chnkDemo",
#      RTYPE = "W", Primary = "GenStock", Secondary = NA, alph = 0.1, B = 100)


#----------------------------------------
# Simulate Lower Granite Data
#----------------------------------------

# set trap rate on weekly basis
my_trap_rate = data.frame(Week = 1:52,
                          trap.rate = 0.07)

# # shut trap down for a few weeks
# my_trap_rate %<>%
#   mutate(trap.rate = ifelse(Week %in% 22:24, 0, trap.rate))

# # change trap rate part-way through season
# my_trap_rate %<>%


iloop <- 10

Nstar <- data.frame(matrix(data=NA, nrow=iloop, ncol = 4))

for(i in 1:iloop){
  tmp <- SimulateLGRdata(N.lgr = 100000,
                       n.pops = c(17,5),
                       run.mu.mu = c(140, 171),
                       run.mu.sd = c(16, 15),
                       run.sd.mu = c(16,25),
                       run.sd.sd = c(8,2),
                       hatch.pop.prob = 0.5,
                       hnc.pop.prob = 0.2,
                       fallback.rate = .05,
                       reascension.rate = 1,
                       night.passage.rate = .05,
                       window.rate = 50/60,
                       marked.rate = 0.05,
                       ladder.det = 0.99,
                       start.date = ymd('20150101'),
                       trap.rate.df = my_trap_rate)

# Format data for use in ISEMP model
model_dat <- SimulateLGRobs(tmp$parameters,
               tmp$sim,
               theta = 100,
               perfect.window = T)

# Format data for use in SCOBI model
scobi_dat <- formatSCOBI_inputs(model_dat$obs,tmp$sim) # uses the same output for ISEMP model run

scobi <- SCOBIv2(adultData = scobi_dat$fish_data, windowData = scobi_dat$window_count, Run = "chnkDemo",
          RTYPE = "W", Primary = "GenStock", Secondary = NA, alph = 0.1, B = 100, writeOutput = FALSE)


Nstar[i,1] <- dim(tmp$sim[tmp$sim$Origin == "NOR" & tmp$sim$Reascent == 0,])[1]  # Truth
Nstar[i,2] <- scobi$Rearing[1,1]
Nstar[i,3] <- scobi$Rearing[1,2] 
Nstar[i,4] <- scobi$Rearing[1,3]
}

names(Nstar) = c("Truth","SCOBI","L_SCOBI","U_SCOBI")
