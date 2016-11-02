# Author: Kevin See
# Purpose: Run branching model scenarios
# Created: 9/29/2016
# Last Modified: 10/25/2016
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
# install.packages(c('lubridate', 'magrittr', 'dplyr', 'tidyr', 'jagsUI', 'MCMCpack'))

#-----------------------------------------------------------------
library(MCMCpack)
library(lubridate)
library(magrittr)
library(plyr)
library(dplyr)
library(tidyr)
library(jagsUI)
# library(RPushbullet)

setwd('PatchOccupancy_Sims')
# setwd('Dropbox/LGR_SimulationStudy/PatchOccupancy_Sims')
# this contains the function to simulate LGR data
source('../SimFnc.R')
# this contains the function to simulate tributary detection data
source('Branch_Sim.R')


# this will send a pushbullet notification about any errors
options(error = function() {
  library(RPushbullet)
  pbPost("note", "Error", geterrmessage(), recipients=1)
})


#-----------------------------------------------------------------
# set mcmc parameters
#-----------------------------------------------------------------
# number of total samples
mcmc.chainLength = 3000

# number of burn-in samples
mcmc.burn = 1500

# thinning interval
mcmc.thin = 10

# number of MCMC chains
mcmc.chains = 4

# how many poseterior samples will we have?
data.frame(Each.Chain = (mcmc.chainLength - mcmc.burn) / mcmc.thin, All.Chains = (mcmc.chainLength - mcmc.burn) / mcmc.thin * mcmc.chains)


#-----------------------------------------------------------------
# set up model location and which parameters to track with JAGS
#-----------------------------------------------------------------
model.loc = 'Simulation_JAGS.txt'

jags.params = c('main_p', 'det_p', 'week_esc', 'tot_esc')

#-----------------------------------------------------------------
# how mnay simulations to do?
n_sim = 100

# set trap rate on weekly basis
# make it consistent and constant throughout the season
my_trap_rate = data.frame(Week = 1:52,
                          trap.rate = 0.08,
                          trap.open = T)

# shut trap down for a few weeks at end of July, beginning of August
my_trap_rate %<>%
  mutate(trap.rate = ifelse(Week %in% 30:32, 0, trap.rate),
         trap.open = ifelse(Week %in% 30:32, F, trap.open))


#-----------------------------------------------------------------
# set up true values to be stored
sim_truth = data.frame(param = c(paste0('det_p[', 1:18, ']'),
                                 paste0('tot_esc[', 1:10, ']')),
                       true_value = c(rep(c(rep(0.5, 3), rep(0.95, 3)), 3),
                                      c(rep(250,3),rep(1000,3),rep(3000,3),12250)))

#-----------------------------------------------------------------
# run simulations
#-----------------------------------------------------------------
# set up list to store results
res_list = vector('list', n_sim)
names(res_list) = 1:n_sim

n_tags = vector('integer', n_sim)
names(n_tags) = 1:n_sim


set.seed(3)
for(i in 1:n_sim) {
  
  sim = SimulateBranchData(trap.rate.df = my_trap_rate)
  
  min_week = min(sim$lgr_truth$Week)
  
  valid_df = sim$valid_tags %>%
    mutate(Week = Week - min_week + 1) %>%
    mutate(id = 1:n())
  
  lgr_truth = sim$lgr_truth %>%
    mutate(Week = Week - min_week + 1) %>%
    mutate(id = 1:n())
  
  sim_true = sim_truth %>%
    bind_rows(lgr_truth %>%
                mutate(Branch = gsub('^Branch-', '', Branch),
                       Branch = revalue(Branch,
                                        c('Black-Box' = 10)),
                       Branch = as.integer(Branch)) %>%
                group_by(branch = Branch, Week) %>%
                summarise(true_value = n_distinct(id)) %>%
                ungroup() %>%
                mutate(param = paste0('week_esc[', Week, ',', branch, ']')) %>%
                select(param, true_value)) %>%
    bind_rows(lgr_truth %>%
                mutate(Branch = gsub('^Branch-', '', Branch),
                       Branch = revalue(Branch,
                                        c('Black-Box' = 10)),
                       Branch = as.integer(Branch)) %>%
                group_by(Week, branch = Branch) %>%
                summarise(branch_fish = n_distinct(id)) %>%
                ungroup() %>%
                left_join(lgr_truth %>%
                            group_by(Week) %>%
                            summarise(tot_fish = n_distinct(id)) %>%
                            ungroup()) %>%
                mutate(true_value = branch_fish / tot_fish,
                       param = paste0('main_p[', Week, ',', branch, ']')) %>%
                select(param, true_value))
  
  
  jags.data = list('n.pops.main' = length(unique(valid_df$Branch)),
                   'n.weeks' = max(valid_df$Week),
                   'n.fish' = nrow(valid_df),
                   'lgr_week' = valid_df$Week,
                   'main_dirch_vec' = rep(1, 10),
                   'obs.mat' = valid_df %>% 
                     select(id, Branch, Lower.obs, Upper.obs) %>% 
                     gather(array, seen, -id, -Branch) %>%
                     mutate(site = paste(Branch, array, sep = '_')) %>%
                     select(id, site, seen) %>%
                     spread(site, seen, fill = 0) %>%
                     select(matches('^Branch')) %>%
                     as.matrix(),
                   'lgr.esc' = lgr_truth %>%
                     group_by(Week) %>%
                     summarise(n_fish = n_distinct(id)) %>%
                     ungroup() %>%
                     left_join(data.frame(Week = 1:max(valid_df$Week))) %>%
                     mutate(n_fish = ifelse(is.na(n_fish), 0, n_fish)) %>%
                     select(n_fish) %>%
                     as.matrix() %>% as.vector())
  
  # set initial values based on observed fish
  jags.inits = jagsInits(valid_df)
  
  # fit model with JAGS
  ptm = proc.time()
  branch_mod = try(jags.basic(data = jags.data, 
                              inits = jags.inits,
                              parameters.to.save = jags.params, 
                              model.file = model.loc, 
                              # n.chains = 1, 
                              # n.burnin = 10, 
                              # n.thin = 2, 
                              # n.iter = 20, 
                              n.chains = mcmc.chains,
                              n.burnin = mcmc.burn,
                              n.thin = mcmc.thin,
                              n.iter = mcmc.chainLength,
                              parallel = T,
                              DIC = FALSE,
                              verbose = T))
  cat(paste('Sim', i, 'took', round(c(proc.time() - ptm)[3] / 60, 2), 'min to run. \n'))
  if(class(branch_mod) == 'try-error') {
    rm(sim, valid_df, lgr_truth, jags.data, branch_mod, sim_true)
    next
  }
  
  # pull out what to save
  res_list[[i]] = summary(branch_mod)$quantiles %>%
    as.data.frame() %>% 
    mutate(param = rownames(.)) %>%
    select(param, everything()) %>%
    tbl_df() %>%
    left_join(sim_true) %>%
    mutate(true_value = ifelse(is.na(true_value), 0, true_value))
  
  n_tags[i] = nrow(valid_df)
  
  rm(sim, valid_df, lgr_truth, jags.data, branch_mod, sim_true)
}

# save results
save(res_list, n_tags, file = 'SimFits/TrapDown_2000tags.rda')

pbPost("note", "Trapdown sims are done", recipients=1)
