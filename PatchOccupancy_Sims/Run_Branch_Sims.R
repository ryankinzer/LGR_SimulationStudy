# Author: Kevin See
# Purpose: Run branching model scenarios
# Created: 9/29/2016
# Last Modified: 9/29/2016
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
# install.packages(c('lubridate', 'magrittr', 'dplyr', 'tidyr', 'jagsUI'))

#-----------------------------------------------------------------
library(MCMCpack)
library(lubridate)
library(magrittr)
library(plyr)
library(dplyr)
library(tidyr)
library(jagsUI)
library(RPushbullet)

setwd('Dropbox/LGR_SimulationStudy/PatchOccupancy_Sims')
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
mcmc.chainLength = 7000

# number of burn-in samples
mcmc.burn = 2000

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

jags.params = c('main_p', 'det_p')

#-----------------------------------------------------------------
# how mnay simulations to do?
n_sim = 2

# set trap rate on weekly basis
# make it consistent and constant throughout the season
my_trap_rate = data.frame(Week = 1:52,
                          trap.rate = 0.08,
                          trap.open = T)

#-----------------------------------------------------------------
# run simulations
#-----------------------------------------------------------------

sim = SimulateBranchData(trap.rate.df = my_trap_rate)

valid_df = sim$valid_tags %>%
  mutate(Week = Week - min(Week) + 1) %>%
  mutate(id = 1:n())

lgr_truth = sim$lgr_truth %>%
  mutate(Week = Week - min(Week) + 1) %>%
  mutate(id = 1:n())

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
                   as.matrix())

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
cat(paste('Took', round(c(proc.time() - ptm)[3] / 60, 2), 'min to run. \n'))
if(class(branch_mod) == 'try-error') next

library(ggmcmc)
# library(stringr)
my_ggs = ggs(branch_mod, family = c('main_p'))
tmp = my_ggs %>%
  filter(value > 0) %>%
  filter(grepl('\\,10\\]', Parameter))
attributes(tmp)[names(attributes(my_ggs))[!names(attributes(my_ggs)) %in% names(attributes(tmp))]] = attributes(my_ggs)[!names(attributes(my_ggs)) %in% names(attributes(tmp))]
attributes(tmp)$class = attributes(my_ggs)$class
my_ggs = tmp

ggs_traceplot(my_ggs) +
  facet_wrap(~ Parameter, scales = 'free')
ggs_density(my_ggs) +
  facet_wrap(~ Parameter, scales = 'free')
ggs_Rhat(my_ggs)
ggs_geweke(my_ggs)


# compare with truth
summ = summary(branch_mod)$quantiles %>%
  as.data.frame() %>% 
  mutate(param = rownames(.)) %>%
  tbl_df() %>%
  filter(grepl('^main_p', param)) %>%
  mutate(week_num = ifelse(grepl('^main_p', param), str_sub(param, 8, 9), NA),
         week_num = gsub('\\,', '', week_num),
         week_num = as.integer(as.character(week_num)),
         branch = sapply(str_split(param, '\\,'), function(y) y[2]),
         branch = gsub(']', '', branch),
         branch = as.integer(branch))

ggplot(summ,
       aes(x = week_num,
           y = `50%`,
           color = as.factor(branch))) +
  geom_ribbon(aes(ymin = `25%`,
                  ymax = `75%`,
                  fill = as.factor(branch)),
              alpha = 0.2,
              color = NA) +
  scale_color_brewer(palette = 'Paired') +
  scale_fill_brewer(palette = 'Paired') +
  geom_line() +
  theme_bw()


lgr_truth %>%
  group_by(Branch) %>%
  summarise(n_fish = n_distinct(id))

lgr_truth %>%
  group_by(Week, Branch) %>%
  summarise(branch_fish = n_distinct(id)) %>%
  left_join(lgr_truth %>%
              group_by(Week) %>%
              summarise(tot_fish = n_distinct(id))) %>%
  mutate(prob = branch_fish / tot_fish,
         prob_se = sqrt(prob * (1 - prob) / tot_fish)) %>%
  ggplot(aes(x = Week,
             y = prob,
             color = Branch)) +
  geom_ribbon(aes(ymin = prob - prob_se,
                  ymax = prob + prob_se,
                  fill = Branch),
              alpha = 0.2,
              color = NA) +
  scale_color_brewer(palette = 'Paired') +
  scale_fill_brewer(palette = 'Paired') +
  geom_line() +
  theme_bw()

lgr_truth %>%
  group_by(Week) %>%
  summarise(tot_fish = n_distinct(id)) %>%
  left_join(summ %>%
              select(Week = week_num,
                     branch,
                     prob = `50%`)) %>%
  mutate(branch_fish = round(tot_fish * prob)) %>%
  group_by(branch) %>%
  summarise(n_fish = sum(branch_fish))

lgr_truth %>%
  group_by(Branch) %>%
  summarise(n_fish = n_distinct(id))

det_est = summary(branch_mod)$quantiles %>%
  as.data.frame() %>% 
  mutate(param = rownames(.)) %>%
  tbl_df() %>%
  filter(grepl('^det_p', param)) %>%
  mutate(site_num = str_sub(param, 7, 8),
         site_num = gsub(']', '', site_num),
         site_num = as.integer(site_num),
         branch = ceiling(site_num / 2),
         site = ifelse(site_num %% 2 == 1, 'Lower', 'Upper')) %>%
  select(branch, site, everything(), -param, -site_num) %>%
  mutate(pop_size = rep(c('small', 'med', 'large'), each = 6),
         true_val = rep(c(0.5, 0.5, 0.5, 0.95, 0.95, 0.95), 3)) %>%
  mutate(inCI = ifelse(`2.5%` <= true_val & `97.5%` >= true_val, T, F))