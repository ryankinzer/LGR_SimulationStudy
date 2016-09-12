# Author: Kevin See
# Purpose: Simulate a number of LGR-like datasets, and estimate escapement using ISEMP model
# Created: 8/31/2016
# Last Modified: 9/2/2016
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
# set mcmc parameters
#-----------------------------------------------------------------
# number of total samples
mcmc.chainLength = 12000

# number of burn-in samples
mcmc.burn = 2000

# thinning interval
mcmc.thin = 20

# number of MCMC chains
mcmc.chains = 4

# how many poseterior samples will we have?
data.frame(Each.Chain = (mcmc.chainLength - mcmc.burn) / mcmc.thin, All.Chains = (mcmc.chainLength - mcmc.burn) / mcmc.thin * mcmc.chains)


#-----------------------------------------------------------------
# set up model location and which parameters to track with JAGS
#-----------------------------------------------------------------
model.loc = 'LGR_TotalEscape_JAGS.txt'

jags.params = c('X.log.all', 'X.all', 'X.day', 'X.night', 'X.reasc', 'X.all.wild.pbt', 'X.new.wild.pbt', 'X.reasc.wild.pbt', 'X.night.wild.pbt', 'X.tot.all', 'X.tot.day', 'X.tot.night', 'X.tot.reasc', 'X.tot.all.wild.pbt', 'X.tot.new.wild.pbt', 'X.tot.night.wild.pbt', 'X.tot.reasc.wild.pbt', 'prop.tagged.pbt', 'X.sigma', 'true.prop', 'win.prop.avg', 'win.prop.true', 'win.prop.sigma', 'hist.prop', 'reasc.avg', 'reasc.true', 'reasc.sigma', 'acf', 'wnc.avg', 'wnc.true', 'wnc.sigma', 'trap.rate.true', 'r', 'k')

# set initial values
jags.inits = function(){
  return(list('X.log.all' = with(jags.data, log(Y.window + 1)),
              'trap.rate.true' = with(jags.data, trap.rate + 0.0001)))
}

#-----------------------------------------------------------------
# how mnay simulations to do?
n_sim = 10

# set trap rate on weekly basis
my_trap_rate = data.frame(Week = 1:52,
                          trap.rate = 0.07)

# # shut trap down for a few weeks
# my_trap_rate %<>%
#   mutate(trap.rate = ifelse(Week %in% 22:24, 0, trap.rate))

# # change trap rate part-way through season
# my_trap_rate %<>%
#   mutate(trap.rate = ifelse(Week > 20, 0.15, trap.rate))


#-----------------------------------------------------------------
# run simulations
#-----------------------------------------------------------------
# set up list to capture comparison results
res = vector('list', n_sim)
names(res) = 1:n_sim
mod_list = obs_list = sim_list = res

set.seed(7)

for(i in 1:n_sim) {
  cat(paste('Starting simulation #', i, '\n'))
  
  # simulate data
  # my_sim = SimulateLGRdata(trap.rate.df = my_trap_rate)

  my_sim = SimulateLGRdata(trap.rate.df = my_trap_rate,
                           fallback.rate = 0.12,
                           reascension.rate = 0.9,
                           night.passage.rate = 0.05,
                           window.rate = 1,
                           marked.rate = 0.07,
                           ladder.det = 0.99)
  
  # pull out data from simulation
  lgr_truth = my_sim$sim
  
  # generate weekly observations
  # theta is used to control how much observation error is put on window counts. Higher theta = less error
  lgr_week = SimulateLGRobs(my_sim$parameters, lgr_truth, theta = 100, perfect.window = F)
  
  # filter for spring/summer Chinook dates
  lgr_truth %<>%
    filter(Start_Date >= ymd('20150301'),
           Start_Date <= ymd('20150817'))
  
  lgr_week$obs %<>%
    inner_join(lgr_truth %>%
                select(Week, Start_Date) %>%
                distinct()) %>%
    select(Start_Date, Week, everything()) %>%
    mutate(Week = Week - min(Week) + 1)
  
  lgr_week$by_origin %<>%
    inner_join(lgr_truth %>%
                 select(Week, Start_Date) %>%
                 distinct()) %>%
    select(Start_Date, Week, Origin, everything()) %>%
    mutate(Week = Week - min(Week) + 1)
  
  # pull data together for JAGS
  jags.data = list('TotLadderWeeks' = nrow(lgr_week$obs),
                   'Y.window' = lgr_week$obs %>% select(win_cnt) %>% as.matrix() %>% as.vector(),
                   'Y.trap' = lgr_week$obs %>% select(Trap) %>% as.matrix() %>% as.vector(),
                   'trap.fish' = lgr_week$obs %>% select(Trap) %>% as.matrix() %>% as.vector(),
                   'wild.pbt' = lgr_week$obs %>% select(wild_fish) %>% as.matrix() %>% as.vector(),
                   'trap.rate' = lgr_week$obs %>% select(trap_rate) %>% as.matrix() %>% as.vector(),
                   'trap.alpha' = lgr_week$obs %>% select(trap_alpha) %>% as.matrix() %>% as.vector(),
                   'trap.beta' = lgr_week$obs %>% select(trap_beta) %>% as.matrix() %>% as.vector(),
                   'Tot.tags' = lgr_week$obs %>% select(Ladder) %>% as.matrix() %>% as.vector(),
                   'ReAsc.tags' = lgr_week$obs %>% select(mark_reasc) %>% as.matrix() %>% as.vector(),
                   'DC.tags' = lgr_week$obs %>% mutate(mark_day = Ladder - mark_night) %>% select(mark_day) %>% as.matrix() %>% as.vector()
  )
  
  # fit model with JAGS
  ptm = proc.time()
  adult.pass.mod = try(jags(data = jags.data, 
                            inits = jags.inits, 
                            parameters.to.save = jags.params, 
                            model.file = model.loc, 
                            n.chains = mcmc.chains, 
                            n.burnin = mcmc.burn, 
                            n.thin = mcmc.thin, 
                            n.iter = mcmc.chainLength, 
                            parallel = T,
                            DIC = FALSE,
                            verbose = F))
  cat(paste('Took', round(c(proc.time() - ptm)[3] / 60, 2), 'min to run. \n'))
  
  # pull out results from posteriors
  attach(adult.pass.mod$sims.list)
  tot_post = ldply(list('All.Fish' = X.tot.all, 
                        'All.Wild.Fish.PBT' = X.tot.all.wild.pbt, 
                        'Unique.Wild.Fish.PBT' = X.tot.new.wild.pbt, 
                        'Daytime.Fish' = X.tot.day, 
                        'Reascent.Fish' = X.tot.reasc, 
                        'Night.Wild.Fish.PBT' = X.tot.night.wild.pbt, 
                        'Wild.Reascents.PBT'= X.tot.reasc.wild.pbt), 
                   .id='Variable') %>% tbl_df() %>%
    gather(iteration, value, -Variable) %>%
    mutate(iteration = gsub('^V', '', iteration),
           iteration = as.integer(iteration))
  
  detach(adult.pass.mod$sims.list)
  
  # summarise posterior
  tot_summ = tot_post %>%
    group_by(Variable) %>%
    summarise(mean = mean(value),
              median = median(value),
              se = sd(value),
              cv = se / mean,
              low_ci = HPDinterval(as.mcmc(value), prob = 0.95)[,1],
              upp_ci = HPDinterval(as.mcmc(value), prob = 0.95)[,2]) %>%
    ungroup()
  
  # compare with "truth" from simulated data
  true_var = lgr_truth %>%
    summarise(All.Fish = length(id),
              Unique.Fish = n_distinct(id),
              Daytime.Fish = sum(Day.passage),
              Reascent.Fish = sum(Reascent)) %>%
    bind_cols(lgr_truth %>%
                filter(Origin == 'NOR') %>%
                summarise(All.Wild.Fish.PBT = length(id),
                          Unique.Wild.Fish.PBT = n_distinct(id),
                          Night.Wild.Fish.PBT = sum(Night.passage),
                          Wild.Reascents.PBT = sum(Reascent))) %>%
    t() %>% as.data.frame()
  true_var = data.frame(Variable = row.names(true_var),
                        Truth = true_var$V1) %>% tbl_df()
  
  
  res[[i]] = true_var %>%
    inner_join(tot_summ) %>%
    mutate(inCI = ifelse(Truth >= low_ci & Truth <= upp_ci, T, F))
  
  sim_list[[i]] = my_sim
  obs_list[[i]] = lgr_week
  mod_list[[i]] = adult.pass.mod
  
  rm(my_sim, lgr_truth, lgr_week, jags.data, adult.pass.mod, true_var, tot_summ, tot_post)
  
  cat(paste('Finished simulation #', i, '\n'))
}

# save results
# save(res, sim_list, obs_list, mod_list, file = 'SimulationFits/SimResults.rda')

#-----------------------------------------------------------------
# Examine results
#-----------------------------------------------------------------
res_df = ldply(res, .id = 'sim') %>% tbl_df

res_df %>%
  group_by(Variable) %>%
  summarise(coverage95 = sum(inCI) / n())

qplot(Truth, median, data = res_df, color = Variable, log = 'xy') + 
  geom_abline()

res_df %>%
  mutate(bias = median - Truth,
         rel_bias = bias / Truth) %>%
  ggplot(aes(x = Variable,
             fill = Variable,
             y = rel_bias)) +
  geom_boxplot() +
  geom_hline(yintercept = 0,
             linetype = 2)

