# Author: Kevin See
# Purpose: Test several versions of ISEMP model of LGR escapement
# Created: 8/30/2016
# Last Modified: 9/1/2016
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

source('../SimFnc.R')

#-----------------------------------------------------------------
# set mcmc parameters
#-----------------------------------------------------------------
# number of total samples
# mcmc.chainLength = 12000
mcmc.chainLength = 40000

# number of burn-in samples
# mcmc.burn = 2000
mcmc.burn = 10000

# thinning interval
# mcmc.thin = 20
mcmc.thin = 30

# number of MCMC chains
mcmc.chains = 4

# how many poseterior samples will we have?
data.frame(Each.Chain = (mcmc.chainLength - mcmc.burn) / mcmc.thin, All.Chains = (mcmc.chainLength - mcmc.burn) / mcmc.thin * mcmc.chains)

#-----------------------------------------------------------------
# generate simulated data
#-----------------------------------------------------------------
# determine trap rate on weekly basis
my_trap_rate = data.frame(Week = 1:52,
                          trap.rate = 0.05)

set.seed(4)
my_sim = SimulateLGRdata(trap.rate.df = my_trap_rate)

my_sim = SimulateLGRdata(N.lgr = 1e5,
                         # n.pops = c(3,1),
                         # run.mu.mu = c(140, 171),
                         # run.mu.sd = c(16, 15),
                         # run.sd.mu = c(16,25),
                         # run.sd.sd = c(8,2),
                         hatch.pop.prob = 0,
                         hnc.pop.prob = 0,
                         fallback.rate = 0,
                         reascension.rate = 1,
                         night.passage.rate = 0,
                         window.rate = 5/6,
                         marked.rate = 0.10,
                         ladder.det = 0.99,
                         trap.rate.df = my_trap_rate)

# pull out data from simulation
lgr_truth = my_sim$sim

# add date for arbitary year, and beginning date of each week
lgr_truth %<>%
  mutate(Date = ymd('20150101') + days(Day)) %>%
  group_by(Week) %>%
  mutate(Start_Date = min(Date)) %>%
  ungroup()

# generate weekly observations
lgr_week = SimulateLGRobs(my_sim$parameters, lgr_truth, theta = 5)

# filter for spring/summer Chinook dates
lgr_week %<>%
  filter(Start_Date >= ymd('20150301'),
         Start_Date <= ymd('20150817')) %>%
  mutate(Week = Week - min(Week) + 1)

lgr_truth %<>%
  filter(Start_Date >= min(lgr_week$Start_Date),
         Start_Date <= max(lgr_week$Start_Date)) %>%
  mutate(Week = Week - min(Week) + 1)


# pull out true value of various parameters
true_var = lgr_week %>%
  summarise(All.Fish = sum(N_tot),
            Daytime.Fish = sum(Day.fish),
            Reascent.Fish = sum(Reascent),
            Unique.Wild.Fish.PBT = sum(N_uniq_wild)) %>%
  bind_cols(lgr_truth %>%
              filter(Origin == 'NOR') %>%
              summarise(All.Wild.Fish.PBT = n_distinct(id) + sum(Reascent),
                        Night.Wild.Fish.PBT = sum(Night.passage),
                        Wild.Reascents.PBT = sum(Reascent))) %>%
  t()
true_var = data.frame(Variable = row.names(true_var),
                      Truth = true_var) %>% tbl_df()


# pull data together for JAGS
jags.data = list('TotLadderWeeks' = nrow(lgr_week),
                 'Y.window' = lgr_week %>% select(win_cnt) %>% as.matrix() %>% as.vector(),
                 'Y.window.log' = lgr_week %>% mutate(log_win = log(win_cnt + 0.5)) %>% select(log_win) %>% as.matrix() %>% as.vector(),
                 'Y.trap' = lgr_week %>% select(wild_fish) %>% as.matrix() %>% as.vector(),
                 'trap.fish' = lgr_week %>% select(Trap) %>% as.matrix() %>% as.vector(),
                 'wild.pbt' = lgr_week %>% select(wild_fish) %>% as.matrix() %>% as.vector(),
                 'trap.rate' = lgr_week %>% select(trap_rate) %>% as.matrix() %>% as.vector(),
                 'trap.alpha' = lgr_week %>% select(trap_alpha) %>% as.matrix() %>% as.vector(),
                 'trap.beta' = lgr_week %>% select(trap_beta) %>% as.matrix() %>% as.vector(),
                 'Tot.tags' = lgr_week %>% select(Ladder) %>% as.matrix() %>% as.vector(),
                 'ReAsc.tags' = lgr_week %>% select(mark_reasc) %>% as.matrix() %>% as.vector(),
                 'DC.tags' = lgr_week %>% mutate(mark_day = Ladder - mark_night) %>% select(mark_day) %>% as.matrix() %>% as.vector()
)


# which parameters to track posteriors of
jags.params = c('X.log.all', 'X.all', 'X.day', 'X.night', 'X.reasc', 'X.all.wild.pbt', 'X.new.wild.pbt', 'X.reasc.wild.pbt', 'X.night.wild.pbt', 'X.tot.all', 'X.tot.day', 'X.tot.night', 'X.tot.reasc', 'X.tot.all.wild.pbt', 'X.tot.new.wild.pbt', 'X.tot.night.wild.pbt', 'X.tot.reasc.wild.pbt', 'prop.tagged.pbt', 'X.sigma', 'true.prop', 'win.prop.avg', 'win.prop.true', 'win.prop.sigma', 'hist.prop', 'reasc.avg', 'reasc.true', 'reasc.sigma', 'acf', 'wnc.avg', 'wnc.true', 'wnc.sigma', 'trap.rate.true', 'r', 'k', 'p', 'alpha')

# set initial values
jags.inits = function(){
  return(list('X.log.all' = with(jags.data, log(Y.window + 1)),
              'trap.rate.true' = with(jags.data, trap.rate + 0.0001)))
}

# versions of model to run
model.locs = c('LGR_TotalEscape_Pois.txt',
               'LGR_TotalEscape_NB1.txt',
               'LGR_TotalEscape_NB2.txt',
               'LGR_TotalEscape_NB3.txt',
               'LGR_TotalEscape_LogNorm.txt')

mod_list = vector('list', length(model.locs))
names(mod_list) = c('Pois', 'NB1', 'NB2', 'NB3', 'Lnorm')

for(i in 1:length(model.locs)) {
  model.loc = model.locs[i]
  cat(paste('Starting model', names(mod_list)[i], '\n'))
  
  ptm = proc.time()
  mod_list[[i]] = try(jags(data = jags.data, 
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
  rm(ptm, model.loc)
}

true_var %>%
  left_join(ldply(mod_list, .id = 'Model',
                  .fun = function(x) {
                    data.frame(All.Fish = x$q50$X.tot.all,
                               Daytime.Fish = x$q50$X.tot.day,
                               Reascent.Fish = x$q50$X.tot.reasc,
                               All.Wild.Fish.PBT = x$q50$X.tot.all.wild.pbt,
                               Unique.Wild.Fish.PBT = x$q50$X.tot.new.wild.pbt,
                               Night.Wild.Fish.PBT = x$q50$X.tot.night.wild.pbt,
                               Wild.Reascents.PBT = x$q50$X.tot.reasc.wild.pbt)
                  }) %>% tbl_df() %>%
              gather(Variable, est, -Model) %>%
              mutate(est = round(est)) %>%
              spread(Model, est))

ldply(mod_list[1:3], .id = 'Model', 
      .fun = function(x) {
        res = x$q50[c('r', 'X.sigma', 'win.prop.sigma', 'reasc.sigma', 'wnc.sigma')] %>%
          ldply(.id = 'param') %>% tbl_df() %>%
          rename(median = V1)
        if(!'r' %in% res$param) res = data.frame(param = 'r',
                                                 median = NA) %>%
            bind_rows(res)
        return(res)
      }) %>%
  spread(Model, median)