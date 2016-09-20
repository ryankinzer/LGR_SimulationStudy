# Author: Kevin See
# Purpose: Simulate a number of LGR-like datasets, and estimate escapement using ISEMP model
# Created: 8/31/2016
# Last Modified: 9/12/2016
# Notes: 

#-----------------------------------------------------------------
# these functions contains the scripts to run the SCOBI estimator
source('SCOBI/SCOBIv2.R')
source('SCOBI/SCSrank.R')

# load needed packages
library(Hmisc) # needed for SCOBI model
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

# this contains the function to simulate data
source('SimFnc.R')

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
model.loc = 'LGR_TotalEscape_AllOrigins.txt'

jags.params = c('X.tot.all', 'X.tot.day', 'X.tot.night', 'X.tot.reasc', 'X.tot.new.wild', 'X.tot.new.hatch', 'X.tot.new.hnc', 'X.sigma', 'true.prop', 'win.prop.avg', 'win.prop.true', 'win.prop.sigma', 'reasc.avg', 'reasc.true', 'reasc.sigma', 'acf', 'org.prop', 'org.sigma', 'trap.rate.true', 'r', 'k')

# set initial values
jags.inits = function(){
  return(list('X.log.all' = with(jags.data, log(Y.window + 1)),
              'trap.rate.true' = with(jags.data, trap.rate + 0.0001)))
}

#-----------------------------------------------------------------
# how mnay simulations to do?
n_sim = 500

# set trap rate on weekly basis
# make it consistent and constant throughout the season
my_trap_rate = data.frame(Week = 1:52,
                          trap.rate = 0.15)

#-----------------------------------------------------------------
# run simulations
#-----------------------------------------------------------------
# set up list to capture comparison results
res = vector('list', n_sim)
names(res) = 1:n_sim
mod_list = obs_list = sim_list = res

set.seed(5)

tot_ptm = proc.time()
i = 1
while(i <= n_sim) {
  cat(paste('Starting simulation #', i, '\n'))
  
  # simulate data
  my_sim = SimulateLGRdata(trap.rate.df = my_trap_rate)
  
  # pull out data from simulation
  lgr_truth = my_sim$sim
  
  # generate weekly observations
  # error rate is the CV of the daily window counts
  lgr_week = SimulateLGRobs(my_sim$parameters, lgr_truth, error_rate = 0.05, perfect.window = T)
  
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
  
  lgr_truth %<>%
    mutate(Week = Week - min(Week) + 1)
  
  #-------------#
  # ISEMP model #
  #-------------#
  
  # pull data together for JAGS
  org_exist = lgr_week$obs %>% select(wild_fish, hatch_fish, HNC_fish) %>% colSums()
  org_exist = ifelse(org_exist > 0, 1, org_exist)
  jags.data = list('TotLadderWeeks' = nrow(lgr_week$obs),
                   'Y.window' = lgr_week$obs %>% select(win_cnt) %>% as.matrix() %>% as.vector(),
                   'Y.trap' = lgr_week$obs %>% select(Trap) %>% as.matrix() %>% as.vector(),
                   'trap.fish' = lgr_week$obs %>% select(Trap) %>% as.matrix() %>% as.vector(),
                   'wild.pbt' = lgr_week$obs %>% select(wild_fish) %>% as.matrix() %>% as.vector(),
                   'hatch.pbt' = lgr_week$obs %>% select(hatch_fish) %>% as.matrix() %>% as.vector(),
                   'hnc.pbt' = lgr_week$obs %>% select(HNC_fish) %>% as.matrix() %>% as.vector(),
                   'org.exist' = org_exist,
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
  if(class(adult.pass.mod) == 'try-error') next
  
  # pull out results from posteriors
  attach(adult.pass.mod$sims.list)
  tot_post = ldply(list('All.Fish' = X.tot.all, 
                        'Daytime.Fish' = X.tot.day, 
                        'Reascent.Fish' = X.tot.reasc, 
                        'Night.Fish' = X.tot.night,
                        'Unique.Wild.Fish' = X.tot.new.wild,
                        'Unique.Hatch.Fish' = X.tot.new.hatch,
                        'Unique.HNC.Fish' = X.tot.new.hnc), 
                   .id='Variable') %>% tbl_df() %>%
    gather(iteration, value, -Variable) %>%
    mutate(iteration = gsub('^V', '', iteration),
           iteration = as.integer(iteration))
  
  detach(adult.pass.mod$sims.list)
  
  # summarise posterior
  tot_summ = tot_post %>%
    group_by(Variable) %>%
    summarise(ISEMP_est = median(value),
              # se = sd(value),
              # cv = se / mean,
              ISEMP_lowCI = HPDinterval(as.mcmc(value), prob = 0.95)[,1],
              ISEMP_uppCI = HPDinterval(as.mcmc(value), prob = 0.95)[,2]) %>%
    ungroup()
  
  #-------------#
  # SCOBI model #
  #-------------#
  
  # Format data for use in SCOBI model
  scobi_dat = formatSCOBI_inputs(lgr_week$obs, lgr_truth) # uses the same output for ISEMP model run
  
  scobi_est = try(SCOBIv2(adultData = scobi_dat$fish_data, 
                          windowData = scobi_dat$window_count, 
                          Run = "chnkDemo",
                          RTYPE = "W", 
                          Primary = "GenStock", 
                          Secondary = NA, 
                          alph = 0.05, 
                          B = 1000, 
                          writeOutput = FALSE))
  
  if(class(scobi_est) == 'try-error') next
  
  scobi_summ = scobi_est$Rearing %>% as.data.frame() %>%
    mutate(RearType = rownames(scobi_est$Rearing)) %>%
    select(Variable = RearType, Estimates:U) %>% tbl_df() %>%
    mutate(Variable = revalue(Variable,
                              c('W' = 'Unique.Wild.Fish',
                                'H' = 'Unique.Hatch.Fish',
                                'HNC' = 'Unique.HNC.Fish'))) %>%
    rename(SCOBI_est = Estimates,
           SCOBI_lowCI = L,
           SCOBI_uppCI = U)
  
  # compare with "truth" from simulated data
  true_var = lgr_truth %>%
    summarise(All.Fish = length(id),
              Unique.Fish = n_distinct(id),
              Daytime.Fish = sum(Day.passage),
              Reascent.Fish = sum(Reascent),
              Night.Fish = sum(Night.passage)) %>%
    bind_cols(lgr_truth %>%
                filter(Origin == 'NOR') %>%
                summarise(Unique.Wild.Fish = n_distinct(id))) %>%
    bind_cols(lgr_truth %>%
                filter(Origin == 'HOR') %>%
                summarise(Unique.Hatch.Fish = n_distinct(id))) %>%
    bind_cols(lgr_truth %>%
                filter(Origin == 'HNC') %>%
                summarise(Unique.HNC.Fish = n_distinct(id))) %>%
    t() %>% as.data.frame()
  true_var = data.frame(Variable = row.names(true_var),
                        Truth = true_var$V1) %>% tbl_df()
  
  
  res[[i]] = true_var %>%
    inner_join(tot_summ) %>%
    left_join(scobi_summ) %>%
    mutate(ISEMP_inCI = ifelse(Truth >= ISEMP_lowCI & Truth <= ISEMP_uppCI, T, F),
           SCOBI_inCI = ifelse(Truth >= SCOBI_lowCI & Truth <= SCOBI_uppCI, T, F))
  
  sim_list[[i]] = my_sim
  obs_list[[i]] = lgr_week
  mod_list[[i]] = adult.pass.mod
  
  rm(my_sim, lgr_truth, lgr_week, jags.data, adult.pass.mod, true_var, tot_summ, tot_post, scobi_dat, scobi_est, scobi_summ)
  
  cat(paste('Finished simulation #', i, '\n'))
  i = i + 1
}
cat(paste('Took', round(c(proc.time() - tot_ptm)[3] / 60, 2), 'min to run all', n_sim, 'sims in total. \n'))

# save results
save(res, sim_list, obs_list, mod_list, file = 'SimulationFits/Sim_Baseline.rda')
