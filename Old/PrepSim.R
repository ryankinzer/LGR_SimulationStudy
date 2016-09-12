# Author: Kevin See
# Purpose: Prep simulated data for JAGS model to estimate total escapement over Lower Granite Dam
# Created: 8/11/2016
# Last Modified: 8/26/2016
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

source('SimFnc.R')
#-----------------------------------------------------------------
# generate simulated data
#-----------------------------------------------------------------
# night.passage.rate = 0.1
# fallback.rate = 0.1
# reascension.rate = 0.2

# determine trap rate on weekly basis
my_trap_rate = data.frame(Week = 1:52,
                          trap.rate = 0.05)

set.seed(4)
my_sim = SimulateLGRdata(trap.rate.df = my_trap_rate)

my_sim = SimulateLGRdata(N.lgr = 1e5,
                         n.pops = c(3,1),
                         run.mu.mu = c(140, 171),
                         run.mu.sd = c(16, 15),
                         run.sd.mu = c(16,25),
                         run.sd.sd = c(8,2),
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

# pull out some parameters from simulation
# these parameters are used elsewhere as well as the simulation function
night.passage.rate = my_sim$parameters$night.passage.rate
window.rate = my_sim$parameters$window.rate


# add date for arbitary year, and beginning date of each week
lgr_truth %<>%
  mutate(Date = ymd('20150101') + days(Day)) %>%
  group_by(Week) %>%
  mutate(Start_Date = min(Date)) %>%
  ungroup()

# True population escapement and proportion at Lower Granite Dam (includes all origins)
N_pop = lgr_truth %>%
  group_by(Population) %>%
  summarise(n_fish = length(Population)) %>%
  ungroup() %>%
  mutate(pop_p = n_fish / nrow(lgr_truth))

#-----------------------------------------------------------------
# mark-recapture estimate of trap rate
trap_mr = lgr_truth %>%
  mutate(id = 1:n()) %>%
  filter(Marked == 1) %>%
  mutate(R = ifelse(Trap == 1 & Ladder == 1, 1, 0)) %>%
  filter(!(Trap == 0 & Ladder == 0)) %>%
  select(id, 
         Week,
         M = Trap,
         C = Ladder,
         R)

mr_n_fish = lgr_truth %>%
  group_by(Week) %>%
  summarise(N = sum(Marked)) %>%
  left_join(trap_mr %>%
              group_by(Week) %>%
              summarise_each(funs(sum), M:R))

# put together list of capture histories
caphist_list = dlply(trap_mr, .(Week), capHistConvert, cols2use = c('M', 'C', 'R'), in.type = 'individual', out.type = 'frequency')

# estimate trap rate based on mark-recapture
trap_rate_mr = ldply(caphist_list, function(x) {
  mod = try(closedp.t(x, dfreq=T), silent = T)
  if(class(mod)[1] == 'try-error') {
    return(data.frame(N_hat = NA, N_se = NA, p = NA, p_se = NA))
  }
  p_mean = inv.logit(coef(mod$glm[['Mt']])[2:3])
  p_se = deltamethod(list(~ 1 / (1 + exp(-x1)), ~ 1 / (1 + exp(-x2))), mean = coef(mod$glm[['Mt']])[2:3], cov = vcov(mod$glm[['Mt']])[2:3,2:3])
  return(data.frame(N_hat = mod$results['Mt', 'abundance'],
                    N_se = mod$results['Mt', 'stderr'],
                    p = p_mean[1], 
                    p_se = p_se[1]))
}, .id = 'Week') %>% tbl_df()


# how is precision of trap rate affected by M?
left_join(mr_n_fish, 
          trap_rate_mr) %>%
  mutate(p_cv = p_se / p) %>%
  filter(R > 0) %>%
  ggplot(aes(x = M,
             y = p_cv)) +
  geom_point() +
  geom_smooth()

# how well do estimates of trap rate match reality?
trap_rate_mr %>%
  left_join(mr_n_fish) %>%
  left_join(my_trap_rate) %>%
  ggplot(aes(x = Week,
             y = p)) +
  geom_errorbar(aes(ymin = p + qnorm(0.025) * p_se,
                    ymax = p + qnorm(0.975) * p_se)) +
  geom_point(aes(size = M)) +
  geom_line(aes(y = trap.rate),
            color = 'red') +
  geom_line(aes(y = M / N),
            color = 'blue') +
  coord_cartesian(ylim = c(0, 0.2)) +
  labs(title = 'Trap Rate')




#-----------------------------------------------------------------
# summarise at weekly scales, including estimate of trap rate
#-----------------------------------------------------------------
lgr_week = lgr_truth %>%
  group_by(Week, Start_Date) %>%
  summarise(N_true = n_distinct(id)) %>%
  ungroup() %>%
  left_join(lgr_truth %>%
              group_by(Week, Origin) %>%
              summarise(n_fish = n_distinct(id)) %>%
              spread(Origin, n_fish, fill = 0)) %>%
  left_join(lgr_truth %>%
              group_by(Week) %>%
              summarise(Day.passage = sum(Night.passage == 0))) %>%
  left_join(lgr_truth %>%
              group_by(Week) %>%
              summarise_each(funs(sum), Window.passage, Trap, Ladder, Fallback, Reascent, Night.passage) %>%
              # assume that reascensions happen within the same day (or at least week), how many fish might be counted at window?
              mutate(Window.fish = Window.passage + (Reascent * rbinom(nrow(.), 1, (1 - night.passage.rate) * window.rate)),
                     # this introduces some error in window counts
                     win_cnt = rnegbin(nrow(.), mu = Window.fish / window.rate, theta = 2))) %>%
                     # win_cnt = round(rnorm(nrow(.), mean = Window.fish / window.rate, sd = 0.05 * (Window.fish / window.rate))),
                     # win_cnt = ifelse(win_cnt < 0, 0, win_cnt))) %>%
                     
                     # this assumes no error in window counts
                     # win_cnt = round(Window.fish / window.rate))) %>%
  left_join(lgr_truth %>%
              filter(Marked == 1,
                     Ladder == 1) %>%
              group_by(Week) %>%
              summarise_each(funs(sum), Reascent, Night.passage) %>%
              rename(mark_reasc = Reascent,
                     mark_night = Night.passage)) %>%
  # # if no hatchery or HNC fish, skip this section
  # left_join(lgr_truth %>%
  #             filter(Trap == 1) %>%
  #             group_by(Week, Origin) %>%
  #             summarise(n_org = length(Origin)) %>%
  #             spread(Origin, n_org, fill = 0) %>%
  #             rename(wild_fish = NOR,
  #                    hatch_fish = HOR,
  #                    HNC_fish = HNC)) %>%
  # mutate_each(funs(ifelse(is.na(.), 0, .)), mark_reasc:HNC_fish) %>%
  left_join(lgr_truth %>%
              filter(Trap == 1) %>%
              group_by(Week, Origin) %>%
              summarise(n_org = length(Origin)) %>%
              spread(Origin, n_org, fill = 0) %>%
              rename(wild_fish = NOR) %>%
              select(Week, wild_fish)) %>%
  mutate_each(funs(ifelse(is.na(.), 0, .)), mark_reasc:wild_fish) %>%
  left_join(trap_rate_mr %>%
              select(Week,
                     trap_rate = p,
                     trap_rate_se = p_se) %>%
              # set up parameters describing trap rate as a beta distribution
              mutate(trap_alpha = ((1 - trap_rate) / trap_rate_se^2 - 1 / trap_rate) * trap_rate^2,
                     trap_beta = trap_alpha * (1 / trap_rate - 1))) %>%
  mutate(trap_alpha = ifelse(is.na(trap_rate_se), 1, trap_alpha),
         trap_beta = ifelse(is.na(trap_rate_se), 1, trap_beta))

# filter for spring/summer Chinook dates
lgr_week %<>%
  filter(Start_Date >= ymd('20150301'),
         Start_Date <= ymd('20150817')) %>%
  mutate(Week = Week - min(Week) + 1)

#-----------------------------------------------------------------
# pull data together for JAGS
jags.data = list('TotLadderWeeks' = nrow(lgr_week),
                 'Y.window' = lgr_week %>% select(win_cnt) %>% as.matrix() %>% as.vector(),
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


#-----------------------------------------------------------------
# set mcmc parameters
#-----------------------------------------------------------------
# number of total samples
mcmc.chainLength = 12000 #40000

# number of burn-in samples
mcmc.burn = 2000 #10000

# thinning interval
mcmc.thin = 20 #30

# number of MCMC chains
mcmc.chains = 4

# how many poseterior samples will we have?
data.frame(Each.Chain = (mcmc.chainLength - mcmc.burn) / mcmc.thin, All.Chains = (mcmc.chainLength - mcmc.burn) / mcmc.thin * mcmc.chains)


#-----------------------------------------------------------------
# run model using JAGS
#-----------------------------------------------------------------

model.loc = 'LGR_TotalEscape_NB1.txt'

jags.params = c('X.log.all', 'X.all', 'X.day', 'X.night', 'X.reasc', 'X.all.wild.pbt', 'X.new.wild.pbt', 'X.reasc.wild.pbt', 'X.night.wild.pbt', 'X.tot.all', 'X.tot.day', 'X.tot.night', 'X.tot.reasc', 'X.tot.all.wild.pbt', 'X.tot.new.wild.pbt', 'X.tot.night.wild.pbt', 'X.tot.reasc.wild.pbt', 'prop.tagged.pbt', 'X.sigma', 'true.prop', 'win.prop.avg', 'win.prop.true', 'win.prop.sigma', 'hist.prop', 'reasc.avg', 'reasc.true', 'reasc.sigma', 'acf', 'wnc.avg', 'wnc.true', 'wnc.sigma', 'trap.rate.true', 'r', 'k', 'p', 'alpha')

# set initial values
jags.inits = function(){
  return(list('X.log.all' = with(jags.data, log(Y.window + 1)),
              'trap.rate.true' = with(jags.data, trap.rate + 0.0001)))
}

start_time = proc.time()
adult.pass.mod = try(jags(data = jags.data, 
                          inits = jags.inits, 
                          parameters.to.save = jags.params, 
                          model.file = model.loc, 
                          n.chains = mcmc.chains, 
                          n.burnin = mcmc.burn, 
                          n.thin = mcmc.thin, 
                          n.iter = mcmc.chainLength, 
                          parallel = T,
                          DIC = FALSE))
end_time = proc.time()
end_time - start_time # returns the CPU time used

#-----------------------------------------------------------------
# model diagnostics
#-----------------------------------------------------------------
library(ggmcmc)
my_ggs = ggs(adult.pass.mod$samples, family = c('X.tot'))
my_ggs = ggs(adult.pass.mod$samples, family = c('sigma'))
my_ggs = ggs(adult.pass.mod$samples, family = c('acf'))
my_ggs = ggs(adult.pass.mod$samples, family = c('avg'))
my_ggs = ggs(adult.pass.mod$samples, family = c('^r$'))
my_ggs = ggs(adult.pass.mod$samples, family = c('^r\\['))
my_ggs = ggs(adult.pass.mod$samples, family = c('k'))
my_ggs = ggs(adult.pass.mod$samples, family = c('^p'))
my_ggs = ggs(adult.pass.mod$samples, family = c('alpha'))

ggs_density(my_ggs) +
  facet_wrap(~ Parameter, scales = 'free')
ggs_traceplot(my_ggs) +
  facet_wrap(~ Parameter, scales = 'free')
ggs_Rhat(my_ggs)
ggs_geweke(my_ggs)
ggs_caterpillar(my_ggs)
ggs_autocorrelation(my_ggs)
ggs_crosscorrelation(my_ggs)
ggs_running(my_ggs)
ggs_compare_partial(my_ggs)

#-----------------------------------------------------------------
# pull out results from posteriors
#-----------------------------------------------------------------
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

other_post = ldply(list('Tot.Rand.Walk.SD' = X.sigma,
                        'Avg.Day' = win.prop.avg,
                        'SD.Day' = win.prop.sigma,
                        'Day.Rate.AutoCorr' = acf[,1] %>% as.matrix(),
                        'Avg.ReAsc' = reasc.avg,
                        'SD.ReAsc' = reasc.sigma,
                        'ReAsc.AutoCorr' = acf[,2] %>% as.matrix(),
                        'Avg.Prop.Wild.PBT' = wnc.avg %>% as.matrix(),
                        'Prop.Wild.PBT.AutoCorr' = acf[,3] %>% as.matrix()),
                        # 'Over.Dispersion.p' = p),
                        # 'Over.Dispersion' = r),
                   .id = 'Variable') %>% tbl_df() %>%
  gather(iteration, value, -Variable) %>%
  mutate(iteration = gsub('^V', '', iteration),
         iteration = as.integer(iteration))

week_post = ldply(list('New.Wild.Fish.PBT' = X.new.wild.pbt, 
                       'All.Wild.Fish.PBT' = X.all.wild.pbt, 
                       'All.Fish' = X.all, 
                       'Daytime.Fish' = X.day,
                       'Night.Wild.Fish' = X.night.wild.pbt,
                       'Wild.Reascents' = X.reasc.wild.pbt, 
                       'Trap.Rate' = trap.rate.true,
                       'Day.Time.Rate' = true.prop, 
                       'Night.Time.Rate' = 1 - true.prop, 
                       'Re-Ascension.Rate' = reasc.true,
                       'Wild.PBT.Rate' = wnc.true), 
                  .id='Variable') %>% tbl_df() %>%
  gather(week_num, value, -Variable)
detach(adult.pass.mod$sims.list)

tot_summ = tot_post %>%
  group_by(Variable) %>%
  summarise(mean = mean(value),
            median = median(value),
            se = sd(value),
            cv = se / mean,
            low_ci = HPDinterval(as.mcmc(value), prob = 0.95)[,1],
            upp_ci = HPDinterval(as.mcmc(value), prob = 0.95)[,2]) %>%
  ungroup()

other_summ = other_post %>%
  group_by(Variable) %>%
  summarise(mean = mean(value),
            median = median(value),
            se = sd(value),
            cv = se / mean,
            low_ci = HPDinterval(as.mcmc(value), prob = 0.95)[,1],
            upp_ci = HPDinterval(as.mcmc(value), prob = 0.95)[,2]) %>%
  ungroup()


week_summ = week_post %>%
  group_by(Variable, week_num) %>%
  summarise(mean = mean(value),
            median = median(value),
            se = sd(value),
            cv = se / mean,
            low_ci = HPDinterval(as.mcmc(value), prob = 0.95)[,1],
            upp_ci = HPDinterval(as.mcmc(value), prob = 0.95)[,2]) %>%
  ungroup() %>%
  mutate(week_num = as.integer(as.character(week_num))) %>%
  select(week_num, everything())

plot_df = left_join(week_summ,
                    lgr_week %>%
                      mutate(week_num = Week)) %>%
  mutate(trap_est = Trap / trap_rate)

#-----------------------------------------------------------------
# make some plots and summary tables
#-----------------------------------------------------------------
true_var = lgr_truth %>%
  summarise(All.Fish = length(id) + sum(Reascent),
            Daytime.Fish = sum(Night.passage == 0),
            Reascent.Fish = sum(Reascent)) %>%
  bind_cols(lgr_truth %>%
              filter(Origin == 'NOR') %>%
              summarise(All.Wild.Fish.PBT = n_distinct(id) + sum(Reascent),
                        Unique.Wild.Fish.PBT = n_distinct(id),
                        Night.Wild.Fish.PBT = sum(Night.passage),
                        Wild.Reascents.PBT = sum(Reascent))) %>%
  t()
true_var = data.frame(Variable = row.names(true_var),
                      Truth = true_var) %>% tbl_df()

true_var %>%
  left_join(tot_summ) %>%
  mutate(inCI = ifelse(Truth >= low_ci & Truth <= upp_ci, T, F))


plot_df %>%
  filter(Variable == 'Daytime.Fish') %>%
  mutate(inCI = ifelse(low_ci <= Day.passage & upp_ci >= Day.passage, T, F),
         bias = median - Day.passage,
         rel_err = ifelse(Day.passage > 0, abs(bias) / Day.passage, NA)) %>%
  summarize(perc_inCI = round(sum(inCI) / n(), 2),
            avg_bias = round(mean(bias)),
            wt_avg_bias = round(weighted.mean(bias, w = Day.passage)),
            avg_abs_bias = round(mean(abs(bias))),
            wt_avg_abs_bias = round(weighted.mean(abs(bias), w = Day.passage)),
            avg_rel_err = round(mean(rel_err, na.rm = T), 3),
            wt_avg_rel_err = round(weighted.mean(rel_err, w = Day.passage, na.rm = T), 3))


# total escapement by week, with window, trap and model estimates
tot_p = plot_df %>%
  filter(Variable == 'All.Fish') %>%
  filter(trap_rate_se / trap_rate < 1) %>%
  ggplot(aes(x = Start_Date, y = median)) +
  geom_ribbon(aes(ymin = low_ci, ymax = upp_ci), fill = 'lightgray') +
  geom_line(aes(y = win_cnt, color = 'Window')) +
  geom_point(aes(y = win_cnt, color = 'Window')) +
  geom_line(aes(y = trap_est, color = 'Trap')) +
  geom_point(aes(y = trap_est, 
                 color = 'Trap')) +
  geom_point(aes(y = N_true,
                 color = 'Truth')) +
  geom_line(aes(y = N_true,
                color = 'Truth')) +
  geom_line(aes(color = 'Model')) +
  geom_point(aes(color = 'Model')) +
  scale_color_brewer(palette = 'Set1') +
  # scale_color_manual(values = c('Model' = 'green', 
  #                               'Window' = 'blue', 
  #                               'Trap' = 'red',
  #                               'Truth' = 'black')) +
  # scale_y_log10() +
  labs(x = 'Date', y = 'Estimate', color = 'Source', title = 'All Fish')
tot_p

# day fish
day_p = plot_df %>%
  filter(Variable == 'Daytime.Fish') %>%
  ggplot(aes(x = Start_Date, y = median)) +
  geom_ribbon(aes(ymin = low_ci, ymax = upp_ci), fill = 'lightgray') +
  geom_line(aes(y = win_cnt, color = 'Window')) +
  geom_point(aes(y = win_cnt, color = 'Window')) +
  geom_point(aes(y = Day.passage,
                 color = 'Truth')) +
  geom_line(aes(y = Day.passage,
                color = 'Truth')) +
  geom_line(aes(color = 'Model')) +
  geom_point(aes(color = 'Model')) +
  scale_color_brewer(palette = 'Set1') +
  # scale_color_manual(values = c('Model' = 'green', 
  #                               'Window' = 'blue', 
  #                               'Trap' = 'red',
  #                               'Truth' = 'black')) +
  # scale_y_log10() +
  labs(x = 'Date', y = 'Estimate', color = 'Source', title = 'All Daytime Fish')
day_p



# wild escapement by week, with window, trap and model estimates
wild_p = plot_df %>%
  filter(Variable == 'All.Wild.Fish.PBT') %>%
  ggplot(aes(x = Start_Date, y = median)) +
  geom_ribbon(aes(ymin = low_ci, ymax = upp_ci), fill = 'lightgray') +
  geom_line(aes(y = win_cnt * (wild_fish / Trap), color = 'Window')) +
  geom_point(aes(y = win_cnt * (wild_fish / Trap), color = 'Window')) +
  geom_line(aes(y = wild_fish / trap_rate, color = 'Trap')) +
  geom_point(aes(y = wild_fish / trap_rate, 
                 color = 'Trap')) +
  geom_point(aes(y = NOR,
                 color = 'Truth')) +
  geom_line(aes(y = NOR,
                color = 'Truth')) +
  geom_line(aes(color = 'Model')) +
  geom_point(aes(color = 'Model')) +
  scale_color_brewer(palette = 'Set1') +
  # scale_y_log10() +
  labs(x = 'Date', y = 'Estimate', color = 'Source', title = 'All Wild Fish')
wild_p

uni_wild_p = plot_df %>%
  filter(Variable == 'New.Wild.Fish.PBT') %>%
  ggplot(aes(x = Start_Date, y = median)) +
  geom_ribbon(aes(ymin = low_ci, ymax = upp_ci), fill = 'lightgray') +
  geom_line(aes(y = win_cnt * (wild_fish / Trap) * (1 - mark_reasc / Ladder), 
                color = 'Window')) +
  geom_point(aes(y = win_cnt * (wild_fish / Trap) * (1 - mark_reasc / Ladder), 
                 color = 'Window')) +
  geom_line(aes(y = wild_fish / trap_rate * (1 - mark_reasc / Ladder),
                color = 'Trap')) +
  geom_point(aes(y = wild_fish / trap_rate* (1 - mark_reasc / Ladder),
                 color = 'Trap')) +
  geom_point(aes(y = NOR,
                 color = 'Truth')) +
  geom_line(aes(y = NOR,
                color = 'Truth')) +
  geom_line(aes(color = 'Model')) +
  geom_point(aes(color = 'Model')) +
  scale_color_brewer(palette = 'Set1') +
  # scale_y_log10() +
  labs(x = 'Date', y = 'Estimate', color = 'Source', title = 'Unique Wild Fish')
uni_wild_p


# look at rates of wild fish
plot_df %>%
  filter(Variable == 'Wild.PBT.Rate') %>%
  ggplot(aes(x = Start_Date, 
             y = median)) +
  geom_ribbon(aes(ymin = low_ci, 
                  ymax = upp_ci, 
                  group = 'Model', 
                  fill = 'Model'), 
              color = NA, 
              alpha = 0.4) +
  geom_hline(data = filter(other_summ,
                           Variable == 'Avg.Prop.Wild.PBT'),
             aes(yintercept =  median),
             lty = 2,
             col = 'darkgreen') +
  geom_line(aes(group = 'Model',
                color = 'Model')) +
  # geom_point(aes(group = Variable, size = win_cnt), position = position_dodge(width = 0.90)) +
  geom_point(aes(color = 'Model', 
                 size = Trap), 
             position = position_dodge(width = 0.90)) +
  geom_line(aes(y = NOR / N_true,
                color = 'Truth')) +
  scale_color_manual(values = c('Truth' = 'blue4', 'Model' = 'red')) +
  scale_fill_manual(values = c('Truth' = 'lightskyblue', 'Model' = 'mistyrose2')) +
  labs(x = 'Week', 
       y = 'Rate', 
       # size = 'Window Count',
       size = 'Fish in Trap',
       title = 'Rate of Wild Fish')

# look at night time passage rate
plot_df %>%
  filter(Variable == 'Night.Time.Rate') %>%
  ggplot(aes(x = Start_Date)) +
  geom_ribbon(aes(ymin = low_ci,
                  ymax = upp_ci,
                  fill = 'Estimate'),
              alpha = 0.4) +
  geom_ribbon(data = data.frame(Start_Date = sort(unique(plot_df$Start_Date)),
                                other_summ %>%
                                  filter(Variable == 'Avg.Day')),
              aes(ymin = 1 - upp_ci,
                  ymax = 1 - low_ci,
                  fill = 'Average'),
              alpha = 0.4) +
  geom_hline(data = filter(other_summ,
                           Variable == 'Avg.Day'),
             aes(yintercept =  1 - median,
                 color = 'Average'),
             lty = 2) +
  geom_line(aes(y = median,
                color = 'Estimate')) +
  geom_point(aes(y = median,
                 color = 'Estimate',
                 size = mark_night)) +
  geom_line(aes(y = Night.passage / N_true,
                color = 'Truth')) +
  geom_point(aes(y = Night.passage / N_true,
                 color = 'Truth')) +
  scale_color_manual(values = c('Truth' = 'blue',
                                'Estimate' = 'red',
                                'Average' = 'darkgreen'),
                     name = 'Source') +
  scale_fill_manual(values = c('Estimate' = 'mistyrose2',
                               'Truth' = 'lightskyblue',
                               'Average' = 'darkseagreen1'),
                    name = 'Source') +
  labs(title = 'Night-Time Passage',
       size = 'Night Marked')

# look at reascension rate
plot_df %>%
  filter(Variable == 'Re-Ascension.Rate') %>%
  ggplot(aes(x = Start_Date)) +
  geom_ribbon(aes(ymin = low_ci,
                  ymax = upp_ci,
                  fill = 'Estimate'),
              alpha = 0.4) +
  geom_ribbon(data = data.frame(Start_Date = sort(unique(plot_df$Start_Date)),
                                other_summ %>%
                                  filter(Variable == 'Avg.ReAsc')),
              aes(ymin = low_ci,
                  ymax = upp_ci,
                  fill = 'Average'),
              alpha = 0.4) +
  geom_hline(data = filter(other_summ,
                           Variable == 'Avg.ReAsc'),
             aes(yintercept =  median,
                 color = 'Average'),
             lty = 2) +
  geom_line(aes(y = median,
                color = 'Estimate')) +
  geom_point(aes(y = median,
                 color = 'Estimate',
                 size = mark_reasc)) +
  geom_line(aes(y = Reascent / N_true,
                color = 'Truth')) +
  geom_point(aes(y = Reascent / N_true,
                 color = 'Truth')) +
  scale_color_manual(values = c('Truth' = 'blue',
                                'Estimate' = 'red',
                                'Average' = 'darkgreen'),
                     name = 'Source') +
  scale_fill_manual(values = c('Estimate' = 'mistyrose2',
                               'Truth' = 'lightskyblue',
                               'Average' = 'darkseagreen1'),
                    name = 'Source') +
  labs(title = 'Re-ascension',
       size = 'Reascending\nMarked')


# look at estimates of trap rate
plot_df %>%
  filter(Variable == 'Trap.Rate') %>%
  ggplot(aes(x = Start_Date)) +
  geom_ribbon(aes(ymin = trap_rate + qnorm(0.025) * trap_rate_se,
                  ymax = trap_rate + qnorm(0.975) * trap_rate_se,
                  fill = 'Prior'),
              alpha = 0.2) +
  geom_line(aes(y = trap_rate,
                color = 'Prior')) +
  geom_ribbon(aes(ymin = low_ci,
                  ymax = upp_ci,
                  fill = 'Estimate'),
              alpha = 0.4) +
  geom_line(aes(y = median,
                color = 'Estimate')) +
  geom_point(aes(y = median,
                 color = 'Estimate',
                 size = Trap)) +
  # geom_line(aes(y = Trap / N_true,
  #               color = 'Truth')) +
  geom_hline(yintercept = unique(my_trap_rate$trap.rate),
             linetype = 2,
             aes(color = 'Truth')) +
  # geom_hline(yintercept = 0.05,
  #            linetype = 2) +
  scale_color_manual(values = c('Truth' = 'blue',
                                'Estimate' = 'red',
                                'Prior' = 'darkgreen'),
                     name = 'Source') +
  scale_fill_manual(values = c('Estimate' = 'mistyrose2',
                               'Truth' = NULL,
                               'Prior' = 'seagreen'),
                    name = 'Source') +
  # coord_cartesian(ylim = c(0, 0.1)) +
  labs(title = 'Trap Rate',
       size = 'Fish in Trap')


plot_df %>%
  filter(Variable == 'Trap.Rate') %>%
  ggplot(aes(x = trap_rate,
             y = median)) +
  geom_point(aes(size = Trap)) +
  geom_abline(color = 'red',
              linetype = 2) +
  geom_smooth(method = lm,
              formula = y ~ -1 + x)

week_post %>%
  filter(Variable == 'Trap.Rate') %>%
  mutate(Week = as.integer(as.character(week_num))) %>%
  left_join(lgr_week) %>%
  mutate(trap_est = Trap / value) %>%
  ggplot(aes(x = N_true,
             y = trap_est)) + 
  geom_point() +
  geom_abline(color = 'red',
              linetype = 2)

lgr_week %>%
  mutate(trap_est = Trap / trap_rate) %>%
  filter(trap_est < 1e12) %>%
  select(Week:N_true, wind_est = win_cnt, trap_est) %>%
  ggplot(aes(x = wind_est, 
             y = trap_est)) +
  geom_point() +
  geom_abline(color = 'red',
              linetype = 2)

lgr_week %>%
  ggplot(aes(x = Window.fish + Reascent,
             y = win_cnt)) +
  geom_point() +
  geom_abline(color = 'red',
              linetype = 2)