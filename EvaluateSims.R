# Author: Kevin See
# Purpose: Make tables and plots based on simulated LGR data and ISEMP and TAC model estimates
# Created: 9/2/2016
# Last Modified: 9/20/2016
# Notes: 

#-----------------------------------------------------------------
# library(MCMCpack)
# library(boot)
# library(msm)
library(lubridate)
library(magrittr)
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(jagsUI)

theme_set(theme_bw())

#-----------------------------------------------------------------
# load results
result_nms = list.files('SimulationFits')
result_nms = result_nms[!grepl('other', result_nms)]
scenario_nms = gsub('^Sim_', '', result_nms)
scenario_nms = gsub('\\.rda$', '', scenario_nms)
length(scenario_nms)

res_list = vector('list', length(scenario_nms))
names(res_list) = scenario_nms
for(i in 1:length(res_list)) {
  load(paste0('SimulationFits/', result_nms[i]))
  res_list[[i]] = res
  rm(res)
}

# pull all results into a data.frame
res_df = ldply(res_list,
               .id = 'Scenario',
               .fun = function(x) {
                 ldply(x,
                       .id = 'sim')
               }) %>% tbl_df()
rm(res_list)

#-----------------------------------------------------------------
# summarise and plot results
#-----------------------------------------------------------------
res_df %>%
  filter(grepl('Unique', Variable)) %>%
  group_by(Scenario, Variable) %>%
  summarise(median_Truth = median(Truth),
            ISEMP_cover95 = sum(ISEMP_inCI) / n(),
            SCOBI_cover95 = sum(SCOBI_inCI) / n(),
            ISEMP_CIwidth = median(ISEMP_uppCI - ISEMP_lowCI),
            SCOBI_CIwidth = median(SCOBI_uppCI - SCOBI_lowCI),
            ISEMP_cv = mean((ISEMP_uppCI - ISEMP_lowCI) / (2 * qnorm(0.975)) / ISEMP_est),
            SCOBI_cv = mean((SCOBI_uppCI - SCOBI_lowCI) / (2 * qnorm(0.975)) / SCOBI_est))

qplot(Truth, SCOBI_est, data = res_df, color = Variable, log = 'xy') + 
  geom_abline() +
  facet_wrap(~ Scenario)

qplot(Truth, ISEMP_est, data = res_df, color = Variable, log = 'xy') + 
  geom_abline()

qplot(SCOBI_est, ISEMP_est, data = res_df, color = Variable, log = 'xy') + 
  geom_abline()


res_df %>%
  select(Scenario:ISEMP_est, SCOBI_est) %>%
  gather(source, est, -(Scenario:Truth)) %>%
  mutate(source = gsub('_est$', '', source)) %>%
  mutate(bias = est - Truth,
         rel_bias = bias / Truth) %>%
  ggplot(aes(x = Variable,
             fill = Variable,
             y = bias)) +
  geom_boxplot() +
  scale_fill_brewer(palette = 'Set1') +
  geom_hline(yintercept = 0,
             linetype = 2) +
  facet_grid(Scenario ~ source, scales = 'free') +
  labs(y = 'Bias') +
  theme(axis.text.x = element_blank()) 

#-----------------------------------------------------------------
# is poor coverage and bias of SCOBI due to night and reascension fish?
#-----------------------------------------------------------------
bias_df = res_df %>%
  filter(Scenario == 'NightReasc',
         Variable %in% c('Night.Fish', 'Reascent.Fish')) %>%
  select(Scenario, sim, Variable, Truth) %>%
  spread(Variable, Truth) %>%
  mutate(pot_bias = Reascent.Fish - Night.Fish) %>%
  left_join(res_df %>%
              filter(Scenario == 'NightReasc',
                     grepl('Unique', Variable)) %>%
              group_by(Scenario, sim) %>%
              summarise(Uni.Fish.Truth = sum(Truth),
                        Tot.Fish.SCOBI = sum(SCOBI_est)) %>%
              mutate(bias = Tot.Fish.SCOBI - Uni.Fish.Truth))

qplot(pot_bias, bias, data = bias_df,
      xlab = 'Reascent - Night',
      ylab = 'Realized Bias') +
  geom_abline(color = 'red',
              linetype = 2) +
  geom_smooth(method = lm)

mod = lm(bias ~ pot_bias, data = bias_df)
summary(mod)
op = par(mfrow = c(2,2))
plot(mod)
par(op)


#------------------------------------------------------
# Re-format res_df into long format of point ests.,
# CI's and CIin for both Models : ISEMP, SCOBI
#------------------------------------------------------
res_long_df <- res_df %>%
  select(Scenario, sim) %>%
  distinct() %>%
  ungroup() %>%
  left_join(res_df) %>%
  gather(key, point, -Scenario, -sim, -Truth, -Variable) %>%
  separate(key, into = c("Model","est"), sep = "_") %>%
  spread(est,point) %>%
  mutate(cv = (uppCI - lowCI) / (2 * qnorm(0.975)) / est,
         bias = est - Truth,
         rel_bias = bias / Truth,
         inCI = ifelse(inCI == 1, T, ifelse(inCI == 0, F, NA)))

res_long_df %>%
  filter(grepl('Unique', Variable)) %>%
  group_by(Scenario, Variable, Model) %>%
  summarise(median_Truth = median(Truth),
            cover95 = sum(inCI) / n(),
            CIwidth = median(uppCI - lowCI),
            median_cv = median(cv),
            rel_bias = median(rel_bias))

res_long_df %>%
  ggplot(aes(x = Model,
             fill = Variable,
             y = rel_bias)) +
  geom_boxplot() +
  scale_fill_brewer(palette = 'Set1') +
  geom_hline(yintercept = 0,
             linetype = 2) +
  facet_wrap(~ Scenario, scales = 'free') +
  labs(y = 'Bias') +
  coord_cartesian(ylim = c(-0.5, 0.5))




pd = .4

res_long_df %>%
  group_by(Scenario) %>%
  sample_n(25) %>%
  filter(grepl('^Unique', Variable)) %>%
  mutate(low_bias = lowCI - Truth,
         upp_bias = uppCI - Truth) %>%
  ggplot(aes(x = sim,
             y = bias,
             color = Model,
             shape = inCI)) +
  scale_shape_manual(values = c('TRUE' = 19,
                                'FALSE' = 1)) +
  geom_errorbar(aes(ymin = low_bias,
                    ymax = upp_bias),
                position = position_dodge(width = pd)) +
  geom_point(position = position_dodge(width = pd)) +
  geom_hline(yintercept = 0,
             linetype = 2) +
  facet_wrap(~Scenario + Variable, scales = 'free') +
  scale_color_brewer(palette = 'Set1')


res_long_df %>%
  group_by(Scenario) %>%
  sample_n(25) %>%
  mutate(low_bias = lowCI - Truth,
         upp_bias = uppCI - Truth) %>%
  ggplot(aes(x = sim,
             y = bias,
             color = Model,
             shape = inCI)) +
  scale_shape_manual(values = c('TRUE' = 19,
                                'FALSE' = 1)) +
  geom_errorbar(aes(ymin = low_bias,
                    ymax = upp_bias),
                position = position_dodge(width = pd)) +
  geom_point(position = position_dodge(width = pd)) +
  geom_hline(yintercept = 0,
             linetype = 2) +
  facet_grid(Variable ~ Scenario, scales = 'free') +
  theme(axis.text.x = element_blank()) +
  scale_color_brewer(palette = 'Set1')



# examine specific model results
i = 8
adult.pass.mod = mod_list[[i]]
lgr_week_obs = obs_list[[i]]$obs
lgr_truth = sim_list[[i]]$sim
sim_params = sim_list[[i]]$parameters
true_var = filter(res_df, sim == i)

#-----------------------------------------------------------------
# model diagnostics
#-----------------------------------------------------------------
library(ggmcmc)
my_ggs = ggs(adult.pass.mod$samples, family = c('X.tot'))
my_ggs = ggs(adult.pass.mod$samples, family = c('sigma'))
my_ggs = ggs(adult.pass.mod$samples, family = c('acf'))
my_ggs = ggs(adult.pass.mod$samples, family = c('avg'))
my_ggs = ggs(adult.pass.mod$samples, family = c('^r$'))
my_ggs = ggs(adult.pass.mod$samples, family = c('k'))
my_ggs = ggs(adult.pass.mod$samples, family = c('theta'))
my_ggs = ggs(adult.pass.mod$samples, family = c('omega'))

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
                        'Prop.Wild.PBT.AutoCorr' = acf[,3] %>% as.matrix(),
                        'Over.Dispersion' = r),
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
                    lgr_week_obs %>%
                      mutate(week_num = Week) %>%
                      select(-(N_tot:N_uniq))) %>%
  mutate(trap_est = Trap / trap_rate) %>%
  left_join(lgr_truth %>%
              group_by(Start_Date) %>%
              summarise(N_tot = length(id),
                        N_uniq = n_distinct(id))) %>%
  left_join(lgr_truth %>%
              group_by(Origin, Start_Date) %>%
              summarise(N_tot = length(id)) %>%
              ungroup() %>%
              mutate(Origin = revalue(Origin,
                                      c('NOR' = 'NOR_tot',
                                        'HOR' = 'HOR_tot',
                                        'HNC' = 'HNC_tot'))) %>%
              spread(Origin, N_tot, fill = 0)) %>%
  left_join(lgr_truth %>%
              group_by(Origin, Start_Date) %>%
              summarise(N_uniq = n_distinct(id)) %>%
              ungroup() %>%
              mutate(Origin = revalue(Origin,
                                      c('NOR' = 'NOR_uniq',
                                        'HOR' = 'HOR_uniq',
                                        'HNC' = 'HNC_uniq'))) %>%
              spread(Origin, N_uniq, fill = 0))

#-----------------------------------------------------------------
# make some plots and summary tables
#-----------------------------------------------------------------
# total escapement by week, with window, trap and model estimates
tot_p = plot_df %>%
  filter(Variable == 'All.Fish') %>%
  filter(trap_rate_se / trap_rate < 1) %>%
  ggplot(aes(x = Start_Date, y = median)) +
  geom_ribbon(aes(ymin = low_ci, 
                  ymax = upp_ci), 
              fill = 'mistyrose2') +
  geom_line(aes(y = win_cnt, 
                color = 'Window')) +
  geom_point(aes(y = win_cnt, 
                 color = 'Window')) +
  geom_line(aes(y = trap_est, 
                color = 'Trap')) +
  geom_point(aes(y = trap_est, 
                 color = 'Trap')) +
  geom_point(aes(y = N_tot,
                 color = 'Truth')) +
  geom_line(aes(y = N_tot,
                color = 'Truth')) +
  geom_line(aes(color = 'Model')) +
  geom_point(aes(color = 'Model')) +
  # scale_color_brewer(palette = 'Set1') +
  scale_color_manual(values = c('Model' = 'red',
                                'Window' = 'darkgreen',
                                'Trap' = 'blue',
                                'Truth' = 'black')) +
  # scale_y_log10() +
  labs(x = 'Date', y = 'Estimate', color = 'Source', title = 'All Fish')
tot_p


# day fish
day_p = plot_df %>%
  filter(Variable == 'Daytime.Fish') %>%
  ggplot(aes(x = Start_Date, y = median)) +
  geom_ribbon(aes(ymin = low_ci, 
                  ymax = upp_ci), 
              fill = 'mistyrose2') +
  geom_line(aes(y = win_cnt, color = 'Window')) +
  geom_point(aes(y = win_cnt, color = 'Window')) +
  geom_point(aes(y = Day.passage,
                 color = 'Truth')) +
  geom_line(aes(y = Day.passage,
                color = 'Truth')) +
  geom_line(aes(color = 'Model')) +
  geom_point(aes(color = 'Model')) +
  # scale_color_brewer(palette = 'Set1') +
  scale_color_manual(values = c('Model' = 'red',
                                'Window' = 'darkgreen',
                                'Trap' = 'blue',
                                'Truth' = 'black')) +
  # scale_y_log10() +
  labs(x = 'Date', y = 'Estimate', color = 'Source', title = 'All Daytime Fish')
day_p



# look at rates of wild fish
plot_df %>%
  filter(Variable == 'Wild.PBT.Rate') %>%
  ggplot(aes(x = Start_Date)) +
  # geom_ribbon(data = data.frame(Start_Date = sort(unique(plot_df$Start_Date)),
  #                               other_summ %>%
  #                                 filter(Variable == 'Avg.Prop.Wild.PBT')),
  #             aes(ymin = upp_ci,
  #                 ymax = low_ci,
  #                 fill = 'Average'),
  #             alpha = 0.4) +
  # geom_hline(data = filter(other_summ,
  #                          Variable == 'Avg.Prop.Wild.PBT'),
  #            aes(yintercept =  median,
  #                color = 'Average'),
#            lty = 2) +
geom_ribbon(aes(ymin = low_ci,
                ymax = upp_ci,
                fill = 'Estimate'),
            alpha = 0.4) +
  geom_line(aes(y = median,
                color = 'Estimate')) +
  geom_point(aes(y = median,
                 color = 'Estimate',
                 size = Trap)) +
  geom_line(aes(y = wild_fish / Trap,
                color = 'Obs')) +
  geom_point(aes(y = wild_fish / Trap,
                 color = 'Obs')) +
  geom_line(aes(y = NOR_tot / N_tot,
                color = 'Truth'),
            linetype = 3) +
  scale_color_manual(values = c('Obs' = 'darkgreen',
                                'Estimate' = 'red',
                                'Average' = 'blue',
                                'Truth' = 'black'),
                     name = 'Source') +
  scale_fill_manual(values = c('Obs' = 'darkseagreen1',
                               'Estimate' = 'mistyrose2',
                               'Average' = 'lightskyblue'),
                    name = 'Source') +
  labs(x = 'Date', 
       y = 'Rate', 
       size = 'Fish in Trap',
       title = 'Rate of Wild Fish')

# look at night time passage rate
plot_df %>%
  filter(Variable == 'Night.Time.Rate') %>%
  ggplot(aes(x = Start_Date)) +
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
  geom_ribbon(aes(ymin = low_ci,
                  ymax = upp_ci,
                  fill = 'Estimate'),
              alpha = 0.4) +
  geom_line(aes(y = median,
                color = 'Estimate')) +
  geom_point(aes(y = median,
                 color = 'Estimate',
                 size = mark_night)) +
  geom_line(aes(y = Night.passage / N_tot,
                color = 'Obs')) +
  geom_point(aes(y = Night.passage / N_tot,
                 color = 'Obs')) +
  geom_line(aes(y = sim_params$night.passage.rate,
                color = 'Truth'),
             linetype = 3) +
  scale_color_manual(values = c('Obs' = 'darkgreen',
                                'Estimate' = 'red',
                                'Average' = 'blue',
                                'Truth' = 'black'),
                     name = 'Source') +
  scale_fill_manual(values = c('Obs' = 'darkseagreen1',
                               'Estimate' = 'mistyrose2',
                               'Average' = 'lightskyblue'),
                    name = 'Source') +
  labs(x = 'Date', 
       y = 'Rate', 
       title = 'Night-Time Passage',
       size = 'Night Marked')

# look at reascension rate
plot_df %>%
  filter(Variable == 'Re-Ascension.Rate') %>%
  ggplot(aes(x = Start_Date)) +
  geom_ribbon(data = data.frame(Start_Date = sort(unique(plot_df$Start_Date)),
                                other_summ %>%
                                  filter(Variable == 'Avg.ReAsc')),
              aes(ymin = upp_ci,
                  ymax = low_ci,
                  fill = 'Average'),
              alpha = 0.4) +
  geom_hline(data = filter(other_summ,
                           Variable == 'Avg.ReAsc'),
             aes(yintercept =  median,
                 color = 'Average'),
             lty = 2) +
  geom_ribbon(aes(ymin = low_ci,
                  ymax = upp_ci,
                  fill = 'Estimate'),
              alpha = 0.4) +
  geom_line(aes(y = median,
                color = 'Estimate')) +
  geom_point(aes(y = median,
                 color = 'Estimate',
                 size = mark_night)) +
  geom_line(aes(y = Reascent / N_tot,
                color = 'Obs')) +
  geom_point(aes(y = Reascent / N_tot,
                 color = 'Obs')) +
  geom_line(aes(y = sim_params$reascension.rate * sim_params$fallback.rate,
                color = 'Truth'),
            linetype = 3) +
  scale_color_manual(values = c('Obs' = 'darkgreen',
                                'Estimate' = 'red',
                                'Average' = 'blue',
                                'Truth' = 'black'),
                     name = 'Source') +
  scale_fill_manual(values = c('Obs' = 'darkseagreen1',
                               'Estimate' = 'mistyrose2',
                               'Average' = 'lightskyblue'),
                    name = 'Source') +
  labs(x = 'Date', 
       y = 'Rate', 
       title = 'Re-ascension',
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
              alpha = 0.2) +
  geom_line(aes(y = median,
                color = 'Estimate')) +
  geom_point(aes(y = median,
                 color = 'Estimate',
                 size = Trap)) +
  geom_line(aes(y = Trap / N_tot,
                color = 'Obs')) +
  geom_line(data = my_trap_rate %>%
              filter(Week > as.integer(difftime(min(plot_df$Start_Date), sim_params$start.date, units = 'weeks'))) %>%
              mutate(Week = Week - min(Week) + 1) %>%
              inner_join(plot_df %>%
                           select(Week, Start_Date) %>%
                           distinct()),
            aes(y = trap.rate,
                color = 'Truth')) +
  scale_color_manual(values = c('Truth' = 'blue',
                                'Estimate' = 'red',
                                'Prior' = 'darkgreen',
                                'Obs' = 'darkblue'),
                     name = 'Source') +
  scale_fill_manual(values = c('Estimate' = 'mistyrose2',
                               'Truth' = NULL,
                               'Prior' = 'seagreen'),
                    name = 'Source') +
  labs(title = 'Trap Rate',
       size = 'Fish in Trap')
