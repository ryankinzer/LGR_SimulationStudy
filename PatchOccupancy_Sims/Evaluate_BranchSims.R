# Author: Kevin See
# Purpose: Evaluate results from simulated branching models
# Created: 10/17/2016
# Last Modified: 10/31/2016
# Notes: 

#-----------------------------------------------------------------
library(lubridate)
library(magrittr)
library(plyr)
library(dplyr)
library(tidyr)
library(jagsUI)
library(stringr)
library(ggplot2)

theme_set(theme_bw())

setwd('PatchOccupancy_Sims')


#-----------------------------------------------------------------
# load results
load('SimFits/Baseline_2000tags.rda')
load('SimFits/TrapDown_2000tags.rda')

res_df = ldply(res_list, .id = 'sim') %>%
  tbl_df() %>%
  mutate(sim = as.integer(as.character(sim))) %>%
  mutate(param_type = NA,
         param_type = ifelse(grepl('det_p', param), 'Detection', param_type),
         param_type = ifelse(grepl('tot_esc', param), 'Tot_Escape', param_type),
         param_type = ifelse(grepl('main_p', param), 'Movement', param_type),
         param_type = ifelse(grepl('week_esc', param), 'Week_Escape', param_type)) %>%
  mutate(branch = NA,
         branch = ifelse(param_type %in% c('Detection', 'Tot_Escape'), gsub(']', '', sapply(str_split(param, '\\['), function(x) x[2])), branch),
         branch = ifelse(grepl('week_esc', param) | grepl('main_p', param), gsub(']', '', sapply(str_split(param, ','), function(x) x[2])), branch),
         site = ifelse(as.integer(branch) %% 2 == 0, 'upstream', 'downstream'),
         site = ifelse(param_type == 'Detection', site, NA),
         branch = ifelse(grepl('det_p', param), ceiling(as.integer(branch) / 2), branch)) %>%
  mutate(week = NA,
         week = ifelse(grepl('week_esc', param) | grepl('main_p', param), gsub('\\[', '', sapply(str_split(param, ','), function(x) str_sub(x[1], -2))), week)) %>%
  mutate(branch = factor(branch, levels = c(1:10)),
         branch = revalue(branch, 
                          c('10' = 'BB')),
         week = as.integer(week)) %>%
  mutate(pop_size = ifelse(branch %in% c(1:3), 'Small',
                           ifelse(branch %in% c(4:6), 'Medium',
                                  ifelse(branch %in% c(7:9), 'Large', 'BB'))),
         pop_size = factor(pop_size, levels = c('Small', 'Medium', 'Large', 'BB'))) %>%
  mutate(det_prob = ifelse(branch %in% c(1, 4, 7), 'Low/Low',
                           ifelse(branch %in% c(2, 5, 8), 'Low/High',
                                  ifelse(branch %in% c(3, 6, 9), 'High/High', 'None'))),
         det_prob = factor(det_prob, levels = c('Low/Low', 'Low/High', 'High/High', 'None'))) %>%
  mutate(bias = `50%` - true_value,
         rel_bias = bias / true_value,
         inCI = ifelse(`2.5%` <= true_value & `97.5%` >= true_value, T, F))

res_df %>%
  filter(true_value > 0) %>%
  group_by(param_type, branch, site) %>%
  mutate_each(funs(times100 = ifelse(param_type %in% c('Detection', 'Movement'), . * 100, .)), `2.5%`:true_value) %>%
  summarise(bias = round(median(bias), 2),
            rel_bias = round(median(rel_bias), 3),
            inCI = round(sum(inCI) / n(), 3),
            RMSE = sqrt(mean(bias^2)),
            RMSE_cv = RMSE / mean(true_value)) %>%
  as.data.frame()


sim_test = 1
sim_test = sample.int(100, 4)

res_df %>%
  filter(sim == sim_test) %>%
  filter(param_type == 'Movement') %>%
  filter(true_value > 0,
         true_value < 1) %>%
ggplot(aes(x = week,
           y = `50%`,
           color = branch)) +
  geom_ribbon(aes(ymin = `25%`,
                  ymax = `75%`,
                  fill = branch),
              alpha = 0.2,
              color = NA) +
  # geom_ribbon(aes(ymin = `2.5%`,
  #                 ymax = `97.5%`,
  #                 fill = branch),
  #             alpha = 0.2,
  #             color = NA) +
  scale_color_brewer(palette = 'Paired') +
  scale_fill_brewer(palette = 'Paired') +
  scale_y_continuous(trans = 'logit') +
  geom_line(aes(y = true_value),
            linetype = 2) +
  geom_line() +
  facet_wrap(~ sim)
  


res_df %>%
  filter(param_type == 'Movement',
         true_value > 0,
         true_value < 1) %>%
  left_join(res_df %>%
              filter(param_type == 'Week_Escape') %>%
              select(sim, branch, week, week_escape = true_value) %>%
              distinct()) %>%
  group_by(branch) %>%
  mutate(rel_escape = week_escape / max(week_escape)) %>%
  ungroup() %>%
  ggplot(aes(x = true_value,
             y = `50%`)) +
  geom_point(aes(alpha = rel_escape)) +
  scale_alpha_continuous(range = c(0.05, 1)) +
  facet_wrap(~ branch, scales = 'fixed') +
  geom_abline(linetype = 2,
              color = 'red') +
  scale_x_continuous(trans = 'logit') +
  scale_y_continuous(trans = 'logit') +
  labs(title = 'Movement Probability',
       x = 'Truth',
       y = 'Estimate',
       alpha = 'Escapement (relative to max)') +
  theme(legend.position = 'bottom')
  
res_df %>%
  filter(param_type == 'Tot_Escape') %>%
  ggplot(aes(x = branch,
             y = rel_bias)) +
  geom_boxplot(aes(fill = pop_size)) +
  geom_hline(yintercept = 0,
             linetype = 2,
             color = 'darkgray') +
  # facet_wrap(~ pop_size, scales = 'free_x') +
  labs(title = 'Total Escapement',
       x = 'Relative Bias',
       y = 'Branch')
