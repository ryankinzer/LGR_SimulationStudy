# Author: Kevin See
# Purpose: Summarise night passage and reascension rates over LGR
# Created: 9/13/2016
# Last Modified: 9/13/2016
# Notes: 

#-----------------------------------------------------------------
library(lubridate)
library(plyr)
library(dplyr)
library(tidyr)
library(WriteXLS)

#-----------------------------------------------------------------
# load data
load('/Users/kevin/Dropbox/ISEMP/Analysis_Projects/LGR_TotalEscapement/DataPrepped.rda')

#-----------------------------------------------------------------
# look at historical night passage data
hist_night = night_passage %>%
  filter(!is.na(tot_fish)) %>%
  group_by(Species, Year) %>%
  summarise_each(funs(sum), day_fish:tot_fish) %>%
  mutate(night_rate = night_fish / tot_fish)

hist_night %>%
  group_by(Species) %>%
  summarise(avg_night_rate = mean(night_rate))

# look at PIT tag data (2010 - 2015)
# night rate
pit_night = lgr_weekly %>%
  group_by(Species, Year = SpawnYear) %>%
  summarise_each(funs(sum(., na.rm = T)), day_tags, tot_tags) %>%
  mutate(night_tags = tot_tags - day_tags,
         night_rate = night_tags / tot_tags)

pit_night %>%
  group_by(Species) %>%
  summarise(avg_night_rate = mean(night_rate))
  
# reascent rate
pit_reasc = lgr_weekly %>%
  group_by(Species, Year = SpawnYear) %>%
  summarise_each(funs(sum(., na.rm = T)), reascent_tags, tot_tags) %>%
  mutate(reasc_rate = reascent_tags / tot_tags)

pit_reasc %>%
  group_by(Species) %>%
  summarise(avg_reasc_rate = mean(reasc_rate))

# save as Excel spreadsheet
WriteXLS(c('hist_night',
           'pit_night',
           'pit_reasc'),
         'Outgoing/Night_Reascenion_Rates.xlsx',
         AdjWidth = T)