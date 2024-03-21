

library(tidyverse)
rm(list = ls())
setwd("~/Documents/_PHD/_BMU/_projects/UKB_CM")
scriptdir<-"1-scripts/2-modelling-core-replace/"



#This script generates a csv file that allows to merge the ROI labelling 
#of glassersub with the usual naming in UKB dfs in this study.
load(paste0(scriptdir,  "0-aux-files/glasser_sub/glassersub_ggseg155_R.bin"))
glassersub$type='cortical'


#Match to new labeling
#You to want to make sure that you join by BRAIN_ROI

glasser_label<- 
  glassersub %>% select(label, hemi, side) %>% 
  mutate(BRAIN_ROI = label) %>%
  mutate(across(BRAIN_ROI, ~ str_replace(.,"^L_", ""))) %>%
  mutate(across(BRAIN_ROI, ~ str_replace(.,"^R_", "")))  %>%
  mutate(across(BRAIN_ROI, ~case_when(
    BRAIN_ROI == "Thalamus" ~ "Thalamus-Proper",
    T==T ~ .
    
  ))) %>% filter(hemi == "right") %>%
  rename("label_BMU" = "label") %>%
  mutate(atlas = "glassersub")

write_csv(glasser_label, paste0(scriptdir,  "0-aux-files/glasser_sub/glassersub_BRAINROI_merge.csv" ))


#A check that you did not exclude ROIs whose names end in L_ or R_
glasser_label$label_BMU[str_detect(glasser_label$label_BMU, "PSL")]
