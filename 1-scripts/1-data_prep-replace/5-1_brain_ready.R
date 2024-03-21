


#BRAIN READY: Make MRI data ready for preprocessing. 


rm(list = ls())
library('tidyverse')
setwd("~/Documents/_PHD/_BMU/_projects/UKB_CM")
scriptdir<-"1-scripts/1-data_prep-replace/"



#################### P1: Data preparation ######################################## 

#1) Input brain data and merge it with generally used nuisance variables
brain<- read_csv(paste0("2-data/3-dataframes/",
                        "S2_MRI_grey_matter_MRI_and_phenotype_exclusion_filtered_MAIN_DF_n21k.csv"))

vars_brain_nuisance<- read_csv(
  "2-data/3-dataframes/S2_NUISANCE_variables_MRI_n_40681.csv"
) 

vars<- read_csv(
  "2-data/3-dataframes/S1_phenotype_data_BASE_VARS_cm_at_bmi_crp_demographics.csv",
  col_select = c("eid","sex_bin" ,"age_MRI" ,'SES' )
)


rois_374<- read_lines(paste0(scriptdir, 
                             "0-prep-aux-files/ROI_to_study_names_374.txt"))

#2) Merge dataframes
df<- 
  left_join(brain, vars, by = "eid") %>% 
  left_join(., vars_brain_nuisance, by = join_by(eid, eid_20904)) 


#3) Add review variables
df.rev<- read_csv("2-data/3-dataframes/S7_REVIEW_VARIABLES.csv") %>%
  filter(Imaging_sample)
df<- df %>% left_join(., df.rev, by = c("eid", "eid_20904"))
rm(df.rev)
df$ethnicity<- as.factor(df$ethnicity)


#################### P2: save ######################################## 
write_csv(df, "2-data/3-dataframes/S5_brain_ready_topreprocess.csv")
