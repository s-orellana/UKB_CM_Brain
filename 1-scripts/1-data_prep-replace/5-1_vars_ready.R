



#GOAL: Make phenotype variables ready to preprocess (nuisance correct)

rm(list = ls())
library('tidyverse')
setwd("~/Documents/_PHD/_BMU/_projects/UKB_CM")

eid<- read_csv("1-scripts/1-data_prep-replace/0-prep-aux-files/EIDs_MRI_sample_21k.csv")


################################## Data prep ##################################
convention<- 
  c("CM" = "cts_sumscore", "AT" = "at_sumscore",
    "CRP" = "CRP_T0", "BMI" = "BMI_T0")


#-----------MRI DATA
vars<- read_csv("2-data/3-dataframes/S2_2_newvars_MRI_n21k.csv")


#Important variable edits
vars<- 
  left_join(eid, vars) %>%  
  select(eid, CRP_T0, BMI_T0, cts_sumscore, at_sumscore,
         SES, age_T0, sex_bin ) %>%
  rename(all_of(convention)) %>% 
  #IMPORTANT: Log transform values
  mutate(across(c(BMI, CRP, AT, CM), ~ log(.))) 


#add review variables
vars_rev<- 
  read_csv("2-data/3-dataframes/S7_REVIEW_VARIABLES.csv") %>%
  filter(Imaging_sample)

vars<- left_join(vars, vars_rev, by = c("eid"))
rm(vars_rev)

vars$ethnicity<- as.factor(vars$ethnicity)

#----------REPLICATION DATA
vars_rep<- read_csv("2-data/3-dataframes/S2_2_newvars_REPLICATION_n116k.csv")


#Important variable edits
vars_rep<- 
  vars_rep %>%
  select(eid, CRP_T0, BMI_T0, cts_sumscore, at_sumscore,
         SES, age_T0, sex_bin ) %>%
  rename(all_of(convention)) %>% 
  #IMPORTANT: Log transform values
  mutate(across(c(BMI, CRP, AT, CM), ~ log(.)))


################################## save ##################################

write_csv(vars, "2-data/3-dataframes/S5_vars_MRIsample_ready_to_process_21k.csv")
write_csv(vars_rep, "2-data/3-dataframes/S5_vars_replication_sample_ready_to_process_116k.csv")
