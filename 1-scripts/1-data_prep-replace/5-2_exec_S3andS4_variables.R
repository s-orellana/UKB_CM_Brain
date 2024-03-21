

#GOAL: derive df's with (1) nuisance regression for phenotypical variables. 
#This script uses functions created for the main analyses, which is why 
#scriptdir is defined as "1-scripts/2-modelling-core-replace/

rm(list = ls())
library('tidyverse')
setwd("~/Documents/_PHD/_BMU/_projects/UKB_CM")



#################### Read in ready to correct data###################################
vars<- read_csv( "2-data/3-dataframes/S5_vars_MRIsample_ready_to_process_21k.csv")
vars_rep<- read_csv("2-data/3-dataframes/S5_vars_replication_sample_ready_to_process_116k.csv")



#----------SCRIPT PREP
#Nuisance regression for phenotype vars function
source(paste0("1-scripts/2-modelling-core-replace/","0-func_pre-imaging-modelling.R"))

#Required for these functions to work:
scriptdir<-"1-scripts/2-modelling-core-replace/"


#################### MRI 1: standard formula ###################################
#vars<- vars %>% mutate(sex_bin = if_else(sex_bin == 1, 0,1))
#vars<- vars %>% mutate(sex_bin = as.factor(sex_bin))

MAIN_formula<- variable ~ sex_bin + age_T0 +  sex_bin * age_T0 + SES 

vars_MAIN<-  nuisance_lm(vars, MAIN_formula) 
nuisance_save(vars_MAIN, ID_var = "MAIN_formula", nature = "MRI")
  
#################### MRI 2: no formula ######################################### 
#No formula 
vars %>% select(eid, CRP, BMI, CM, AT) %>% 
write_csv(paste0( scriptdir, "1-model_outputs/nuisance_regressed/", 
  "NO_nuisance_CORRECTION_MRI_log_transformed_variables.csv"))

#################### MRI 3: standard formula NO SES ############################
#NO SES 
NOSES_formula<- variable ~ sex_bin + age_T0 +  sex_bin * age_T0

vars_NOSES<-  nuisance_lm(vars, NOSES_formula) 
nuisance_save(vars_MAIN, ID_var = "noSES_formula", nature = "MRI")


#################### MRI 4: standard formula + MDD + anxiety #############################
#PLUS MDD (PLUS ANXIETY!)
#1. Compute new df
df.MH<- 
  vars %>%
  #Keep only complete cases
  filter(Imaging_sample & !is.na(sums_MDD) & !is.na(sums_anxiety)) 


#2. Nuisance regress
MDD_formula<- variable ~ sex_bin + age_T0 +  sex_bin * age_T0 + SES + sums_MDD + sums_anxiety
lm_MDD<-  nuisance_lm(df.MH, MDD_formula) 
nuisance_save(lm_MDD, ID_var = "MAIN_plus_MDD", nature = "MRI")

#################### MRI - BONUS: standard formula + COMBI #############################
COMBI_formula<- variable ~ sex_bin + age_T0 +  sex_bin * age_T0 + SES + 
  sums_MDD + sums_anxiety + ethnicity

#Use vars MH!
combi_lm<- df.MH %>% 
  filter(!is.na(ethnicity)) %>%
  nuisance_lm(., COMBI_formula) 

nuisance_save(combi_lm, ID_var = "MAIN_COMBI", nature = "MRI")

#################### MRI 5: standard formula + ethnicity #############################
et_formula<- variable ~ sex_bin + age_T0 +  sex_bin * age_T0 + SES + ethnicity

ET_lm<-  
  vars  %>% #keep complete cases
  filter(!is.na(ethnicity)) %>% nuisance_lm(., et_formula) 

nuisance_save(ET_lm, ID_var = "MAIN_ET", nature = "MRI")

#################### MRI 5: standard formula + ever diagnosed #############################
mh_formula<- variable ~ sex_bin + age_T0 +  sex_bin * age_T0 + SES + mh_ever_diagnosed
df.MH$mh_ever_diagnosed<- as.factor(df.MH$mh_ever_diagnosed)

mh_lm<-  nuisance_lm(vars, mh_formula) 
nuisance_save(mh_lm, ID_var = "MAIN_EVER", nature = "MRI")





#####################***** REPLICATION SAMPLE ***** ############################
#Standard formula 
rep_vars_MAIN<-  nuisance_lm(vars_rep, MAIN_formula) 
nuisance_save(rep_vars_MAIN, ID_var = "MAIN_formula_AUG23", nature = "REPLICATION")






