





#GOAL Obtain filtered phenotypical data for  non-imaging (replication)
#and MRI subjects and


library('tidyverse')
rm(list = ls())
setwd("~/Documents/_PHD/_BMU/_projects/UKB_CM")
scriptdir<-"1-scripts/2-modelling-core-replace/"


#Derive MRI subjects and non-MRI subjects

MRI<- read_csv("2-data/3-dataframes/S2_MRI_grey_matter_MRI_and_phenotype_exclusion_filtered_MAIN_DF_n21k.csv")
vars_unfiltered<- read_csv("2-data/3-dataframes/S1_phenotype_data_BASE_VARS_cm_at_bmi_crp_demographics.csv")


eid_MR<- MRI$eid

#newvars for replication sample:
vars_replication <-
  vars_unfiltered %>% 
  filter(excl_incomplete_phenotype)  %>%
  select(!excl_incomplete_phenotype) %>% 
  #eliminate subjects who are being used in MRI
  filter(!eid %in% eid_MR)

#newvars for imaging sample:
vars_MRI<- 
  vars_unfiltered %>% 
  filter(eid %in% eid_MR)




#save
write_csv(vars_MRI, 
          "2-data/3-dataframes/S2_2_newvars_MRI_n21k.csv")

write_csv(vars_replication, 
          "2-data/3-dataframes/S2_2_newvars_REPLICATION_n116k.csv")






