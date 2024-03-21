
#-----------------------------------------------------------------------------

#GOAL: Prepare model working dataframe for non-imaging variables.

#INPUT: should be a FILTERED UKB data release. 

#NOTE: A lot of this file is necessarily hardcoded.

#EXCLUSION: The file also creates a vector for excluding subjects that do not 
#meet inclusion criteria according to the main study. 

#-----------------------------------------------------------------------------

rm(list = ls())
library('tidyverse')
setwd("~/Documents/_PHD/_BMU/_projects/UKB_CM")
scriptdir<-"1-scripts/1-data_prep-replace/"


prefix<- "filtered_release_ukb41173_BY_"
path<- "2-data/0-filtered_data/"

#Filtered UKB data release
df<- read_csv(paste0(path, prefix, "cts_items_n153633", ".csv"))

#Input CSV file listing the fields you want to extract 
baseline_fields<- 
  read_csv(paste0(scriptdir,"0-prep-aux-files/",
                  "UKB_fields_of_interest_baseline_analyses.csv")) %>% select(1:4)

#load function for extracting UKB fields of interest:
source(paste0(scriptdir,'/1-func_UKB_column_extract.R'))



########### Create baseline dataframes (newvars) #################################


var_names<-
  baseline_fields %>% 
  #do not include MH diagnosis for now
  filter(field_ID != "20544") %>%
  pull(field_ID)
  

#Extract data for my baseline variables of interest
df.base<- extract_columns(df, var_names, visit_num = c(0, 2))


#Extract "ever diagnosed with a mental disorder"
df.mh<- extract_columns(df, "20544") 
  
  #Turn df into categorical variable for "ever diagnosed with a
  #mental disorder
  df.mh<-
    df.mh %>%
    mutate(ever_diagnosed_mentally_ill =
           rowSums(!is.na(select_if(. ,
                                    !startsWith(names(.), "eid")))) > 0)


########### DO: Recode df.base (newvars) ######################################


  
#A. Renaming basic model variables: ------------------------------------
  df.base %>% 
    #1) Recoding of  CRP variables 
    rename(CRP_T0 = "30710-0.0") %>% 
    #2) Age
    rename("age_T0" = "21003-0.0",
           "age_MRI" = "21003-2.0") %>% 
    #3) BMI
    rename(
      "BMI_T0" = "21001-0.0",
      "BMI_MRI" = "21001-2.0" ) %>% 
    #4) Sex 
    rename(
      "sex_bin" = "31-0.0") %>% 
    mutate(sex = if_else(sex_bin == 1, "male", "female")) %>% 
    relocate(starts_with("sex"), .after = starts_with("eid")) %>% 
    
    #5) SES
    rename("SES" = "189-0.0") -> df.base
    


#B. Recoding AT and cts  --------------------------------------------
  #CM computed from the following fields:
  # - 20489: Felt loved as a child (Reverse)
  # - 20488: Physically abused by family as a child
  # - 20487: Felt hated by family as a child
  # - 20490: sexually molested as a child
  # - 20491: someone to take me to the doctor when needed (Reverse)  
  
  
  ##---Adult trauma---
  #- I have been in a confiding relationship (20522) (reversed)
  #-	A partner of ex-partner deliberately hit me or used violence in any other way (20523)
  #-	A partner of ex-partner repeatedly belittled me to the extent that I felt worthless (20521)
  #-	A partner of ex-partner sexually interfered with me, or forced me to have sex against my wishes (20524)
  #-	There was money to pay the rent or mortgage when I needed it (20525) (reversed)
  
  
  
cts_items<- paste(c('20489','20488', '20487', '20490', '20491'), collapse = "|")
  
at_items<-  paste(c('20522', '20523', '20521', '20524', '20525'), collapse = "|")

  df.base<- 
      df.base %>% 
    
      
      #CM--------- 
      #Reverse code "reverse items"
      mutate(across(c("20489-0.0", "20491-0.0"), ~(4 -.)))  %>% 
      #Make sumscores
      mutate(cts_sumscore = rowSums(select(., matches(cts_items)), na.rm = F)) %>% 
      relocate(cts_sumscore, .after = "20491-0.0") %>% 
    
      #AT--------- 
      #select(contains("eid"), matches(at_items)) %>% 
      #Turn "prefer not to answer" (-818) into NA
      mutate(across( matches(at_items),  
                ~replace(., . == -818.00, NA_integer_))) %>% 
      
      #Reverse code "reverse items"
      mutate(across(c("20522-0.0", "20525-0.0"), ~(4 -.))) %>% 
    
      #create sumscore
      mutate(at_sumscore = rowSums(select(., matches(at_items)), na.rm = F)) %>% 
      relocate(at_sumscore, .after = "20525-0.0") %>% 
    
      #create exclusion column
      mutate(excl_at = rowSums(is.na(select(., matches(at_items)))) > 0)
      

      #Sum 1 to cts and at to allow for later logtransforms!
    df.base<- 
        df.base %>% 
        mutate(across(c(at_sumscore, cts_sumscore), ~.+1))
 
#clean df for extraction     
df.base<-  df.base %>%  select(!matches("-0.|-2."))
     
#Quick: perform an exclusion tally
df.base %>%
  filter(!(is.na(CRP_T0) | is.na(BMI_T0) |
             is.na(SES) | excl_at))


#generate exclusion vector:
df.base<- 
  df.base %>% 
  mutate(excl_incomplete_phenotype = 
           if_else(
             !(is.na(CRP_T0) | is.na(BMI_T0) |
                 is.na(SES) | excl_at), T, F)) %>%
  filter(excl_incomplete_phenotype) 


#Minor correction
df.base$sex_bin<- as.factor(df.base$sex_bin)

################### Extract subjects who withdrew by the start of study #########

withdrawed<- read_csv("1-scripts/1-data_prep-replace/0-prep-aux-files/eids_without_withdrawed_subs.csv")
withdrawed<- setdiff(df.base$eid, withdrawed$eid)

df.base <- df.base %>% filter(!eid %in% withdrawed)
    
########### DO: Clean and save ######################################

#Prefix for all dfs created in this script
prefix<- "S1_phenotype_data_"


write_csv(df.base, paste0("2-data/3-dataframes/", prefix,
          "BASE_VARS_cm_at_bmi_crp_demographics.csv"))
    

write_csv(df.base, paste0("2-data/3-dataframes/", prefix,
          "MH_diagnosis_field20544_raw.csv"))
    

    
    
    
  

