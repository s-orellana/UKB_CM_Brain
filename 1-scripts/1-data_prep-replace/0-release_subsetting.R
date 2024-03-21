


#GOAL: Function that subsets your UKB release of interest. It will 
#also automatically bridges between applications 54358 and 20904.
#Ultimately the goal is to reduce the number of rows to work with. 


rm(list = ls())
library('tidyverse')
setwd("~/Documents/_PHD/_BMU/_projects/UKB_CM")
scriptdir<-"1-scripts/1-data_prep-replace/"



################# Basic input: ################################################

#Input full dataset associated with the innitial data request: UKB 53358
df<- read_csv('2-data/0-full_data/ukb41173.csv')

#Bridge dataset to the Psychiatry data request: UKB 20904
bridge<- read_csv("2-data/0-full_data/bridgeFile_54358_20904.csv") 
df <- 
  df %>%
  left_join(., bridge, by =  c("eid" = "eid_54358"))


################ FUNC: release filtering ###############################################

UKB_release_subset<-  function(ukb_release, filtration.items) {
  
  #ukb_release: df of UKB release to filtrate
  
  #filtration.items: vector with the field ID's of the columns whose
  #absence of data you will use to subset (i.e. if NA in any of these
  #columns then we filter the release according to those)
  
  
  #1.----- Obtain columns to base exclusion on
  #1.2 If no data is given about the instances, we assume we take 
  #the baseline visit and first measure of the variable
  filtration.items<- 
    filtration.items %>% 
    as.data.frame() %>% 
    mutate_all(~if_else(str_detect(.,"-"), ., str_c(., "-0.0"))) %>% 
    pull()
  
  #Include eids in call
  items<- paste(c('eid',filtration.items), collapse = "|")
  
  #Create df with input items
  excl.columns <- 
    ukb_release %>% 
    select(matches(items)) %>% 
    relocate(starts_with("eid"), everything())
  
  #2.-- Exclude "prefer not to answer" - set to missing
  excl.columns<- 
    excl.columns %>%
    #"Prefer not to answer" is coded as (-818) 
    mutate(across(all_of(filtration.items), 
                  ~replace(., . == -818.00, NA_integer_))) 
  
  
  #3.-- Create exclusion vector:
  #Tells you whether a subject has an NA on any of the
  #columns selected by filtration.items.
  excl.columns<- 
    excl.columns %>% 
    mutate(has_na = 
             ifelse(
               rowSums(is.na(
                 select(.,contains(filtration.items)))) > 0, T, F)) %>%
    select(starts_with("eid"),has_na)
  
  
  #4. Use the information in excl.columns to filter your data
  
  ukb_release %>% 
    left_join(., excl.columns, join_by(eid, eid_20904)) %>% 
    filter(!has_na) %>% return()
  
  
  
}


################ DATA GEN: generate filtered dfs ######################

#1). All subjects with CM data
filtration.cts<- c('20489','20488', '20487', '20490', '20491') 
UKB_cts<- UKB_release_subset(df, filtration.cts)

#2) All subjects with T1 structural images
UKB_MRIT1<- UKB_release_subset(df, "20252-2.0")


############### Organized data saving ##############################
prefix<- "filtered_release_ukb41173_BY_"
path<- "2-data/0-filtered_data/"


write_csv(UKB_cts, paste0(path, prefix, "cts_items_n153633", ".csv" ))
write_csv(UKB_MRIT1, paste0(path, prefix, "T1_acquired_n40681", ".csv" ))







