


#GOAL: This script (1) collects all imaging nuisance regressors and other MRI 
#relevant information into a single dataframe and (2) should also merge all imaging 
#data together into a single df.



rm(list = ls())
library('tidyverse')
setwd("~/Documents/_PHD/_BMU/_projects/UKB_CM")
scriptdir<-"1-scripts/2-modelling-core-replace/"




############################# NUISANCE REGRESSORS DF COLLECTION #############################
#1. Read in eids of UKB subjects with T1 data
MRI<- 
  read_csv("2-data/0-filtered_data/filtered_release_ukb41173_BY_T1_acquired_n40681.csv",
         col_types = cols_only("eid" = col_guess()))

#Read full dataset for later variable extraction:
ukb_release<- #Release filtered for presence of T1
  read_csv("2-data/0-filtered_data/filtered_release_ukb41173_BY_T1_acquired_n40681.csv")



#2) ------- Applications bridge:
bridge<- 
  read.csv('2-data/2-imaging/1-bridgeFile_54358_20904.csv') %>%
  rename("eid" = "eid_54358")



#3) ------- Acquisition:
acq<- 
  read_csv("2-data/2-imaging/Acquisition.csv") %>% select(-1) %>%
  rename("acquisition" = "Acq",
         "eid_20904_long" = "eid") %>%
  mutate(eid_20904 = stringr::str_replace(eid_20904_long, "UKB", ""),
         across(eid_20904, ~as.double(.)))


#------- Attach

MRI <- 
  MRI %>% 
  left_join(., bridge, by = "eid") %>%
  left_join(., acq, by = "eid_20904")

#4) ------------------- Euler and motion
euler1<- read_csv('2-data/2-imaging/QC_R1_Euler.csv')
euler2<- read_csv('2-data/2-imaging/QC_R2_Euler.csv')
motion1<- read_csv('2-data/2-imaging/QC_R1_MotionParameters.csv')
motion2<- read_csv('2-data/2-imaging/QC_R2_MotionParameters.csv')

euler<- bind_rows(euler1, euler2)
euler<- euler %>% 
  mutate(euler_total = rowSums(.[c('eulerLeft', 'eulerRight')])) %>% 
  rename("eid_20904" = "oldID", 
         "eid_20904_long" = "eid")

MRI<- 
  MRI %>%
  left_join(.,euler, by = join_by(eid_20904, eid_20904_long))

#NOTE: fd is for "frame displacement"
motion1<- motion1 %>% select(-1, -eid) %>% rename("eid_20904_long" = 'ID')
motion2<- motion2 %>% select(-1) %>% rename("eid_20904_long" = "eid")

motion<- bind_rows(motion1, motion2)  
rm(motion1, motion2)
rm(euler1, euler2)


#------- Attach
MRI <- 
  MRI %>% left_join(.,motion, by = join_by(eid_20904_long))


#5) -------- Extra: intracraneal volume
TIV<- read_csv("2-data/2-imaging/UKB_extended_MRI_data.csv",
               col_select = c("participant","eTIV")) %>%
  mutate(eid_20904 = stringr::str_replace(participant, "UKB", ""),
         across(eid_20904, ~as.double(.))) %>% 
  select(eid_20904, eTIV)


MRI<- MRI %>% left_join(., TIV, by = join_by(eid_20904))
#--------- SCANNING PARAMETERS UKB: ----------------------------------
#Head and coil position: suggested by UKB as nuisance regressors
# 25756: Lateral X position
# 25757: Lateral Y position
# 25758: Lateral Z position
# 25759: scanner table position - (i.e. Z-coordinate of the coil)

# 54: scan site

var.names<- c("eid","^25756-2.0", "^25757-2.0", "^25758-2.0", "^25759-2.0",
              "^54-2.0")
ukb_param<- 
  ukb_release %>%
  select(matches(paste(var.names, collapse = ))) 

ukb_param %>% 
  rename("head.position_X" = "25756-2.0",
         "head.position_Y" = "25757-2.0",
         "head.position_Z" = "25758-2.0",
         "coil.position"   = "25759-2.0",
         "MRI_scanning.site" = "54-2.0")  %>%
  
  #Across all variables: recode the value -999999. i.e. "measure non-
  #recoverable from the data. 
  
  mutate(across(!contains("eid"),  
                ~ if_else(. == -999999, NA, . ))) %>% 
  
  #-------- MERGE
  left_join(MRI, ., by = c("eid", "eid_20904")) -> MRI

#--------------------------- SAVE ----------------------------------

write_csv(MRI, paste0("2-data/3-dataframes/", "S2_NUISANCE_variables_MRI_n_40681.csv"))

############################# ARRANGING IMAGING DATAFRAMES #############################
rm(list = ls())
MRI<- read_csv( paste0("2-data/3-dataframes/", "S2_NUISANCE_variables_MRI_n_40681.csv"))



#1)---- subcortical structures
subctx1<- read_csv('2-data/2-imaging/Struct_R1_SubcorticalVolume.csv')
subctx2<- read_csv('2-data/2-imaging/Struct_R2_SubcorticalVolume.csv')
identical(names(subctx1), names(subctx2)) 
subctx <- bind_rows(subctx1, subctx2)
rm(subctx1,subctx2)

#Format as needed:
subctx<- 
  subctx %>% select(-"...1") %>% 
  rename("eid_20904_long" = "NewID") 
  #select(contains("eid"), all_of(temp_study_rois))


#2)----- Cortical structures:
cortex_2<- read_csv("2-data/2-imaging/Struct_CT.csv")
cortex_1<- read_csv("2-data/2-imaging/Struct_CT_R1.csv")
identical(names(cortex_2), names(cortex_2))
cortex<- 
  bind_rows(cortex_1, cortex_2) %>% select(-"...1") %>%
  select(-SA, -GM_Vol) %>%
  rename("eid_20904_long" = "NewID") 
rm(cortex_1, cortex_2)

#--------------- MERGE:
#input list of ROIs used in study
study_rois<- 
  read_csv(paste0("1-scripts/1-data_prep-replace/",
                  "0-prep-aux-files/ROI_to_study_names_374.txt"), col_names = FALSE ) %>% pull()

 MRI_grey_matter<-
   MRI %>% select(starts_with("eid")) %>%
    left_join(., cortex, by = "eid_20904_long") %>% 
    left_join(., subctx, by = "eid_20904_long") %>% 
   # Keep in only study-relevant data
    select( all_of(study_rois),eid, eid_20904)
    
 
 
 
 #--------------------------- SAVE ----------------------------------
 
 prefix<- "S2_MRI_grey_matter_"
 
 
 #.1 Save unfiltered
 write_csv(MRI_grey_matter, 
           paste0("2-data/3-dataframes/", prefix, 
                       "unfiltered_n_40681.csv"))
 
 #2. Save MRI  - filtered 
 MRI_grey_matter %>%  
   filter(complete.cases(.)) %>% 
   write_csv(., paste0("2-data/3-dataframes/", prefix, 
                       "COMPLETE_cases_n_37457.csv"))
 
 #3. Save MRI + nuisance filtered -------> BASIS OF STUDIES
 MRI_grey_matter %>%
   left_join(., MRI, by = join_by(eid, eid_20904)) %>%
   #filter by complete cases & presence of T1&T2
   filter(complete.cases(.) &
            acquisition == "T1T2") %>% 
   #keep only study relevant data, again,
   select(all_of(study_rois), eid, eid_20904) %>% 
   write_csv(., paste0("2-data/3-dataframes/", prefix, 
                       "AND_nuisance-vars_QC-FILTERED_n35k.csv"))
 
 
 #4. Save with filtering of phenotype variables (for base analysis variables)
 # and with filtering for MRI nuisance
 newvars<- 
   read_csv("2-data/3-dataframes/S1_phenotype_data_BASE_VARS_cm_at_bmi_crp_demographics.csv",
            col_types = cols_only("eid" = col_guess(),
                                  "excl_incomplete_phenotype" = col_guess()))
 
 MRI_grey_matter %>%
   #filter by complete cases & presence of T1&T2
   left_join(., MRI, by = join_by(eid, eid_20904)) %>%
   filter(complete.cases(.) &
            acquisition == "T1T2") %>% 
  #Add filter: complete phenotypical var cases
   left_join(., newvars, by = "eid") %>% 
   filter(excl_incomplete_phenotype) %>%
   #keep only study relevant data, again,
   select(all_of(study_rois), eid, eid_20904)  %>% 
   write_csv(., paste0("2-data/3-dataframes/", prefix, 
                       "MRI_and_phenotype_exclusion_filtered_MAIN_DF_n21k.csv"))
  
 

 
