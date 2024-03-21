#GOAL: derive df's with (1) nuisance regression and (2) outlier extractions
# and (3) lateralization - if requested, for MRI data




##### env def:
rm(list = ls())
library('tidyverse')
setwd("~/Documents/_PHD/_BMU/_projects/UKB_CM")
scriptdir<-"1-scripts/1-data_prep-replace/"



#################### P1: Data preparation ######################################## 

#1. Input brain data ready to preprocess
df<- read_csv("2-data/3-dataframes/S5_brain_ready_topreprocess.csv")

#2. Input ROI list 
rois_374<- read_lines(paste0(scriptdir, 
                             "0-prep-aux-files/ROI_to_study_names_374.txt"))


#3) Input working functions
source(paste0(scriptdir,"3-func_MRI-nuisance-regression.R"))
source(paste0(scriptdir,"4-func_MRI_outlier_n_lateralization.R"))



#################### P1: standard formula ######################################## 
MAIN_formula<- 
  grey_matter ~  
  sex_bin + age_MRI +  sex_bin * age_MRI + SES +
  euler_total +
  MRI_scanning.site + head.position_X + head.position_Y + head.position_Z +
  coil.position + fd + fd_max 

#1)Correct for nuisance
df_MAIN<- MRI_nuisance(df, MAIN_formula,  rois_374)

#2)Correct for outliers
df_MAIN_out<- MRI_outlier_correction(df_MAIN, rois_374, ID_nuisance = "ORIGINAL",  save = T)

#3)
df_MAIN_out_180<- MRI_lateralized(df_MAIN_out, ID_nuisance = "ORIGINAL",  save = T)


#################### P2: Imaging nuisance only formula ######################################## 
IM_ONLY_formula<- 
  grey_matter ~  
  euler_total +
  MRI_scanning.site + head.position_X + head.position_Y + head.position_Z +
  coil.position + fd + fd_max 


#1)Correct for niuisance
df_IM_ONLY<- MRI_nuisance(df, IM_ONLY_formula,  rois_374)

#2)Correct for outliers
df_IM_ONLY_out<- MRI_outlier_correction(df_IM_ONLY, rois_374)

#3)
df_IM_ONLY_out_180<- MRI_lateralized(df_IM_ONLY_out, ID_nuisance = "IM_ONLY",  save = T)


#################### P3: standard formula + ICV ######################################## 

ICV_MAIN_formula<- 
  grey_matter ~  
  euler_total +
  MRI_scanning.site + head.position_X + head.position_Y + head.position_Z +
  coil.position + fd + fd_max + eTIV


#1)Correct for niuisance
df_ICV_MAIN<- MRI_nuisance(df, ICV_MAIN_formula,  rois_374)

#2)Correct for outliers
df_ICV_MAIN_out<- MRI_outlier_correction(df_ICV_MAIN, rois_374)

#3)
df_ICV_MAINout_180<- MRI_lateralized(df_ICV_MAIN_out, ID_nuisance = "ORIGINAL_plus_ICV",  save = T)



#################### P4: standard formula -MINUS- SES ######################################## 


MAIN_no_SES<- 
  grey_matter ~  
  sex_bin + age_MRI +  sex_bin * age_MRI + 
  euler_total +
  MRI_scanning.site + head.position_X + head.position_Y + head.position_Z +
  coil.position + fd + fd_max 

#1)Correct for niuisance
df_noSES<- MRI_nuisance(df, MAIN_no_SES,  rois_374)

#2)Correct for outliers
df_noSES_out<- MRI_outlier_correction(df_noSES, rois_374)

#3)
df_noSES_out_180<- MRI_lateralized(df_noSES_out, ID_nuisance = "ORIGINAL_no_SES",  save = T)



#################### P5: standard formula + MDD ######################################## 
#Number of subjects will be equal to complete MDD cases (17302) + PLUS ANXIETY NOW
#Input MDD scores
df.rev<- read_csv("2-data/3-dataframes/S7_REVIEW_VARIABLES.csv") %>%
  filter(Imaging_sample)

df.MH<- 
  df.rev %>% filter(!is.na(sums_MDD) & !is.na(sums_anxiety)) %>%
  left_join(.,df)

MAIN_formula_mdd<- 
  grey_matter ~  
  sex_bin + age_MRI +  sex_bin * age_MRI + SES +
  euler_total +
  MRI_scanning.site + head.position_X + head.position_Y + head.position_Z +
  coil.position + fd + fd_max + sums_MDD + sums_anxiety

#1)Correct for nuisance
df_mdd<- MRI_nuisance(df.MH, MAIN_formula_mdd,  rois_374)

#2)Correct for outliers
df_mdd_out<- MRI_outlier_correction(df_mdd, rois_374, ID_nuisance = "ORIGINAL_MDD",  save = T)

#3)
df_MAIN_out_180<- MRI_lateralized(df_mdd_out, ID_nuisance = "ORIGINAL_MDD",  save = T)

##################### P7 WORKFLOW FORMULA TO SIMPLIFY  ####################################
run_workflow <- function(data, formula = NULL, roi_list, id_nuisance) {
  
  # Your function for nuisance regression can go here
  df_nuisance <- MRI_nuisance(data, formula, roi_list)
  
  # Your function for outlier correction can go here
  df_out <- MRI_outlier_correction(df_nuisance, roi_list, ID_nuisance = id_nuisance, save = F)
  
  # Your function for lateralization can go here
  df_final <- MRI_lateralized(df_out, ID_nuisance = id_nuisance, save=T)
  
  return(df_final)
}

#################### P6: standard formula + ethnicity ######################################## 

MAIN_et<- 
  grey_matter ~  
  sex_bin + age_MRI +  sex_bin * age_MRI + SES +
  euler_total +
  MRI_scanning.site + head.position_X + head.position_Y + head.position_Z +
  coil.position + fd + fd_max + ethnicity 


run_workflow(df, formula = MAIN_et, rois_374, "ORIGINAL_plusETNI" )


#################### P7 Ever diagnosed ####################################################
df$mh_ever_diagnosed<- as.factor(df$mh_ever_diagnosed)

MAIN_ever<- 
  grey_matter ~  
  sex_bin + age_MRI +  sex_bin * age_MRI + SES +
  euler_total +
  MRI_scanning.site + head.position_X + head.position_Y + head.position_Z +
  coil.position + fd + fd_max + ethnicity + mh_ever_diagnosed

run_workflow(df, formula = MAIN_ever, rois_374, "ORIGINAL_ever" )


#################### P8: (COMBI) standard formula + ethnicity + ANXIETY ######################################## 

MAIN_combi<- 
  grey_matter ~  
  sex_bin + age_MRI +  sex_bin * age_MRI + SES +
  euler_total +
  MRI_scanning.site + head.position_X + head.position_Y + head.position_Z +
  coil.position + fd + fd_max + ethnicity + sums_MDD + sums_anxiety


run_workflow(df, formula = MAIN_combi, rois_374, "ORIGINAL_combi" )

