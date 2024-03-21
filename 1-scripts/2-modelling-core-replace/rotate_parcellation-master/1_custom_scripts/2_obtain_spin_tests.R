


#SPIN TESTS - Organized 


rm(list = ls())
setwd("~/Documents/_PHD/_BMU/_projects/UKB_CM")
scriptdir<-"1-scripts/2-modelling-core-replace/"
library('tidyverse')

#Current data iteration
script_ID<- "ID_MAIN"


################# 0. Input and prep data to compute spin tests for:##################################


#---------------A. Linear models of BMI, CRP, AT.
LM<- read_csv(paste0(scriptdir, 
                     "1-model_outputs/data_frames/", 
                     script_ID, "_180_LM_Models.csv"))

LM.p<- LM %>% filter(layer == "ctx") %>%
  select(BRAIN_ROI, variable, t.statistic) %>%
  tidyr::pivot_wider(names_from = variable, 
                     values_from = "t.statistic", 
                     names_prefix = "tstat.")


#--------------B. MDD case-control contrasts
#Input MDD data 
MDD_contrast<- data.table::fread(paste0(scriptdir,'0-aux-files/MDD_contrast.txt'))

#Adaptation code:
MDD_contrast<-
  MDD_contrast %>% 
  as_tibble() %>% 
  #1. Identify by hemisphere
  mutate(hemi = case_when(
    str_detect(ROI, "^L_") == T ~ "left",
    str_detect(ROI, "^R_") == T ~ "right")) %>% 
  #2. Adapt naming
  mutate(across(ROI, ~ str_replace(.,"^L_", ""))) %>%
  mutate(across(ROI, ~ str_replace(.,"^R_", ""))) %>%
  mutate(across(ROI, ~str_c(., "_ROI"))) %>%
  #3. Separate left and right
  pivot_wider(
    names_from = "hemi",
    values_from = "tstat",
    names_prefix =  "mdd_tstat_") %>%
  #4. Compute mean of two hemispheres
  mutate(mdd_laterlized_tstat = 
           rowMeans(select(., starts_with("mdd_tstat")))) %>%
  #5. Adapt to my scripts:
  rename("BRAIN_ROI" = "ROI")

#Merge with dfs
LM.p<- LM.p %>% left_join(., MDD_contrast)


################# 1. spin tests :##################################
#NOTE: The function perm.sphere assumes that my ROIs are ordered as in SI Table 1
#of Glasser et al, Nature, 2016. (They should be)


#A. Input relevant functions and data for permutation testing 
#Function available in https://github.com/frantisekvasa/rotate_parcellation
source(paste0(scriptdir,"rotate_parcellation-master/R/perm.sphere.p.R" )) 

#Input matrix of random permutations (unilateral)
load(paste0(scriptdir,"rotate_parcellation-master/2_output_matrices/",
            'perm.id.unilateral.RData'))

#B. Prep 
combi<- c("tstat.CRP", "tstat.BMI" , "tstat.AT",  "mdd_laterlized_tstat"  )
pairs<- as.tibble(t(combn(combi,2)))
pairs$p.permutation<- NA

#C. Perform permutation testing
for (num in 1:nrow(pairs)) {
  x<- LM.p[pairs[[num,1]]] %>% pull()
  y<- LM.p[pairs[[num,2]]] %>% pull()
  
  pairs[num,3] <-  perm.sphere.p(x,y,perm.id.unilateral,corr.type='spearman')
  
}

pairs$significant<- pairs$p.permutation<0.05

print(pairs)


