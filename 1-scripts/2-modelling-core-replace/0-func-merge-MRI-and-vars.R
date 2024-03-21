

merger<- function(brain, script_ID){
  
  #nuisance regressed according to current ID 
  vars<- read_csv(paste0(
    "1-scripts/2-modelling-core-replace/1-model_outputs/nuisance_regressed/",
    script_ID,
    "_MR_nuisance_REGRESSED_variables.csv"))
  
  #bridging eid:
  eid_bridge<- read_csv(
    "2-data/3-dataframes/S1_phenotype_data_BASE_VARS_cm_at_bmi_crp_demographics.csv",
    col_types = cols_only("eid" = col_double(), 
                          "eid_20904" = col_character())
  )

  brain<- 
    brain %>% 
    left_join(., vars) %>%
    left_join(.,eid_bridge, by = "eid") %>% 
    relocate( eid_20904, .after = everything())

    return(brain)
}
