



################### BASIC MODELLING FUNCTIONS ####################################


#NOTE: 4 functions do the following:
  #1.Fit lavaan models per ROI
  #2. Prep lavaan model results for data writing and plotting; add FDR
  #corrections considering results of all ROIS.
  #3.Computes a lavaan::anova for any two models.
  #4. Wraps around the functions above to fit and contrast any two models
  #of interest. 




########################### 1.A Model fitting #######################################

#SEPARATE so that this is run for a given amount of models
#Make a different object for model comparison.

df_lavaan<- function(dfo, model,  prefix = NULL ) {
  
  #DF to use
  #model = lavaan model syntax to evaluate
  #Prefix = prefix added to column names so that you know which model these
  #results refer to when merging with other df's
  
  #FIT MEASURES TO EXTRACT!: Only robust measures are taken as we use a "robust"
  #MLR estimator. 
  measures<- c( "cfi.robust", "srmr",
                "rmsea.robust", "rmsea.ci.lower.robust", "rmsea.ci.upper.robust",
                #Just in case these are necessary
                "aic", "bic",
                "chisq", "pvalue",
                "chisq.scaled", "pvalue.scaled"
                
                )
    
  
  columnames<-
    read_csv(paste0(scriptdir,"0-aux-files/BRAIN_ROIS_187.csv")) %>% pull()
  
    
  #DF DEFINITION
  dfo %>% select(-eid, -eid_20904) %>%
  pivot_longer(all_of(columnames), names_to = "BRAIN_ROI", 
               values_to = "grey_matter") %>% 
  dplyr::select(BRAIN_ROI, grey_matter, everything())  %>%
  group_by(BRAIN_ROI) %>%
  nest() %>%
    
  #1. fit lavaan model
  mutate(sem_fit =
           map(data, ~ lavaan::sem(model, data = scale(.x), estimator = "MLR")
           )) %>%
  
  #2.Extract results
  mutate(res_parEstimates = 
          map(sem_fit, ~ .x %>%  lavaan::parameterEstimates() )) %>%
  mutate(res_stdSolution = 
           map(sem_fit, ~ .x %>%  lavaan::standardizedSolution() %>% #edit
                 rename("est" = "est.std")  %>% 
                 #Edits to allow for plotting
                 mutate(NA_in_Z = z ) %>%
                 mutate(across(z, ~if_else(is.na(.), 0,.)))
               
               
               )) %>%
  
 
  
  #3.Create df of GOF values
  mutate(GOF_measures = 
            map(sem_fit, ~ .x %>% lavaan::fitmeasures(., measures))) %>%
    #Make GOF results tidy
  mutate(across(GOF_measures, ~map(.x, ~ .x %>% enframe() %>%
                             pivot_wider(names_from = name, values_from = value) %>%
                             mutate(across(.cols= everything(), as.numeric)) 
    )))
  
  }

#NOTE: you want results derived with lavaan::standardizedSolution() when you
#scale your data beforehand!



###################### 1.B.  UNNEST AND GGSEG PREP LAVAAN RESULTS ############################

#LAVAAN RESULTS UNNESTING: Function gives you un-nested results for the 
#lavaan models, plus by-path FDR-corrected p-values.

df_toplot <- function(fdo, to_unnest) {
  
  # df: Df with nested lavaan model results.
  # to_unnest: column of lavaan results that you want to 
  # unnest. i.e. whether fits were extracted with lavaan's:
  #parameterEstimates( )
  #standardizedSolution( )
  
  #NOTE: You need a column of anatomical labels called "BRAIN_ROI"
  
  
  
  #Redundant paths not necessary for plotting
  paths_exclude<- 
    c("c","a", "b", "a2", "b2")
  
  #Order of levels we want when plotting
  # levels<- c("BMI%->%ROI" , "CRP%->%ROI", "AT%->%ROI"  , 
  #            "CM%->%BMI%->%ROI", "CM%->%CRP%->%ROI", "CM%->%AT%->%ROI", 
  #            "CM%->%BMI%->%CRP%->%ROI" , "CM%->%AT%->%CRP%->%ROI" )
  
  levels<- c("BMI%->%MRI" , "CRP%->%MRI", "AT%->%MRI"  , 
             "CM%->%BMI%->%MRI", "CM%->%CRP%->%MRI", "CM%->%AT%->%MRI", 
             "CM%->%BMI%->%CRP%->%MRI" , "CM%->%AT%->%CRP%->%MRI" )
  
  #----- Begin df editing
  fdo %>% 
    select(BRAIN_ROI, 
           {{to_unnest}}) %>%
    unnest({{to_unnest}}) %>%
    
    #1.Make df friendlier to downstream processes
    filter(op != "~~") %>% 
    mutate(across(.cols = c("lhs","label"), 
                  ~case_when(
                    . == "med1" ~ "a*e",
                    . == "med2" ~ "a*b*f",
                    . == "med3" ~ "a2*g",
                    . == "med4" ~ "a2*b2*f",
                    . == "med5" ~ "c*f",
                    TRUE ~ . ))) %>% 
    
    
    #1.2 Filter out non-brain paths:
    filter(
      !label %in% all_of(paths_exclude)) %>%
    
    #2. Avoid clash with ggplot. Rename "label" to "PATH"
    rename("PATH" = "label") %>%
    
    
    #3. FDR correct by ROI inside each path (i.e. ROI 180)
    ungroup() %>%
    group_by(PATH ) %>% #there E as many versions of a path as ROIs
    mutate(pvalue_FDR_180 = p.adjust(pvalue, method = "fdr")) %>%
    ungroup() %>% group_by(BRAIN_ROI) %>%
    
    
    #------- GGSEG agreemet alterations:
    
    #1.code for cerebral layer:
    mutate(layer = 
             case_when(
               stringr::str_detect(BRAIN_ROI, "_ROI") == T ~ "ctx",
               T == T ~ "sub" ))  %>% 
    
    
    #2. code hemisphere and reorganize 
    mutate(hemisphere = "right") %>% 
    
    
    #3. make a label column for ggplot:
    mutate(label = BRAIN_ROI) %>%
    
    #5. eliminate regions not recognized by ggplot (in theory)
    mutate(ggout = case_when( 
      label == "Accumbens-area" ~ T,
      T ~ F)) %>%
    
    #6. reorganize
    relocate("BRAIN_ROI", "layer", "hemisphere", "label") %>%
    
    
    #---------------- Specific labelling edits for easthetics
    #and to enable use of BMU plotting of glasser:
    
    
    # mutate(variable = 
    #          case_when(PATH == "f" ~ "CRP%->%ROI",
    #                    PATH == "e" ~ "BMI%->%ROI",
    #                    PATH == "g" ~ "AT%->%ROI",
    #                    PATH == "a*e" ~ "CM%->%BMI%->%ROI",
    #                    PATH == "a*b*f" ~ "CM%->%BMI%->%CRP%->%ROI",
    #                    PATH == "a2*b2*f" ~ "CM%->%AT%->%CRP%->%ROI",
    #                    PATH == "a2*g" ~ "CM%->%AT%->%ROI",
    #                    PATH == "c*f" ~ "CM%->%CRP%->%ROI" )) %>%
  
  
  mutate(variable = 
           case_when(PATH == "f" ~ "CRP%->%MRI",
                     PATH == "e" ~ "BMI%->%MRI",
                     PATH == "g" ~ "AT%->%MRI",
                     PATH == "a*e" ~ "CM%->%BMI%->%MRI",
                     PATH == "a*b*f" ~ "CM%->%BMI%->%CRP%->%MRI",
                     PATH == "a2*b2*f" ~ "CM%->%AT%->%CRP%->%MRI",
                     PATH == "a2*g" ~ "CM%->%AT%->%MRI",
                     PATH == "c*f" ~ "CM%->%CRP%->%MRI" )) %>%
    
    #gives the order we want in plotting
    mutate(across(variable, ~factor(., levels = levels))) %>%
    
    #7. BMU-brain maps edit
    rename("label_conventional" = "label") %>%
    mutate(across(label_conventional, ~case_when(
      layer == "ctx" ~ str_replace(., "_ROI", ""),
      TRUE ~ .
    ))) %>% 
    
    mutate(across(label_conventional, ~case_when(
      layer == "sub"  ~str_c("Right-",.),
      layer == "ctx"  ~str_c("rh_R_",.)
      
    ))) %>% 
    mutate(across(label_conventional, as.factor)) -> output
  
  
  #Allow plotting of near null effects (PVAL-edit)
  #This merely allows plotting, does not change results.
  if (to_unnest == "res_stdSolution" ) {
    
    output %>% 
      mutate(NA_in_pval_FDR_180 = pvalue_FDR_180) %>% 
      mutate(across(NA_in_pval_FDR_180, 
                    ~if_else(is.na(.), 1, .))) -> output
    
  }
  
  return(output)
    
  
}


########################### 2. Model COMPARISON #######################################


#Must require the input of two different df's (for two different models) from the above
#segment. 

df_xy_anova <- function(fit_x, fit_y, naming = NULL) {
  
  #INPUTS:
  #fit_x: first model to compare
  #fit_y: second model to compare
  #naming: string vector of the names we can identify the models
  #with. IMPORTANT: order must be the same as that of fit_x and fit_y,
  #otherwise you will mislabel the outputs. 
  
  
  #OUTPUT: a df with lavaan::anova results for comparing the two models
  #inputted into the function.
  
  #0.Check
  if(is.null(naming) | missing(naming) ) { 
    stop("no names provided for the models")}
  
  
  #1.Joining input models
  fit_x<- fit_x %>% select(BRAIN_ROI, sem_fit)
  fit_y<- fit_y %>% select(BRAIN_ROI, sem_fit) 
  
  
  fit_df <- left_join(fit_x, fit_y, by = c("BRAIN_ROI"),
                      suffix = c(".x", ".y"))
  
  
  #2. model comparison & output prep
  fit_df %>%
    mutate(anova_res  = 
             map2(sem_fit.x, sem_fit.y, ~ lavaan::anova(.x,.y))) %>% 
    mutate(across(anova_res, ~map(.x, as_tibble))) %>% #turn into tibble
    mutate(across(anova_res, ~map(.x,                 #add column to dfs
                                  ~.x %>% 
                                    mutate(model =
                                             if_else(Df == 2, 
                                                     naming[1], 
                                                     naming[2]))))) %>%
    select(BRAIN_ROI, anova_res) %>%
    unnest(anova_res) %>% relocate(BRAIN_ROI, model) %>%
    rename(
      'chis_diff' = "Chisq diff",
      "DF_diff"   = "Df diff",
      "p.value_anova" = "Pr(>Chisq)") %>%
  
  
  #------------------ Output from this function should also be GGSEG friendly
  #------- GGSEG agreement alterations:
  
  #1.code for cerebral layer:
  mutate(layer = 
           case_when(
             stringr::str_detect(BRAIN_ROI, "_ROI") == T ~ "ctx",
             T == T ~ "sub" ))  %>% 
    
    
    #2. code hemisphere and reorganize 
    mutate(hemisphere = "right") %>% 
    
    
    #3. make a label column for ggplot:
    mutate(label = BRAIN_ROI) %>%
    
    #5. eliminate regions not recognized by ggplot (in theory)
    mutate(ggout = case_when( 
      label == "Accumbens-area" ~ T,
      T ~ F)) %>%
    
    #6. reorganize
    relocate("BRAIN_ROI", "layer", "hemisphere", "label") %>%
    
    #---------------- Specific labelling edits for easthetics
    #and to enable use of BMU plotting of glasser:
    
    #preserve structure
    group_by(BRAIN_ROI) %>%
    
    rename("label_conventional" = "label") %>%
    mutate(across(label_conventional, ~case_when(
      layer == "ctx" ~ str_replace(., "_ROI", ""),
      TRUE ~ .
    ))) %>% 
    
    mutate(across(label_conventional, ~case_when(
      layer == "sub"  ~str_c("Right-",.),
      layer == "ctx"  ~str_c("rh_R_",.)
      
    ))) %>% 
    mutate(across(label_conventional, as.factor))
    
  
}



########################### 3. RESULTS IN ONE GO (2 MODELS) #######################################

  

completeModelling <- function(data, modelx, modely) {
  
  #INPUT:
  #data = data with all ROIs (there mustn't be redundant ROIs)
  #modelx = lavaan syntax of first model
  #modely = lavaan syntax of second model
  
  
  
  #OUTPUT: For two different lavaan models, fit to a given amount of ROIS
  #this function fits those models at all ROIs (independently) and outputs:
  
  #A list with two items (results within each have results across all ROIs)
  
  #"lavaan_fits": ggseg ready lavaan fit results grouped by model
  #"model_comparisons": ggseg ready output of per-region ANOVA model comparison.  
  
  #0.output hold
  output<- list()
  
  #1.obtain model names for labelling 
  modelx_name<- deparse(substitute(modelx))  %>% str_replace("model.", "")
  modely_name<- deparse(substitute(modely))  %>% str_replace("model.", "")
  
  
  #2.Independently fit lavaan models
  fit_x<- df_lavaan(data, modelx) 
  fit_y<- df_lavaan(data, modely)
  
  
  
  #---------------Plotting output
  
  tempx<- 
    fit_x %>% mutate(model = modelx_name ) %>% relocate(model) %>%
    group_by(model) %>% nest("lavaan_results" = !group_cols()) %>%
    
    mutate("results_parEstimates" =
             map(lavaan_results,  ~ .x  %>%  df_toplot(., "res_parEstimates" ))) %>%
    
    mutate("res_stdSolution" =
             map(lavaan_results,  ~ .x  %>%  df_toplot(., "res_stdSolution" )))
  
  
  tempy<- 
    fit_y %>% mutate(model = modely_name ) %>% relocate(model) %>%
    group_by(model) %>% nest("lavaan_results" = !group_cols()) %>%
    
    mutate("results_parEstimates" =
             map(lavaan_results,  ~ .x  %>%  df_toplot(., "res_parEstimates" ))) %>%
    
    mutate("res_stdSolution" =
             map(lavaan_results,  ~ .x  %>%  df_toplot(., "res_stdSolution" )))
  
  
  
  output$lavaan_fits<- tempx %>% bind_rows(.,tempy)
  
  #---------------Comparison output
  
  model_names<- c(modelx_name, modely_name)
  output$model_comparisons<- df_xy_anova(fit_x, fit_y, naming = model_names)
  
  return(output)
  
  
}








