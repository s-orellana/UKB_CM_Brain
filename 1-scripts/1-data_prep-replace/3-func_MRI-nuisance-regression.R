


#Function: Regresses nuisance variables as requested by the formula
#input. This is for imaging data only.

#CAUTIONS: Note how we add the mean back to the residuals to preserve the
#biological interpretability of ct/vol values.

MRI_nuisance <- function(data, formula, columns_to_regress) {
  
  #This only works for nuisance regression of brain data!
  #Thus, all pivot_longer values will be "grey_matter", and 
  # "columns_to_regress" should be ROI names
  

  
  #0. Set function inputs
  formula <- as.formula(formula)
  
  #1. Nest data to allow for modelling
  data %>% 
    pivot_longer(all_of(columns_to_regress), names_to = "BRAIN_ROI",
                 values_to = "grey_matter")  %>% 
    group_by(BRAIN_ROI) %>% 
    nest() %>% 
    
    #2. Compute model of interest
    mutate(model = map(data, ~lm(formula, data = .,  na.action="na.exclude")),
           residuals = map(model, resid),
           eid = map(data, ~ pluck(.x, "eid")),
           #!!We add the mean for biological plausibility (not obtaining negative CT/VOL values)
           mean_gm = map_dbl(data, ~mean(.x$grey_matter, na.rm = TRUE)),
           residuals = map2(residuals, mean_gm, ~ .x + .y)) %>%
    
    #3. Unnest data showing residuals
    select(BRAIN_ROI, residuals, eid) %>%
    unnest( c(residuals, eid)) %>%
    pivot_wider(names_from  = BRAIN_ROI,
                id_cols = eid,
                values_from = residuals) %>% 
    relocate("eid", .after = last_col())
}



