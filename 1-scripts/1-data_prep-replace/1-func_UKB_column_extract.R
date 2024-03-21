

#Function: Extract fields of interest --------------------------------------
#This function allows you to specify the field ID, the visit/instance(s) to extract,
# and the array entries. All but the field ID are optional. If nothing optional is 
# specified you get everything encoded by a given field ID. 

# Define function to extract columns of interest from a dataset

extract_columns <- function(df, var_names, visit_num = NULL, array_idx = NULL) {
  
 
  # Create regex pattern for matching variable names and visit numbers
  if (!is.null(visit_num)) {
    visit_num <- as.character(visit_num)
    if (length(visit_num) > 1) {
      visit_pattern <- paste0("((", paste(visit_num, collapse = "|"), "))")
      var_pattern <- paste(var_names, visit_pattern, sep = "-")
    } else {
      var_pattern <- paste(var_names, visit_num, sep = "-")
    }
  } else {
    var_pattern <- var_names
  }
  
  # If array_idx is specified, create regex pattern for matching array columns
  if (!is.null(array_idx)) {
    array_idx <- as.character(array_idx)
    array_pattern <- paste0("(", paste(array_idx, collapse = "|"), ")$")
    var_pattern <- paste(var_pattern, array_pattern, sep = "-")
  }
  
  # Create the final regex pattern
  final_pattern <- paste(var_pattern, collapse = "|")
  
  # Include special columns whose names start with "eid"
  eid_cols <- df %>% select(starts_with("eid"))
  
  # Return the selected columns along with 'eid' columns
  df %>% 
    select(matches(final_pattern)) %>%
    bind_cols(eid_cols, .) %>% 
    relocate(starts_with("eid"), everything())
}
