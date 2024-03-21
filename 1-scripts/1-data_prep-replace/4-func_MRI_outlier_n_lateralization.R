



#Input: dataframes corrected for nuisance variables
#1) Correct for ROI outliers (set to NA)
#2) Lateralize across both hemispheres

############### Eliminate regional outliers ################################ 

MRI_outlier_correction<- function(df, roi_names, ID_nuisance = NULL,  save = NULL){
  
  #For every region check where a subject's value is ± 5 MAD from median
  #if so, this is an outlier, set to NA
  
  # define a function to replace outliers with NA
  replace_outliers <- function(x) {
    b <- mad(x, na.rm = TRUE, constant = 1.4826)
    a <- median(x, na.rm = TRUE)
    ifelse(x > (5 * b + a) | x < (a - 5 * b), NA, x)
  }
  
  # apply the function to every column in CT
  df <- df %>% mutate(across(all_of(roi_names), replace_outliers))
  
  #save if requested
  if (!is.null(save)) {
    write_csv(df , paste0("2-data/3-dataframes/", 
                          "S4_MRI_NO_OUTLIERS_", 
                          "nuisance_",
                          ID_nuisance, "_",  
                          length(roi_names), 
                          "_rois_n_",tally(df)[[1]],".csv"))
    
    return(df)
  } else{
    return(df)}
  
  
}




############### Generate "lateralized" dataframe ############################## 

MRI_lateralized<- function(df, ID_nuisance = NULL, save = NULL){
  
  #You only lateralize a df with 374 ROIs
  
  if (ncol(df) > 200) {
    ########### Separate hemispheres
    left<- df %>% select(starts_with(c("Left", "lh", "eid"))) 
    right<- df %>% select(starts_with(c("Right", "rh", "eid"))) 
    
    # Function to standardize column names
    standardize_colnames <- function(cols) {
      cols %>%
        str_replace("^(lh_L_|Left-)", "") %>%
        str_replace("^(rh_R_|Right-)", "")}
      
      
      # Standardize column names for left and right hemisphere data
      left_cols <- left %>% select(starts_with(c("Left", "lh", "eid"))) %>% 
        names() %>% standardize_colnames()
      right_cols <- right %>% select(starts_with(c("Right", "rh", "eid"))) %>% 
        names() %>% standardize_colnames()
      
      names(left)<- left_cols
      names(right)<- right_cols
      
      
      # Check if standardized column names are identical
      if (!identical(left_cols, right_cols)) {
        stop("Standardized column names are not identical")
      }
      
      
      nhemi<-187
      temp<- data.frame(matrix(ncol = nhemi, nrow = tally(df)[[1]] ))
      brainmean<-  left %>% select(starts_with("eid")) %>% bind_cols(.,temp)
      rm(temp)
      names(brainmean)[2:length(left)]<- left_cols[!str_detect(left_cols, "eid")]
      
      
      for (x in  left_cols[!str_detect(left_cols, "eid")] ) {
        if (names(right[, x]) ==  names(left[, x])) {
          
          a<- bind_cols(right[, x], left[,x])
          #a<- bind_cols(right[, x], left[,x], .warn_conflicts = FALSE)
          brainmean[,x]<- rowMeans(a, na.rm = T)
          rm(a)
        }
        else(print("columns not in the same order"))
      }
      
      brainmean<-
        brainmean %>% 
        relocate(starts_with("eid"), .after = last_col()) 
      
      
      #save if requested
      if (!is.null(save)) {
        write_csv(brainmean , paste0("2-data/3-dataframes/", 
                                     "S4_MRI_NO_OUTLIERS_LATERALIZED_",
                                     "nuisance_",
                                     ID_nuisance,
                                     "_187_rois_n_",tally(df)[[1]],".csv"))
        
        return(brainmean)
      } else{return(brainmean)}
    
  } else {
    print("This does not look like a whole brain df")
    break } 
  
}
    
  
  
  