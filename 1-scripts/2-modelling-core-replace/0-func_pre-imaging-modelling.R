

####### non-imaging modelling & output functions

     


##################### Non-MRI modelling / analyses####################################

    #NOTE: 
    #A) For simplicity we refer to BMI, CRP, CM, AT as 'phenotype' variables.
    #B) These functions would be useful making a new (if necessary) version
    #   of the phen-variables ready for analysis. 


    #1. UKB_to_df: Log transform raw phenotype values + add qqploots
    #2. nuisance_lm: Corrects phen-variables for nuisance variables.
    #3. nuisance_save: saves the nuisance-corrected variables

    #4. table_nuisance_lm: Makes a table of nuisance variable effects & 
    #   results
    #5. Matrix_scatter: generates a matrix with correlations, PDFs, and 
    #  scatterplots of the relationships between phe-vars
    #6. lavaan_fit_n_table: generates the lavaan fit of the phen-vars of the
    #  study's principal (non-MRI) path model & saves results to a table.
    #7. preimaging_wrapper: Generates all results relying on these functions
    #   for the imaging and non-imaging sample (with the basic data of the
    #   sutdy)

######################## 0. Dataframe prep ######################## 
#Function performs two things:

#1)Take the collection of UKB-filtered variables and turn them into the required 
#format of a df log-tranformed "BMI", "CRP", "cts" and "at".

#2)Provides a q-q plot of the data raw and log-transformed.

UKB_to_df <- function(dfo, ID_var = NULL, qqplot = T) {
  
  require('ggplot2')
  
  dfo<- dfo %>% select(eid, CRP_T0, BMI_T0, 
                       cts_sumscore, at_sumscore,
                       age_T0, sex, SES)
  
  names(dfo)<- c("eid", "CRP", "BMI", "CM", "AT",  
                 "age", "sex", "SES")
  
  
  #Also create a log dfo!
  dfo_log<- dfo %>% 
    mutate_at(vars(CM, AT), ~log(.+1)) %>%
    mutate_at(vars(CRP, BMI), ~log(.))
  
  
  #If plot is requested (False for NON-imaging data)
  if (qqplot == T) {
    
  
  #Store Q-Q plots of model variables nonetheless:
  plot_qq<- function(input, log = F){
    
    pldfo<- input
    pldfo<-pldfo %>% rename(
      "CM" = "CM",
      "AT" = "AT"
    ) %>% 
      select("CRP", "BMI", "CM", "AT") %>%
      pivot_longer( cols = everything())
    
    
    plt_title<- "raw values"
    nature<- "raw_values"
    
    if (log == T) {
      plt_title<- "log trasformed values"
      nature<- "log_transf_values"}
    
    
    
    my_qqplot <- function(.data, .title) {
      ggplot(data = .data, mapping = aes(sample = value)) +
        qqplotr::stat_qq_band(alpha=0.5) +
        #stat_qq_line() +
        qqplotr::stat_qq_point() +
        facet_wrap(~name , scales = "free") +
        labs(x = "Theoretical Quantiles", y = "Sample Quantiles", title = .title)
    }
    
    ggsave(paste0(scriptdir, "3-plots/P0-no-MRI-modelling/",
                  ID_var, "_QQ_PLOTS_", nature, "_", ".png"),
           my_qqplot(pldfo, plt_title),
           width = 6, height = 6)
  }
 
  plot_qq(dfo, log = F)
  plot_qq(dfo_log, log = T)
  }
  
  
  #Always return the log-transformed data
  return(dfo_log) 
  
}


########################## 1. Nuisance regression & saving ouput ##########################
nuisance_lm<- function(df, formula) {
  

  if (all(c("at", "cts") %in% names(df))) {
    columnames <- c("CRP", "BMI", "at", "cts")
  } else if (all(c("AT", "CM") %in% names(df))) {
    columnames <- c("CRP", "BMI", "AT", "CM")
  } else {
    print("Column names don't match expected values")
  }

  #0. Set function inputs
  rhs<- str_split_fixed(as.character(formula), " ~ ", 1)[[3]]
  new_formula<- as.formula(paste("variable_value ~", rhs))
  
  df %>% 
  #1.---------Organize data, compute models, and residuals
  pivot_longer(all_of(columnames), names_to = "variable",
               values_to = "variable_value") %>%
  group_by(variable) %>%
  nest() %>% 
  mutate(model = 
           map(data, 
               ~lm(new_formula,
                   data = .))) %>%
  mutate(model_summary = 
           map(model, ~broom::tidy(.))) %>%
   
   #Obtain residual values 
   mutate(residuals = 
    map(model, ~ pluck(.x, "residuals"))     
   )  %>%
   mutate(eid = 
            map(data, ~ pluck(.x, "eid"))     
   )  %>% return()

} 

nuisance_save <- function(df, ID_var = NULL, nature = NULL) {
  
    df.regressed<- 
      df %>%  
      select(variable, residuals, eid) %>%
      unnest(cols = c("eid", "residuals"))  %>%
      pivot_wider(names_from = "variable", values_from = "residuals")
    
    
    write_csv(df.regressed, paste0(
      scriptdir, "1-model_outputs/nuisance_regressed/", 
      ID_var, "_", nature, "_nuisance_REGRESSED_variables.csv"
    ))
}


########################## 2. Table of nuisance variable effects ##########################
table_nuisance_lm <- function(df, ID_var =NULL, label = NULL) {
  
  #Note: input data must have been modelled already
  #If you want the variables such as age, ses and sex to look nice 
  # (e.g. Age instead of age_t0) this should be formatted in the input
  #formula and df in the preceding function.
  
  require("kableExtra")
  
  df %>% 
    select(variable, model_summary) %>%
    unnest(cols = c(model_summary)) %>%
    ungroup() %>%
    
    #2.-------Change numeric displays
    mutate_at(vars(p.value), ~case_when(
      . == 0 ~ "0.0000",
      . > 0 & . < 0.0001 ~ "<0.0001",
      . > 0 ~ sprintf("%.4f", round(.,4)) ,
      TRUE ~ "NA")) %>%
    mutate_if(is.numeric, ~sprintf("%.3f", round(.,3))) %>%
    
    #3.-------Aesthetic edits:
    filter(term != "(Intercept)") %>%
    mutate_at(vars(term), ~case_when(
      . == "sexmale" ~ "sex (male)",
      . == "age:sexmale" ~ "age*sex (male)",
      TRUE ~ .
    )) %>%
    rename(
      "betaplace" = "estimate",
      "SE" = "std.error",
      "pplace" = "p.value",
      "t" = "statistic",
      "nuisance variable" = "term"
    ) %>%
    #4.--------Format into printable table
    select(-variable) %>%
    kableExtra::kbl(booktabs = T, format = "latex") %>%
    kableExtra::kable_styling(latex_options = c("hold_position")) %>%
    kableExtra::row_spec(0,bold=TRUE) %>%
    pack_rows("CRP", 1,4) %>%
    pack_rows("BMI", 5,8) %>%
    pack_rows("AT", 9,12) %>%
    pack_rows("CM", 13,16)  -> a
  
  
  
  ######### Latex table post-processing:
  #1. Latex beta // CI
  a<- sub("betaplace", "$\\beta$", a, fixed=TRUE)
  a<- sub("pplace", "$P-$value", a,fixed = T)
  
  #2.Add label
  temp<- sprintf("\\label{tbl:%s}\n\\end{table}\n", label)
  a<- sub("\\end{table}",  temp, a, fixed=TRUE)
  
  #3.customization
  a<-sub("\\begin{table}[!h]\n", "\\begin{table}[hb!t]\n", a, fixed=TRUE)
  
  
  
  #4. Store
  file_name<-  paste0(scriptdir, "2-tables/", ID_var, "-", 
                      "Latex_table", "_", label, ".txt") 
  cat(as.character(a), sep = "\n", file = file_name, append = TRUE)

  
  
}

########################## 3. Matrix of variable descriptive plots ##########################
matrix_scatter <- function(dfo, ID_var = NULL,  nature = NULL) {
  
  require("GGally")
  require("ggplot2")
  
  ################## 1. Dataframe prep ############################################################
  dfo<- 
    dfo %>% select(variable, eid, residuals) %>% 
    unnest(cols = c("residuals", "eid")) %>%
    pivot_wider( names_from = "variable" ,
                 values_from = "residuals") %>%
    select(-eid) 
  
  
  ##################### 2. GGALLY setup for scatterplot + density matrices #################### 
  
  #THEME:
  plot_theme<- theme(
    text=element_text(size=30, family="Helvetica"),
    axis.text.x = element_text(angle = 50, size = 20),
    axis.text.y = element_text(angle = 0, size = 20))
  
  
  
  #UPPER Triangular: should give correlation values in plain text
  upperFn <- function(data, mapping, method, ...){
    map_xx <- eval_data_col(data, mapping$x)
    map_yy <- eval_data_col(data, mapping$y)
    corr <- cor(map_xx, map_yy, method=method, use='complete.obs')
    
    #Empty plot: # Draw ggplot2 plot with text only
    p<- ggplot() +                      
      annotate("text", x = 1, y = 1, size = 8,
               label = as.character(round(corr, 3))) + 
      theme_void() 
    p
  }
  
  #LOWER Triangular: should give customized scatterplots &
  #fit an "lm" line. 
  
  lowerFn <- function(data, mapping, method = "lm", color = "some_color", ...) {
    
    map_xx <- eval_data_col(data, mapping$x) 
    map_yy <- eval_data_col(data, mapping$y)
    
    p <- 
      ggplot(data = data, mapping = mapping) +
      geom_point(colour = "black", alpha = 0.3, size = 0.1, shape = 21, fill = "azure4") +
      geom_smooth(method = method, color = "black", ...) +
      geom_density2d(colour = color , size = 0.3) +
      #Code for having only two ticks:
      scale_y_continuous(limits =round(range(map_yy)),
                         breaks =round(range(map_yy))) +     
      scale_x_continuous(limits =round(range(map_xx)),
                         breaks =round(range(map_xx)))
    
    
    p
  }
  
  #DIAGONAL: Allows to customize the ticks of the diagoinal as well
  diagFn <- function(data, mapping,  ...) {
    
    map_xx <- eval_data_col(data, mapping$x) 
    map_yy_den <- density(map_xx)$y
    
    p <- 
      ggplot(data = data, mapping = mapping) +
      geom_density() +
      scale_x_continuous(limits =round(range(map_xx)),
                         breaks =round(range(map_xx))) +
      
      #KEY ofr only two ticks on Y axis
      scale_y_continuous(limits =round(range(map_yy_den),3), #so that the ends are not cut-off
                         breaks =round(range(map_yy_den),1))
    
    p
  }
  
  
  ##################### 3. Generate plot ######################################################## 
  
  ggsave(paste0(scriptdir, "3-plots/P0-no-MRI-modelling/",
                ID_var, "_MATRIX_", nature, ".png"),
         
         ggpairs(
           dfo, 
           lower = list(continuous = wrap(lowerFn, method = "lm", color = "#64dc35")),
           diag = list(continuous = wrap(diagFn)),
           upper = list(continuous = wrap(upperFn, method = 'spearman')) 
         )  + plot_theme,
         width = 6,
         height = 6
         
  )
}


########################## 3. Lavaan fit & table of model outputs ##########################
#Fit data to BASELINE lavaan model and print output


lavaan_fit_n_table<-  function(df.regressed.big, ID_var = NULL, label){
  
  
  ##################################A. Produce lavaan fit
  #-------------0. Input lavaan model syntax.
  source(paste0(scriptdir,"0-aux-files/lavaan_models/", "0-lavaan_model_first_fig.R" ))
  
  
  
  #-----------1. Obtain lavaan fit model and output table 
  data.ready<- 
    df.regressed.big  %>%  
    select(variable, residuals, eid) %>%
    unnest(cols = c("eid", "residuals"))  %>%
    pivot_wider(names_from = "variable", values_from = "residuals")
  
  fit<- sem(model, data = scale(data.ready), estimator = "MLR")
  
  #save fit
  save(fit, file =
         paste0(scriptdir, "1-model_outputs/",
                "BASELINE_MODEL_", 
                ID_var,
                "_lavaan_fit_", label, ".rds"))
  
  
  ##################################B. Output table of results:
  
  
  bold <- function(x) {paste('\\textbf{',x,'}', sep ='')}
  #NOTE: results are reported with the Standardized Solution because
  #we normalized the data beforehand.
  a<- 
    #parameterestimates(fit) %>% 
    standardizedSolution(fit) %>% 
    as_tibble() %>% relocate(op, rhs, lhs) %>%
    rename( "Relationship" = "op", 
            "Outcome" = "lhs",  
            "Predictor" = "rhs",
            "Path" = "label",
            "SE" = "se" ) %>% 
    mutate_at(vars(Path), ~ if_else(Relationship == ":=", Predictor, .)) %>%
    filter(Outcome != Predictor & Outcome != "totalCRP" )  %>%
    mutate_at(vars(Outcome, Predictor), ~case_when(
      . == "AT" ~"AT",
      . == "CM" ~ "CM", 
      Relationship == ":=" ~ "",
      TRUE ~ . ))  %>%
    mutate_at(vars(Relationship), ~case_when(
      Relationship == "~" ~ "Direct",
      Relationship == ":=" ~ "Indirect")) %>%
    
    mutate_at(vars(Predictor), ~if_else(Relationship == "Indirect", "CM", .)) %>%
    mutate_at(vars(Outcome), ~case_when(
      Path == "a*b" ~ "BMICRP",
      Path == "a2*b2" ~ "AACRP",
      TRUE ~.
    )) %>% 
    
    mutate_at(vars(Path), ~case_when(
      Path == "c" ~ "a2",
      Path == "a" ~ "a1",
      Path ==  "a2" ~ "a3",
      Path == "b" ~ "b1",
      Path == "a*b" ~"a1*b1",
      Path == "a2*b2" ~"a3*b2",
      TRUE ~ .
    )) %>% 
    arrange(factor(.$Path, 
                   levels = c("a1", "a2", "a3", "b1", "b2",
                              "", "a1*b1", "a3*b2"))) %>%
    #Change how the p-value is displayed
    mutate_at(vars(pvalue), ~case_when(
      . == 0 ~ "0.0000",
      . > 0 & . < 0.0001 ~ "<0.0001",
      . > 0 ~ sprintf("%.4f", round(.,4)) ,
      TRUE ~ "NA")) %>%
    #Round to 4 digits and keep trailing zeroes
    mutate_if(is.numeric, ~sprintf("%.4f", round(.,4))) %>%
    mutate( "95\\% CI" = stringr::str_c(ci.lower, ", ", ci.upper), 
            #Drop columns used
            ci.upper=NULL, ci.lower=NULL ) %>%
    rename("betaplaceholder" = "est.std",
           "$P-$value" = "pvalue") %>%  
    xtable::xtable(caption = "write caption here") %>% 
    xtable::print.xtable(., include.rownames=FALSE, sanitize.colnames.function=bold, 
                         print.results = getOption("xtable.print.results", FALSE))
  
  #LATEX TABLE POST-PROCESSING
  #1. Latex beta // & ->
  a<- sub("betaplaceholder", "$\\beta$", a, fixed=TRUE)
  a<- sub("BMICRP", "BMI$->$CRP", a, fixed=TRUE)
  a<- sub("AACRP", "AT$->$CRP", a, fixed=TRUE)
  
  #2. Make table fit dimensions of latex page
  a<- sub("\\centering\n", "\\centering\n%Force table to fit within page\n\\resizebox{\\linewidth}{!}{%\n", a, fixed=TRUE)
  a<- sub("\\end{tabular}\n", "\\end{tabular}}\n", a, fixed=TRUE)
  
  #3.Add label
  temp<- sprintf("\\label{tbl:%s}\n\\end{table}\n", label)
  a<- sub("\\end{table}\n",  temp, a, fixed=TRUE)
  
  #4.customization
  a<-sub("\\begin{table}[ht]\n", "\\begin{table}[hb!t]\n", a, fixed=TRUE)
  
  #5. Make caption fit with my custom changes 
  #with then latex "caption" package
  a<- sub("\\caption{write caption here} \n",  "\\caption[]% \n{write caption here} \n", 
          a, fixed=TRUE)
  
  
  
  #6.Store
  file_name<-  paste0(scriptdir, "2-tables/", ID_var, "-", 
                      "Latex_table_", "LAVAAN_FIT_", label, ".txt") 
  cat(as.character(a), sep = "\n", file = file_name, append = TRUE)
  
  
  ##################################C. Fit indices
    fit_indexes<- fitMeasures(fit, c("srmr", 
                     "rmsea.robust", "rmsea.ci.lower.robust", "rmsea.ci.upper.robust",
                     "cfi.robust",
                     "chisq.scaled", "pvalue.scaled")) %>% 
    as.matrix() %>% 
    as.data.frame()  %>% tibble::rownames_to_column() %>%
    as.tibble() %>%
    rename("fit_index" = "rowname", "value" = "V1" ) %>%
    mutate_if(is.numeric, signif, digits = 6) %>%
    mutate_if(is.numeric, round, digits = 4)
    
    file_name_2<- paste0(scriptdir, "2-tables/", ID_var, "_",
                         "FIT_INDICES_", label, ".csv")
    
    write_csv(fit_indexes, file_name_2)

  
}


######################### 4. Wrapper function for all first analyses ########################
preimaging_wrapper <- function(dfBIG, dfMR, ID_variable = NULL, formula = NULL) {
  
  
  
  ################## MEGA-FUNCTION SKELETON ##################
  
  
  #2. Obtain df's in the desired format + print qqplots for the imaging sample.
  dfBIG<- UKB_to_df(dfBIG, ID_var = ID_variable, qqplot = F ) #no qqplts
  dfMR<- UKB_to_df(dfMR, ID_var = ID_variable, qqplot = T )
  
  
  #3.---------- Model out nuisance variables & store df of regressed vars
  #formula is preset for main analyses
  formula_nui<- formula
  
  dfBIG<-nuisance_lm(dfBIG, formula_nui)
  dfMR<-nuisance_lm(dfMR, formula_nui) 
  
  nuisance_save(dfBIG, ID_var = ID_variable, "noImaging")
  nuisance_save(dfMR, ID_var = ID_variable, "MR")
  
  #3.2 Print nuisance variable table results
  table_nuisance_lm(dfBIG, ID_var = ID_variable, label = "nonImaging_nuisance")
  table_nuisance_lm(dfMR, ID_var = ID_variable,  label = "MR_nuisance")
  
  #3.3 Obtain matrix figure
  matrix_scatter(dfBIG, ID_var = ID_variable, nature = "NoImaging")
  matrix_scatter(dfMR, ID_var = ID_variable, nature = "MR")
  
  
  #4.----------- Print lavaan model table results + save fit
  lavaan_fit_n_table(dfBIG, ID_var = ID_variable, label = "nonImaging_fit")
  lavaan_fit_n_table(dfMR, ID_var = ID_variable, label = "MR_fit")
  
}
