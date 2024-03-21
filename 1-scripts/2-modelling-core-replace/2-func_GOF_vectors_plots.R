  

#############
    #1. Generates PDFs of a goodness of fit statistics for any two models
    #   that were contrasted (and exist in the same df)
    #2. Generates brain plot for the contrast between any two models


########################## GOF STATISTICS DISTRIBUTIONS ############################################
#NOTE: Supplementary figure - runs  


plt_lavaan_gof <- function(big.df, 
                           model.type,
                           ID_var = NULL) {
  
  #GOAL: Gives supplementary figure of GOF statistics for a given model.
  #The figure being the distributions of the corresponding statistics 
  # across all ROIs. 
  
  require('ggplot2')
  
  #A.---------- Define model to work with
  if (model.type == "simple") {
    model_ID<- "sparse"
    color_pal<- "#4b7d7b"
  }
  
  if (model.type == "complex") {
    model_ID<- "full"
    color_pal<- "#e4960f"
  }
  
  
  #B.----------Prep data
  #1.EXTACT Lavaan model fit indices
  df_gof<- 
    big.df  %>% select(model, lavaan_results) %>%
    unnest(lavaan_results)  %>% 
    select(model, BRAIN_ROI, GOF_measures) %>%
    unnest(GOF_measures) %>%
    filter(grepl(model_ID, model))
  
  mdl.name<- df_gof %>% slice(1) %>% pull(model)
  
  #2. Apply desirable naming conventions
  df_gof<- df_gof %>% select(model, BRAIN_ROI, 
                             rmsea.robust, 
                             srmr,  
                             cfi.robust,
                             chisq.scaled, pvalue.scaled) %>% 
    rename("rmsea" = "rmsea.robust",
           "cfi" = "cfi.robust", 
           "chisq_pvalue" = "pvalue.scaled",
           "chisq" = "chisq.scaled")
  
  
  #C. --------- Chi^2 plot 
  p_chi<- 
    ggpubr::ggdensity(df_gof, 
                      x = "chisq",
                      bin = 5,  rug = TRUE,
                      color = color_pal, 
                      fill = color_pal) + 
    guides(colour = "none", fill = "none") 
  
  if (
    #See if there is S=Z according to a chi^2 test
    #(i.e. there is no significant difference in obs and modeled data)
    any(df_gof$chisq_pvalue>0.05) 
  ){
    
    lims<- df_gof %>% filter(chisq_pvalue<0.05) %>%
      slice_max(chisq_pvalue) %>%
      select(model, chisq, chisq_pvalue) %>%
      pull(chisq) %>% signif(., digits = 3)
    
    p_chi<- p_chi + geom_vline(xintercept = lims ,
                               color = "black",
                               linetype="dashed", size = 1.3)
  }
  
  
  
  
  
  #D. --------- Other fit indices plots
  
  
  p1<- ggpubr::ggdensity(df_gof, 
                         x = "rmsea",
                         bin = 5,  rug = TRUE,
                         color = color_pal, 
                         fill = color_pal
  ) + 
    guides(colour = "none", fill = "none") +
    geom_vline(xintercept = 0.05, color = "black",
               linetype="dashed", size = 1.3) 
  
  
  p2<- 
    ggpubr::ggdensity(df_gof , 
                      x = "cfi",
                      bin = 5,  rug = TRUE,
                      color = color_pal, 
                      fill = color_pal
    ) +
    guides(colour = "none", fill = "none") +
    geom_vline(xintercept = 0.97, color = "black",
               linetype="dashed", size = 1.3) 
  
  
  
  p3<- 
    ggpubr::ggdensity(df_gof,
                      x = "srmr",
                      bin = 5,  rug = TRUE,
                      color = color_pal, 
                      fill = color_pal
    ) +
    guides(colour = "none", fill = "none") +
    geom_vline(xintercept = 0.05, color = "black", linetype="dashed", size = 1.3) 
  
  
  #E. ------------- Save
  ggsave(
    file = paste0(scriptdir, "3-plots/P2/", ID_var,
                  "_", mdl.name, "_", "GOF_dists.png"),
    cowplot::plot_grid(p_chi, p1, p2, p3, ncol = 2, align = "h", axis = "b"),
    width = 6,
    height = 6
  )
  
}




########################## ANOVA CONTRAST MAP & VECTOR #########################
plt_vct_anova<- function(anova.df, ID_var = NULL, return = NULL) {
  
  
  
  #GOAL build discrete mask of ROIs where one model performs better than the other
  #(chi^2 values are ignored here. Only p-values are used).
  
  #1.Prep the dataframe
  anova.df<- anova.df %>% ungroup() %>% tidyr::drop_na()
  anova.df$p.value_anova_FDR<- p.adjust(anova.df$p.value_anova)
  
  
  #2.Set of ROIS where the more complex model fits better than the less
  #complex model
  anova.df <- anova.df%>% mutate(performance_complex = 
                                   if_else(p.value_anova_FDR<0.05, T,F))
  
  #3.Set of ROIS where the complex model is not an improvement on the more
  #complex one.
  anova.df$performance_simpler <- !anova.df$performance_complex
  
  #4. Single vector representation
  anova.df$contrast_color<-  "simple" #green (simple)
  anova.df$contrast_color[anova.df$performance_complex]<- "complex"
  anova.df$contrast_color<- as.factor(anova.df$contrast_color)
  
  #5. Save this vector for later use
  #make it specific for this data
  temp <-  anova.df %>% select(BRAIN_ROI,
                               performance_simpler, performance_complex, contrast_color)
  
  write_csv(temp,
            paste0(scriptdir, "0-aux-files/performance_vectors/vect_relative_performance_", 
                   ID_var, ".csv"))
  
  
  #--------------------------------------- execute plots
  require('ggseg')
  require('ggplot2')
  require(ggsegGlasser)
  
  
  load(paste0(scriptdir, "0-aux-files/glasser_sub/glassersub_ggseg155_R.bin"),
       envir = globalenv())
  glassersub$type<<- 'cortical'
  lab.merge<- read_csv(paste0(scriptdir,
                              "0-aux-files/glasser_sub/glassersub_BRAINROI_merge.csv"))
  
  
  plot_dir<-  paste0(scriptdir, "3-plots/P2/", ID_var, "_", "brain_performance_contrast.png")
  
  
  plot.anova<-
    anova.df %>% filter( ggout == F ) %>%
    left_join(., lab.merge) %>%
    rename("label" = "label_BMU") %>%
    ggseg(., atlas= 'glassersub', 
          mapping=aes(fill = contrast_color), 
          #show.legend = FALSE,
          color = "gray3", size =0.1,
          hemisphere = "right", 
          position='stacked') +
    scale_fill_manual(values = c( "#e4960f", "#4b7d7b"), 
                      na.translate = F ) +
    theme_void() +
    labs(fill = "Performace") +
    theme(legend.text = element_text( size = rel(1.2)),
          legend.title = element_text( size = rel(1.5)))
  
  
  #save
  ggsave(plot_dir,plot.anova)
  
  #return if requested
  if (!is.null(return) && return == T) {
    return(plot.anova)
  }
  
  
}
