



#Three plotting functions for path models!



#1) Plot generating ####################################################

plt_obj<- function(dfo, path, filling, 
                   ID_var = NULL, model.type = NULL,
                   threshold = NULL) {
  
  #Threshold options: T or NULL
  
  #1. Input variables allowing merger with custom gglasser atlas 
  lab.merge<- 
    read_csv(paste0(scriptdir, "0-aux-files/glasser_sub/glassersub_BRAINROI_merge.csv"))
  
  #2. The following needs to be read outside the function
  load(paste0(scriptdir, "0-aux-files/glasser_sub/glassersub_ggseg155_R.bin"),
       envir = globalenv())
  glassersub$type<<- 'cortical'
  
  
  
  
  
  #3. Add performance mask! 
  mask<- read_csv(paste0(scriptdir, "0-aux-files/performance_vectors/", 
                         "vect_relative_performance_", ID_var, ".csv")) 
  
  #3.2 Threshold by mask and (if requested) p-value
  
  if ( threshold == F ) {
    dfo<-
      dfo %>% left_join(., mask, by = "BRAIN_ROI") %>%
      mutate(fill_mask = 
               if_else(contrast_color == model.type,
                       {{filling}}, NA_real_))
  }
  
  if(threshold == T){
    dfo<-
      dfo %>% left_join(., mask, by = "BRAIN_ROI") %>%
      mutate(fill_mask = 
               if_else(contrast_color == model.type &
                         pvalue_FDR_180<0.05,
                       {{filling}}, NA_real_))
    
  }
  
  
  
  #3. Generate plot
  dfo  %>% left_join(., lab.merge) %>% #merge by BRAIN_ROI!
    rename("label" = "label_BMU") %>%
    filter( ggout == F) %>% 
    filter(PATH == path) %>%
    ggseg(., atlas= 'glassersub', 
          mapping=aes(fill = fill_mask), 
          color = "gray3", size =0.3,
          hemisphere = "right", 
          position='stacked')
  
}


#----------- Aesthetics and limits! 

#2) Aesthetics and limits generating functions ############################
plt_aes<- function(dfo, path, filling, fill_title = NULL) {
  
  
  #1.----------------- Define legend limits
  npaths<-  dfo %>% pull(PATH) %>% n_distinct()
  
  path_limits<- dfo %>% group_by(PATH) %>%
    mutate(lim_lower = min({{filling}})) %>% 
    mutate(lim_upper = max({{filling}})) %>%
    #store 
    rowwise() %>%
    #keep only the max
    mutate(limit_value = pmax(abs(lim_lower), abs(lim_upper)))  %>%
    #round
    mutate(across( limit_value, ~ signif(.,2))) %>%
    select(PATH, variable, limit_value) %>%
    ungroup() %>%
    slice_head(n=npaths) 
  
  limval<- path_limits$limit_value[path_limits$PATH == path]

  
  #2.----------------- Define plot title
  titl<- as.character(path_limits$variable[path_limits$PATH == path])
  
  #3.-----------------  Define fill name
  if (fill_title == "z") {
    fill_name<- "z"
    
  }
  if (fill_title == "beta") {
    fill_name<- "beta"
    
  }
  
  #4.-----------------  COLOR PALLETE
  fake_scico <- RColorBrewer::brewer.pal(n = 9, name = "RdBu")

  
  #--------------------------- List of aesthetic elements to add:
  ggproto<- list(
  scale_fill_gradientn(colors = c(low=fake_scico[1], 
                                  mid=fake_scico[5], 
                                  high=fake_scico[9]),
                       na.value = "transparent",
                       limits = c(-limval, limval),
                       breaks = c(-limval, limval),
                       
                       #Lena's solution
                       #limits = c(-5, 5),
                       #breaks = c(-5, 5),
                       #labels = c(expression("">-5), expression(5<"")),
                       oob= scales::squish),
  
  
  
    theme_void(),
    theme(
          legend.position="bottom",
          legend.key.height = unit(rel(0.5), "line"),
          legend.key.width =   unit(rel(4), "line"),
          legend.text = element_text(size = rel(3.5)),
          legend.title = element_text(hjust = 0.5, size = rel(4)),
          plot.title = element_text(hjust = 0, size = rel(2.5), 
                                    color = "black",
                                    face = "bold")
          
          ),
    guides(fill=guide_colourbar(title.position="top"),
         title.hjust =0.5),
    labs(title = rlang::parse_expr(titl)),
    labs(fill = rlang::parse_expr(fill_name))
  
  )
  
}





#3) Single brain wrapper ####################################################
#This function wraps all the functions to plot a single path 

plot_wrap <- function(dfo, path = NULL, 
                      filling, fill_title= NULL, 
                      ID_var = NULL,
                      model.type = NULL,
                      threshold = NULL,
                      diagram = NULL,
                      save.todir = NULL,
                      suffix = NULL,
                      results_col = NULL
) {
  
  #Inputs: 
  # --dfo: plotting-ready df
  # --path: path to plot with our defined notation
  # --filling: CHARACTER or "z" or "beta"
  # --ID_var:  dataset identifyer 
  # --model.type: are we plotting the "simple" or the
  #"complex" model? - if complex, output will be green
  # --threshold: should ROIs shown on the plot be only those
  # that pass significance with an FDR=%5 for 180 ROIs? (T/F)
  # --fill_title: character. Specify between "z" and "beta"
  # --diagram: (T/F) plot with diagram or not - OPTION TO OMIT
  # --save.todir" (T/F) saving output as individual plot in special
  # directory. OPTION TO OMIT. 
  
  ############################ 1. Examine preset conditions
  if (missing(diagram) || is.null(diagram) == T) { diagram = F }
  if (missing(save.todir) || is.null(save.todir) == T) { save.todir = F }
  
  
  
  ############################ 2. Plotting
  plot<- 
    
    #1. Make plot object 
    plt_obj(dfo, path, {{filling}}, 
            ID_var = ID_var,
            model.type = model.type , 
            threshold = threshold) + 
    #2. Add common desired easthetic formatting
    plt_aes(dfo, path, {{filling}}, fill_title = fill_title)
  
  
  #Add corner diagram if wanted!
  if (diagram == T) {
    
    plot<- 
      cowplot::ggdraw() +
      cowplot::draw_image(
        file.path(paste0(scriptdir, "0-aux-files/aux_plots/",
                         model.type, "_", path, ".png" )),
        scale = rel(0.14),
        x = 0.35, y = 0.31) + #y0.27
      cowplot::draw_plot(plot)
  }
  
  ############################ 3. Saving
  #---------- saving conditions: type of fill
  if (fill_title == "z") {fill_dir<- "zetas"}
  
  if (fill_title == "beta") {fill_dir<- "betas"}
  
  #---------- saving conditions: threshold
  if (threshold == T) { thresh_ID<- "THRESH" }
  if (threshold == F) { thresh_ID<- "UNthresh" }
  
  
  #------Determine wether to save or print
  
  ############## Plot directories ############################
  
  #1. Set directory for plots
  # Set the desired directory path
  desired_dir <- paste0(scriptdir, "3-plots/P3/",
                        results_col, "/",
                        ID_var, "_", suffix, 
                        "/", fill_dir)
  
  # Check if the directory exists
  if (!file.exists(desired_dir)) {
    # Create the directory if it does not exist
    dir.create(file.path(desired_dir), recursive = TRUE)
  }

  
  if (save.todir==T) {
    ggsave(
      filename = paste0(desired_dir, "/",
                        ID_var, "_", model.type, "_",
                        thresh_ID, "_", path, ".png" ),
      plot
    )
    
    return(plot)
  }
  
  if(save.todir==F){ print(plot)}
 
  
  
}


#4) Plot ALL brains wrapper ####################################################
#Pass only df of results and print ALL study brain plots (of path analyses)
#at the corresponding directory

allplots <- function(big.df, 
                     results_col = NULL,
                     model.name = NULL, 
                     model.type = NULL,
                     ID_var = NULL,
                     suffix = NULL
) {
  
  require(ggplot2)
  require(ggseg)
  require(ggsegGlasser)
  
  ############## Plot directories ############################
  
  
  dir_zetas<- paste0(scriptdir, "3-plots/P3/",
                     results_col, "/",
                     ID_var, "_", suffix, 
                    "/","zetas")
  
  dir_betas<- paste0(scriptdir, "3-plots/P3/",
                     results_col, "/",
                     ID_var, "_", suffix, 
                     "/","betas")

  
  ############# Predefined variables ##########################
  if (model.type == "complex") {
    paths_loop<- c( "a*e" , "a*b*f", "a2*g" , "a2*b2*f" ,"c*f" )
    opt_ncol<- 2
    opt_nrow<- 3
    op_width <- 15
    op_height <- 25
    
  
  }
  
  if(model.type == "simple"){
    paths_loop<- c( "a*b*f", "a2*b2*f" ,"c*f" )
    opt_ncol<- 3
    opt_nrow<- 1
    op_width <- 20
    op_height <- 12
  }
  
  
  ############# 1. Extract results to plot 
  dfi<- 
    big.df %>% filter(model == model.name) %>% 
    select(model, {{results_col}}) %>%
    unnest(cols = {{results_col}} ) %>% group_by(PATH)
  
  
  
  ############# 2. Plot zetas ########################
  #---A. Untresholded ------ principal 
  plist<- list()
  for (i in 1:length(paths_loop)) {
    plist[[i]]<-
      plot_wrap(dfi, path = paths_loop[i],
                z, fill_title = "z",
                ID_var = ID_var,
                model.type = model.type,
                threshold = F, #---------Key
                diagram = T,
                save.todir = T,
                suffix = suffix,
                results_col = results_col )
  }
  
  ggsave( paste0(dir_zetas, "/",
                 ID_var,"_", model.type,"_UNTH", "_allresults.png" ),
          ggpubr::ggarrange(plotlist = plist, 
                            ncol = opt_ncol, nrow = opt_nrow),
          width = op_width, 
          height = op_height)
  rm(plist)
  
  
  #---B. Thresholded (Principal results)
  plist<- list()
  for (i in 1:length(paths_loop)) {
    plist[[i]]<-
      plot_wrap(dfi, path = paths_loop[i],
                z, fill_title = "z",
                ID_var = ID_var,
                model.type = model.type,
                threshold = T, #---------Key
                diagram = T,
                save.todir = T,
                suffix = suffix,
                results_col = results_col )
  }
  ggsave( paste0( dir_zetas, "/",
                 ID_var,"_", model.type,"_THRESH", "_allresults.png" ),
          ggpubr::ggarrange(plotlist = plist, 
                            ncol = opt_ncol, nrow = opt_nrow),
          width = op_width,
          height = op_height)
  rm(plist)
  
  
  ############# 2. Plot betas ########################
  #---A. Untresholded
  plist<- list()
  for (i in 1:length(paths_loop)) {
    plist[[i]]<-
      plot_wrap(dfi, path = paths_loop[i],
                est, fill_title = "beta",
                ID_var = ID_var,
                model.type = model.type,
                threshold = F, #---------Key
                diagram = T,
                save.todir = T,
                suffix = suffix,
                results_col = results_col )
  }
  ggsave( paste0( dir_betas,"/",
                 ID_var,"_", model.type,"_UNTH", "_allresults.png" ),
          ggpubr::ggarrange(plotlist = plist, ncol= opt_ncol, nrow = opt_nrow),
          width = op_width,
          height = op_height)
  rm(plist)
  
  
  #---B. Thresholded
  plist<- list()
  for (i in 1:length(paths_loop)) {
    plist[[i]]<- 
      plot_wrap(dfi, path = paths_loop[i],
                est, fill_title = "beta",
                ID_var = ID_var,
                model.type = model.type, 
                threshold = T, #---------Key
                diagram = T,
                save.todir = T,
                suffix = suffix,
                results_col = results_col )
  }
  ggsave( paste0(dir_betas, "/",
                 ID_var, "_", model.type,"_THRESH", "_allresults.png" ),
          ggpubr::ggarrange(plotlist = plist, ncol= opt_ncol, nrow = opt_nrow),
          width = op_width, 
          height = op_height)
  rm(plist)
  
  
}
