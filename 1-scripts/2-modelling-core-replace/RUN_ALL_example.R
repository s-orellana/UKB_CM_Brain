
library('tidyverse')
library('lavaan')
rm(list = ls())
setwd("~/Documents/_PHD/_BMU/_projects/UKB_CM")
scriptdir<-"1-scripts/2-modelling-core-replace/"

#-------------------------------0 BASELINE RUN --------------------------
#A. Specify an ID for the run of this script & nuisance formula
#This is the ID of we are rerunning everything to check for consistency 
ID_CURRENT<- "ID_2024" 

#-------------------------------1. RUN H1 RESULTS --------------------------
#This section generates the H1 results and runs the whole nuisance correction 
#from scratch. 
source(paste0(scriptdir, "0-func_pre-imaging-modelling.R"))

#B. Read newvars for MRI and phen - raw and untouched
dfBIG<- read_csv("2-data/3-dataframes/S2_2_newvars_REPLICATION_n116k.csv")
dfMR<- read_csv("2-data/3-dataframes/S2_2_newvars_MRI_n21k.csv")


#CORRECTION: Code adds the 1 for log transform. But we added it ealier too, so 
#we take it out now.
dfBIG$cts_sumscore<- dfBIG$cts_sumscore - 1
dfBIG$at_sumscore<- dfBIG$at_sumscore - 1
dfMR$cts_sumscore<- dfMR$cts_sumscore - 1
dfMR$at_sumscore<- dfMR$at_sumscore - 1


#C. Run function to obtain all tables & results
origional_formula<- variable_value ~ SES + age + sex + age*sex
preimaging_wrapper(dfBIG, dfMR, ID_variable = ID_CURRENT, formula = origional_formula )


#-----------------------------------2. RUN LAVAAN MODELS --------------------------

#BASIC. Prep brain data
#Input: (1) outlier corrected, (2) nuisance corrected, (3) lateralized brain
source(paste0(scriptdir, "0-func-merge-MRI-and-vars.R"))
brain<- read_csv("2-data/3-dataframes/S4_MRI_NO_OUTLIERS_LATERALIZED_nuisance_ORIGINAL_187_rois_n_21738.csv")

#BASIC. Merge with nuisance corrected phenotypical data
brain.ready<- merger(brain, ID_CURRENT)

#A.Input: lavaan modelling function and lavaan models
source(paste0(scriptdir, "0-aux-files/lavaan_models/","1-lavaan-models-MRI.r"))
source(paste0(scriptdir, "1-func_modelling_lavaan.r"))

#B. Run modelling function and save results
result<- completeModelling(brain.ready, model.full, model.sparse)

saveRDS(result[[1]], file = paste0(scriptdir,"1-model_outputs/",
                                   ID_CURRENT, "_", "results_full_sparse_lavaan.rds"))
saveRDS(result[[2]], file = paste0(scriptdir,"1-model_outputs/",
                                   ID_CURRENT, "_", "results_full_sparse_anova.rds"))


#---------------------------2.1. RUN MODEL CONTRAST RESULTS--------------------------


#A. Read GOF analysis functions
source(paste0(scriptdir, "2-func_GOF_vectors_plots.r"))

#If models have been derived already:
store<- readRDS(file = paste0(scriptdir,"1-model_outputs/",
                              ID_CURRENT, "_" ,"results_full_sparse_lavaan.rds"))
df.anova<- readRDS(file = paste0(scriptdir,"1-model_outputs/",
                                 ID_CURRENT,"_" ,"results_full_sparse_anova.rds"))


#A. Execute model-individual GOF analyses
plt_lavaan_gof(store, model.type = "simple", ID_var = ID_CURRENT)
plt_lavaan_gof(store, model.type = "complex", ID_var = ID_CURRENT)

#2. Execute model contrasts
plt_vct_anova(df.anova, ID_var = ID_CURRENT)

#---------------------------2.2. RUN MODEL CONTRAST PLOTS-------------------------- 

source(paste0(scriptdir, "3-func_plotting_lavaan_mlds.R")) 


allplots(store, 
         results_col = "res_stdSolution",
         model.name = "sparse", 
         model.type = "simple", 
         ID_var = ID_CURRENT,
         suffix = "res_stdSolution")

allplots(store, 
         results_col = "res_stdSolution",
         model.name = "full", #sparse
         model.type = "complex", #simple
         ID_var = ID_CURRENT,
         suffix = "res_stdSolution")


#---------------------------3. RUN H2 --------------------------------
source(paste0(scriptdir, "4-func_imaging_lms.R"))

#A.OPTIONAL: Prep linear models ----------------------------
#Execute if not done already. Otherwise ignore. You have the option of rerunning
#but not saving the output to an external file. 
#LM<- lm_creator(brain.ready, roi_n = 180, save=T, script_ID = ID_CURRENT, lateral=T)


#-------------------2. Obtain linear model plots ----------------------------
#Read in linear models if not running the lineas above
LM<- read_csv(paste0(scriptdir, 
                     "1-model_outputs/data_frames/", 
                     ID_CURRENT, "_180_LM_Models.csv"))

#generate SUPP table of results if necessary:
lm_table_output(LM,script_ID =ID_CURRENT)

#A. Brain
plt_brain_wrapper(LM, script_ID = ID_CURRENT)

#B. Top ROI scatterplots
plt_scatt_lm(script_ID = ID_CURRENT, brain.ready)


#=================== Obtain session info ==========================
#All packages and their versions used in running this code:
session_info <- capture.output(sessionInfo())

# Write the captured information to a file
write(session_info, file = paste0(scriptdir,"sessionInfo.txt"))

