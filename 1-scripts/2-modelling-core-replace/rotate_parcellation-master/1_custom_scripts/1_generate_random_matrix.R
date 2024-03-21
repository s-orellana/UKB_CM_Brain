

rm(list = ls())
setwd("~/Documents/_PHD/_BMU/_projects/UKB_CM")
scriptdir<-"1-scripts/2-modelling-core-replace/"

#Code adapted from L. Dorfschmidt

#1st create a matrix of random permutations using F. Vasa's code. 


#Run this code once for generating a single matrix
#File & Functions available from https://github.com/frantisekvasa/rotate_parcellation
hcp.centroids = 
  as.matrix(read.table(
    paste0(scriptdir,'rotate_parcellation-master/sphere_HCP.txt')))

source(paste0(scriptdir,"rotate_parcellation-master/R/rotate.parcellation.R" ))
library("matrixStats")

# obtain coordinates of regions on L and R hemispheres
coord.l = hcp.centroids[1:180,];              # coordinates of L hemisphere
coord.r = hcp.centroids[181:360,];            # coordinates of R hemisphere

# run spherical permutation script
perm.id = rotate.parcellation(coord.l,coord.r,nrot=10000)
perm.id.unilateral = perm.id[1:180,]
save(perm.id, file = paste0(scriptdir, "rotate_parcellation-master/2_output_matrices/",
                            'perm.id.RData'))
save(perm.id.unilateral, file = paste0(scriptdir, "rotate_parcellation-master/2_output_matrices/",
                                       'perm.id.unilateral.RData'))


