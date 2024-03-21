# Childhood maltreatment influences adult brain structure through its effects on immune, metabolic and psychosocial factors

Author: Sofia C. Orellana

This repository contains code to generate the main results of the manuscript "Childhood maltreatment influences adult brain structure through its effects on immune, metabolic and psychosocial factors" by Sofia C. Orellana, Richard A.I. Bethlehem,  Ivan Simpson-Kent, Anne-Laura van Harmelen,  Petra E. Vértes and Edward Bullmore. *Proceedings of the National Academy of Sciences* (In press)


## Data
All data used in the construction of the manuscript were obtained via the UK BIOBANK and downloaded from its databases. Researchers wishing to access these data may do so by submitting a request at https://www.ukbiobank.ac.uk/.


## Code

### Directories
These scripts expect to be embedded within a directory called **1-scripts** and within a larger directory structure as bellow, meant to store some of the data inputs and ouputs:

```
.
├── 1-scripts
│   ├── 1-data_prep-replace
│   └── 2-modelling-core-replace
└── 2-data
    ├── 0-filtered_data
    ├── 0-full_data
    ├── 2-imaging
    └── 3-dataframes
```

Other empty directories in the repository are expected by the scripts as output directories, or as the locations where processed data is expected to be. 

### Scripts

- Scripts in 1-data_prep-replace prepare UKB MRI and phenotypical data for main analyses. This involves daata cleanage, merging, application of exclusion criteria and nuisance correction of MRI data. These scripts run sequentially by their numbering. Those starting with _func_ are non executable scripts containing the functions called by other scripts (usually those labelled with the same number). 

- Scripts in 2-modelling-core-replace contain the code necessary to run the main analyses in the study. Each script starting with _func_  is a cluster of functions necessary to run a given set of the analyses. The script **RUN_ALL_example.R** is an example usage that runs the entire study.


**Caveats:**
            - Please note that some of the code in the directory *1-data_prep-replace* is designed for supplementary analyses, but not fully supported by this repository.  
            - Code for spin tests was cloned and adapted from [this repository](https://github.com/frantisekvasa/rotate_parcellation).


