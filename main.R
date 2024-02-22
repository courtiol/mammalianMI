# Load custom functions ---------------------------------------------------

source("functions.R")

# Checking dependencies ---------------------------------------------------
check_dependencies_all()

# Data preparation --------------------------------------------------------

## Import the phylogenetic tree
tree_full <- ape::read.tree("data/RAxML_bipartitions.result_FIN4_raw_rooted_wBoots_4098mam1out_OK.newick")


## Preparation of the Maternal Investment dataset

### Load the raw full dataset with Maternal Investment data
MI_raw <- read.csv2("data/MI.csv", dec = ".", na.strings = "")

### Format the full dataset (see functions.R for details)
MI_full <- prepare_df_MIfull(MI_raw)

### Prepare subsample with no missing data for modelling

### Prepare subsample with for comparison between orders

### Prepare subsample with for comparison between subclasses

### Prepare subsample for the 20 indicator species
