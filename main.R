# Load custom functions ---------------------------------------------------

source("functions.R")

# Checking dependencies ---------------------------------------------------
check_dependencies_all(c("ape", "nlme", "smatr", "spaMM"))


# Load dependencies -------------------------------------------------------
library(spaMM)
library(smatr)

# Data preparation --------------------------------------------------------

## Import the phylogenetic tree
tree <- ape::read.tree("data/RAxML_bipartitions.result_FIN4_raw_rooted_wBoots_4098mam1out_OK.newick")


## Preparation of the Maternal Investment dataset

### Load the raw full dataset with Maternal Investment data
MI_raw <- read.csv2("data/MI.csv", dec = ".", na.strings = "")

### Format the full dataset (see functions.R for details)
MI_full <- prepare_df_MIfull(MI_raw)
nrow(MI_full) # 1349

### Prepare subsample with for comparison between subclasses
MI_subclasses <- droplevels(MI_full[!is.na(MI_full$Adult_mass) &
                            !is.na(MI_full$Litter_mass) &
                            MI_full$Key %in% tree[["tip.label"]], ])
nrow(MI_subclasses) # 370
str(MI_subclasses)

### Prepare subsample with for comparison between orders
orders_vs_N    <- aggregate(MI_subclasses[, "Key", drop = FALSE], list(Order = MI_subclasses$Order), length)
orders_to_keep <- as.character(orders_vs_N$Order[orders_vs_N$Key >= 15])
MI_orders <- droplevels(MI_subclasses[MI_subclasses$Order %in% orders_to_keep, ])
unique(MI_orders$Order) # [1] Cetartiodactyla Carnivora Dasyuromorphia Diprotodontia Eulipotyphla Primates Rodentia    
nrow(MI_orders) # 327
str(MI_orders)

### Prepare subsample with no missing data for modelling
MI_models <- droplevels(MI_subclasses[!is.na(MI_subclasses$Investment_duration), ])
nrow(MI_models) # 348
str(MI_models)

### Prepare subsample for the 20 indicator species
indicator_species <- c("Tailless tenrec", "Red fox", "Blue whale", "Eurasian shrew",
                       "American bison", "Rat", "African bush elephant", "Impala", "Red panda",
                       "Tammar wallaby", "Tiger", "Grey seal", "Chimpanzee",
                       "Geoffroy's spider monkey", "Southern hairy-nosed wombat", "European hare",
                       "Red kangaroo", "Greater short-nosed fruit bat", "Short-beaked echidna", "Tasmanian devil")
MI_indicators <- droplevels(MI_subclasses[MI_subclasses$Name %in% indicator_species, ])
nrow(MI_indicators) # 20
str(MI_indicators)


# Fitting regression models -----------------------------------------------

## Fitting SLR model for method comparison

fit_SLR_models <- fitme(log(Litter_mass) ~ log(Adult_mass), data = MI_models)
plot(fit_SLR_models, ask = FALSE, which = "mean")    ## diagnostics (good!)
plot(fit_SLR_models, ask = FALSE, which = "predict") ## diagnostics (good!)
extract_fit_summary(fit_SLR_models)
#       estimate  lower  upper
# elev    -0.389 -0.458 -0.320
# slope    0.762  0.741  0.783

## Fitting PLMM model for method comparison

fit_PLMM_models <- fitme_phylo_lambdafree(formula = log(Litter_mass) ~ log(Adult_mass) + corrMatrix(1|Key),
                                          data = MI_models, tree = tree)
plot(fit_PLMM_models, ask = FALSE, which = "mean")  ## diagnostics (good!)
plot(fit_PLMM_models, ask = FALSE, which = "ranef") ## diagnostics (good!)
plot(fit_PLMM_models, ask = FALSE, which = "predict", re.form = NA) ## diagnostics (bad: residual variance captured by random variance)
plot(log(MI_models$Litter_mass), predict(fit_PLMM_models, re.form = NA, type = "link")[, 1]) ## diagnostics, excluding ranef (good!)
extract_fit_summary(fit_PLMM_models)
#           estimate  lower upper
# elevation   -0.453 -1.379 0.473
# slope        0.758  0.725 0.792
# lambda       0.777  0.630 0.873

## Fitting SMA model for method comparison

fit_SMA_models <- sma(Litter_mass ~ Adult_mass, data = MI_models, log = "xy", method = "SMA")
plot(fit_SMA_models,which = "default") ## diagnostics (good!)
plot(fit_SMA_models,which = "residual") ## diagnostics (good!)
plot(fit_SMA_models,which = "qq") ## diagnostics (good!)
extract_fit_summary(fit_SMA_models)
#           estimate  lower  upper
# elevation   -0.171 -0.202 -0.141
# slope        0.788  0.767  0.810

## Fitting MA model for method comparison

fit_MA_models <- sma(Litter_mass ~ Adult_mass, data = MI_models, log = "xy", method = "MA")
plot(fit_MA_models,which = "default") ## diagnostics (good!)
plot(fit_MA_models,which = "residual") ## diagnostics (good!)
plot(fit_MA_models,which = "qq") ## diagnostics (good!)
extract_fit_summary(fit_MA_models)
#           estimate  lower  upper
# elevation   -0.171 -0.201 -0.141
# slope        0.782  0.760  0.804

## Fitting MSLR model for method comparison

fit_MSLR_models <- fitme(log(Litter_mass) ~ log(Adult_mass) + log(Investment_duration), data = MI_models)
plot(fit_MSLR_models, ask = FALSE, which = "mean")    ## diagnostics (good!)
plot(fit_MSLR_models, ask = FALSE, which = "predict") ## diagnostics (good!)
extract_fit_summary(fit_MSLR_models)
#              estimate  lower  upper
# elevation       0.472 -0.146  1.090
# slope           0.801  0.766  0.837
# slope_InvDur   -0.176 -0.301 -0.051
