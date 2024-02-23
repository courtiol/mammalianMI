# Load custom functions ---------------------------------------------------

source("functions.R")

# Checking dependencies ---------------------------------------------------
check_dependencies_all(c("ape", "ggplot2", "nlme", "scales", "smatr", "spaMM", "tidyr"))


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
nrow(MI_full) # 382

### Prepare subsample with for comparison between subclasses
MI_subclasses <- droplevels(MI_full[MI_full$Key %in% tree[["tip.label"]], ])
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


# Descriptive statistics --------------------------------------------------

cor_global <- cor.test(log(MI_full$Adult_mass), log(MI_full$Litter_mass))
round(cor_global$estimate, digits = 3)[[1]] # correlation estimate
# [1] 0.966
signif(cor_global$p.value, digits = 3) # pvalue
# [1] 2.17e-226


# Fitting regression models -----------------------------------------------

## Fitting SLR model for method comparison

fit_SLR_models <- fitme(log(Litter_mass, 10) ~ log(Adult_mass, 10), data = MI_models)
plot(fit_SLR_models, ask = FALSE, which = "mean")    ## diagnostics (good!)
plot(fit_SLR_models, ask = FALSE, which = "predict") ## diagnostics (good!)
extract_fit_summary(fit_SLR_models)
#              estimate  lower  upper
# intercept      -0.169 -0.199 -0.139
# 10^intercept    0.678  0.632  0.726
# slope           0.762  0.741  0.783
compure_r2(fit_SLR_models)
#    estimate lower upper         p
# r2    0.934 0.959 0.973 1.04e-206

## Fitting PLMM model for method comparison

fit_PLMM_models <- fitme_phylo_lambdafree(formula = log(Litter_mass, 10) ~ log(Adult_mass, 10) + corrMatrix(1|Key),
                                          data = MI_models, tree = tree)
plot(fit_PLMM_models, ask = FALSE, which = "mean")  ## diagnostics (good!)
plot(fit_PLMM_models, ask = FALSE, which = "ranef") ## diagnostics (good!)
plot(fit_PLMM_models, ask = FALSE, which = "predict") ## diagnostics (bad: residual variance captured by random variance)
plot(log(MI_models$Litter_mass, 10), predict(fit_PLMM_models, re.form = NA, type = "link")[, 1]) ## diagnostics, excluding ranef (good!)
extract_fit_summary(fit_PLMM_models)
#              estimate  lower upper
# intercept      -0.197 -0.599 0.206
# 10^intercept    0.636  0.252 1.610
# slope           0.758  0.725 0.792
# lambda          0.777  0.630 0.873
compure_r2(fit_PLMM_models) ## same as above!
#    estimate lower upper         p
# r2    0.934 0.959 0.973 1.04e-206

## Fitting SMA model for method comparison

fit_SMA_models <- sma(Litter_mass ~ Adult_mass, data = MI_models, log = "xy", method = "SMA")
plot(fit_SMA_models,which = "default") ## diagnostics (good!)
plot(fit_SMA_models,which = "residual") ## diagnostics (good!)
plot(fit_SMA_models,which = "qq") ## diagnostics (good!)
extract_fit_summary(fit_SMA_models)
#              estimate  lower  upper
# intercept      -0.171 -0.202 -0.141
# 10^intercept    0.674  0.628  0.723
# slope           0.788  0.767  0.810
compure_r2(fit_SMA_models)
#    estimate lower upper         p
# r2    0.983  0.99 0.993 1.29e-309

## Fitting MA model for method comparison

fit_MA_models <- sma(Litter_mass ~ Adult_mass, data = MI_models, log = "xy", method = "MA")
plot(fit_MA_models,which = "default") ## diagnostics (good!)
plot(fit_MA_models,which = "residual") ## diagnostics (good!)
plot(fit_MA_models,which = "qq") ## diagnostics (good!)
extract_fit_summary(fit_MA_models)
#              estimate  lower  upper
# intercept      -0.171 -0.201 -0.141
# 10^intercept    0.675  0.629  0.724
# slope           0.782  0.760  0.804
compure_r2(fit_MA_models)
#    estimate lower upper         p
# r2    0.974 0.984  0.99 1.09e-277

## Fitting MSLR model for method comparison

fit_MSLR_models <- fitme(log(Litter_mass, 10) ~ log(Adult_mass, 10) + log(Investment_duration, 10), data = MI_models)
plot(fit_MSLR_models, ask = FALSE, which = "mean")    ## diagnostics (good!)
plot(fit_MSLR_models, ask = FALSE, which = "predict") ## diagnostics (good!)
extract_fit_summary(fit_MSLR_models)
#              estimate   lower   upper
# intercept       0.205 -0.0633  0.4730
# 10^intercept    1.600  0.8640  2.9700
# slope           0.801  0.7660  0.8370
# slope_InvDur   -0.176 -0.3010 -0.0505
compure_r2(fit_MSLR_models)
#    estimate lower upper         p
# r2    0.936  0.96 0.973 2.46e-208

## Fitting MPLMM model for method comparison
fit_MPLMM_models <- fitme_phylo_lambdafree(formula = log(Litter_mass, 10) ~ log(Adult_mass, 10) + log(Investment_duration, 10) + corrMatrix(1|Key),
                                           data = MI_models, tree = tree)
plot(fit_MPLMM_models, ask = FALSE, which = "mean")  ## diagnostics (good!)
plot(fit_MPLMM_models, ask = FALSE, which = "ranef") ## diagnostics (good!)
plot(fit_MPLMM_models, ask = FALSE, which = "predict", re.form = NA) ## diagnostics (bad: residual variance captured by random variance)
plot(log(MI_models$Litter_mass, 10), predict(fit_MPLMM_models, re.form = NA, type = "link")[, 1]) ## diagnostics, excluding ranef (good!)
extract_fit_summary(fit_MPLMM_models)
#              estimate   lower upper
# intercept      -0.453 -1.0200 0.110
# 10^intercept    0.353  0.0965 1.290
# slope           0.737  0.6900 0.783
# slope_InvDur    0.127 -0.0600 0.314
# lambda          0.796  0.6540 0.885
compure_r2(fit_MPLMM_models)
#    estimate lower upper         p
# r2    0.932 0.957 0.972 1.22e-203


# Figure 1 ----------------------------------------------------------------

draw_figure_1(data_models = MI_models,
              fit_SLR = fit_SLR_models, fit_PLMM = fit_PLMM_models,
              fit_SMA = fit_SMA_models, fit_MA = fit_MA_models,
              fit_MSLR = fit_MSLR_models, fit_MPLMM = fit_MPLMM_models)
ggplot2::ggsave(filename = "figures/Fig1.pdf", scale = 0.6)
ggplot2::ggsave(filename = "figures/Fig1.png", scale = 0.6)
