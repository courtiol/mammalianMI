# Load custom functions ---------------------------------------------------

source("functions.R")

# Checking dependencies ---------------------------------------------------
check_dependencies_all(c("ape", "coin", "doSNOW", "ggdist", "ggplot2", "nlme",
                         "patchwork", "rphylopic", "smatr", "spaMM", "tidyr"))


# Load dependencies used here (and not just in functions imported above) ----
library(spaMM)
library(smatr)

# General settings --------------------------------------------------------

## By default, the script will not run intensive computational steps, replace FALSE by TRUE if you want too
run_slow <- FALSE


# Data preparation --------------------------------------------------------

## Import the phylogenetic tree (downloaded from https://datadryad.org/stash/dataset/doi:10.5061/dryad.q2bvq83r2)
tree <- ape::read.tree("data/RAxML_bipartitions.result_FIN4_raw_rooted_wBoots_4098mam1out_OK.newick")


## Preparation of the Maternal Investment dataset

### Load the raw full dataset with Maternal Investment data
MI_raw <- read.csv2("data/MI.csv", dec = ".", na.strings = "")

### Format the full dataset (see functions.R for details)
MI_full <- prepare_df_MIfull(MI_raw)
nrow(MI_full) # 1041

### Prepare subsample with for comparison between subclasses
MI_subclasses <- droplevels(MI_full[!is.na(MI_full$Investment_duration) & MI_full$Key %in% tree[["tip.label"]], ])
nrow(MI_subclasses) # 738
str(MI_subclasses)

### Prepare subsample without dropping species without investment duration
MI_subclasses_noD <- droplevels(MI_full[MI_full$Key %in% tree[["tip.label"]], ])
nrow(MI_subclasses_noD) # 801
str(MI_subclasses_noD)

### Prepare subsample for the comparison between orders
MI_orders <- prepare_df_MIfull.orders(MI_subclasses)
nrow(MI_orders) # 699
str(MI_orders)

#### Prepare subsample for the comparison between orders solely including Eutheria
MI_orders_euth <- droplevels(MI_orders[MI_orders$Subclass == "Eutheria", ])
nrow(MI_orders_euth)
# [1] 632

#### Prepare subsample for the comparison between orders solely including Metatheria
MI_orders_meta <- droplevels(MI_orders[MI_orders$Subclass == "Metatheria", ])
nrow(MI_orders_meta)
# [1] 67

### Prepare subsample with no missing data for modelling 
MI_models <- MI_subclasses
nrow(MI_models) # 738
str(MI_models)

### Prepare subsample with no missing data for mass proxy comparison
MI_mass <- droplevels(MI_subclasses[!is.na(MI_subclasses$Female_adult_mass), ])
nrow(MI_mass) # 105
str(MI_mass)
sum(MI_mass$Family %in% c("Cercopithecidae")) ## number of Cercopithecidae
# [1] 44
sum(MI_mass$Family %in% c("Odobenidae", "Otariidae", "Phocidae")) ## number of Pinnipedia
# [1] 6


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

cor_global <- cor.test(MI_full$Adult_mass_log10, MI_full$Litter_mass_log10)
round(cor_global$estimate, digits = 3)[[1]] # correlation estimate
# [1] 0.967
cor_global$p.value # pvalue
# [1] 0


# Fitting regression models -----------------------------------------------

## Fitting SLR model for method comparison

### Fitting SLR model
fit_SLR_models <- fitme(Litter_mass_log10 ~ Adult_mass_log10, data = MI_models)

### Checking SLR model assumptions
plot(fit_SLR_models, ask = FALSE, which = "mean")    ## diagnostics (good!)
plot(fit_SLR_models, ask = FALSE, which = "predict") ## diagnostics (good!)

### Computing CI for estimates in SLR model
pretty(extract_fit_summary(fit_SLR_models))
#                  estimate lower_asymptotic upper_asymptotic
# (Intercept)        -0.196           -0.216           -0.176
# Adult_mass_log10    0.778            0.765            0.791
# 10^(Intercept)      0.637            0.608            0.666

### Computing R2 in SLR model
compure_r2(fit_SLR_models)
#    estimate lower upper    p
# r2    0.947 0.939 0.954 0.00


## Fitting PLMM model for method comparison

### Fitting PLMM model with estimation of best Pagel's Lambda (version 1: very slow)
if (run_slow) {
  fit_PLMM_models <- fitme_phylo_lambdafree(
    tree = tree, data = MI_models,
    args_spaMM = list(formula = Litter_mass_log10 ~ Adult_mass_log10 + corrMatrix(1|Key),
                      # fixed = list(phi = 1e-5) ## fix phi if no resid.model!
                      resid.model =  ~ Adult_mass_log10 + (1|Key)))
}

### Fitting PLMM model without estimation of best Pagel's Lambda (version 2: much faster)
### Note: since best Pagel's Lambda is 1 on these data, the same output can be quickly obtained as follows
fit_PLMM_models <- fitme_phylo_lambdafixed(
  lambda = 1, data = MI_models, tree = tree,
  args_spaMM = list(formula = Litter_mass_log10 ~ Adult_mass_log10 + corrMatrix(1|Key),
                    resid.model =  ~ Adult_mass_log10 + (1|Key)))

### Profile for Pagel's lamba
if (run_slow) {
  profile_lambda_PLMM <- profile_lambda(fit_PLMM_models)
  plot(logLik ~ Pagel_lambda, data = profile_lambda_PLMM, type = "o")
}

### Checking PLMM model assumptions
plot(fit_PLMM_models, ask = FALSE, which = "mean")  ## diagnostics (heteroscedastic, but this is accounted for)
plot(fit_PLMM_models, ask = FALSE, which = "ranef") ## diagnostics (ok)
plot(fit_PLMM_models, ask = FALSE, which = "predict") ## diagnostics (bad: residual variance partially captured by random variance)
plot(MI_models$Litter_mass_log10, predict(fit_PLMM_models, re.form = NA, type = "link")[, 1]) ## diagnostics, excluding ranef (good!)

### Computing CI for estimates in PLMM model (very computationally intensive)
if (run_slow) {
  PLMM_summary <- extract_fit_summary(fit_PLMM_models)
  pretty(PLMM_summary$fixef) # Note: we report the basic intervals in the MS
#                   estimate lower_normal upper_normal lower_percent upper_percent lower_basic upper_basic
#  (Intercept)        -0.192       -0.544        0.157        -0.556         0.150      -0.534       0.173
#  Adult_mass_log10    0.806        0.781        0.832         0.780         0.831       0.782       0.833
#  10^(Intercept)      0.643        0.286         1.44         0.278          1.41       0.293        1.49
  pretty(PLMM_summary$Pagel_Lambda)
# estimate    lower    upper 
#   "1.00"  "0.999"   "1.00" 
}

### Computing R2 in PLMM model
compure_r2(fit_PLMM_models) ## same as above
#    estimate lower upper    p
# r2    0.947 0.939 0.954 0.00

## Fitting SMA model for method comparison

### Fitting SMA model
fit_SMA_models <- sma(Litter_mass_log10 ~ Adult_mass_log10, data = MI_models, method = "SMA")

### Checking SMA model assumptions
plot(fit_SMA_models, which = "default") ## diagnostics (good!)
plot(fit_SMA_models, which = "residual") ## diagnostics (good!)
plot(fit_SMA_models, which = "qq") ## diagnostics (ok)

### Computing CI for estimates in SMA model
pretty(extract_fit_summary(fit_SMA_models))
#              estimate  lower  upper
# intercept      -0.190 -0.210 -0.170
# 10^intercept    0.645  0.616  0.675
# slope           0.799  0.786  0.813

### Computing R2 in SMA model
compure_r2(fit_SMA_models)
#    estimate lower upper    p
# r2    0.987 0.984 0.988 0.00

## Fitting MA model for method comparison

### Fitting MA model
fit_MA_models <- sma(Litter_mass_log10 ~ Adult_mass_log10, data = MI_models, method = "MA")

### Checking MA model assumptions
plot(fit_MA_models, which = "default") ## diagnostics (good!)
plot(fit_MA_models, which = "residual") ## diagnostics (good!)
plot(fit_MA_models, which = "qq") ## diagnostics (ok)

### Computing CI for estimates in MA model
pretty(extract_fit_summary(fit_MA_models))
#              estimate  lower  upper
# intercept      -0.192 -0.212 -0.172
# 10^intercept    0.643  0.614  0.673
# slope           0.795  0.781  0.808

### Computing R2 in MA model
compure_r2(fit_MA_models)
#    estimate lower upper    p
# r2    0.980 0.977 0.983 0.00


## Fitting MSLR model for method comparison

### Fitting MSLR model
fit_MSLR_models <- fitme(Litter_mass_log10 ~ Adult_mass_log10 + Investment_duration_log10, data = MI_models)

### Checking MSLR model assumptions
plot(fit_MSLR_models, ask = FALSE, which = "mean")    ## diagnostics (good!)
plot(fit_MSLR_models, ask = FALSE, which = "predict") ## diagnostics (good!)

### Computing CI for estimates in MSLR model
pretty(extract_fit_summary(fit_MSLR_models))
#                           estimate lower_asymptotic upper_asymptotic
# (Intercept)                  0.664            0.508            0.819
# Adult_mass_log10             0.867            0.847            0.887
# Investment_duration_log10   -0.398           -0.470           -0.327
# 10^(Intercept)                4.61             3.22             6.60

### Computing R2 in MSLR model
compure_r2(fit_MSLR_models)
#    estimate lower upper    p
# r2    0.954 0.947 0.960 0.00


## Fitting MPLMM model for method comparison

### Fitting MPLMM model with estimation of best Pagel's Lambda (version 1: very slow)
if (run_slow) {
  fit_MPLMM_models <- fitme_phylo_lambdafree(
    data = MI_models, tree = tree,
    args_spaMM = list(formula = Litter_mass_log10 ~ Adult_mass_log10 + Investment_duration_log10 + corrMatrix(1|Key),
                      resid.model =  ~ Adult_mass_log10 + (1|Key)))
}

### Fitting MPLMM model without estimation of best Pagel's Lambda (version 2: much faster)
### Note: since best Pagel's Lambda is 1 on these data, the same output can be quickly obtained as follows
fit_MPLMM_models <- fitme_phylo_lambdafixed(
  lambda = 1, data = MI_models, tree = tree,
  args_spaMM = list(formula = Litter_mass_log10 ~ Adult_mass_log10 + Investment_duration_log10 + corrMatrix(1|Key),
                    resid.model =  ~ Adult_mass_log10 + (1|Key)))

### Profile for Pagel's lamba
if (run_slow) {
  profile_lambda_MPLMM <- profile_lambda(fit_MPLMM_models)
  plot(logLik ~ Pagel_lambda, data = profile_lambda_MPLMM, type = "o")
}

### Checking MPLMM model assumptions
plot(fit_MPLMM_models, ask = FALSE, which = "mean")  ## diagnostics (heteroscedastic, but this is accounted for)
plot(fit_MPLMM_models, ask = FALSE, which = "ranef") ## diagnostics (ok)
plot(fit_MPLMM_models, ask = FALSE, which = "predict") ## diagnostics (bad: residual variance captured by random variance)
plot(MI_models$Litter_mass_log10, predict(fit_MPLMM_models, re.form = NA, type = "link")[, 1]) ## diagnostics, excluding ranef (good!)

### Computing CI for estimates in MPLMM model (very computationally intensive)
if (run_slow) {
  MPLMM_summary <- extract_fit_summary(fit_MPLMM_models)
  pretty(MPLMM_summary$fixef) # Note: we report the basic intervals in the MS
#                            estimate lower_normal upper_normal lower_percent upper_percent lower_basic upper_basic
#  (Intercept)                  0.204       -0.225        0.625        -0.221         0.663      -0.255       0.628
#  Adult_mass_log10             0.836        0.805        0.868         0.802         0.868       0.804       0.870
#  Investment_duration_log10   -0.199       -0.326      -0.0701        -0.328       -0.0638      -0.334     -0.0701
#  10^(Intercept)                1.60        0.596         4.21         0.602          4.60       0.556        4.25
  
  pretty(MPLMM_summary$Pagel_Lambda) # Pagel's lambda
# estimate    lower    upper 
#   "1.00"  "0.999"   "1.00" 
}

### Computing R2 in MPLMM model
compure_r2(fit_MPLMM_models)
#    estimate lower upper    p
# r2    0.952 0.945 0.959 0.00

# Figure 1 ----------------------------------------------------------------

draw_figure_1(data_models = MI_models,
              fit_SLR = fit_SLR_models, fit_PLMM = fit_PLMM_models,
              fit_SMA = fit_SMA_models, fit_MA = fit_MA_models,
              fit_MSLR = fit_MSLR_models, fit_MPLMM = fit_MPLMM_models)
ggplot2::ggsave(filename = "figures/Fig1.pdf", scale = 1.2, width = 15, height = 10, units = "cm")
ggplot2::ggsave(filename = "figures/Fig1.png", scale = 1.2, width = 15, height = 10, units = "cm")


# Models comparison -------------------------------------------------------

## Computation of MI metrics
MI_models$MI_SLR   <- residuals(fit_SLR_models)
                      # same as: log(MI_models$Litter_mass / ((10^fixef(fit_SLR_models)["(Intercept)"]) * MI_models$Adult_mass^fixef(fit_SLR_models)["log(Adult_mass, 10)"]), base = 10)
MI_models$MI_PLMM  <- MI_models$Litter_mass_log10 - predict(fit_PLMM_models, re.form = NA, type = "link")[, 1]
                      # same as: log(MI_models$Litter_mass / ((10^fixef(fit_PLMM_models)["(Intercept)"]) * MI_models$Adult_mass^fixef(fit_PLMM_models)["log(Adult_mass, 10)"]), base = 10)
MI_models$MI_SMA   <- residuals(fit_SMA_models)
                      # same as: log(MI_models$Litter_mass / ((10^coef(fit_SMA_models) ["elevation"]) * MI_models$Adult_mass^coef(fit_SMA_models)["slope"]), base = 10)
MI_models$MI_MA    <- residuals(fit_MA_models)
                      # same as: log(MI_models$Litter_mass / ((10^coef(fit_MA_models) ["elevation"]) * MI_models$Adult_mass^coef(fit_MA_models)["slope"]), base = 10)
MI_models$MI_MSLR  <- residuals(fit_MSLR_models)
                      # same as: log(MI_models$Litter_mass / ((10^fixef(fit_MSLR_models)["(Intercept)"]) * MI_models$Adult_mass^fixef(fit_MSLR_models)["log(Adult_mass, 10)"] * MI_models$Investment_duration^fixef(fit_MSLR_models)["log(Investment_duration, 10)"]), base = 10)
MI_models$MI_MPLMM <- MI_models$Litter_mass_log10 - predict(fit_MPLMM_models, re.form = NA, type = "link")[, 1]
                      # same as: log(MI_models$Litter_mass / ((10^fixef(fit_MPLMM_models)["(Intercept)"]) * MI_models$Adult_mass^fixef(fit_MPLMM_models)["log(Adult_mass, 10)"] * MI_models$Investment_duration^fixef(fit_MPLMM_models)["log(Investment_duration, 10)"]), base = 10)

## Comparison of non phylogenetic methods
corMI <- cor(MI_models[, c("MI_SLR", "MI_SMA", "MI_MA", "MI_MSLR")])
diag(corMI) <- NA
corMI
pretty(range(corMI, na.rm = TRUE)) # Range of correlation coefficients between models
# [1] "0.921" "1.00" 

quade.test(as.matrix(MI_models[, c("MI_SLR", "MI_SMA", "MI_MA", "MI_MSLR")]))
# data:  as.matrix(MI_models[, c("MI_SLR", "MI_SMA", "MI_MA", "MI_MSLR")])
# Quade F = 0.64218, num df = 3, denom df = 2211, p-value = 0.5879

## Comparison between phylogenetic and non-phylogenetic counterpart
corMI <- cor(MI_models[, c("MI_PLMM", "MI_MPLMM", "MI_SLR", "MI_SMA", "MI_MA", "MI_MSLR")])
diag(corMI) <- NA
corMI
pretty(range(corMI, na.rm = TRUE)) # Range of correlation coefficients between models
# [1] "0.917" "1.00" 

## Comparison of models by LRT
if (run_slow) {
  univariate_phylo_test <- compute_LRT(fit = fit_PLMM_models, fit_null = fit_SLR_models)
  # ======== Bootstrap: ========
  #   Raw simulated p-value: 0.000999
  univariate_phylo_test$basicLRT
  # chi2_LR df p_value
  # p_v 932.6449 NA      NA
  1 + compute_df(fit_PLMM_models) - compute_df(fit_SLR_models) # dfs include + 1 for Pagel
  # [1] 3
}

if (run_slow) {
  multivariate_phylo_test <- compute_LRT(fit = fit_MPLMM_models, fit_null = fit_MSLR_models)
  # ======== Bootstrap: ========
  #   Raw simulated p-value: 0.000999
  multivariate_phylo_test$basicLRT
  # chi2_LR df p_value
  # p_v   823.7 NA      NA
  1 + compute_df(fit_MPLMM_models) - compute_df(fit_MSLR_models) # dfs include + 1 for Pagel
  # [1] 3
}


# Comparison of mass proxies ----------------------------------------------

fit_MPLMM_mass_default <- fitme_phylo_lambdafixed(
  lambda = 1, tree = tree, data = MI_mass, 
  args_spaMM = list(formula = Litter_mass_log10 ~ Adult_mass_log10 + Investment_duration_log10 + corrMatrix(1|Key),
                    resid.model =  ~ Adult_mass_log10 + (1|Key)))

fit_MPLMM_mass_females <- fitme_phylo_lambdafixed(
  lambda = 1, tree = tree, data = MI_mass,
  args_spaMM = list(formula = Litter_mass_log10 ~ Female_adult_mass_log10 + Investment_duration_log10 + corrMatrix(1|Key),
                    resid.model =  ~ Female_adult_mass_log10 + (1|Key))) ## resid model does not quite converge...

MI_mass$MI_default  <- MI_mass$Litter_mass_log10 - predict(fit_MPLMM_mass_default, re.form = NA, type = "link")[, 1]
MI_mass$MI_females  <- MI_mass$Litter_mass_log10 - predict(fit_MPLMM_mass_females, re.form = NA, type = "link")[, 1]

coin::wilcoxsign_test(MI_default ~ MI_females, data = MI_mass, distribution = "exact")
# data:  y by x (pos, neg) 
# stratified by block
# Z = 1.6608, p-value = 0.09704
# alternative hypothesis: true mu is not equal to 0

## Percentage of species for which difference between estimates is > 0.1
pretty(100*mean(abs(MI_mass$MI_females - MI_mass$MI_default) > 0.1))
# [1] "15.2"

## 0.1 expressed in SD of MI
pretty(mean(c(0.1/sd(MI_mass$MI_default),
              0.1/sd(MI_mass$MI_females))), digits = 2)
# [1] "0.40"


# Figure 2 ----------------------------------------------------------------

draw_figure_2(data_mass = MI_mass, fit_default = fit_MPLMM_mass_default, fit_females = fit_MPLMM_mass_females)
ggplot2::ggsave(filename = "figures/Fig2.pdf", scale = 1.2, width = 15, height = 10, units = "cm")
ggplot2::ggsave(filename = "figures/Fig2.png", scale = 1.2, width = 15, height = 10, units = "cm")


# Figure 3 ----------------------------------------------------------------

draw_figure_3(data_mass = MI_mass, fit_default = fit_MPLMM_mass_default, fit_females = fit_MPLMM_mass_females)
ggplot2::ggsave(filename = "figures/Fig3.pdf", scale = 1.2, width = 15, height = 10, units = "cm")
ggplot2::ggsave(filename = "figures/Fig3.png", scale = 1.2, width = 15, height = 10, units = "cm")


# Comparison of Subclasses ------------------------------------------------

MI_subclasses$MI  <- MI_subclasses$Litter_mass_log10 - predict(fit_MPLMM_models, newdata = MI_subclasses, re.form = NA, type = "link")[, 1]

## Figure 4 
draw_figure_4(MI_subclasses)
ggplot2::ggsave(filename = "figures/Fig4.pdf", scale = 1.2, width = 15, height = 10, units = "cm")
ggplot2::ggsave(filename = "figures/Fig4.png", scale = 1.2, width = 15, height = 10, units = "cm")

## Descriptive statistics
pretty(aggregate(MI ~ Subclass, data = MI_subclasses, \(x) c(mean = mean(x), sd = sd(x), N = length(x))) )
#      Subclass MI.mean MI.sd MI.N
# 1    Eutheria  0.0483 0.242 654.
# 2  Metatheria -0.0708 0.313 81.0
# 3 Monotremata  -0.356 0.697 3.00

## MI comparisons
kruskal.test(MI ~ Subclass, data = MI_subclasses)
# Asymptotic Kruskal-Wallis Test
# 
# data:  MI by Subclass (Eutheria, Metatheria, Monotremata)
# chi-squared = 16.218, df = 2, p-value = 0.0003008

coin::wilcox_test(MI ~ Subclass, data = droplevels(MI_subclasses[MI_subclasses$Subclass != "Monotremata", ]), distribution = "asymptotic")
# Asymptotic Wilcoxon-Mann-Whitney Test
# 
# data:  MI by Subclass (Eutheria, Metatheria)
# Z = 3.874, p-value = 0.0001071
# alternative hypothesis: true mu is not equal to 0

## Test using MPLMM

fit_MPLMM_subclass <- fitme_phylo_lambdafixed(
  lambda = 1, tree = tree, data = MI_subclasses[MI_subclasses$Subclass != "Monotremata", ], 
  args_spaMM = list(formula = Litter_mass_log10 ~ Adult_mass_log10 + Investment_duration_log10 + corrMatrix(1|Key),
                    resid.model =  ~ Adult_mass_log10 + (1|Key)))

fit_MPLMM_subclass_S <- fitme_phylo_lambdafixed(
  lambda = 1, tree = tree, data = MI_subclasses[MI_subclasses$Subclass != "Monotremata", ], 
  args_spaMM = list(formula = Litter_mass_log10 ~ Subclass + Adult_mass_log10 + Investment_duration_log10 + corrMatrix(1|Key),
                    resid.model =  ~ Adult_mass_log10 + (1|Key)))

fit_MPLMM_subclass_SM <- fitme_phylo_lambdafixed(
  lambda = 1, tree = tree, data = MI_subclasses[MI_subclasses$Subclass != "Monotremata", ], 
  args_spaMM = list(formula = Litter_mass_log10 ~ Subclass*Adult_mass_log10 + Investment_duration_log10 + corrMatrix(1|Key),
                    resid.model =  ~ Adult_mass_log10 + (1|Key),
                    control.HLfit = list(NbThreads = 2)))

fit_MPLMM_subclass_SD <- fitme_phylo_lambdafixed(
  lambda = 1, tree = tree, data = MI_subclasses[MI_subclasses$Subclass != "Monotremata", ], 
  args_spaMM = list(formula = Litter_mass_log10 ~ Adult_mass_log10 + Subclass*Investment_duration_log10 + corrMatrix(1|Key),
                    resid.model =  ~ Adult_mass_log10 + (1|Key),
                    control.HLfit = list(NbThreads = 2)))

fit_MPLMM_subclass_SMD <- fitme_phylo_lambdafixed(
  lambda = 1, tree = tree, data = MI_subclasses[MI_subclasses$Subclass != "Monotremata", ], 
  args_spaMM = list(formula = Litter_mass_log10 ~ Subclass*(Adult_mass_log10 + Investment_duration_log10) + corrMatrix(1|Key),
                    resid.model =  ~ Adult_mass_log10 + (1|Key),
                    control.HLfit = list(NbThreads = 2)))

if (run_slow) {
  subclass_test_SM.vs.S <- compute_LRT(fit_MPLMM_subclass_SM, fit_MPLMM_subclass_S)
  subclass_test_SM.vs.S
  # ======== Bootstrap: ========
  #  Raw simulated p-value: 0.023
  subclass_test_SM.vs.S$basicLRT
  # chi2_LR df p_value
  # p_v 7.216617 NA      NA
  compute_df(fit_MPLMM_subclass_SM) - compute_df(fit_MPLMM_subclass_S) # dfs
  # [1] 1
}

if (run_slow) {
  subclass_test_SD.vs.S <- compute_LRT(fit_MPLMM_subclass_SD, fit_MPLMM_subclass_S)
  subclass_test_SD.vs.S
  # ======== Bootstrap: ========
  #   Raw simulated p-value: 0.014
  subclass_test_SD.vs.S$basicLRT
  # chi2_LR df p_value
  # p_v 7.500569 NA      NA
  compute_df(fit_MPLMM_subclass_SD) - compute_df(fit_MPLMM_subclass_S) # dfs
  # [1] 1
}


fit_MPLMM_euth <- fitme_phylo_lambdafixed(
  lambda = 1, tree = tree, data = MI_subclasses[MI_subclasses$Subclass == "Eutheria", ], 
  args_spaMM = list(formula = Litter_mass_log10 ~ Adult_mass_log10 + Investment_duration_log10 + corrMatrix(1|Key),
                    resid.model =  ~ Adult_mass_log10 + (1|Key)))

if (run_slow) {
  MPLMM_euth_summary <- extract_fit_summary(fit_MPLMM_euth, lambdaCI = FALSE)
  pretty(MPLMM_euth_summary)
#                            estimate lower_normal upper_normal lower_percent upper_percent lower_basic upper_basic
#  (Intercept)                  0.382       -0.104        0.875       -0.0977         0.889      -0.124       0.862
#  Adult_mass_log10             0.848        0.816        0.881         0.816         0.880       0.816       0.880
#  Investment_duration_log10   -0.190       -0.314      -0.0653        -0.322       -0.0662      -0.313     -0.0575
#  10^(Intercept)                2.41        0.788         7.50         0.799          7.74       0.751        7.28
}

fit_MPLMM_meta <- fitme_phylo_lambdafixed(
  lambda = 1, tree = tree, data = MI_subclasses[MI_subclasses$Subclass == "Metatheria", ], 
  args_spaMM = list(formula = Litter_mass_log10 ~ Adult_mass_log10 + Investment_duration_log10 + corrMatrix(1|Key),
                    resid.model =  ~ Adult_mass_log10 + (1|Key)))

if (run_slow) {
  MPLMM_meta_summary <- extract_fit_summary(fit_MPLMM_meta, lambdaCI = FALSE)
  pretty(MPLMM_meta_summary)
#                           estimate lower_normal upper_normal lower_percent upper_percent lower_basic upper_basic
# (Intercept)                 -0.109        -1.33         1.03         -1.33          1.07       -1.29        1.11
# Adult_mass_log10             0.762        0.647        0.871         0.643         0.873       0.651       0.880
# Investment_duration_log10  -0.0806       -0.552        0.411        -0.555         0.411      -0.572       0.393
# 10^(Intercept)               0.778       0.0470         10.8        0.0470          11.9      0.0509        12.9
}

fit_MPLMM_euth_noD <- fitme_phylo_lambdafixed( ## noD as for no investment duration
  lambda = 1, tree = tree, data = MI_subclasses[MI_subclasses$Subclass == "Eutheria", ], 
  args_spaMM = list(formula = Litter_mass_log10 ~ Adult_mass_log10 + corrMatrix(1|Key),
                    resid.model =  ~ Adult_mass_log10 + (1|Key)))

if (run_slow) {
  MPLMM_euth_noD_summary <- extract_fit_summary(fit_MPLMM_euth_noD, lambdaCI = FALSE)
  pretty(MPLMM_euth_noD_summary)
#                  estimate lower_normal upper_normal lower_percent upper_percent lower_basic upper_basic
# (Intercept)        0.0430       -0.409        0.501        -0.418         0.500      -0.414       0.504
# Adult_mass_log10    0.820        0.793        0.848         0.792         0.846       0.793       0.847
# 10^(Intercept)       1.10        0.390         3.17         0.382          3.16       0.385        3.19
  
}

fit_MPLMM_meta_noD <- fitme_phylo_lambdafixed(
  lambda = 1, tree = tree, data = MI_subclasses[MI_subclasses$Subclass == "Metatheria", ], 
  args_spaMM = list(formula = Litter_mass_log10 ~ Adult_mass_log10 + corrMatrix(1|Key),
                    resid.model =  ~ Adult_mass_log10 + (1|Key)))

if (run_slow) {
  MPLMM_meta_noD_summary <- extract_fit_summary(fit_MPLMM_meta_noD, lambdaCI = FALSE)
  pretty(MPLMM_meta_noD_summary)
#                  estimate lower_normal upper_normal lower_percent upper_percent lower_basic upper_basic
# (Intercept)        -0.286       -0.870        0.269        -0.838         0.283      -0.855       0.265
# Adult_mass_log10    0.748        0.671        0.823         0.669         0.824       0.672       0.826
# 10^(Intercept)      0.517        0.135         1.86         0.145          1.92       0.140        1.84
}


# Comparison of Orders ---------------------------------------------------

MI_orders_euth$MI  <- MI_orders_euth$Litter_mass_log10 - predict(fit_MPLMM_euth, newdata = MI_orders_euth, re.form = NA, type = "link")[, 1]
MI_orders_meta$MI  <- MI_orders_meta$Litter_mass_log10 - predict(fit_MPLMM_meta, newdata = MI_orders_meta, re.form = NA, type = "link")[, 1]

## Figure 5
fig5A <- draw_figure_5(MI_orders_euth, tag = "A.", col_begin = 0, col_end = 0.7, scale = 0.9)
fig5B <- draw_figure_5(MI_orders_meta, dotsize = 0.8, tag = "B.", col_begin = 0.8, col_end = 0.9, scale = 0.2)
patchwork::wrap_plots(fig5A, fig5B, ncol = 2, widths = c(0.8, 0.3))
ggplot2::ggsave(filename = "figures/Fig5.pdf", scale = 1.2, width = 22, height = 10, units = "cm")
ggplot2::ggsave(filename = "figures/Fig5.png", scale = 1.2, width = 22, height = 10, units = "cm")

## Descriptive statistics
pretty(aggregate(MI ~ Order, data = MI_orders_euth, \(x) c(mean = mean(x), sd = sd(x), N = length(x))))
#             Order MI.mean MI.sd MI.N
# 1       Carnivora  -0.134 0.258 85.0
# 2 Cetartiodactyla -0.0823 0.203 61.0
# 3      Chiroptera  -0.318 0.126 71.0
# 4    Eulipotyphla   0.183 0.152 58.0
# 5      Lagomorpha  -0.128 0.244 22.0
# 6        Primates  -0.318 0.191 83.0
# 7        Rodentia  -0.144 0.201 252.

pretty(aggregate(MI ~ Order, data = MI_orders_meta, \(x) c(mean = mean(x), sd = sd(x), N = length(x))))
#            Order MI.mean MI.sd MI.N
# 1 Dasyuromorphia   0.247 0.281 23.0
# 2  Diprotodontia  -0.198 0.183 44.0

## MI comparisons
kruskal.test(MI ~ Order, data = MI_orders_euth)
# Asymptotic Kruskal-Wallis Test
# 
# data:  MI by
# Order (Carnivora, Cetartiodactyla, Chiroptera, Eulipotyphla, Lagomorpha, Primates, Rodentia)
# chi-squared = 196.35, df = 6, p-value < 2.2e-16

coin::wilcox_test(MI ~ Order, data = MI_orders_meta, distribution = "exact")
# Exact Wilcoxon-Mann-Whitney Test
# 
# data:  MI by Order (Dasyuromorphia, Diprotodontia)
# Z = 5.4803, p-value = 1.481e-09

## fit models for LRT
fit_MPLMM_subclass_O_euth <- fitme_phylo_lambdafixed(
  lambda = 1, tree = tree, data = MI_orders_euth, 
  args_spaMM = list(formula = Litter_mass_log10 ~ Order + Adult_mass_log10 + Investment_duration_log10 + corrMatrix(1|Key),
                    resid.model =  ~ Adult_mass_log10 + (1|Key)))

fit_MPLMM_subclass_OM_euth <- fitme_phylo_lambdafixed(
  lambda = 1, tree = tree, data = MI_orders_euth, 
  args_spaMM = list(formula = Litter_mass_log10 ~ Order*Adult_mass_log10 + Investment_duration_log10 + corrMatrix(1|Key),
                    resid.model =  ~ Adult_mass_log10 + (1|Key),
                    control.HLfit = list(NbThreads = 2)))

fit_MPLMM_subclass_OD_euth <- fitme_phylo_lambdafixed(
  lambda = 1, tree = tree, data = MI_orders_euth, 
  args_spaMM = list(formula = Litter_mass_log10 ~ Adult_mass_log10 + Order*Investment_duration_log10 + corrMatrix(1|Key),
                    resid.model =  ~ Adult_mass_log10 + (1|Key),
                    control.HLfit = list(NbThreads = 2)))

fit_MPLMM_subclass_O_meta <- fitme_phylo_lambdafixed(
  lambda = 1, tree = tree, data = MI_orders_meta, 
  args_spaMM = list(formula = Litter_mass_log10 ~ Order + Adult_mass_log10 + Investment_duration_log10 + corrMatrix(1|Key),
                    resid.model =  ~ Adult_mass_log10 + (1|Key)))

fit_MPLMM_subclass_OM_meta <- fitme_phylo_lambdafixed(
  lambda = 1, tree = tree, data = MI_orders_meta, 
  args_spaMM = list(formula = Litter_mass_log10 ~ Order*Adult_mass_log10 + Investment_duration_log10 + corrMatrix(1|Key),
                    resid.model =  ~ Adult_mass_log10 + (1|Key),
                    control.HLfit = list(NbThreads = 2)))

fit_MPLMM_subclass_OD_meta <- fitme_phylo_lambdafixed(
  lambda = 1, tree = tree, data = MI_orders_meta, 
  args_spaMM = list(formula = Litter_mass_log10 ~ Adult_mass_log10 + Order*Investment_duration_log10 + corrMatrix(1|Key),
                    resid.model =  ~ Adult_mass_log10 + (1|Key),
                    control.HLfit = list(NbThreads = 2)))

## LRTs
if (run_slow) {
  subclass_test_OM_euth <- compute_LRT(fit_MPLMM_subclass_OM_euth, fit_MPLMM_subclass_O_euth)
  subclass_test_OM_euth
  # ======== Bootstrap: ========
  # Raw simulated p-value: 0.0571
  subclass_test_OM_euth$basicLRT
  # chi2_LR df p_value
  # p_v 21.64892 NA      NA
  compute_df(fit_MPLMM_subclass_OM_euth) - compute_df(fit_MPLMM_subclass_O_euth) # dfs
  # [1] 6
}

if (run_slow) {
  subclass_test_OD_euth <- compute_LRT(fit_MPLMM_subclass_OD_euth, fit_MPLMM_subclass_O_euth)
  subclass_test_OD_euth
  # ======== Bootstrap: ========
  # Raw simulated p-value: 0.134
  subclass_test_OD_euth$basicLRT
  # chi2_LR df p_value
  # p_v 12.51678 NA      NA
  compute_df(fit_MPLMM_subclass_OD_euth) - compute_df(fit_MPLMM_subclass_O_euth) # dfs
  # [1] 6
}

if (run_slow) {
  subclass_test_OM_meta <- compute_LRT(fit_MPLMM_subclass_OM_meta, fit_MPLMM_subclass_O_meta)
  subclass_test_OM_meta
  # ======== Bootstrap: ========
  # Raw simulated p-value: 0.036
  subclass_test_OM_meta$basicLRT
  # chi2_LR df p_value
  # p_v 16.64098 NA      NA
  compute_df(fit_MPLMM_subclass_OM_meta) - compute_df(fit_MPLMM_subclass_O_meta) # dfs
  # [1] 1
}

if (run_slow) {
  subclass_test_OD_meta <- compute_LRT(fit_MPLMM_subclass_OD_meta, fit_MPLMM_subclass_O_meta)
  subclass_test_OD_meta
  # ======== Bootstrap: ========
  # Raw simulated p-value: 0.314
  subclass_test_OD_meta$basicLRT
  # chi2_LR df p_value
  # p_v 1.61142 NA      NA
  compute_df(fit_MPLMM_subclass_OD_meta) - compute_df(fit_MPLMM_subclass_O_meta) # dfs
  # [1] 1
}


# 20 Indicator Species --------------------------------------------------

MI_indicators$MI  <- MI_indicators$Litter_mass_log10 - predict(fit_MPLMM_models, newdata = MI_indicators, re.form = NA, type = "link")[, 1]

pretty(MI_indicators[order(-MI_indicators$MI), c("Name", "Species", "Order", "Subclass", "MI")])
#                              Name                Species           Order    Subclass      MI
# 2                 Tailless tenrec       Tenrec ecaudatus    Afrosoricida    Eutheria   0.803
# 998                Eurasian shrew          Sorex araneus    Eulipotyphla    Eutheria   0.625
# 61                        Red fox          Vulpes vulpes       Carnivora    Eutheria   0.458
# 150                    Blue whale  Balaenoptera musculus Cetartiodactyla    Eutheria   0.216
# 7                  American bison            Bison bison Cetartiodactyla    Eutheria   0.180
# 563         African bush elephant     Loxodonta africana     Proboscidea    Eutheria   0.134
# 896                           Rat          Rattus rattus        Rodentia    Eutheria  0.0915
# 54                      Red panda        Ailurus fulgens       Carnivora    Eutheria  0.0576
# 4                          Impala     Aepyceros melampus Cetartiodactyla    Eutheria  0.0523
# 339                Tammar wallaby       Macropus eugenii   Diprotodontia  Metatheria  0.0101
# 297               Tasmanian devil   Sarcophilus harrisii  Dasyuromorphia  Metatheria -0.0143
# 124                     Grey seal     Halichoerus grypus       Carnivora    Eutheria -0.0200
# 545                    Chimpanzee        Pan troglodytes        Primates    Eutheria -0.0270
# 449      Geoffroy's spider monkey       Ateles geoffroyi        Primates    Eutheria -0.0425
# 77                          Tiger        Panthera tigris       Carnivora    Eutheria -0.0628
# 398                 European hare        Lepus europaeus      Lagomorpha    Eutheria  -0.110
# 191 Greater short-nosed fruit bat      Cynopterus sphinx      Chiroptera    Eutheria  -0.134
# 345                  Red kangaroo         Macropus rufus   Diprotodontia  Metatheria  -0.304
# 385   Southern hairy-nosed wombat  Lasiorhinus latifrons   Diprotodontia  Metatheria  -0.319
# 421          Short-beaked echidna Tachyglossus aculeatus     Monotremata Monotremata  -0.464


# Figure 6 ----------------------------------------------------------------

if (run_slow) { ## slow to run as it downloads the silhouettes
  fig6 <- draw_figure_6(MI_indicators)
  ggplot2::ggsave(filename = "figures/Fig6.pdf", scale = 1.2, width = 15, height = 10, units = "cm")
  ggplot2::ggsave(filename = "figures/Fig6.png", scale = 1.2, width = 15, height = 10, units = "cm")
  fig6$Phylopic_artist
#  [1] "Becky Barnes"          "Sam Arman"             "T. Michael Keesey"     "Baheerathan Murugavel" "Ferran Sayol"         
#  [6] "Gabriela Palomo-Munoz" "Andy Wilson"           "Kai Caspar"            "Margot Michaud"        "Margot Michaud"       
#  [11] "Geoff Shaw"            "T. Michael Keesey"     "T. Michael Keesey"     "Ferran Sayol"          "T. Michael Keesey"    
#  [16] "Gabriela Palomo-Munoz" "T. Michael Keesey"     "Rebecca Groom"         "Becky Barnes"          "Yan Wong"  
}