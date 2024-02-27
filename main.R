# Load custom functions ---------------------------------------------------

source("functions.R")

# Checking dependencies ---------------------------------------------------
check_dependencies_all(c("ape", "coin", "doSNOW", "ggplot2", "nlme", "scales", "smatr", "spaMM", "tidyr"))


# Load dependencies -------------------------------------------------------
library(spaMM)
library(smatr)
library(tidyr)


# Data preparation --------------------------------------------------------

## Import the phylogenetic tree
tree <- ape::read.tree("data/RAxML_bipartitions.result_FIN4_raw_rooted_wBoots_4098mam1out_OK.newick")


## Preparation of the Maternal Investment dataset

### Load the raw full dataset with Maternal Investment data
MI_raw <- read.csv2("data/MI.csv", dec = ".", na.strings = "")

### Format the full dataset (see functions.R for details)
MI_full <- prepare_df_MIfull(MI_raw)
nrow(MI_full) # 1040

### Prepare subsample with for comparison between subclasses
MI_subclasses <- droplevels(MI_full[MI_full$Key %in% tree[["tip.label"]], ])
nrow(MI_subclasses) # 799
str(MI_subclasses)

### Prepare subsample with for comparison between orders
orders_vs_N    <- aggregate(MI_subclasses[, "Key", drop = FALSE], list(Order = MI_subclasses$Order), length)
orders_to_keep <- as.character(orders_vs_N$Order[orders_vs_N$Key >= 15])
MI_orders <- droplevels(MI_subclasses[MI_subclasses$Order %in% orders_to_keep, ])
sort(unique(MI_orders$Order))
# [1] Carnivora       Cetartiodactyla Chiroptera      Dasyuromorphia  Diprotodontia   Eulipotyphla   
# [7] Lagomorpha      Primates        Rodentia  
nrow(MI_orders) # 756
str(MI_orders)

### Prepare subsample with no missing data for modelling
MI_models <- droplevels(MI_subclasses[!is.na(MI_subclasses$Investment_duration), ])
nrow(MI_models) # 736
str(MI_models)

### Prepare subsample with no missing data for mass proxy comparison
MI_mass <- droplevels(MI_subclasses[!is.na(MI_subclasses$Male_adult_mass) & !is.na(MI_subclasses$Female_adult_mass), ])
nrow(MI_mass) # 108
str(MI_mass)

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
extract_fit_summary(fit_SLR_models)
#              estimate  lower  upper
# intercept      -0.195 -0.215 -0.175
# 10^intercept    0.638  0.610  0.668
# slope           0.779  0.765  0.792

### Computing R2 in SLR model
compure_r2(fit_SLR_models)
#    estimate lower upper    p
# r2    0.948 0.970 0.977 0.00


## Fitting PLMM model for method comparison

### Fitting PLMM model with estimation of best Pagel's Lambda (version 1: very slow)
if (FALSE) { # switch FALSE to TRUE to run
  fit_PLMM_models <- fitme_phylo_lambdafree(tree = tree, data = MI_models,
                                            args_spaMM = list(formula = Litter_mass_log10 ~ Adult_mass_log10 + corrMatrix(1|Key),
                                                              # fixed = list(phi = 1e-5) ## if no resid.model!
                                                              resid.model =  ~ Adult_mass_log10 + (1|Key)))
}

### Fitting PLMM model without estimation of best Pagel's Lambda (version 2: much faster)
### Note: since best Pagel's Lambda is 1 on these data, the same output can be quickly obtained as follows
fit_PLMM_models <- fitme_phylo_lambdafixed(lambda = 1, data = MI_models, tree = tree,
                                           args_spaMM = list(formula = Litter_mass_log10 ~ Adult_mass_log10 + corrMatrix(1|Key),
                                                             resid.model =  ~ Adult_mass_log10 + (1|Key)))

### Checking PLMM model assumptions
plot(fit_PLMM_models, ask = FALSE, which = "mean")  ## diagnostics (heteroscedastic, but this is accounted for)
plot(fit_PLMM_models, ask = FALSE, which = "ranef") ## diagnostics (ok)
plot(fit_PLMM_models, ask = FALSE, which = "predict") ## diagnostics (bad: residual variance partially captured by random variance)
plot(MI_models$Litter_mass_log10, predict(fit_PLMM_models, re.form = NA, type = "link")[, 1]) ## diagnostics, excluding ranef (good!)

### Computing CI for estimates in PLMM model (very computationally intensive)
PLMM_summary <- extract_fit_summary(fit_PLMM_models)

pretty(PLMM_summary$fixef)
#                  estimate lower_normal upper_normal lower_percent upper_percent lower_basic upper_basic
# (Intercept)        -0.215       -0.313       -0.106        -0.326        -0.122      -0.307      -0.103
# Adult_mass_log10    0.807        0.792        0.828         0.788         0.824       0.791       0.827
# 10^(Intercept)      0.610        0.486        0.783         0.472         0.755       0.493       0.789

pretty(PLMM_summary$Pagel_Lambda)
# estimate    lower    upper 
#   "1.00"  "0.999"   "1.00" 

### Computing R2 in PLMM model
compure_r2(fit_PLMM_models) ## same as above!
#    estimate lower upper    p
# r2    0.948 0.970 0.977 0.00


## Fitting SMA model for method comparison

### Fitting SMA model
fit_SMA_models <- sma(Litter_mass_log10 ~ Adult_mass_log10, data = MI_models, method = "SMA")

### Checking SMA model assumptions
plot(fit_SMA_models, which = "default") ## diagnostics (good!)
plot(fit_SMA_models, which = "residual") ## diagnostics (good!)
plot(fit_SMA_models, which = "qq") ## diagnostics (ok)

### Computing CI for estimates in SMA model
extract_fit_summary(fit_SMA_models)
#              estimate  lower  upper
# intercept      -0.189 -0.209 -0.170
# 10^intercept    0.646  0.618  0.677
# slope           0.800  0.787  0.813

### Computing R2 in SMA model
compure_r2(fit_SMA_models)
#    estimate lower upper    p
# r2    0.987 0.992 0.994 0.00


## Fitting MA model for method comparison

### Fitting MA model
fit_MA_models <- sma(Litter_mass_log10 ~ Adult_mass_log10, data = MI_models, method = "MA")

### Checking MA model assumptions
plot(fit_MA_models, which = "default") ## diagnostics (good!)
plot(fit_MA_models, which = "residual") ## diagnostics (good!)
plot(fit_MA_models, which = "qq") ## diagnostics (ok)

### Computing CI for estimates in MA model
extract_fit_summary(fit_MA_models)
#              estimate  lower  upper
# intercept      -0.191 -0.211 -0.171
# 10^intercept    0.645  0.616  0.675
# slope           0.795  0.781  0.809

### Computing R2 in MA model
compure_r2(fit_MA_models)
#    estimate lower upper    p
# r2    0.980 0.989 0.991 0.00


## Fitting MSLR model for method comparison

### Fitting MSLR model
fit_MSLR_models <- fitme(Litter_mass_log10 ~ Adult_mass_log10 + Investment_duration_log10, data = MI_models)

### Checking MSLR model assumptions
plot(fit_MSLR_models, ask = FALSE, which = "mean")    ## diagnostics (good!)
plot(fit_MSLR_models, ask = FALSE, which = "predict") ## diagnostics (good!)

### Computing CI for estimates in MSLR model
extract_fit_summary(fit_MSLR_models)
#              estimate   lower   upper
# intercept       0.672  0.519  0.826
# 10^intercept     4.70   3.30   6.70
# slope           0.869  0.849  0.889
# slope_InvDur   -0.401 -0.472 -0.331

### Computing R2 in MSLR model
compure_r2(fit_MSLR_models)
#    estimate lower upper    p
# r2    0.956 0.974 0.981 0.00


## Fitting MPLMM model for method comparison

### Fitting MPLMM model with estimation of best Pagel's Lambda (version 1: very slow)
if (FALSE) { # switch FALSE to TRUE to run
  fit_MPLMM_models <- fitme_phylo_lambdafree(data = MI_models, tree = tree,
                                             args_spaMM = list(formula = Litter_mass_log10 ~ Adult_mass_log10 + Investment_duration_log10 + corrMatrix(1|Key),
                                                               resid.model =  ~ Adult_mass_log10 + (1|Key)))
}

### Fitting MPLMM model without estimation of best Pagel's Lambda (version 2: much faster)
### Note: since best Pagel's Lambda is 1 on these data, the same output can be quickly obtained as follows
fit_MPLMM_models <- fitme_phylo_lambdafixed(lambda = 1, data = MI_models, tree = tree,
                                            args_spaMM = list(formula = Litter_mass_log10 ~ Adult_mass_log10 + Investment_duration_log10 + corrMatrix(1|Key),
                                                              resid.model =  ~ Adult_mass_log10 + (1|Key)))

### Checking MPLMM model assumptions
plot(fit_MPLMM_models, ask = FALSE, which = "mean")  ## diagnostics (heteroscedastic, but this is accounted for)
plot(fit_MPLMM_models, ask = FALSE, which = "ranef") ## diagnostics (ok)
plot(fit_MPLMM_models, ask = FALSE, which = "predict") ## diagnostics (bad: residual variance captured by random variance)
plot(MI_models$Litter_mass_log10, predict(fit_MPLMM_models, re.form = NA, type = "link")[, 1]) ## diagnostics, excluding ranef (good!)

### Computing CI for estimates in MPLMM model (very computationally intensive)
MPLMM_summary <- extract_fit_summary(fit_MPLMM_models)

pretty(MPLMM_summary$fixef)
#                           estimate lower_normal upper_normal lower_percent upper_percent lower_basic upper_basic
# (Intercept)                  0.173       -0.141        0.249        0.0958         0.484      -0.139       0.249
# Adult_mass_log10             0.837        0.809        0.852         0.821         0.865       0.808       0.852
# Investment_duration_log10   -0.195       -0.219      -0.0501        -0.340        -0.170      -0.220     -0.0500
# 10^(Intercept)                1.49        0.723         1.77          1.25          3.05       0.725        1.78

pretty(MPLMM_summary$Pagel_Lambda)
# estimate    lower    upper 
#   "1.00"  "0.999"   "1.00" 

### Computing R2 in MPLMM model
compure_r2(fit_MPLMM_models)
#    estimate lower upper    p
# r2    0.953 0.973 0.980 0.00

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
range(corMI, na.rm = TRUE)
# [1] 0.9187744 0.9996684

quade.test(as.matrix(MI_models[, c("MI_SLR", "MI_SMA", "MI_MA", "MI_MSLR")]))
# data:  as.matrix(MI_models[, c("MI_SLR", "MI_SMA", "MI_MA", "MI_MSLR")])
# Quade F = 0.66567, num df = 3, denom df = 2205, p-value = 0.5731

## Comparison between phylogenetic and non-phylogenetic counterpart
corMI <- cor(MI_models[, c("MI_PLMM", "MI_MPLMM", "MI_SLR", "MI_SMA", "MI_MA", "MI_MSLR")])
diag(corMI) <- NA
corMI
range(corMI, na.rm = TRUE)
# [1] 0.9134880 0.9996684

univariate_phylo_test <- data.frame(LRT = unname(-2*(logLik(fit_SLR_models) - logLik(fit_PLMM_models))))
univariate_phylo_test$df <- 1 ## that should be more than 1 since there are df for the residual model...
univariate_phylo_test$p <- with(univariate_phylo_test, pchisq(LRT, df, lower.tail = FALSE))
univariate_phylo_test
#        LRT df           p
# 1 932.2787  1 9.439609e-205 # OUTDATED and WRONG

multivariate_phylo_test <- data.frame(LRT = unname(-2*(logLik(fit_MSLR_models) - logLik(fit_MPLMM_models))))
multivariate_phylo_test$df <- 1 ## that should be more than 1 since there are df for the residual model...
multivariate_phylo_test$p <- with(multivariate_phylo_test, pchisq(LRT, df, lower.tail = FALSE))
multivariate_phylo_test
#       LRT df             p
# 1 810.2442  1 3.19747e-178 # OUTDATED and WRONG

## Comparison between the 2 PLMMs
multivariatePLMM_phylo_test <- data.frame(LRT = unname(-2*(logLik(fit_PLMM_models) - logLik(fit_MPLMM_models))))
multivariatePLMM_phylo_test$df <- 1
multivariatePLMM_phylo_test$p <- with(multivariatePLMM_phylo_test, pchisq(LRT, df, lower.tail = FALSE))
multivariatePLMM_phylo_test
#        LRT df p
# -0.6251432  1 1 ## OUTDATED and WRONG


# Comparison of mass proxies ----------------------------------------------

fit_PLMM_mass_default <- fitme_phylo_lambdafixed(lambda = 1, tree = tree, data = MI_mass, 
                                                 args_spaMM = list(formula = Litter_mass_log10 ~ Adult_mass_log10 + corrMatrix(1|Key),
                                                                   resid.model =  ~ Adult_mass_log10 + (1|Key)))

fit_PLMM_mass_males <- fitme_phylo_lambdafixed(lambda = 1, tree = tree, data = MI_mass,
                                               args_spaMM = list(formula = Litter_mass_log10 ~ Male_adult_mass_log10 + corrMatrix(1|Key),
                                                                 resid.model =  ~ Male_adult_mass_log10 + (1|Key)))

fit_PLMM_mass_females <- fitme_phylo_lambdafixed(lambda = 1, tree = tree, data = MI_mass,
                                                 args_spaMM = list(formula = Litter_mass_log10 ~ Female_adult_mass_log10 + corrMatrix(1|Key),
                                                                   resid.model =  ~ Female_adult_mass_log10 + (1|Key)))

MI_mass$MI_default_full  <- MI_mass$Litter_mass_log10 - predict(fit_PLMM_models, newdata = MI_mass, re.form = NA, type = "link")[, 1]

MI_mass$MI_default  <- MI_mass$Litter_mass_log10 - predict(fit_PLMM_mass_default, re.form = NA, type = "link")[, 1]
MI_mass$MI_females  <- MI_mass$Litter_mass_log10 - predict(fit_PLMM_mass_females, re.form = NA, type = "link")[, 1]
MI_mass$MI_males    <- MI_mass$Litter_mass_log10 - predict(fit_PLMM_mass_males, re.form = NA, type = "link")[, 1]

MI_mass_long <- as.data.frame(pivot_longer(MI_mass, cols = c("MI_default", "MI_females", "MI_males"), names_to = "Sex", values_to = "MI"))
MI_mass_long$Sex <- as.factor(MI_mass_long$Sex)
MI_mass_long$Name <- as.factor(MI_mass_long$Name)

quade.test(as.matrix(MI_mass[, c("MI_default", "MI_females", "MI_males")]))
# Quade F = 1.6138, num df = 2, denom df = 214, p-value = 0.2015

## Percentage of species for which difference between estimates is > 0.1
pretty(100*mean(abs(MI_mass$MI_females - MI_mass$MI_default) > 0.1))
# [1] "13.0"

## 0.1 expressed in SD of MI
pretty(mean(c(0.1/sd(MI_mass$MI_default),
              0.1/sd(MI_mass$MI_females),
              0.1/sd(MI_models$MI_PLMM))), digits = 2)
# [1] "0.35"


# Figure 2 ----------------------------------------------------------------

draw_figure_2(data_mass = MI_mass, fit_default = fit_PLMM_mass_default, fit_males = fit_PLMM_mass_males, fit_females = fit_PLMM_mass_females)
ggplot2::ggsave(filename = "figures/Fig2.pdf", scale = 1.2, width = 15, height = 10, units = "cm")
ggplot2::ggsave(filename = "figures/Fig2.png", scale = 1.2, width = 15, height = 10, units = "cm")


# Figure 3 ----------------------------------------------------------------

draw_figure_3(data_mass = MI_mass, fit_default = fit_PLMM_mass_default, fit_females = fit_PLMM_mass_females)
ggplot2::ggsave(filename = "figures/Fig3.pdf", scale = 1.2, width = 15, height = 10, units = "cm")
ggplot2::ggsave(filename = "figures/Fig3.png", scale = 1.2, width = 15, height = 10, units = "cm")

draw_figure_x(data_mass = MI_mass, fit_default = fit_PLMM_mass_default, fit_males = fit_PLMM_mass_males, fit_females = fit_PLMM_mass_females)
ggplot2::ggsave(filename = "figures/FigS1.pdf", scale = 1.2, width = 15, height = 10, units = "cm")
ggplot2::ggsave(filename = "figures/FigS1.png", scale = 1.2, width = 15, height = 10, units = "cm")

draw_figure_xx(data_mass = MI_mass, fit_default = fit_PLMM_mass_default, fit_females = fit_PLMM_mass_females)
ggplot2::ggsave(filename = "figures/FigS2.pdf", scale = 1.2, width = 15, height = 10, units = "cm")
ggplot2::ggsave(filename = "figures/FigS2.png", scale = 1.2, width = 15, height = 10, units = "cm")
