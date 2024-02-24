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
nrow(MI_full) # 1056

### Prepare subsample with for comparison between subclasses
MI_subclasses <- droplevels(MI_full[MI_full$Key %in% tree[["tip.label"]], ])
nrow(MI_subclasses) # 814
str(MI_subclasses)

### Prepare subsample with for comparison between orders
orders_vs_N    <- aggregate(MI_subclasses[, "Key", drop = FALSE], list(Order = MI_subclasses$Order), length)
orders_to_keep <- as.character(orders_vs_N$Order[orders_vs_N$Key >= 15])
MI_orders <- droplevels(MI_subclasses[MI_subclasses$Order %in% orders_to_keep, ])
sort(unique(MI_orders$Order))
# [1] Carnivora       Cetartiodactyla Chiroptera      Dasyuromorphia  Didelphimorphia Diprotodontia   Eulipotyphla    Lagomorpha     
# [9] Primates        Rodentia 
nrow(MI_orders) # 785
str(MI_orders)

### Prepare subsample with no missing data for modelling
MI_models <- droplevels(MI_subclasses[!is.na(MI_subclasses$Investment_duration), ])
nrow(MI_models) # 750
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

cor_global <- cor.test(MI_full$Adult_mass_log10, MI_full$Litter_mass_log10)
round(cor_global$estimate, digits = 3)[[1]] # correlation estimate
# [1] 0.966
cor_global$p.value # pvalue
# [1] 0


# Fitting regression models -----------------------------------------------

## Fitting SLR model for method comparison

fit_SLR_models <- fitme(Litter_mass_log10 ~ Adult_mass_log10, data = MI_models)
plot(fit_SLR_models, ask = FALSE, which = "mean")    ## diagnostics (good!)
plot(fit_SLR_models, ask = FALSE, which = "predict") ## diagnostics (good!)
extract_fit_summary(fit_SLR_models)
#              estimate  lower  upper
# intercept      -0.194 -0.214 -0.174
# 10^intercept    0.639  0.610  0.669
# slope           0.780  0.766  0.793
compure_r2(fit_SLR_models)
#    estimate lower upper p
# r2    0.945 0.968 0.976 0

## Fitting PLMM model for method comparison

fit_PLMM_models <- fitme_phylo_lambdafree(formula = Litter_mass_log10 ~ Adult_mass_log10 + corrMatrix(1|Key),
                                          #resid.model =  ~ Adult_mass_log10,
                                          data = MI_models, tree = tree)
plot(fit_PLMM_models, ask = FALSE, which = "mean")  ## diagnostics (okish)
plot(fit_PLMM_models, ask = FALSE, which = "ranef") ## diagnostics (okish)
plot(fit_PLMM_models, ask = FALSE, which = "predict") ## diagnostics (bad: residual variance captured by random variance)
plot(MI_models$Litter_mass_log10, predict(fit_PLMM_models, re.form = NA, type = "link")[, 1]) ## diagnostics, excluding ranef (good!)
extract_fit_summary(fit_PLMM_models)
#              estimate  lower upper
# intercept      -0.173 -0.710 0.364
# 10^intercept    0.671  0.195  2.31
# slope           0.791  0.760 0.822
# lambda          0.945  0.913 0.966
compure_r2(fit_PLMM_models) ## same as above!
#    estimate lower upper p
# r2    0.945 0.968 0.976 0

## Fitting SMA model for method comparison

fit_SMA_models <- sma(Litter_mass_log10 ~ Adult_mass_log10, data = MI_models, method = "SMA")
plot(fit_SMA_models, which = "default") ## diagnostics (good!)
plot(fit_SMA_models, which = "residual") ## diagnostics (good!)
plot(fit_SMA_models, which = "qq") ## diagnostics (ok)
extract_fit_summary(fit_SMA_models)
#              estimate  lower  upper
# intercept      -0.188 -0.208 -0.168
# 10^intercept    0.648  0.619  0.679
# slope           0.802  0.789  0.816
compure_r2(fit_SMA_models)
#    estimate lower upper p
# r2    0.986 0.992 0.994 0

## Fitting MA model for method comparison

fit_MA_models <- sma(Litter_mass_log10 ~ Adult_mass_log10, data = MI_models, method = "MA")
plot(fit_MA_models, which = "default") ## diagnostics (good!)
plot(fit_MA_models, which = "residual") ## diagnostics (good!)
plot(fit_MA_models, which = "qq") ## diagnostics (ok)
extract_fit_summary(fit_MA_models)
#              estimate  lower  upper
# intercept      -0.190 -0.210 -0.169
# 10^intercept    0.646  0.617  0.677
# slope           0.797  0.783  0.811
compure_r2(fit_MA_models)
#    estimate lower upper p
# r2    0.979 0.988 0.991 0

## Fitting MSLR model for method comparison

fit_MSLR_models <- fitme(Litter_mass_log10 ~ Adult_mass_log10 + Investment_duration_log10, data = MI_models)
plot(fit_MSLR_models, ask = FALSE, which = "mean")    ## diagnostics (good!)
plot(fit_MSLR_models, ask = FALSE, which = "predict") ## diagnostics (good!)
extract_fit_summary(fit_MSLR_models)
#              estimate   lower   upper
# intercept       0.642  0.485  0.800
# 10^intercept     4.39   3.06   6.30
# slope           0.865  0.845  0.886
# slope_InvDur   -0.387 -0.459 -0.315
compure_r2(fit_MSLR_models)
#    estimate lower upper p
# r2    0.952 0.972 0.979 0

## Fitting MPLMM model for method comparison
fit_MPLMM_models <- fitme_phylo_lambdafree(formula = Litter_mass_log10 ~ Adult_mass_log10 + Investment_duration_log10 + corrMatrix(1|Key),
                                           data = MI_models, tree = tree)
plot(fit_MPLMM_models, ask = FALSE, which = "mean")  ## diagnostics (okish)
plot(fit_MPLMM_models, ask = FALSE, which = "ranef") ## diagnostics (okish)
plot(fit_MPLMM_models, ask = FALSE, which = "predict") ## diagnostics (bad: residual variance captured by random variance)
plot(MI_models$Litter_mass_log10, predict(fit_MPLMM_models, re.form = NA, type = "link")[, 1]) ## diagnostics, excluding ranef (good!)
extract_fit_summary(fit_MPLMM_models)
#              estimate   lower upper
# intercept      -0.553  -1.18 0.0756
# 10^intercept    0.280 0.0657   1.19
# slope           0.763  0.726  0.801
# slope_InvDur    0.197 0.0439  0.351
# lambda          0.952  0.924  0.971
compure_r2(fit_MPLMM_models)
#    estimate lower upper p
# r2    0.936 0.963 0.972 0


# Figure 1 ----------------------------------------------------------------

draw_figure_1(data_models = MI_models,
              fit_SLR = fit_SLR_models, fit_PLMM = fit_PLMM_models,
              fit_SMA = fit_SMA_models, fit_MA = fit_MA_models,
              fit_MSLR = fit_MSLR_models, fit_MPLMM = fit_MPLMM_models)
ggplot2::ggsave(filename = "figures/Fig1.pdf", scale = 1)
ggplot2::ggsave(filename = "figures/Fig1.png", scale = 1)


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
# [1] 0.9268376 0.9996569

quade.test(as.matrix(MI_models[, c("MI_SLR", "MI_SMA", "MI_MA", "MI_MSLR")]))
# data:  as.matrix(MI_models[, c("MI_SLR", "MI_SMA", "MI_MA", "MI_MSLR")])
# Quade F = 0.87625, num df = 3, denom df = 2247, p-value = 0.4526

## Comparison between phylogenetic and non-phylogenetic counterpart
corMI <- cor(MI_models[, c("MI_PLMM", "MI_MPLMM", "MI_SLR", "MI_SMA", "MI_MA", "MI_MSLR")])
diag(corMI) <- NA
corMI
range(corMI, na.rm = TRUE)
# [1] 0.8568404 0.9996569

univariate_phylo_test <- data.frame(LRT = unname(-2*(logLik(fit_SLR_models) - logLik(fit_PLMM_models))))
univariate_phylo_test$df <- 1
univariate_phylo_test$p <- with(univariate_phylo_test, pchisq(LRT, df, lower.tail = FALSE))
univariate_phylo_test
#        LRT df             p
# 1 547.0137  1 5.618796e-121

multivariate_phylo_test <- data.frame(LRT = unname(-2*(logLik(fit_MSLR_models) - logLik(fit_MPLMM_models))))
multivariate_phylo_test$df <- 1
multivariate_phylo_test$p <- with(multivariate_phylo_test, pchisq(LRT, df, lower.tail = FALSE))
multivariate_phylo_test
#        LRT df             p
# 1 449.6248  1 8.704853e-100

## Comparison between the 2 PLMMs
anova(fit_PLMM_models, fit_MPLMM_models) ## asymptotic test
#      chi2_LR df    p_value
# p_v 6.027052  1 0.01408824

# spaMM.options(nb_cores = parallel::detectCores() - 1) # optional: for quicker output
# library(doSNOW) # optional: for quicker output
# anova(fit_PLMM_models, fit_MPLMM_models, boot.repl = 999, fit_env = list(corM = fit_PLMM_models$phylo$corM)) # very slow



