
# Nice display ------------------------------------------------------------

pretty <- function(x, digits = 3) {
  names <- names(x)
  format <- paste0("%#.", digits, "g")
  res <- sprintf(format, signif(x, digits = digits))
  names(res) <- names
  res
}


# Check dependencies ------------------------------------------------------

## This function checks that all the package dependencies are met

check_dependencies_all <- function(pkgs) {
  
  check_dependencies_pkg <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste0("You need to install the package ", pkg, " using install.packages('", pkg, "') and rerun 'check_dependencies_all()'"))
    }
  }
  
  sapply(pkgs, \(pkg) {
    check_dependencies_pkg(pkg)
    message(paste("Package", pkg, "available. Version =", packageVersion(pkg), "found on your system."))
  })
  print("All dependencies have been met, you can continue 🎉")
}


# Prepare data frame with full data for Maternal Investment ---------------

## This function prepares the data containing information required to compute the Maternal Investment metrics

## It produces a dataset with 15 columns:
## - Species: the species name
## - Key: a key based on the taxonomy which is used to match the tips in the phylogenetic tree
## - Subclass: the mammalian subclass
## - Order: the mammalian order
## - Name: the vernacular name for the species
## - Adult_mass: the adult mass for the species, in kg
## - Adult_mass_log10: the log (base 10) adult mass for the species
## - Male_adult_mass: the feamle adult mass for the species, in kg
## - Male_adult_mass_log10: the log (base 10) female adult mass for the species
## - Female_adult_mass: the feamle adult mass for the species, in kg
## - Female_adult_mass_log10: the log (base 10) female adult mass for the species
## - Litter_mass: the mass of the litter at weaning age, in kg
## - Litter_mass_log10: the log (base 10) mass of the litter at weaning age
## - Investment_duration: the duration of maternal investment (gestation + lactation), in days
## - Investment_duration_log10: the log (base 10) duration of maternal investment (gestation + lactation)

prepare_df_MIfull <- function(raw_df) {
  
  ## Rename column
  raw_df$Key  <- raw_df$Name
  raw_df$Name <- raw_df$Common.name
  
  ## Compute derived columns
  ### Note: all masses are ultimately expressed in kg and duration in days
  raw_df$Species <- paste(raw_df$Genus, raw_df$Species)
  raw_df$Adult_mass   <- raw_df$Adult_mass_g/1000
  raw_df$Male_adult_mass   <- raw_df$Male_adult_mass_g/1000
  raw_df$Female_adult_mass   <- raw_df$Female_adult_mass_g/1000
  raw_df$Weaning_mass <- raw_df$Weaning_mass_g/1000
  raw_df$Litter_mass <- raw_df$Weaning_mass*raw_df$Litter.Clutch.size
  raw_df$Investment_duration <- raw_df$Gestation_days + raw_df$Lactation_days
  raw_df$Adult_mass_log10 <- log(raw_df$Adult_mass, base = 10)
  raw_df$Male_adult_mass_log10 <- log(raw_df$Male_adult_mass, base = 10)
  raw_df$Female_adult_mass_log10 <- log(raw_df$Female_adult_mass, base = 10)
  raw_df$Litter_mass_log10 <- log(raw_df$Litter_mass, base = 10)
  raw_df$Investment_duration_log10 <-  log(raw_df$Investment_duration, base = 10)
  
  ## Format character-columns into factors
  raw_df$Subclass <- as.factor(raw_df$Subclass)
  raw_df$Order    <- as.factor(raw_df$Order)
  raw_df$Key      <- as.factor(raw_df$Key)
  raw_df$Name     <- as.factor(raw_df$Name)
  
  ## Drop row for which critical information is missing
  raw_df <- raw_df[!is.na(raw_df$Litter_mass) & !is.na(raw_df$Adult_mass), ]
  
  ## Remove species for which the offspring are more than 15% heavier than the adult mass
  too_big <- raw_df$Weaning_mass_g > 1.15*raw_df$Adult_mass_g
    if (length(too_big) > 0) {
      message(paste("The following", sum(too_big),  "species have been discarded since the offspring are, at weaning age, more than 15% heavier than adults:"))
      message(paste(raw_df$Name[too_big], collapse = "\n"))
    }
  
  raw_df <- raw_df[!too_big, ]
  
  ## Remove species for which the adult default mass is more than 15% heavier than the heaviest adult sex
  too_big2 <- !is.na(raw_df$Male_adult_mass) & !is.na(raw_df$Female_adult_mass) & 1.15*pmax(raw_df$Male_adult_mass, raw_df$Female_adult_mass) < raw_df$Adult_mass
  
  if (length(too_big2) > 0) {
    message(paste("\nThe following", sum(too_big2),  "species have been discarded since the default adult mass was more than 15% heavier than the heaviest sex in adults:"))
    message(paste(raw_df$Name[too_big2], collapse = "\n"))
  }
  
  raw_df <- raw_df[!too_big2, ]
  
  ## Remove species for which the adult default mass is more than 15% lighter than the lighter adult sex
  too_small <- !is.na(raw_df$Male_adult_mass) & !is.na(raw_df$Female_adult_mass) & 0.85*pmin(raw_df$Male_adult_mass, raw_df$Female_adult_mass) > raw_df$Adult_mass
  
  if (length(too_small) > 0) {
    message(paste("\nThe following", sum(too_small),  "species have been discarded since the default adult mass was more than 15% lighter than the lightest sex in adults:"))
    message(paste(raw_df$Name[too_small], collapse = "\n"))
  }
  
  raw_df <- raw_df[!too_small, ]

  ## Reorder and select columns
  raw_df <- raw_df[, c("Species", "Key", "Subclass", "Order", "Name",
                       "Adult_mass", "Adult_mass_log10",
                       "Male_adult_mass", "Male_adult_mass_log10",
                       "Female_adult_mass", "Female_adult_mass_log10",
                       "Litter_mass", "Litter_mass_log10",
                       "Investment_duration", "Investment_duration_log10")]
  
  ## Final cleaning
  raw_df <- droplevels(raw_df)
  rownames(raw_df) <- NULL
  
  raw_df
}


# Fit P(G)LMM models ------------------------------------------------------

## This function fits a phylogenetic GLMM given a lambda value using spaMM

fitme_phylo_lambdafixed_old <- function(lambda = 1, tree, data, cor_fn = ape::corPagel, return.fit = TRUE, ...) {
  corM <- nlme::corMatrix(nlme::Initialize(cor_fn(lambda, form = ~ Key, phy = tree), data = data)) 
  fit <- spaMM::fitme(corrMatrix = corM, data = data, ...)
  fit$phylo <- list(lambda = lambda, tree = tree, cor_fn = cor_fn, corM = corM)
  if (return.fit) return(fit)
  logLik(fit)
}


fitme_phylo_lambdafixed <- function(lambda = 1, tree, data, cor_fn = ape::corPagel, return.fit = TRUE, args_spaMM = list()) {
  corM <- nlme::corMatrix(nlme::Initialize(cor_fn(lambda, form = ~ Key, phy = tree), data = data)) 
  args_spaMM <- c(args_spaMM, list(corrMatrix = corM, data = data))
  fit <- do.call(spaMM::fitme, args_spaMM)
  fit$phylo <- list(lambda = lambda, tree = tree, cor_fn = cor_fn, corM = corM)
  if (return.fit) return(fit)
  logLik(fit)
}

## This function fits a phylogenetic GLMM using spaMM and estimate the value of the parameter of the correlation structure (lambda)
## by outer estimation

fitme_phylo_lambdafree <- function(tree, data, cor_fn = ape::corPagel, args_spaMM = list()) {
  message("Fitting the P(G)LMM... be patient")
  best_lambda <- optimize(fitme_phylo_lambdafixed, interval = c(0, 1), maximum = TRUE,
                          tree = tree, data = data, return.fit = FALSE, args_spaMM = args_spaMM)
  message(paste("estimated lambda =", round(best_lambda$maximum, digits = 3)))
  fit <- fitme_phylo_lambdafixed(lambda = best_lambda$maximum, tree = tree, data = data,
                                 cor_fn = cor_fn, return.fit = TRUE, args_spaMM = args_spaMM)
  fit
}

## This function estimates the 95% CI of the parameter of the correlation structure (lambda)

confint_lambda <- function(bestfit) {
  lambda_ref <- bestfit$phylo$lambda
  CI_fit <- function(x, ...) {
    message(paste("optimise for Pagel's Lambda =", x))
    abs(fitme_phylo_lambdafixed(lambda = x, ...) - (logLik(bestfit) - qchisq(0.95, df = 1)/2)) ## asymptotic
  }
  
  if (is.null(bestfit$residModel)) {
      fixed <- list(phi = 1e-5)
    } else {
      fixed <- list()
      bestfit$residModel$family <- NULL # to avoid bug in spaMM 
    }

  message("Estimating upper boundary for lambda... be patient")
  if (lambda_ref == 1) {
    lambda_upr <- 1 ## no estimation needed since estimate already at parameter boundary
  } else {
  lambda_upr <- optimize(CI_fit, interval = c(lambda_ref, 1), tree = bestfit$phylo$tree, data = bestfit$data, cor_fn = bestfit$phylo$cor_fn,
                         args_spaMM = list(formula = formula(bestfit), resid.model = bestfit$residModel, fixed = fixed),
                         return.fit = FALSE)$minimum
  }
    
  message("Estimating lower boundary for lambda... be patient")
  if (lambda_ref == 0) {
    lambda_lwr <- 0 ## no estimation needed since estimate already at parameter boundary
  } else {
    lambda_lwr <- optimize(CI_fit, interval = c(0, lambda_ref), tree = bestfit$phylo$tree, data = bestfit$data, cor_fn = bestfit$phylo$cor_fn,
                           args_spaMM = list(formula = formula(bestfit), resid.model = bestfit$residModel, fixed = fixed),
                           return.fit = FALSE)$minimum
  }
  
  c(estimate = lambda_ref, lower = lambda_lwr, upper = lambda_upr)
}

## This function estimates the 95% CI for the fixed effects of models fitted with spaMM

confint_fixef <- function(fit, boot_args = NULL) {
  
  if (!is.null(boot_args) && boot_args$nb_cores > parallel::detectCores()) {
    stop("CI not computed under default settings; this is computationally challenging and therefore is best done on a large computer.")
  }
  
  message("Estimating CI for fixed effects... be patient")
  
  requireNamespace("doSNOW", quietly = TRUE)
  params <- names(spaMM::fixef(fit))
  res <- confint(fit, parm = params, verbose = FALSE, boot_args = boot_args, format = "stats")
  res <- lapply(res, \(block) {
    colnames(block) <- c("lower", "upper")
    block
    })
  res_wide <- do.call("cbind", res)
  colnames(res_wide) <- paste0(rep(names(res), each = 2), "_", colnames(res_wide))
  cbind(estimate = spaMM::fixef(fit), res_wide)
}


# Extract information from fits -------------------------------------------

## This functions extract parameters and their 95% confidence intervals

extract_fit_summary <- function(fit, digits = 3, boot_args = list(nb_cores = 50, nsim = 1000, seed = 123), lambdaCI = TRUE) {
  
  if (inherits(fit, what = "HLfit")) {
    
    if (!is.null(fit$phylo)) {
      corM <<- fit$phylo$corM ## to circumvent spaMM scoping issue in confint()
    } else {
      message("Simple models, so parametric bootstrap not used.")
      boot_args <- NULL
    }
    
    fixef_stats <- confint_fixef(fit, boot_args = boot_args)
    fixef_stats <- rbind(fixef_stats, "10^(Intercept)" = 10^fixef_stats["(Intercept)", ] )

    if (!is.null(fit$phylo) & lambdaCI) {
      lambda_stats <- confint_lambda(fit)
      stats <- rbind(fixef_stats, lambda_stats)
    } else {
      stats <- fixef_stats
    }
    
  } else if (inherits(fit, what = "sma")) {
    stats <- fit$coef[[1]]
    rownames(stats)[1] <- "intercept"
    stats <- rbind("10^intercept" = 10^stats[1, ], stats)
    stats <- stats[c("intercept", "10^intercept", "slope"),]
    colnames(stats) <- c("estimate", "lower", "upper")
  } else {
    stop("object class not recognized")
  }
  
  res <- do.call("data.frame", lapply(as.data.frame(stats), \(x) pretty(x, digits = digits)))
  rownames(res) <- rownames(stats)
  res
}


## This functions extract the r-squared value associated with the model
## Note: We compute it as the squared coefficient for the Pearson correlation between prediction and observations at the log scale
##       For mixed modes, we do not account for the random effect when computing the predictions

compure_r2 <- function(fit, digits = 3) {
  
  obs  <- fit$data$Litter_mass

  if (inherits(fit, what = "HLfit")) {
    obs <- log(obs, 10)
    pred <- predict(fit, re.form = NA)[, 1]
  } else if (inherits(fit, what = "sma")) {
    pred <- fitted(fit)
  } else {
    stop("object class not recognized")
  }
  
  cor_res <- cor.test(pred, obs)
  stats <- data.frame(estimate = pretty(cor_res$estimate^2, digits = digits),
                      lower = pretty(cor_res$conf.int[1][[1]], digits = digits),
                      upper = pretty(cor_res$conf.int[2][[1]], digits = digits),
                      p = pretty(cor_res$p.value, digits = digits))
  rownames(stats) <- "r2"
  stats
}



# Figures -----------------------------------------------------------------

## This function draws figure 1

draw_figure_1 <- function(data_models, fit_SLR, fit_PLMM, fit_SMA, fit_MA, fit_MSLR, fit_MPLMM) {
  
  data_pred <- data.frame(Adult_mass = c(0.001, 1e6),
                          Adult_mass_log10 = log(c(0.001, 1e6), base = 10),
                          Investment_duration_log10 = mean(data_models$Investment_duration_log10))
  data_pred$SLR   <- 10^(predict(fit_SLR, newdata = data_pred)[, 1])
  data_pred$PLMM  <- 10^(predict(fit_PLMM, newdata = data_pred, re.form = NA)[, 1])
  data_pred$SMA   <- 10^as.numeric(as.matrix(cbind(1, data_pred[, "Adult_mass_log10", drop = FALSE])) %*% matrix(coef(fit_SMA)))
  data_pred$MA    <- 10^as.numeric(as.matrix(cbind(1, data_pred[, "Adult_mass_log10", drop = FALSE])) %*% matrix(coef(fit_MA)))
  data_pred$MSLR  <- 10^(predict(fit_MSLR, newdata = data_pred)[, 1])
  data_pred$MPLMM <- 10^(predict(fit_MPLMM, newdata = data_pred, re.form = NA)[, 1])
  
  data_pred <- tidyr::pivot_longer(data_pred, cols = SLR:MPLMM, names_to = "Model", values_to = "Predict")
  data_pred$Model <- factor(data_pred$Model, levels = c("SLR", "PLMM", "SMA", "MA", "MSLR", "MPLMM"))
  
  fig <- ggplot2::ggplot(data = data_models, ggplot2::aes(Adult_mass, Litter_mass)) + 
    ggplot2::scale_x_continuous(trans = "log10", breaks = c(0.1, 1, 10, 100, 1000, 10000, 100000), 
                                labels = scales::number_format(accuracy = 0.1), expand = c(0, 0)) + 
    ggplot2::scale_y_continuous(trans = "log10",
                                breaks = c(0.1, 1, 10, 100, 1000, 10000),
                                labels = scales::number_format(accuracy = 0.1), expand = c(0, 0)) +
    ggplot2::scale_shape_manual(values = 21:23) +
    ggplot2::scale_fill_manual(values = c("steelblue", "darkred", "#FCC501")) +
    ggplot2::scale_color_viridis_d() +
    ggplot2::geom_point(ggplot2::aes(shape = Subclass, fill = Subclass), alpha = 0.3, size = 2) +
    ggplot2::geom_line(ggplot2::aes(y = Predict, x = Adult_mass, colour = Model), data = data_pred,
                       linewidth = 0.7, alpha = 0.8, inherit.aes = FALSE) +
    ggplot2::labs(x = 'Adult mass (kg)', y = 'Litter mass at weaning age (kg)') +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "right")
    
    print(fig)
}

## This function draws figure 2

draw_figure_2 <- function(data_mass, fit_default, fit_males, fit_females) {
  data_pred <- data.frame(Adult_mass = c(0.001, 1e4),
                          Adult_mass_log10 = log(c(0.001, 1e4), base = 10),
                          Male_adult_mass = c(0.001, 1e4),
                          Male_adult_mass_log10 = log(c(0.001, 1e4), base = 10),
                          Female_adult_mass = c(0.001, 1e4),
                          Female_adult_mass_log10 = log(c(0.001, 1e4), base = 10))
  
  data_pred$default <- 10^(predict(fit_default, newdata = data_pred, re.form = NA)[, 1])
  data_pred$males   <- 10^(predict(fit_males, newdata = data_pred, re.form = NA)[, 1])
  data_pred$females <- 10^(predict(fit_females, newdata = data_pred, re.form = NA)[, 1])
  
  data_pred <- tidyr::pivot_longer(data_pred, cols = default:females, names_to = "Sex", values_to = "Predict")
  data_pred$Sex <- factor(data_pred$Sex, levels = c("default", "males", "females"))
  
  fig <- ggplot2::ggplot(data = data_mass, ggplot2::aes(x = Adult_mass, y = Litter_mass)) + 
    ggplot2::scale_x_continuous(trans = "log10", breaks = c(0.1, 1, 10, 100, 1000, 10000, 100000), 
                                labels = scales::number_format(accuracy = 0.1), expand = c(0, 0)) + 
    ggplot2::scale_y_continuous(trans = "log10",
                                breaks = c(0.1, 1, 10, 100, 1000, 10000),
                                labels = scales::number_format(accuracy = 0.1), expand = c(0, 0)) +
    ggplot2::scale_shape_manual(values = 21:23) +
    ggplot2::scale_fill_manual(values = c("steelblue", "darkred", "#FCC501")) +
    ggplot2::scale_color_viridis_d() +
    ggplot2::geom_point(ggplot2::aes(shape = Subclass, fill = Subclass), alpha = 0.3, size = 2) +
    ggplot2::geom_line(ggplot2::aes(y = Predict, x = Adult_mass, colour = Sex), data = data_pred,
                       linewidth = 0.7, alpha = 0.8, inherit.aes = FALSE) +
    ggplot2::labs(x = 'Adult mass (kg)', y = 'Litter mass at weaning age (kg)', colour = 'Source for body mass') +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "right")
  
  print(fig)
}

draw_figure_2bis <- function(data_mass, fit_default, fit_default_full, fit_males, fit_females) {
  data_pred <- data.frame(Adult_mass = c(0.001, 1e4),
                          Adult_mass_log10 = log(c(0.001, 1e4), base = 10),
                          Male_adult_mass = c(0.001, 1e4),
                          Male_adult_mass_log10 = log(c(0.001, 1e4), base = 10),
                          Female_adult_mass = c(0.001, 1e4),
                          Female_adult_mass_log10 = log(c(0.001, 1e4), base = 10))
  
  data_pred$default <- 10^(predict(fit_default, newdata = data_pred, re.form = NA)[, 1])
  data_pred$`default (larger dataset)` <- 10^(predict(fit_default_full, newdata = data_pred, re.form = NA)[, 1])
  data_pred$males   <- 10^(predict(fit_males, newdata = data_pred, re.form = NA)[, 1])
  data_pred$females <- 10^(predict(fit_females, newdata = data_pred, re.form = NA)[, 1])

  data_pred <- tidyr::pivot_longer(data_pred, cols = default:females, names_to = "Sex", values_to = "Predict")
  data_pred$Sex <- factor(data_pred$Sex, levels = c("default (larger dataset)", "default", "males", "females"))
  
  fig <- ggplot2::ggplot(data = data_mass, ggplot2::aes(x = Adult_mass, y = Litter_mass)) + 
    ggplot2::scale_x_continuous(trans = "log10", breaks = c(0.1, 1, 10, 100, 1000, 10000, 100000), 
                                labels = scales::number_format(accuracy = 0.1), expand = c(0, 0)) + 
    ggplot2::scale_y_continuous(trans = "log10",
                                breaks = c(0.1, 1, 10, 100, 1000, 10000),
                                labels = scales::number_format(accuracy = 0.1), expand = c(0, 0)) +
    ggplot2::scale_shape_manual(values = 21:23) +
    ggplot2::scale_fill_manual(values = c("steelblue", "darkred", "#FCC501")) +
    ggplot2::scale_color_viridis_d() +
    ggplot2::geom_point(ggplot2::aes(shape = Subclass, fill = Subclass), alpha = 0.3, size = 2) +
    ggplot2::geom_line(ggplot2::aes(y = Predict, x = Adult_mass, colour = Sex), data = data_pred,
                       linewidth = 0.7, alpha = 0.8, inherit.aes = FALSE) +
    ggplot2::labs(x = 'Adult mass (kg)', y = 'Litter mass at weaning age (kg)', colour = 'Source for body mass') +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "right")
  
  print(fig)
}

## This function draws figure x (not used)

draw_figure_x <- function(data_mass, fit_default, fit_males, fit_females) {

  data_mass$default <- data_mass$Litter_mass_log10 - predict(fit_default, re.form = NA, type = "link")[, 1]
  data_mass$males   <- data_mass$Litter_mass_log10 - predict(fit_males, re.form = NA, type = "link")[, 1]
  data_mass$females <- data_mass$Litter_mass_log10 - predict(fit_females, re.form = NA, type = "link")[, 1]
  
  data_pred <- tidyr::pivot_longer(data_mass, cols = default:females, names_to = "Sex", values_to = "Predict")
  data_pred$Sex <- factor(data_pred$Sex, levels = c("default", "males", "females"))
  
  fig <- ggplot2::ggplot(data = data_pred, ggplot2::aes(x = Adult_mass, y = Predict, shape = Subclass, fill = Sex)) + 
    ggplot2::geom_point(alpha = 0.8, size = 2) +
    ggplot2::geom_text(ggplot2::aes(y = Predict + 0.2, label = Name), data = data_pred[data_pred$Predict > 1, ]) +
    ggplot2::scale_x_continuous(trans = "log10", breaks = c(0.1, 1, 10, 100, 1000, 10000, 100000), 
                                labels = scales::number_format(accuracy = 0.1), expand = c(0.1, 0.1)) + 
    ggplot2::scale_y_continuous(breaks = seq(-0.5, 4, by = 0.5),
                                labels = scales::number_format(accuracy = 0.1)) +
    ggplot2::scale_shape_manual(values = 21:23) +
    ggplot2::scale_fill_manual(values = c("grey", "#1fc3aa", "#8624f5")) +
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(shape = 21))) +
    ggplot2::labs(x = 'Adult mass (kg)', y = 'Maternal investment metric (MI)', fill = 'Source for body mass') +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "right")
  
  print(fig)
}

## This function draws figure 3

draw_figure_3 <- function(data_mass, fit_default, fit_females) {
  
  data_mass$default <- data_mass$Litter_mass_log10 - predict(fit_default, re.form = NA, type = "link")[, 1]
  data_mass$females <- data_mass$Litter_mass_log10 - predict(fit_females, re.form = NA, type = "link")[, 1]
  
  fig <- ggplot2::ggplot(data = data_mass, ggplot2::aes(y = females, x = default, shape = Subclass, fill = Subclass)) + 
    ggplot2::geom_point(alpha = 0.8, size = 2) +
    ggplot2::geom_text(ggplot2::aes(y = females, x = default - 0.02, label = Name), hjust = 1, data = data_mass[(data_mass$females - data_mass$default) > 0.1, ], size = 2) +
    ggplot2::geom_text(ggplot2::aes(y = females, x = default + 0.02, label = Name), hjust = 0, data = data_mass[(data_mass$default - data_mass$females) > 0.1, ], size = 2) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggplot2::scale_x_continuous(breaks = seq(-0.5, 4, by = 0.5), 
                                labels = scales::number_format(accuracy = 0.1)) + 
    ggplot2::scale_y_continuous(breaks = seq(-0.5, 4, by = 0.5),
                                labels = scales::number_format(accuracy = 0.1)) +
    ggplot2::scale_shape_manual(values = 21:23) +
    ggplot2::scale_fill_manual(values = c("steelblue", "darkred", "#FCC501")) +
    ggplot2::coord_fixed(xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75)) +
    ggplot2::labs(y = 'Maternal investment predicted using female adult mass', x = 'Maternal investment predicted using default adult mass') +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "right")
  
  print(fig)
}


## This function draws figure xx

draw_figure_xx <- function(data_mass, fit_default, fit_females) {
  
  data_mass$default <- data_mass$Litter_mass_log10 - predict(fit_default, re.form = NA, type = "link")[, 1]
  data_mass$females <- data_mass$Litter_mass_log10 - predict(fit_females, re.form = NA, type = "link")[, 1]
  data_mass$average_mass <- 0.5*(data_mass$Male_adult_mass + data_mass$Female_adult_mass)
  data_mass$departure_avg <- abs(data_mass$Adult_mass - data_mass$average_mass)
  
  fig <- ggplot2::ggplot(data = data_mass, ggplot2::aes(y = females - default, x = departure_avg, shape = Subclass, fill = Subclass)) + 
    ggplot2::geom_point(alpha = 0.8, size = 2) +
    ggplot2::geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
    ggplot2::geom_text(ggplot2::aes(y = females - default + 0.05, label = Name), data = data_mass[(data_mass$females - data_mass$default) > 0.1, ], size = 2) +
    ggplot2::geom_text(ggplot2::aes(y = females - default - 0.05, label = Name), data = data_mass[(data_mass$default - data_mass$females) > 0.1, ], size = 2) +
    ggplot2::scale_x_continuous(trans = "log1p", breaks = c(0.1, 1, 10, 100, 1000, 10000, 100000), 
                                labels = scales::number_format(accuracy = 0.1)) + 
    ggplot2::scale_y_continuous(breaks = seq(-0.5, 4, by = 0.5),
                                labels = scales::number_format(accuracy = 0.1)) +
    ggplot2::scale_shape_manual(values = 21:23) +
    ggplot2::scale_fill_manual(values = c("steelblue", "darkred", "#FCC501")) +
    ggplot2::coord_fixed(ylim = c(-0.75, 0.75)) +
    ggplot2::labs(y = 'MI (females) - MI (default)', x = 'abs(Adult mass - Average adult mass across sexes)') +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "right")
  
  print(fig)
}



