
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

## It produces a dataset with 7 columns:
## - Key: a key based on the taxonomy which is used to match the tips in the phylogenetic tree
## - Subclass: the mammalian subclass
## - Order: the mammalian order
## - Name: the vernacular name for the species
## - Adult_mass: the adult mass for the species, in kg
## - Litter_mass: the mass of the litter at weaning age, in kg
## - Investment_duration: the duration of maternal investment (gestation + lactation), in days

prepare_df_MIfull <- function(raw_df) {
  
  ## Rename column
  raw_df$Key  <- raw_df$Name
  raw_df$Name <- raw_df$Common.name
  
  ## Compute derived columns
  ### Note: all masses are ultimately expressed in kg and duration in days
  raw_df$Adult_mass   <- raw_df$Adult_mass_g/1000
  raw_df$Weaning_mass <- raw_df$Weaning_mass_g/1000
  raw_df$Litter_mass <- raw_df$Weaning_mass*raw_df$Litter.Clutch.size
  raw_df$Investment_duration <- raw_df$Gestation_days + raw_df$Lactation_days
  
  ## Format character-columns into factors
  raw_df$Subclass <- as.factor(raw_df$Subclass)
  raw_df$Order    <- as.factor(raw_df$Order)
  raw_df$Key      <- as.factor(raw_df$Key)
  raw_df$Name     <- as.factor(raw_df$Name)
  
  ## Remove columns not needed
  raw_df$Kingdom            <- NULL
  raw_df$Phylum             <- NULL
  raw_df$Class              <- NULL
  raw_df$Family             <- NULL
  raw_df$Genus              <- NULL
  raw_df$Species            <- NULL
  raw_df$OrderC             <- NULL
  raw_df$FamilyC            <- NULL
  raw_df$Common.name        <- NULL
  raw_df$Adult_mass_g       <- NULL
  raw_df$Weaning_mass_g     <- NULL
  raw_df$Weaning_mass       <- NULL
  raw_df$Gestation_days     <- NULL
  raw_df$Lactation_days     <- NULL
  raw_df$Litter.Clutch.size <- NULL
  
  ## Reorder columns
  raw_df <- raw_df[, c("Key", "Subclass", "Order", "Name", "Adult_mass", "Litter_mass", "Investment_duration")]
  
  ## Drop row for which critical information is missing
  raw_df[!is.na(MI_full$Adult_mass) & !is.na(MI_full$Litter_mass), ]
}


# Fit P(G)LMM models ------------------------------------------------------

## This function fits a phylogenetic GLMM given a lambda value using spaMM

fitme_phylo_lambdafixed <- function(lambda = 1, tree, data, cor_fn = ape::corPagel, return.fit = TRUE, ...) {
  corM <- nlme::corMatrix(nlme::Initialize(cor_fn(lambda, form = ~ Key, phy = tree), data = data)) 
  fit <- spaMM::fitme(corrMatrix = corM, data = data, ...)
  fit$phylo <- list(lambda = lambda, tree = tree, cor_fn = cor_fn, corM = corM)
  if (return.fit) return(fit)
  logLik(fit)
}

## This function fits a phylogenetic GLMM using spaMM and estimate the value of the parameter of the correlation structure (lambda)
## by outer estimation

fitme_phylo_lambdafree <- function(tree, data, cor_fn = ape::corPagel, ...) {
  message("Fitting the P(G)LMM... be patient")
  best_lambda <- optimize(fitme_phylo_lambdafixed, interval = c(0, 1), maximum = TRUE,
                          fixed = list(phi = 1e-5), tree = tree, data = data, return.fit = FALSE, ...)
  message(paste("estimated lambda =", round(best_lambda$maximum, digits = 3)))
  fit <- fitme_phylo_lambdafixed(lambda = best_lambda$maximum, tree = tree, data = data,
                                 cor_fn = cor_fn, return.fit = TRUE, fixed = list(phi = 1e-5), ...)
  fit
}

## This function estimates the 95% CI of the parameter of the correlation structure (lambda)

confint_lambda <- function(bestfit) {
  lambda_ref <- bestfit$phylo$lambda
  CI_fit <- function(x, ...) abs(fitme_phylo_lambdafixed(lambda = x, ...) - (logLik(bestfit) - qchisq(0.95, df = 1)/2)) ## asymptotic
  
  message("Estimating lower boundary for lambda... be patient")
  lambda_lwr <- optimize(CI_fit, interval = c(0, lambda_ref),
                         formula = formula(bestfit), tree = bestfit$phylo$tree, data = bestfit$data, cor_fn = bestfit$phylo$cor_fn,
                         fixed = list(phi = 1e-5), return.fit = FALSE)$minimum
  
  message("Estimating upper boundary for lambda... be patient")
  lambda_upr <- optimize(CI_fit, interval = c(lambda_ref, 1),
                         formula = formula(bestfit), tree = bestfit$phylo$tree, data = bestfit$data, cor_fn = bestfit$phylo$cor_fn,
                         fixed = list(phi = 1e-5), return.fit = FALSE)$minimum
  
  c(estimate = lambda_ref, lower = lambda_lwr, upper = lambda_upr)
}



# Extract information from fits -------------------------------------------

## This functions extract parameters and their 95% confidence intervals

extract_fit_summary <- function(fit, digits = 3) {

  if (inherits(fit, what = "HLfit")) {
    if (!is.null(fit$phylo)) corM <<- fit$phylo$corM ## to circumvent spaMM scoping issue in confint()
    intercept <- c(estimate = fixef(fit)["(Intercept)"][[1]],
                   confint(fit, parm = "(Intercept)", verbose = FALSE)$interval)
    intercept_transformed <- c(estimate = 10^fixef(fit)["(Intercept)"][[1]],
                               10^confint(fit, parm = "(Intercept)", verbose = FALSE)$interval)
    slope <- c(estimate = fixef(fit)["log(Adult_mass, 10)"][[1]],
              confint(fit, parm = "log(Adult_mass, 10)", verbose = FALSE)$interval)
    stats <- rbind(intercept, intercept_transformed, slope)
    if (grepl(pattern = ".*Investment_duration", x = as.character(formula(fit))[3])) {
      slope_InvDur <- c(estimate = fixef(fit)["log(Investment_duration, 10)"][[1]],
                       confint(fit, parm = "log(Investment_duration, 10)", verbose = FALSE)$interval)
      stats <- rbind(stats, slope_InvDur)
    }
    colnames(stats) <- c("estimate", "lower", "upper")
    if (!is.null(fit$phylo)) {
      lambda <- confint_lambda(fit)
      stats <- rbind(stats, lambda)
    }
  } else if (inherits(fit, what = "sma")) {
    stats <- fit$coef[[1]]
    rownames(stats)[1] <- "intercept"
    stats <- rbind(intercept_transformed = 10^stats[1, ], stats)
    stats <- stats[c("intercept", "intercept_transformed", "slope"),]
    colnames(stats) <- c("estimate", "lower", "upper")
  } else {
    stop("object class not recognized")
  }
  
  rownames(stats)[rownames(stats) == "intercept_transformed"] <- "10^intercept"
  rows <- rownames(stats)
  res <- do.call("data.frame", lapply(as.data.frame(stats), \(x) sprintf('%#.3g', signif(x, digits = digits))))
  rownames(res) <- rows
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
  stats <- data.frame(estimate = round(cor_res$estimate^2, digits = digits),
                      lower = round(cor_res$conf.int[1][[1]], digits = digits),
                      upper = round(cor_res$conf.int[2][[1]], digits = digits),
                      p = signif(cor_res$p.value, digits = digits))
  rownames(stats) <- "r2"
  stats
}



# Figures -----------------------------------------------------------------

## This function draws figure 1

draw_figure_1 <- function(data_models, fit_SLR, fit_PLMM, fit_SMA, fit_MA, fit_MSLR, fit_MPLMM) {
  
  data_pred <- data.frame(Adult_mass = c(0.001, 1e6),
                          Investment_duration = median(data_models$Investment_duration))
  data_pred$SLR   <- 10^(predict(fit_SLR, newdata = data_pred)[, 1])
  data_pred$PLMM  <- 10^(predict(fit_PLMM, newdata = data_pred, re.form = NA)[, 1])
  data_pred$SMA   <- 10^as.numeric(as.matrix(cbind(1, log(data_pred[, "Adult_mass", drop = FALSE], 10))) %*% matrix(coef(fit_SMA)))
  data_pred$MA    <- 10^as.numeric(as.matrix(cbind(1, log(data_pred[, "Adult_mass", drop = FALSE], 10))) %*% matrix(coef(fit_MA)))
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
    ggplot2::xlab('Adult mass (kg)') +
    ggplot2::ylab('Litter mass at weaning age (kg)') +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "right")
    
    print(fig)
  }
