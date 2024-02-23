
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
  raw_df[, c("Key", "Subclass", "Order", "Name", "Adult_mass", "Litter_mass", "Investment_duration")]
}


# Extract information from fits -------------------------------------------

extract_fit_summary <- function(fit, digits = 3) {

  if ("HLfit" %in% class(fit)) {
    elev_stats <- c(estimate = fixef(fit)["(Intercept)"][[1]],
                    confint(fit, parm = "(Intercept)", verbose = FALSE)$interval)
    scale_stats <- c(estimate = fixef(fit)["log(Adult_mass)"][[1]],
                     confint(fit, parm = "log(Adult_mass)", verbose = FALSE)$interval)
    stats <- rbind(elev_stats, scale_stats)
    colnames(stats) <- c("estimate", "lower", "upper")
  }
  
  round(stats, digits = digits) 
}
