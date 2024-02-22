##############################################################################
################# Maternal investment metric introduction ##################
##############################################################################

## Loading packages used here 
library(stringr)
library(dunn.test)
library(readxl)
library(ggstatsplot)
library(rphylopic)
library(tidyverse)
library(scales)
library(data.table)
library(gridExtra)
library(smatr)
library(ape)
library(phytools)
library(nlme)
library(visreg)

## Session info
sessionInfo() 
# R version 4.2.1 (2022-06-23)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Monterey 12.6.5

# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib

# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
# [1] data.table_1.14.8  scales_1.2.1       lubridate_1.9.2    forcats_1.0.0     
# [5] dplyr_1.1.2        purrr_1.0.1        readr_2.1.4        tidyr_1.3.0       
# [9] tibble_3.2.1       ggplot2_3.4.2      tidyverse_2.0.0    rphylopic_1.1.1   
# [13] ggstatsplot_0.11.1 dunn.test_1.3.5    stringr_1.5.0      readxl_1.4.2          

##############################################################################

## Making certain data numeric
dataset <- read_excel("~/Desktop/Research/Introduction MI final/MI allometry dataset.xlsx",
                      col_types = c("text", "text", "text", "text", "text", "text", "text", "text", "text", "numeric",
                                    "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric",
                                    "numeric", "numeric", "numeric", "text", "text", "text", "text", "numeric",
                                    "numeric", "numeric", "numeric", "numeric", "text", "text", "text"))

## Add calculated values to the dataset
# Masses in kg
dataset$Adult_mass_kg <- dataset$Adult_mass_g/1000
dataset$Weaning_mass_kg <- dataset$Weaning_mass_g/1000
dataset$Litter_weaning_mass <- dataset$Weaning_mass_kg*dataset$`Litter/Clutch size`

#Total investment duration
dataset$Investment_duration_days <- (dataset$Gestation_days+dataset$Lactation_days)

##############################################################################
### Allometric scaling of the maternal investment equation
## Mass correlation
#Correlation between Adult mass weaning litter mass
cor.test(dataset$Litter_weaning_mass, dataset$Adult_mass_kg, use = "complete.obs", method="pearson")

MA_fit <- sma(Litter_weaning_mass~Adult_mass_kg, dataset, log="xy", method="MA")
summary(MA_fit)
MA_coefs <- coef(MA_fit)

SMA_fit <- sma(Litter_weaning_mass~Adult_mass_kg, dataset, log="xy", method="SMA")
summary(SMA_fit)
SMA_coefs <- coef(SMA_fit)

ML_fit <- sma(Litter_weaning_mass~Adult_mass_kg, dataset, log="xy", method="OLS")
summary(ML_fit)
ML_coefs <- coef(ML_fit)

MLR_fit <- lm(log(Litter_weaning_mass)~log(Adult_mass_kg) + log(Investment_duration_days), data=dataset)
summary(MLR_fit)
MLR_coefs <- coef(MLR_fit)

# f <- 1.65*0.360*(dataset$Adult_mass_kg^0.804)
# f <- exp(coef(MLR_fit)["(Intercept)"])*(dataset$Adult_mass_kg^coef(MLR_fit)["log(Adult_mass_kg)"])**(251.08^coef(MLR_fit)["log(Investment_duration_days)"])

data_MLR_plot <- data.frame(Adult_mass_kg = c(0.001, 1e6),
                            Investment_duration_days = median(dataset$Investment_duration_days, na.rm = TRUE))

data_MLR_plot$Litter_weaning_mass_expected <- exp(predict(MLR_fit, newdata = data_MLR_plot))

tree_full <- ape::read.tree("/Users/timhuijsmans/RAxML_bipartitions.result_FIN4_raw_rooted_wBoots_4098mam1out_OK.newick")
data_clean <- as.data.frame(dataset[!is.na(dataset$Litter_weaning_mass) & !is.na(dataset$Adult_mass_kg) & dataset$Name %in% tree_full[["tip.label"]], ])
tree <- ape::drop.tip(tree_full, setdiff(tree_full[["tip.label"]], data_clean$Name))

phylofit <- function(lambda = 1, form, data, return.fit = FALSE, ...) {
  corM <- nlme::corMatrix(nlme::Initialize(ape::corPagel(lambda, form = ~ Name, phy = tree), data = data)) 
  fit <- spaMM::fitme(as.formula(form), corrMatrix = corM, fixed = list(phi=1e-5), data = data, ...)
  if (return.fit) {
    print(paste("best lambda = ", round(lambda, 4)))
    return(fit)
  }
  logLik(fit)
}

best_lambda <- optimize(phylofit, interval = c(0, 1), maximum = TRUE,
                        form = log(Litter_weaning_mass) ~ log(Adult_mass_kg) + corrMatrix(1|Name),
                        data = data_clean)

best_fit_phylo <- phylofit(best_lambda$maximum,
                           form = log(Litter_weaning_mass) ~ log(Adult_mass_kg) + corrMatrix(1|Name),
                           data = data_clean, return.fit = TRUE)

data_phylo_plot <- data_MLR_plot
data_phylo_plot$Litter_weaning_mass_expected <- exp(predict(best_fit_phylo, newdata = data_phylo_plot, re.form = NA)[,1])

#Plot of mass ratio
ggplot(data = data_clean, aes((Adult_mass_kg), Litter_weaning_mass, col=Subclass)) + 
  geom_point(shape=1) + 
  scale_x_continuous(trans = "log10", labels = number_format(accuracy=0.1), expand = c(0, 0)) + 
  scale_y_continuous(trans = "log10", labels = number_format(accuracy=0.1), , expand = c(0, 0)) +
  xlab('Adult mass (kg)') +
  ylab('Litter weaning mass (kg)') +
  theme_classic() +
  geom_abline(intercept = ML_coefs[1], slope = ML_coefs[2], color="blue4",  linewidth=0.7) +
  geom_abline(intercept = MA_coefs[1], slope = MA_coefs[2], color="black",  linewidth=0.7) +
  geom_abline(intercept = SMA_coefs[1], slope = SMA_coefs[2], color="orange3",  linewidth=0.7) +
  geom_line(aes(y = Litter_weaning_mass_expected, x = Adult_mass_kg), data = data_MLR_plot, color = "red",  linewidth=0.7) +
  geom_line(aes(y = Litter_weaning_mass_expected, x = Adult_mass_kg), data = data_phylo_plot, color = "turquoise",  linewidth=0.7) +
  theme(legend.position = "NULL") +
  ggtitle('Adult mass vs litter weaning mass') +
  theme(plot.title = element_text(size=12, face="bold", hjust=0.5)) +
  scale_color_manual(values = c("steelblue", "darkred", "#FCC501")) 

## Calculate litter weaning mass expected ##NOTE TIM: sorry, stil the old version
#Equation litter weaning mass expected
dataset$Litter_weaning_mass_expected <- exp(coef(ML_fit)["(Intercept)"])*(dataset$Adult_mass_kg^coef(ML_fit)["log(Adult_mass_kg)"])

dataset$Litter_weaning_mass_expected_2 <- 1.65*(dataset$Adult_mass_kg^0.804)*(dataset$Investment_duration_days^(-0.185))

dataset$Litter_weaning_mass_expected_3 <- 0.659*(dataset$Adult_mass_kg^0.786)

dataset$Litter_weaning_mass_expected_4 <- 0.659*(dataset$Adult_mass_kg^0.792)

#Add maternal investment to dataset
dataset$MI <- log(dataset$Litter_weaning_mass/dataset$Litter_weaning_mass_expected)

dataset$MI2 <- log(dataset$Litter_weaning_mass/dataset$Litter_weaning_mass_expected_2)

dataset$MI3 <- log(dataset$Litter_weaning_mass/dataset$Litter_weaning_mass_expected_3)

dataset$MI4 <- log(dataset$Litter_weaning_mass/dataset$Litter_weaning_mass_expected_4)

#Calculate  + sd
mean((dataset$MI), na.rm=TRUE)
#  = 8.913065e-05

sd((dataset$MI), na.rm=TRUE)
# sd = 0.6605445


### Testing for outliers in Maternal investment
ggplot(dataset, aes(x = factor(1), y = MI)) +
  geom_boxplot(width = 0.4, fill = "white") +
  labs(x = NULL) 

boxplot.stats(dataset$MI)$out
out <- boxplot.stats(dataset$MI)$out
out_ind <- which(dataset$MI %in% c(out))
out_ind
dataset[out_ind, ]
#All values seem valid. 

#Plot of adult mass vs maternal investment ## NOTE TIM: this is four times the same plot now, axis names, line of best fit and plot title still have to be adjusted
plot1 <- ggplot(data =dataset, aes((Adult_mass_kg), MI, col=Subclass)) + 
  geom_point(shape=1) + 
  scale_x_continuous(trans = "log10", labels = number_format(accuracy=0.1)) + 
  xlab('Adult mass (kg)') +
  ylab('Maternal investment') +
  theme_classic() +
  geom_smooth(method="lm", linewidth=0.4, color="black") +
  theme(legend.position = "NULL") +
  ggtitle('Adult mass vs maternal investment') +
  theme(plot.title = element_text(size=12, face="bold", hjust=0.5)) +
  scale_color_manual(values = c("steelblue", "darkred", "#FCC501")) 

plot2 <- ggplot(data =dataset, aes((Adult_mass_kg), MI2, col=Subclass)) + 
  geom_point(shape=1) + 
  scale_x_continuous(trans = "log10", labels = number_format(accuracy=0.1)) + 
  xlab('Adult mass (kg)') +
  ylab('Maternal investment') +
  theme_classic() +
  geom_smooth(method="lm", linewidth=0.4, color="black") +
  theme(legend.position = "NULL") +
  ggtitle('Adult mass vs maternal investment') +
  theme(plot.title = element_text(size=12, face="bold", hjust=0.5)) +
  scale_color_manual(values = c("steelblue", "darkred", "#FCC501")) 

plot3 <- ggplot(data =dataset, aes((Adult_mass_kg), MI3, col=Subclass)) + 
  geom_point(shape=1) + 
  scale_x_continuous(trans = "log10", labels = number_format(accuracy=0.1)) + 
  xlab('Adult mass (kg)') +
  ylab('Maternal investment') +
  theme_classic() +
  geom_smooth(method="lm", linewidth=0.4, color="black") +
  theme(legend.position = "NULL") +
  ggtitle('Adult mass vs maternal investment') +
  theme(plot.title = element_text(size=12, face="bold", hjust=0.5)) +
  scale_color_manual(values = c("steelblue", "darkred", "#FCC501")) 

plot4 <- ggplot(data =dataset, aes((Adult_mass_kg), MI4, col=Subclass)) + 
  geom_point(shape=1) + 
  scale_x_continuous(trans = "log10", labels = number_format(accuracy=0.1)) + 
  xlab('Adult mass (kg)') +
  ylab('Maternal investment') +
  theme_classic() +
  geom_smooth(method="lm", linewidth=0.4, color="black") +
  theme(legend.position = "NULL") +
  ggtitle('Adult mass vs maternal investment') +
  theme(plot.title = element_text(size=12, face="bold", hjust=0.5)) +
  scale_color_manual(values = c("steelblue", "darkred", "#FCC501"))  

plot5 <- ggplot(data =dataset, aes((Adult_mass_kg), MI5, col=Subclass)) + 
  geom_point(shape=1) + 
  scale_x_continuous(trans = "log10", labels = number_format(accuracy=0.1)) + 
  xlab('Adult mass (kg)') +
  ylab('Maternal investment') +
  theme_classic() +
  geom_smooth(method="lm", linewidth=0.4, color="black") +
  theme(legend.position = "NULL") +
  ggtitle('Adult mass vs maternal investment') +
  theme(plot.title = element_text(size=12, face="bold", hjust=0.5)) +
  scale_color_manual(values = c("steelblue", "darkred", "#FCC501")) 

plot6 <- ggplot(data =dataset, aes((Adult_mass_kg), MI6, col=Subclass)) + 
  geom_point(shape=1) + 
  scale_x_continuous(trans = "log10", labels = number_format(accuracy=0.1)) + 
  xlab('Adult mass (kg)') +
  ylab('Maternal investment') +
  theme_classic() +
  geom_smooth(method="lm", linewidth=0.4, color="black") +
  theme(legend.position = "NULL") +
  ggtitle('Adult mass vs maternal investment') +
  theme(plot.title = element_text(size=12, face="bold", hjust=0.5)) +
  scale_color_manual(values = c("steelblue", "darkred", "#FCC501"))  

# Combine 6 plots
(plot1 + plot2)/(plot3 + plot4)/(plot5 + plot6)


#Test for correlation between adult mass and maternal investment
cor.test(dataset$MI, dataset$Adult_mass_kg, use = "complete.obs", method="pearson")

##############################################################################
### Female mass
ds <- Amniote_Database_Aug_2015

ds$Weaning_mass_kg <- ds$weaning_weight_g/1000
ds$Female_adult_mass_kg <- ds$female_body_mass_g/1000
ds$Litter_weaning_mass_kg <- ds$Weaning_mass_kg*ds$litter_or_clutch_size_n

allometry5 <- lm(log(Litter_weaning_mass_kg)~log(Female_adult_mass_kg), data=ds)
summary(allometry5)
confint(allometry5)

ds$Litter_weaning_mass_expected_7 <- 0.606*(ds$Female_adult_mass_kg^0.752)

ds$MI7 <- log(ds$Litter_weaning_mass_kg/ds$Litter_weaning_mass_expected_7)

##############################################################################
### Applying metric on emperical data
#Plot of subclass vs Maternal investment
datasubclass <- dataset
datasubclass$Subclass <- factor(datasubclass$Subclass,     # Reorder factor levels
                                c("Monotremata", "Metatheria", "Eutheria"))
ggplot(data = datasubclass, aes(Subclass, MI, col=Subclass)) + 
  geom_boxplot() + 
  geom_jitter (shape=1) +
  ylab('Maternal investment') +
  xlab('Subclass') +
  theme_classic() +
  geom_smooth(method="lm", linewidth=0.5) +
  theme(legend.position = "NULL") +
  ggtitle('Maternal investment per subclass') +
  theme(plot.title = element_text(size=12, face="bold", hjust=0.5)) +
  scale_color_manual(values = c("NA", "darkred", "steelblue")) +
  geom_point(data = datasubclass[datasubclass$Subclass == "Monotremata",],
             mapping = aes(y = MI, x = 1), colour = "#FCC501", shape=1)

#Test for signigicant differences in maternal investment between subclasses
dunn.test(dataset$MI, dataset$Subclass, method = "bonferroni")

setDT(dataset)[, list(Subclass_mean = mean(MI, na.rm=TRUE), Subclass_sd = sd(MI, na.rm=TRUE)), 
               by = c("Subclass")]

wilcox.test(dataset$MI, dataset$subclass, paired=F)

#Test for signigicant differences in maternal investment between subclasses
dunn.test(dataset$MI2, dataset$Subclass, method = "bonferroni")

setDT(dataset)[, list(Subclass_mean = mean(MI2, na.rm=TRUE), Subclass_sd = sd(MI2, na.rm=TRUE)), 
               by = c("Subclass")]

wilcox.test(dataset$MI2, dataset$subclass, paired=F)

#Test for signigicant differences in maternal investment between subclasses
dunn.test(dataset$MI3, dataset$Subclass, method = "bonferroni")

setDT(dataset)[, list(Subclass_mean = mean(MI3, na.rm=TRUE), Subclass_sd = sd(MI3, na.rm=TRUE)), 
               by = c("Subclass")]

wilcox.test(dataset$MI3, dataset$subclass, paired=F)

#Test for signigicant differences in maternal investment between subclasses
dunn.test(dataset$MI4, dataset$Subclass, method = "bonferroni")

setDT(dataset)[, list(Subclass_mean = mean(MI4, na.rm=TRUE), Subclass_sd = sd(MI4, na.rm=TRUE)), 
               by = c("Subclass")]

wilcox.test(dataset$MI, dataset$subclass, paired=F)


#Allometric coefficient per order > 15 sample species
Artio <- subset(dataset, Order == "Artiodactyla")

Artio1 <- lm(log(Litter_weaning_mass)~log(Adult_mass_kg), data=Artio)
confint(Artio1)
summary(Artio1)

# log(Adult_mass_kg)  0.80046    0.04510  17.748 2.64e-15 ***
# log(Adult_mass_kg)  0.7073725 0.8935416
# Multiple R-squared:  0.9292,	Adjusted R-squared:  0.9263 
# F-statistic:   315 on 1 and 24 DF,  p-value: 2.644e-15

Carni <- subset(dataset, Order == "Carnivora")

Carni1 <- lm(log(Litter_weaning_mass)~log(Adult_mass_kg), data=Carni)
confint(Carni1)
summary(Carni1)

# log(Adult_mass_kg)  0.72940    0.02978   24.49   <2e-16 ***
# log(Adult_mass_kg)  0.6700131 0.78878850
# Multiple R-squared:  0.8941,	Adjusted R-squared:  0.8927 
# F-statistic: 599.7 on 1 and 71 DF,  p-value: < 2.2e-16

Rode <- subset(dataset, Order == "Rodentia")

Rode1 <- lm(log(Litter_weaning_mass)~log(Adult_mass_kg), data=Rode)
confint(Rode1)
summary(Rode1)

# log(Adult_mass_kg)  0.80100    0.02824  28.369  < 2e-16 ***
# log(Adult_mass_kg)  0.7450735  0.8569301
# Multiple R-squared:  0.875,	Adjusted R-squared:  0.8739 
# F-statistic: 804.8 on 1 and 115 DF,  p-value: < 2.2e-16

Euli <- subset(dataset, Order == "Eulipotyphla")

Euli1 <- lm(log(Litter_weaning_mass)~log(Adult_mass_kg), data=Euli)
confint(Euli1)
summary(Euli1)

# log(Adult_mass_kg)   0.6430     0.0612  10.506 1.38e-08 ***
# log(Adult_mass_kg)  0.5132718 0.7727650
# Multiple R-squared:  0.8734,	Adjusted R-squared:  0.8655 
# F-statistic: 110.4 on 1 and 16 DF,  p-value: 1.376e-08

Dasy <- subset(dataset, Order == "Dasyuromorphia")

Dasy1 <- lm(log(Litter_weaning_mass)~log(Adult_mass_kg), data=Dasy)
confint(Dasy1)
summary(Dasy1)

# log(Adult_mass_kg)  0.64843    0.07643   8.483 1.05e-07 ***
# log(Adult_mass_kg)  0.4878434 0.8090083
# Multiple R-squared:  0.7999,	Adjusted R-squared:  0.7888 
# F-statistic: 71.97 on 1 and 18 DF,  p-value: 1.052e-07

Dipro <- subset(dataset, Order == "Diprotodontia")

Dipro1 <- lm(log(Litter_weaning_mass)~log(Adult_mass_kg), data=Dipro)
confint(Dipro1)
summary(Dipro1)

# log(Adult_mass_kg)  0.71603    0.04315  16.596  < 2e-16 ***
# log(Adult_mass_kg)  0.628686  0.8033741
# Multiple R-squared:  0.8788,	Adjusted R-squared:  0.8756 
# F-statistic: 275.4 on 1 and 38 DF,  p-value: < 2.2e-16

Prima <- subset(dataset, Order == "Primates")

Prima1 <- lm(log(Litter_weaning_mass)~log(Adult_mass_kg), data=Prima)
confint(Prima1)
summary(Prima1)

# log(Adult_mass_kg)  0.69898    0.02922   23.92  < 2e-16 ***
# log(Adult_mass_kg)  0.6398763  0.7580801
# Multiple R-squared:  0.9362,	Adjusted R-squared:  0.9346 
# F-statistic: 572.2 on 1 and 39 DF,  p-value: < 2.2e-16

Eutheria <- subset(dataset, Subclass == "Eutheria")

Euthe <- lm(log(Litter_weaning_mass)~log(Adult_mass_kg), data=Eutheria)
confint(Euthe)
summary(Euthe)

# log(Adult_mass_kg)  0.78245    0.01014   77.13   <2e-16 ***
# log(Adult_mass_kg)  0.7624838  0.8024081
# Multiple R-squared:  0.9511,	Adjusted R-squared:  0.9509 
# F-statistic:  5949 on 1 and 306 DF,  p-value: < 2.2e-16

Metatheria <- subset(dataset, Subclass == "Metatheria")

Meta <- lm(log(Litter_weaning_mass)~log(Adult_mass_kg), data=Metatheria)
confint(Meta)
summary(Meta)

# log(Adult_mass_kg)  0.61550    0.03147  19.558  < 2e-16 ***
# log(Adult_mass_kg)  0.5527175  0.6782820
# Multiple R-squared:  0.8472,	Adjusted R-squared:  0.845 
# F-statistic: 382.5 on 1 and 69 DF,  p-value: < 2.2e-16

#########################
Artio2 <- lm(log(Litter_weaning_mass) ~ log(Adult_mass_kg) + log(Investment_duration_days), data=Artio)
confint(Artio2)
summary(Artio2)

Carni2 <- lm(log(Litter_weaning_mass) ~ log(Adult_mass_kg) + log(Investment_duration_days), data=Carni)
confint(Carni2)
summary(Carni2)

Rode2 <- lm(log(Litter_weaning_mass) ~ log(Adult_mass_kg) + log(Investment_duration_days), data=Rode)
confint(Rode2)
summary(Rode2)

Euli2 <- lm(log(Litter_weaning_mass) ~ log(Adult_mass_kg) + log(Investment_duration_days), data=Euli)
confint(Euli2)
summary(Euli2)

Dasy2 <- lm(log(Litter_weaning_mass) ~ log(Adult_mass_kg) + log(Investment_duration_days), data=Dasy)
confint(Dasy2)
summary(Dasy2)

Dipro2 <- lm(log(Litter_weaning_mass) ~ log(Adult_mass_kg) + log(Investment_duration_days), data=Dipro)
confint(Dipro2)
summary(Dipro2)

Prima2 <- lm(log(Litter_weaning_mass) ~ log(Adult_mass_kg) + log(Investment_duration_days), data=Prima)
confint(Prima2)
summary(Prima2)

#################

## Showing maternal investment of 20 example animals
#Set up dataset
part <- dataset[c(7, 946, 242, 181, 334, 212, 898, 461, 1191, 707, 784, 1329, 17, 10, 392, 686, 585, 638, 645, 732),]

part$name <- str_c(part$Genus, ' ', part$Species)

#Plot actual litter weaning mass vs expected litter weaning mass
plotMI1 <- ggplot(dataset) + 
  geom_point(aes(x = Litter_weaning_mass_expected, y = Litter_weaning_mass, col=Subclass), shape = 1) +
  geom_phylopic(aes(x = Litter_weaning_mass_expected, y = Litter_weaning_mass, name = name, col=Subclass), 
                size = 0.3, alpha = 0.7, data=part) +
  coord_cartesian(xlim = c(0.01,100000), ylim = c(0.01, 10000)) +
  scale_x_continuous(trans = "log10", labels = number_format(accuracy=0.1)) + 
  scale_y_continuous(trans = "log10", labels = number_format(accuracy=0.1)) +
  ylab('Actual litter weaning mass (kg)') +
  xlab('Predicted litter weaning mass (kg)') +
  theme_bw() +
  theme(legend.position = "NULL") +
  ggtitle('Predicted litter weaning mass vs actual litter weaning mass') +
  theme(plot.title = element_text(size=12, face="bold", hjust=0.4)) +
  geom_abline(size=0.8, color="black") +
  scale_color_manual(values = c("steelblue", "darkred", "#FCC501")) +
  geom_point(data = part, aes(x=Litter_weaning_mass_expected, y = Litter_weaning_mass), color="black", shape = 1) 

#Plot maternal investment for 20 indicator species 
plotMI2 <- 
  ggplot(part) + 
  geom_phylopic(aes(x = MI, y = reorder(`Common name`, -MI), name = name, col=Subclass), 
                size = 0.8) +
  coord_cartesian(xlim = c(-1.6,2)) +
  xlab('Maternal investment') +
  ylab('Species') +
  theme_bw() +
  theme(legend.position = "NULL") +
  ggtitle('Maternal investment for 20 indicator species') +
  theme(plot.title = element_text(size=12, face="bold")) +
  geom_vline(xintercept=0, size=0.8, color="black") +
  geom_point(aes(x=MI, y = `Common name`), color="black", shape = 1) +
  scale_color_manual(values = c("steelblue", "darkred", "#FCC501")) 

grid.arrange(plotMI1, plotMI2, ncol = 1)

##############################################################################
##Citation

RStudio.Version()

spp <- list(sp1='Loxodonta africana',
            sp2='Panthera tigris', 
            sp3='Tenrec ecaudatus', 
            sp4='Ailurus fulgens',
            sp5='Halichoerus grypus',
            sp6='Vulpes vulpes',  
            sp7='Pan troglodytes',
            sp8='Cynopterus sphinx',
            sp9='Rattus rattus',
            sp10='Lepus europaeus',
            sp11='Ateles geoffroyi',
            sp12='Sorex araneus',
            sp13='Bison bison',
            sp14='Aepyceros melampus',
            sp15='Balaenoptera musculus',
            sp16='Lasiorhinus latifrons',
            sp17='Sarcophilus harrisii',
            sp18='Macropus eugenii',
            sp19='Macropus rufus',
            sp20='Tachyglossus aculeatus')
uuid_list <- vector('list', length = length(spp))
att_cont <- vector('list', length = length(spp))

## Let's get the uuid first for each species in our list
for(i in 1:length(att_cont)){
  uuid_list[[i]] <- get_uuid(spp[[i]])
}

# Now the contributor. You might select to include more fields, but I only used contributor. Check them using the $
for (i in 1:length(spp)){
  att_cont[[i]] <- get_attribution(uuid = uuid_list[[i]])
  att_cont[[i]] <- paste0(spp[[i]], ', ', att_cont[[i]]$contributor)
}

#Now let's collapse species and contributors into one single string
sp_att <- paste(unlist(att_cont), collapse='; ')

ggplot()+
  labs(caption = paste0('Attribution: ', sp_att))
sp_att

citation("rphylopic")