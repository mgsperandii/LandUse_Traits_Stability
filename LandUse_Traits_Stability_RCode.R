# R code to reproduce analyses in "Functional traits mediate the effect of land use on drivers of community stability within and across trophic levels"

# R version 4.1.2
# Package versions:
# tidyverse: 2.0.0
# MuMIn: 1.47.5
# piecewiseSEM: 2.2.0
# semEff: 0.6.1


# load packages-----
library(tidyverse)
library(piecewiseSEM)
library(semEff)
library(MuMIn)

# load data-----
# numeric variables are already scaled and centered
load("YourDirectory/LandUse_Traits_Stability_Datasets.RData")

# useful functions-------
# Var_to_incl: Function for model selection; returns names of variables to retain based on a specified evidence ratio threshold.
# mod_obj: The model object on which model selection is performed.
# Fixed: A one-sided formula for the fixed component of the model.
Var_to_incl <- function(mod_obj, thr = 2.72, fixed = NULL) {
  require(MuMIn)
  sw_obj <- sw(dredge(mod_obj, fixed = fixed))
  names_sw <- names(sw_obj)
  er <- setNames(as.numeric(sw_obj)/((1 - as.numeric(sw_obj))), nm = names_sw)
  er <- er[er > thr]
  resp <- formula(mod_obj)[[2]]
  res <- paste(resp, paste("~", paste(names(er), collapse = " + ")))
  return(res)
}

# ---------------------------------------------------------------------------------
# Analyses within individual taxonomic and trophic levels (Hypothesis 1)-----------
# ---------------------------------------------------------------------------------

# 1. Plants------------------------------------------------------------------------
# ---------------------------------------------------------------------------------
# legend----
# A legend linking variable names from the PLANTS dataframes (plants_gr; plants_for)
# to the terminology used in the manuscript:
# NAME USED IN THE CODE -> NAME USED IN THE MS

# Useful_EP_PlotID  -> Plot
# meancov           -> Mean Total Abundance (MeanTotAbu in  figures)
# stability_cov     -> Plant Community Stability (computed on cover)
# eta_w             -> Synchrony index
# weighted_avg_var  -> Weighted average population variability
# LUI               -> Land-use intensity index (grasslands)
# Formi             -> Land-use intensity index (forests)
# rao               -> Multi-Trait Rao's Q index for functional diversity (FD in figures)
# PC1               -> First axis of the PCA (representing the dominant ecological strategy)
# PC2               -> Second axis of the PCA (representing the dominant ecological strategy)
# meanrich          -> Mean Species Richness
# meantreecov       -> Mean Plot Tree Cover in forests

# 1.1 grasslands------------------------------------------------------------------
# Fitting hypothesised and reduced SEMs -------------------------------------
# Fitting hypothesised model
SEM_plants_gr_hypo <- psem(
  lm(PC1 ~ LUI, data = plants_gr),
  lm(PC2 ~ LUI, data = plants_gr),
  lm(meanrich ~ rao + PC1 + PC2 + LUI, data = plants_gr), 
  lm(rao ~ LUI, data = plants_gr),
  lm(eta_w ~  LUI + meanrich + rao + PC2, data = plants_gr), # PC2 included based on d-separation test
  lm(meancov ~ meanrich + rao + PC1 + PC2 + LUI + eta_w, data = plants_gr), 
  lm(weighted_avg_var ~ meanrich + meancov + PC1 + PC2 + LUI + rao, data = plants_gr), # rao included based on d-separation test
  lm(stability_cov ~ LUI + meanrich + meancov + eta_w + rao + PC1 + PC2 + weighted_avg_var, data = plants_gr), 
  rao %~~% PC1,
  rao %~~% PC2,
  data = plants_gr)

# Check model fit
summary(SEM_plants_gr_hypo)

# Model selection for each of the models in the SEM
lapply(list(
  PC1 = lm(PC1 ~ LUI, na.action = "na.fail", data = plants_gr), 
  PC2 = lm(PC2 ~ LUI, na.action = "na.fail", data = plants_gr), 
  RICH = lm(meanrich ~ rao + PC1 + PC2 + LUI, na.action = "na.fail", data = plants_gr),
  RAO = lm(rao ~ LUI, na.action = "na.fail", data = plants_gr), 
  SYNC = lm(eta_w ~  LUI + meanrich + rao + PC2, na.action = "na.fail", data = plants_gr), 
  COV = lm(meancov ~ meanrich + rao +  PC1 + PC2 + LUI + eta_w, na.action = "na.fail", data = plants_gr),
  AVG_POP_VAR = lm(weighted_avg_var ~ meanrich + meancov + PC1 + PC2 + LUI + rao, na.action = "na.fail", data = plants_gr),
  STAB = lm(stability_cov ~ LUI + meanrich + meancov + eta_w + rao + PC1 + PC2 + weighted_avg_var, na.action = "na.fail", data = plants_gr)), Var_to_incl)

# Fitting reduced model
SEM_plants_gr_rd <-psem(
  lm(PC1 ~ LUI, data = plants_gr),
  #lm(PC2 ~ 1, data = plants_gr), # Uncomment to avoid fitting intercept-only model
  lm(meanrich ~ PC1 + PC2 + rao + LUI, data = plants_gr), 
  lm(rao ~ LUI, data = plants_gr),
  lm(eta_w ~  PC2 + LUI, data = plants_gr), 
  lm(meancov ~ meanrich + eta_w + PC1, data = plants_gr), 
  lm(weighted_avg_var ~ PC1 + meancov + PC2 + meanrich + rao, data = plants_gr), # rao included based on d-separation test
  lm(stability_cov ~ eta_w + weighted_avg_var + meanrich + rao, data = plants_gr), 
  rao %~~% PC1,
  rao %~~% PC2,
  data = plants_gr)

# Check model fit
summary(SEM_plants_gr_rd) 

# comparing AICc values
AIC_psem(SEM_plants_gr_hypo) 
AIC_psem(SEM_plants_gr_rd) 

# Bootstrapping and calculating standardised effects------ 
# Bootstrapping model effects with 10,000 resamples for the reduced SEM
SEM_plants_gr_boot <- bootEff(SEM_plants_gr_rd, R = 10000, catch.err = FALSE, parallel = "snow", seed = 23)
# Calculating standardised effects (direct, indirect, total) for endogenous variables n the bootstrapped model
SEM_plants_gr_boot_eff <- semEff(SEM_plants_gr_boot)

# Summary of effects and confidence intervals for endogenous variables i
summary(SEM_plants_gr_boot_eff) 
# Same as above, for the stability model only
summary(SEM_plants_gr_boot_eff, responses = "stability.cov") 
# Effect of LUI on stability and mediating effect of functional features 
summary(semEff(SEM_plants_gr_boot, predictor = "LUI")) 


# 1.2 forests------------------
# Fitting hypothesised and reduced SEMs -------------------------------------
# Fitting hypothesised model
SEM_plants_for_hypo <- psem(
  lm(PC1 ~ Formi + meantreecov, data = plants_for),
  lm(PC2 ~ Formi + meantreecov, data = plants_for),
  lm(meanrich ~ Formi + meantreecov + PC1 + PC2 + rao, data = plants_for),
  lm(rao ~ Formi + meantreecov, data = plants_for),
  lm(eta_w ~ Formi + meanrich + meantreecov + rao, data = plants_for),
  lm(meancov ~ meanrich + rao + PC1 + PC2 + Formi + eta_w + meantreecov, data = plants_for),
  lm(weighted_avg_var ~ meanrich + meancov + PC1 + PC2 + Formi + meantreecov + rao, data = plants_for), # rao included based on d-separation test
  lm(stability_cov ~ Formi + meanrich + meancov + eta_w + rao + PC1 + PC2 + weighted_avg_var + meantreecov, data = plants_for),
  rao %~~% PC1,
  rao %~~% PC2,
  data = plants_for)

# Check model fit
summary(SEM_plants_for_hypo) 

# Model selection for each of the models in the SEM
lapply(list(
  PC1 = lm(PC1 ~ Formi + meantreecov, na.action = "na.fail", data = plants_for), 
  PC2 = lm(PC2 ~ Formi + meantreecov, na.action = "na.fail", data = plants_for),
  RICH = lm(meanrich ~ Formi + meantreecov + PC1 + PC2 + rao, na.action = "na.fail", data = plants_for), 
  RAO = lm(rao ~ Formi + meantreecov, na.action = "na.fail", data = plants_for), 
  SYNC = lm(eta_w ~  Formi + meanrich + meantreecov + rao, na.action = "na.fail", data = plants_for), 
  COV = lm(meancov ~ meanrich + rao +  PC1 + PC2 + Formi + eta_w + meantreecov, na.action = "na.fail", data = plants_for), 
  AVG_POP_VAR = lm(weighted_avg_var ~ meanrich + meancov + PC1 + PC2 + Formi + meantreecov + rao, na.action = "na.fail", data = plants_for),
  STAB = lm(stability_cov ~ Formi + meanrich + meancov + eta_w + rao + PC1 + PC2 + weighted_avg_var + meantreecov, na.action = "na.fail", data = plants_for)), Var_to_incl)

# Fitting reduced model
SEM_plants_for_rd <- psem(
  lm(PC1 ~ meantreecov, data = plants_for), 
  #lm(PC2 ~ 1, data = plants_for) # Uncomment to avoid fitting an intercept-only model
  lm(meanrich ~ PC2 + Formi + meantreecov + PC1, data = plants_for), 
  lm(rao ~ Formi + meantreecov, data = plants_for),  
  #lm(eta_w ~ 1, data = plants_for) # Uncomment to avoid fitting an intercept-only model
  lm(meancov ~ meanrich + PC1 + meantreecov + rao + PC2, data = plants_for),
  lm(weighted_avg_var ~ rao + meancov + PC1, data = plants_for),
  lm(stability_cov ~ eta_w + weighted_avg_var + rao + meancov, data = plants_for), 
  rao %~~% PC1,
  rao %~~% PC2,
  data = plants_for)

# Check model fit
summary(SEM_plants_for_rd) 

# comparing AICc
AIC_psem(SEM_plants_for_hypo) 
AIC_psem(SEM_plants_for_rd) 

# Bootstrapping and calculating standardised effects------ 
# Bootstrapping model effects with 10,000 resamples for the reduced SEM
SEM_plants_for_boot <- bootEff(SEM_plants_for_rd, R = 10000, catch.err = FALSE, parallel = "snow", seed = 23)
# Calculating standardised effects (direct, indirect, total) for endogenous variables in the bootstrapped model
SEM_plants_for_boot_eff <- semEff(SEM_plants_for_boot)

# Summary of effects and confidence intervals for endogenous variables 
summary(SEM_plants_for_boot_eff) # summary of the SEM
# Same as above, for the stability model only
summary(SEM_plants_for_boot_eff, responses = "stability.cov") # results for the stability model only
# Effect of Formi on stability and mediating effect of functional features 
summary(semEff(SEM_plants_for_boot, predictor = "Formi")) 

# ---------------------------------------------------------------------------------
# 2. Total arthropods-------
# ---------------------------------------------------------------------------------
# legend--------
# A legend linking variable names from the ARTHROPODS dataframes (arthro_gr; arthro_for)
# to the terminology used in the manuscript:
# NAME USED IN THE CODE -> NAME USED IN THE MS

# PlotID            -> Plot
# MeanTotAbu        -> Arthropod Mean Total Abundance (MeanTotAbu in figures)
# Stab_Abu          -> Arthropod community stability (computed on arthropod abundance)
# eta_w             -> Synchrony index
# weighted_avg_var  -> Weighted average population variability
# LUI               -> Land-use intensity index (grasslands)
# Formi             -> Land-use intensity index (forests)
# rao               -> Multi-Trait Rao's Q index for arthropod functional diversity (FD in figures)
# PC1               -> First axis of the PCA (representing the dominant ecological strategy)
# PC2               -> Second axis of the PCA (representing the dominant ecological strategy)
# mean_richness     -> Mean species richness
# MeanTreeCov       -> Mean plot tree cover in forests
# PlantFD           -> Multi-Trait Rao's Q index for plant functional diversity 
# PlantPC1          -> First axis of the PCA computed on plant functional traits (representing the dominant ecological strategy in plants)
# PlantPC2          -> Second axis of the PCA computed on plant functional traits (representing the dominant ecological strategy in plants)

# 2.1 grasslands-------------------------------------------------------------
# Fitting hypothesised and reduced SEMs -------------------------------------
# Fitting hypothesised model
SEM_arthro_gr_hypo <- psem(
  lm(PC1 ~ LUI, data = arthro_gr), 
  lm(PC2 ~ LUI, data = arthro_gr),
  lm(mean_richness ~ LUI + PC1 + PC2, data = arthro_gr), 
  lm(rao ~ PC1 + PC2 + LUI, data = arthro_gr),
  lm(eta_w ~ mean_richness + LUI + rao + PC1 + PC2, data = arthro_gr),  
  lm(MeanTotAbu ~ mean_richness + rao + LUI + eta_w + PC1 + PC2, data = arthro_gr), 
  lm(weighted_avg_var ~ LUI + mean_richness + MeanTotAbu + PC1 + PC2 + rao, data = arthro_gr),
  lm(Stab_Abu ~ eta_w + mean_richness + MeanTotAbu + PC1 + PC2 + rao + LUI + weighted_avg_var, data = arthro_gr), 
  rao %~~% mean_richness,
  data = arthro_gr)

# Check model fit
summary(SEM_arthro_gr_hypo) 

# Model selection for each of the models in the SEM
lapply(list(
  PC1 = lm(PC1 ~ LUI, na.action = "na.fail", data = arthro_gr), 
  PC2 = lm(PC2 ~ LUI, na.action = "na.fail", data = arthro_gr), 
  RICH = lm(mean_richness ~ LUI + PC1 + PC2, na.action = "na.fail", data = arthro_gr), 
  RAO = lm(rao ~ PC1 + PC2 + LUI, na.action = "na.fail", data = arthro_gr),
  SYNC = lm(eta_w ~ mean_richness + LUI + rao + PC1 + PC2, na.action = "na.fail", data = arthro_gr), 
  MEANTOTABUND = lm(MeanTotAbu ~ mean_richness + rao + LUI + eta_w + PC1 + PC2, na.action = "na.fail", data = arthro_gr),  
  AVG_POP_VAR = lm(weighted_avg_var ~ LUI + mean_richness + MeanTotAbu + PC1 + PC2 + rao, na.action = "na.fail", data = arthro_gr),
  STAB = lm(Stab_Abu ~ eta_w + mean_richness + MeanTotAbu + PC1 + PC2 + rao + LUI + weighted_avg_var, na.action = "na.fail", data = arthro_gr)),
  Var_to_incl)

# Fitting reduced model
SEM_arthro_gr_rd <- psem(
  #lm(PC1 ~ 1, data = arthro_gr), # Uncomment to avoid fitting intercept-only model for PC1
  lm(PC2 ~ LUI, data = arthro_gr),
  lm(mean_richness ~ PC1 + LUI, data = arthro_gr), 
  lm(rao ~ PC1, data = arthro_gr),
  lm(eta_w ~ PC2 + LUI, data = arthro_gr),  
  lm(MeanTotAbu ~ mean_richness + eta_w + PC1 + rao, data = arthro_gr), 
  lm(weighted_avg_var ~ mean_richness + PC1 + rao + MeanTotAbu, data = arthro_gr),
  lm(Stab_Abu ~ eta_w + weighted_avg_var + MeanTotAbu, data = arthro_gr),
  rao %~~% mean_richness,
  data = arthro_gr)

# Check model fit
summary(SEM_arthro_gr_rd) 

# comparing AICc values
AIC_psem(SEM_arthro_gr_hypo) 
AIC_psem(SEM_arthro_gr_rd) 

# Bootstrapping and calculating standardised effects------ 
# Bootstrapping model effects with 10,000 resamples for the reduced SEM
SEM_arthro_gr_boot <- bootEff(SEM_arthro_gr_rd, R = 10000, catch.err = FALSE, parallel = "snow", seed = 23)
# Calculating standardised effects (direct, indirect, total) for endogenous variables in the bootstrapped model
SEM_arthro_gr_boot_eff <- semEff(SEM_arthro_gr_boot)

# Summary of effects and confidence intervals for endogenous variables
summary(SEM_arthro_gr_boot_eff) 
# Same as above, for the stability model only
summary(SEM_arthro_gr_boot_eff, responses = "Stab.Abu") 
# Effect of LUI on stability and mediating effect of functional features
summary(semEff(SEM_arthro_gr_boot, predictor = "LUI")) 

# 2.2 forests--------------
# Fitting hypothesised and reduced SEMs -------------------------------------
# Fitting hypothesised model
SEM_arthro_for_hypo <- psem(
  lm(PC1 ~ Formi + MeanTreeCov, data = arthro_for), 
  lm(PC2 ~ Formi + MeanTreeCov, data = arthro_for), 
  lm(mean_richness ~ Formi + PC1 + PC2 + MeanTreeCov, data = arthro_for),  
  lm(rao ~ PC1 + PC2 + Formi + MeanTreeCov, data = arthro_for),  
  lm(eta_w ~  mean_richness + Formi + rao + PC1 + PC2 + MeanTreeCov, data = arthro_for), 
  lm(MeanTotAbu ~ mean_richness + rao + Formi + eta_w + PC1 + PC2 + MeanTreeCov, data = arthro_for), 
  lm(weighted_avg_var ~ Formi + mean_richness + MeanTotAbu + PC1 + PC2 + rao + MeanTreeCov, data = arthro_for), 
  lm(Stab_Abu ~ eta_w + mean_richness + MeanTotAbu + PC1 + PC2 + rao + Formi + weighted_avg_var + MeanTreeCov, data = arthro_for), 
  rao %~~% mean_richness,
  data = arthro_for)

# Check model fit
summary(SEM_arthro_for_hypo) 

# Model selection for each of the models in the SEM
lapply(list(
  PC1 = lm(PC1 ~ Formi + MeanTreeCov, na.action = "na.fail", data = arthro_for),
  PC2 = lm(PC2 ~ Formi + MeanTreeCov, na.action = "na.fail", data = arthro_for),
  RICH = lm(mean_richness ~ Formi + PC1 + PC2 + MeanTreeCov, na.action = "na.fail", data = arthro_for),
  RAO = lm(rao ~ PC1 + PC2 + Formi + MeanTreeCov, na.action = "na.fail", data = arthro_for),
  SYNC = lm(eta_w ~  mean_richness + Formi + rao + PC1 + PC2 + MeanTreeCov, na.action = "na.fail", data = arthro_for),
  MEANTOTABU = lm(MeanTotAbu ~ mean_richness + rao + Formi + eta_w + PC1 + PC2 + MeanTreeCov, na.action = "na.fail", data = arthro_for),
  AVGPOPVAR = lm(weighted_avg_var ~ Formi + mean_richness + MeanTotAbu + PC1 + PC2 + rao + MeanTreeCov, na.action = "na.fail", data = arthro_for),
  STAB = lm(Stab_Abu ~ eta_w + mean_richness + MeanTotAbu + PC1 + PC2 + rao + Formi + weighted_avg_var + MeanTreeCov, na.action = "na.fail", data = arthro_for)),
  Var_to_incl)

# Fitting reduced model
SEM_arthro_for_rd <- psem(
  lm(PC1 ~ Formi, data = arthro_for), 
  lm(PC2 ~ Formi, data = arthro_for), 
  lm(mean_richness ~ PC1 + Formi, data = arthro_for),  
  lm(rao ~ PC2 + PC1, data = arthro_for),  
  lm(eta_w ~ MeanTreeCov + rao, data = arthro_for), 
  lm(MeanTotAbu ~ mean_richness, data = arthro_for), 
  lm(weighted_avg_var ~ rao, data = arthro_for), 
  lm(Stab_Abu ~ eta_w + weighted_avg_var + mean_richness + MeanTotAbu, data = arthro_for), 
  rao %~~% mean_richness,
  data = arthro_for)

# Check model fit
summary(SEM_arthro_for_rd) 

# comparing AICc values
AIC_psem(SEM_arthro_for_hypo) 
AIC_psem(SEM_arthro_for_rd) 

# Bootstrapping and calculating standardised effects------ 
# Bootstrapping model effects with 10,000 resamples for the reduced SEM
SEM_arthro_for_boot <- bootEff(SEM_arthro_for_rd, R = 10000, catch.err = FALSE, parallel = "snow", seed = 23)
# Calculating standardised effects (direct, indirect, total) for endogenous variables in the bootstrapped model
SEM_arthro_for_boot_eff <- semEff(SEM_arthro_for_boot)

# Summary of effects and confidence intervals for endogenous variables
summary(SEM_arthro_for_boot_eff) 
# Same as above, for the stability model only
summary(SEM_arthro_for_boot_eff, responses = "Stab.Abu") 
# Effect of Formi on stability and mediating effect of functional features 
summary(semEff(SEM_arthro_for_boot, predictor = "Formi")) 

# ---------------------------------------------------------------------------------
# 3. Herbivores-------
# --------------------------------------------------------------------------------
# legend--------
# A legend linking variable names from the HERBIVORES dataframes (herbiv_gr; herbiv_for)
# to the terminology used in the manuscript:
# NAME USED IN THE CODE -> NAME USED IN THE MS

# PlotID            -> Plot
# MeanTotAbu_H      -> Herbivore Mean Total Abundance (MeanTotAbu in figures)
# Stab_Abu_H        -> Herbivore community stability (computed on arthropod abundance)
# eta_w             -> Synchrony index
# weighted_avg_var  -> Weighted average population variability
# LUI               -> Land-use intensity index (grasslands)
# Formi             -> Land-use intensity index (forests)
# rao               -> Multi-Trait Rao's Q index for herbivore functional diversity (FD in figures)
# PC1               -> First axis of the PCA (representing the dominant ecological strategy)
# PC2               -> Second axis of the PCA (representing the dominant ecological strategy)
# MeanRich_H        -> Herbivore mean species richness
# MeanTreeCov       -> Mean plot tree cover in forests
# PlantFD           -> Multi-Trait Rao's Q index for plant functional diversity 
# PlantPC1          -> First axis of the PCA computed on plant functional traits (representing the dominant ecological strategy in plants)
# PlantPC2          -> Second axis of the PCA computed on plant functional traits (representing the dominant ecological strategy in plants)

# 3.1 grasslands--------------
# Fitting hypothesised and reduced SEMs -------------------------------------
# Fitting hypothesised model
SEM_herb_gr_hypo <- psem(
  lm(PC1 ~ LUI, data = herbiv_gr), 
  lm(PC2 ~ LUI, data = herbiv_gr),
  lm(MeanRich_H ~ LUI + PC1 + PC2, data = herbiv_gr),  
  lm(rao ~ PC1 + PC2 + LUI, data = herbiv_gr),  
  lm(eta_w ~ MeanRich_H + LUI + rao + PC1 + PC2, data = herbiv_gr),  
  lm(MeanTotAbu_H ~ MeanRich_H + rao + LUI + eta_w + PC1 + PC2, data = herbiv_gr), 
  lm(weighted_avg_var ~ LUI + MeanRich_H + MeanTotAbu_H + PC1 + PC2 + rao, data = herbiv_gr),
  lm(Stab_Abu_H ~ eta_w + MeanRich_H + MeanTotAbu_H + PC1 + PC2 + rao + LUI + weighted_avg_var, data = herbiv_gr), 
  rao %~~% MeanRich_H,
  data = herbiv_gr)

summary(SEM_herb_gr_hypo) 

# Model selection for each of the models in the SEM
lapply(list(
  PC1 = lm(PC1 ~ LUI, na.action = "na.fail", data = herbiv_gr),
  PC2 = lm(PC2 ~ LUI, na.action = "na.fail", data = herbiv_gr),
  RICH = lm(MeanRich_H ~ LUI + PC1 + PC2, na.action = "na.fail", data = herbiv_gr),
  RAO = lm(rao ~ PC1 + PC2 + LUI, na.action = "na.fail", data = herbiv_gr),
  SYNC = lm(eta_w ~ MeanRich_H + LUI + rao + PC1 + PC2, na.action = "na.fail", data = herbiv_gr),
  MEANTOTABU = lm(MeanTotAbu_H ~ MeanRich_H + rao + LUI + eta_w + PC1 + PC2, na.action = "na.fail", data = herbiv_gr),
  AVGPOPVAR = lm(weighted_avg_var ~ LUI + MeanRich_H + MeanTotAbu_H + PC1 + PC2 + rao, na.action = "na.fail", data = herbiv_gr),
  STAB = lm(Stab_Abu_H ~ eta_w + MeanRich_H + MeanTotAbu_H + PC1 + PC2 + rao + LUI + weighted_avg_var, na.action = "na.fail", data = herbiv_gr)),
  Var_to_incl)

# Fitting reduced model
SEM_herb_gr_rd <- psem(
  #lm(PC1 ~ 1, data = herbiv_gr), # Uncomment to avoid fitting intercept-only model 
  lm(PC2 ~ LUI, data = herbiv_gr),
  lm(MeanRich_H ~ PC2 + LUI, data = herbiv_gr),  
  lm(rao ~ PC2 + PC1, data = herbiv_gr),  
  lm(eta_w ~ rao, data = herbiv_gr), 
  lm(MeanTotAbu_H ~ MeanRich_H + eta_w + rao + PC1 + LUI, data = herbiv_gr), # MeanRich_H included based on d-separation test
  lm(weighted_avg_var ~ MeanRich_H + rao + PC1 + MeanTotAbu_H + PC2, data = herbiv_gr),
  lm(Stab_Abu_H ~ eta_w + weighted_avg_var + MeanTotAbu_H, data = herbiv_gr),
  rao %~~% MeanRich_H,
  data = herbiv_gr)

# Check model fit
summary(SEM_herb_gr_rd) 

# comparing AICc values
AIC_psem(SEM_herb_gr_hypo) 
AIC_psem(SEM_herb_gr_rd) 

# Bootstrapping and calculating standardised effects------ 
# Bootstrapping model effects with 10,000 resamples for the reduced SEM
SEM_herb_gr_boot <- bootEff(SEM_herb_gr_rd, R = 10000, catch.err = FALSE, parallel = "snow", seed = 23)
# Calculating standardised effects (direct, indirect, total) for endogenous variables in the bootstrapped model
SEM_herb_gr_boot_eff <- semEff(SEM_herb_gr_boot)

# Summary of effects and confidence intervals for endogenous variables
summary(SEM_herb_gr_boot_eff) 
# Same as above, for the stability model only
summary(SEM_herb_gr_boot_eff, responses = "Stab.Abu.H") 
# Effect of LUI on stability and mediating effect of functional features 
summary(semEff(SEM_herb_gr_boot, predictor = "LUI")) 

# 3.2 forests--------------
# Fitting hypothesised and reduced SEMs -------------------------------------
# Fitting hypothesised model
SEM_herb_for_hypo <- psem(
  lm(PC1 ~ Formi + MeanTreeCov, data = herbiv_for), 
  lm(PC2 ~ Formi + MeanTreeCov, data = herbiv_for),
  lm(MeanRich_H ~ Formi + PC1 + PC2 + MeanTreeCov, data = herbiv_for),
  lm(rao ~ PC1 + PC2 + Formi + MeanTreeCov, data = herbiv_for), 
  lm(eta_w ~  MeanRich_H + Formi + rao + PC1 + PC2 + MeanTreeCov, data = herbiv_for),
  lm(MeanTotAbu_H ~ MeanRich_H + rao + Formi + eta_w + PC1 + PC2 + MeanTreeCov, data = herbiv_for), 
  lm(weighted_avg_var ~ Formi + MeanRich_H + MeanTotAbu_H + PC1 + PC2 + rao + MeanTreeCov, data = herbiv_for), 
  lm(Stab_Abu_H ~ eta_w + MeanRich_H + MeanTotAbu_H + PC1 + PC2 + rao + Formi + weighted_avg_var + MeanTreeCov, data = herbiv_for), 
  rao %~~% MeanRich_H,
  data = herbiv_for)

# Check model fit
summary(SEM_herb_for_hypo) 

# Model selection for each of the models in the SEM
lapply(list(
  PC1 = lm(PC1 ~ Formi + MeanTreeCov, na.action = "na.fail", data = herbiv_for), 
  PC2 = lm(PC2 ~ Formi + MeanTreeCov, na.action = "na.fail", data = herbiv_for), 
  RICH = lm(MeanRich_H ~ Formi + PC1 + PC2 + MeanTreeCov, na.action = "na.fail", data = herbiv_for), 
  RAO = lm(rao ~ PC1 + PC2 + Formi + MeanTreeCov, na.action = "na.fail", data = herbiv_for),  
  SYNC = lm(eta_w ~  MeanRich_H + Formi + rao + PC1 + PC2 + MeanTreeCov, na.action = "na.fail", data = herbiv_for), 
  MEANTOTABUND = lm(MeanTotAbu_H ~ MeanRich_H + rao + Formi + eta_w + PC1 + PC2 + MeanTreeCov, na.action = "na.fail", data = herbiv_for), 
  AVG_POP_VAR = lm(weighted_avg_var ~ Formi + MeanRich_H + MeanTotAbu_H + PC1 + PC2 + rao + MeanTreeCov, na.action = "na.fail", data = herbiv_for),
  STAB = lm(Stab_Abu_H ~ eta_w + MeanRich_H + MeanTotAbu_H + PC1 + PC2 + rao + Formi + weighted_avg_var + MeanTreeCov, na.action = "na.fail", data = herbiv_for)),
  Var_to_incl)

# Fitting reduced model
SEM_herb_for_rd <- psem(
  #lm(PC1 ~ 1, data = herbiv_for), # Uncomment to avoid fitting intercept-only model 
  lm(PC2 ~ Formi, data = herbiv_for),
  lm(MeanRich_H ~ PC1, data = herbiv_for), 
  lm(rao ~ PC1, data = herbiv_for), 
  lm(eta_w ~  PC1, data = herbiv_for), 
  lm(MeanTotAbu_H ~ PC2 + MeanRich_H + rao, data = herbiv_for), 
  lm(weighted_avg_var ~ Formi, data = herbiv_for), 
  lm(Stab_Abu_H ~ weighted_avg_var + eta_w + PC2 + PC1, data = herbiv_for), 
  rao %~~% MeanRich_H,
  data = herbiv_for)

# Check model fit
summary(SEM_herb_for_rd) 

# comparing AICc values
AIC_psem(SEM_herb_for_hypo) 
AIC_psem(SEM_herb_for_rd) 

# Bootstrapping and calculating standardised effects------ 
# Bootstrapping model effects with 10,000 resamples for the reduced SEM
SEM_herb_for_boot <- bootEff(SEM_herb_for_rd, R = 10000, catch.err = FALSE, parallel = "snow", seed = 23)
# Calculating standardised effects (direct, indirect, total) for endogenous variables in the bootstrapped model
SEM_herb_for_boot_eff <- semEff(SEM_herb_for_boot)

# Summary of effects and confidence intervals for endogenous variables
summary(SEM_herb_for_boot_eff) 
# Same as above, for the stability model only
summary(SEM_herb_for_boot_eff, responses = "Stab.Abu.H") 
# Effect of Formi on stability and mediating effect of functional features 
summary(semEff(SEM_herb_for_boot, predictor = "Formi")) 

# ---------------------------------------------------------------------------------
# 4. Carnivores----
# ---------------------------------------------------------------------------------
# legend--------
# A legend linking variable names from the CARNIVORES dataframes (carniv_gr; carniv_for)
# to the terminology used in the manuscript:
# NAME USED IN THE CODE -> NAME USED IN THE MS

# PlotID            -> Plot
# MeanTotAbu_Carn   -> Carnivore Mean total Abundance (MeanTotAbu in figures)
# Stab_Abu_Carn     -> Carnivore community stability (computed on carnivore abundance)
# eta_w             -> Synchrony index
# weighted_avg_var  -> Weighted average population variability
# LUI               -> Land-use intensity index (grasslands)
# Formi             -> Land-use intensity index (forests)
# rao               -> Multi-Trait Rao's Q index for carnivore functional diversity (FD in figures)
# PC1               -> First axis of the PCA (representing the dominant ecological strategy)
# PC2               -> Second axis of the PCA (representing the dominant ecological strategy)
# MeanRich_Carn     -> Carnivore mean species richness
# MeanTreeCov       -> Mean plot tree cover in forests
# HerbPC1           -> First axis of the PCA computed on herbivores functional traits (representing the dominant ecological strategy in herbivores)
# HerbPC2           -> Second axis of the PCA computed on herbivores functional traits (representing the dominant ecological strategy in herbivores)
# Herbrao           -> Multi-Trait Rao's Q index for herbivore functional diversity 
# PlantFD           -> Multi-Trait Rao's Q index for plant functional diversity 
# PlantPC1          -> First axis of the PCA computed on plant functional traits (representing the dominant ecological strategy in plants)
# PlantPC2          -> Second axis of the PCA computed on plant functional traits (representing the dominant ecological strategy in plants)

# 4.1 grasslands--------------
# Fitting hypothesised and reduced SEMs -------------------------------------
# Fitting hypothesised model
SEM_carn_gr_hypo <- psem(
  lm(PC1 ~ LUI, data = carniv_gr),  
  lm(PC2 ~ LUI, data = carniv_gr),
  lm(MeanRich_Carn ~ LUI + PC1 + PC2, data = carniv_gr), 
  lm(rao ~ PC1 + PC2 + LUI, data = carniv_gr),  
  lm(eta_w ~ MeanRich_Carn + LUI + rao + PC1 + PC2, data = carniv_gr),  
  lm(MeanTotAbu_Carn ~ MeanRich_Carn + rao + LUI + eta_w + PC1 + PC2, data = carniv_gr),
  lm(weighted_avg_var ~ LUI + MeanRich_Carn + MeanTotAbu_Carn + PC1 + PC2 + rao, data = carniv_gr),
  lm(Stab_Abu_Carn ~ eta_w + MeanRich_Carn + MeanTotAbu_Carn + PC1 + PC2 + rao + LUI + weighted_avg_var, data = carniv_gr), 
  rao %~~% MeanRich_Carn,
  data = carniv_gr)

# Check model fit
summary(SEM_carn_gr_hypo) 

# Model selection for each of the models in the SEM
lapply(list(
  PC1 = lm(PC1 ~ LUI, na.action = "na.fail", data = carniv_gr),
  PC2 = lm(PC2 ~ LUI, na.action = "na.fail", data = carniv_gr),
  RICH = lm(MeanRich_Carn ~ LUI + PC1 + PC2, na.action = "na.fail", data = carniv_gr), 
  RAO = lm(rao ~ PC1 + PC2 + LUI, na.action = "na.fail", data = carniv_gr),
  SYNC = lm(eta_w ~ MeanRich_Carn + LUI + rao + PC1 + PC2, na.action = "na.fail", data = carniv_gr), 
  MEANTOTABUND = lm(MeanTotAbu_Carn ~ MeanRich_Carn + rao + LUI + eta_w + PC1 + PC2, na.action = "na.fail", data = carniv_gr),
  AVG_POP_VAR = lm(weighted_avg_var ~ LUI + MeanRich_Carn + MeanTotAbu_Carn + PC1 + PC2 + rao, na.action = "na.fail", data = carniv_gr),
  STAB = lm(Stab_Abu_Carn ~ eta_w + MeanRich_Carn + MeanTotAbu_Carn + PC1 + PC2 + rao + LUI + weighted_avg_var, na.action = "na.fail", data = carniv_gr)), 
  Var_to_incl)

# fitting reduced model
SEM_carn_gr_rd <- psem(
  #lm(PC1 ~ 1, data = carniv_gr), # Uncomment to avoid fitting intercept-only model 
  lm(PC2 ~ LUI, data = carniv_gr), 
  lm(MeanRich_Carn ~ LUI, data = carniv_gr), 
  lm(rao ~ PC1, data = carniv_gr),  
  lm(eta_w ~ MeanRich_Carn, data = carniv_gr),  
  lm(MeanTotAbu_Carn ~ MeanRich_Carn + PC1 + rao, data = carniv_gr),
  lm(weighted_avg_var ~ MeanRich_Carn + LUI + rao + PC1, data = carniv_gr),
  lm(Stab_Abu_Carn ~ eta_w + MeanRich_Carn + MeanTotAbu_Carn + weighted_avg_var, data = carniv_gr), # eta_w included based on d-separation test
  rao %~~% MeanRich_Carn,
  data = carniv_gr)

# Check model fit
summary(SEM_carn_gr_rd)  

# comparing AICc values
AIC_psem(SEM_carn_gr_hypo) 
AIC_psem(SEM_carn_gr_rd) 

# Bootstrapping and calculating standardised effects------ 
# Bootstrapping model effects with 10,000 resamples for the reduced SEM
SEM_carn_gr_boot <- bootEff(SEM_carn_gr_rd, R = 10000, catch.err = FALSE, parallel = "snow", seed = 23)
# Calculating standardised effects (direct, indirect, total) for endogenous variables in the bootstrapped model
SEM_carn_gr_boot_eff <- semEff(SEM_carn_gr_boot)

# Summary of effects and confidence intervals for endogenous variables 
summary(SEM_carn_gr_boot_eff) 
# Same as above, for the stability model only
summary(SEM_carn_gr_boot_eff, responses = "Stab.Abu.Carn") 
# Effect of LUI on stability and mediating effect of functional features 
summary(semEff(SEM_carn_gr_boot, predictor = "LUI")) 
# ********************************************************************************************************

# 4.2 forests--------------
# Fitting hypothesised and reduced SEMs -------------------------------------
# Fitting hypothesised model
SEM_carn_for_hypo <- psem(
  lm(PC1 ~ Formi + MeanTreeCov, data = carniv_for),   
  lm(PC2 ~ Formi + MeanTreeCov, data = carniv_for),
  lm(MeanRich_Carn ~ Formi + PC1 + PC2 + MeanTreeCov, data = carniv_for), 
  lm(rao ~ PC1 + PC2 + Formi + MeanTreeCov, data = carniv_for),  
  lm(eta_w ~ MeanRich_Carn + Formi + rao + PC1 + PC2 + MeanTreeCov , data = carniv_for),  
  lm(MeanTotAbu_Carn ~ MeanRich_Carn + rao + Formi + eta_w + PC1 + PC2 + MeanTreeCov, data = carniv_for),
  lm(weighted_avg_var ~ Formi + MeanRich_Carn + MeanTotAbu_Carn + PC1 + PC2 + rao + MeanTreeCov, data = carniv_for),
  lm(Stab_Abu_Carn ~ eta_w + MeanRich_Carn + MeanTotAbu_Carn + PC1 + PC2 + rao + Formi + weighted_avg_var + MeanTreeCov, data = carniv_for), 
  rao %~~% MeanRich_Carn,
  data = carniv_for)

# Check model fit
summary(SEM_carn_for_hypo) 

# Model selection for each of the models in the SEM
lapply(list(
  PC1 = lm(PC1 ~ Formi + MeanTreeCov, na.action = "na.fail", data = carniv_for),
  PC2 = lm(PC2 ~ Formi + MeanTreeCov, na.action = "na.fail", data = carniv_for),
  RICH = lm(MeanRich_Carn ~ Formi + PC1 + PC2 + MeanTreeCov, na.action = "na.fail", data = carniv_for), 
  RAO = lm(rao ~ PC1 + PC2 + Formi + MeanTreeCov, na.action = "na.fail", data = carniv_for),
  SYNC = lm(eta_w ~ MeanRich_Carn + Formi + rao + PC1 + PC2 + MeanTreeCov, na.action = "na.fail", data = carniv_for), 
  MEANTOTABUND = lm(MeanTotAbu_Carn ~ MeanRich_Carn + rao + Formi + eta_w + PC1 + PC2 + MeanTreeCov, na.action = "na.fail", data = carniv_for), 
  AVG_POP_VAR = lm(weighted_avg_var ~ Formi + MeanRich_Carn + MeanTotAbu_Carn + PC1 + PC2 + rao + MeanTreeCov, na.action = "na.fail", data = carniv_for),
  STAB = lm(Stab_Abu_Carn ~ eta_w + MeanRich_Carn + MeanTotAbu_Carn + PC1 + PC2 + rao + Formi + weighted_avg_var + MeanTreeCov, na.action = "na.fail", data = carniv_for)),
  Var_to_incl)

# Fitting reduced model
SEM_carn_for_rd <- psem(
  # lm(PC1 ~ 1, data = carniv_for), # Uncomment to avoid fitting intercept-only model 
  # lm(PC2 ~ 1, data = carniv_for), # Uncomment to avoid fitting intercept-only model 
  # lm(MeanRich_Carn ~ 1, data = carniv_for), # Uncomment to avoid fitting intercept-only model 
  # lm(rao ~ 1, data = carniv_for), # Uncomment to avoid fitting intercept-only model 
  # lm(eta_w ~ 1, data = carniv_for), # Uncomment to avoid fitting intercept-only model 
  lm(MeanTotAbu_Carn ~ MeanRich_Carn + Formi, data = carniv_for), 
  #lm(weighted_avg_var ~ 1, data = carniv_for), # Uncomment to avoid fitting intercept-only model 
  lm(Stab_Abu_Carn ~ eta_w + weighted_avg_var + MeanTotAbu_Carn + MeanRich_Carn, data = carniv_for),
  rao %~~% MeanRich_Carn,
  data = carniv_for)

# Check model fit
summary(SEM_carn_for_rd) 

# comparing AICc values
AIC_psem(SEM_carn_for_hypo) 
AIC_psem(SEM_carn_for_rd) 

# Bootstrapping and calculating standardised effects------ 
# Bootstrapping model effects with 10,000 resamples for the reduced SEM
SEM_carn_for_boot <- bootEff(SEM_carn_for_rd, R = 10000, catch.err = FALSE, parallel = "snow", seed = 23)
# Calculating standardised effects (direct, indirect, total) for endogenous variables in the bootstrapped model
SEM_carn_for_boot_eff <- semEff(SEM_carn_for_boot)

# Summary of effects and confidence intervals for endogenous variables
summary(SEM_carn_for_boot_eff) 
# Same as above, for the stability model only
summary(SEM_carn_for_boot_eff, responses = "Stab.Abu.Carn") 
# Effect of Formi on stability and mediating effect of functional features 
summary(semEff(SEM_carn_for_boot, predictor = "Formi")) # effect of Formi on stability + mediating eff of functional features


# ---------------------------------------------------------------------------------
# Analyses across individual taxonomic and trophic levels (Hypothesis 2)-----------
# ---------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------
# 5. Total arthropods-----
# ---------------------------------------------------------------------------------
# legend--------
# A legend linking variable names from the ARTHROPODS dataframes (arthro_gr; arthro_for)
# to the terminology used in the manuscript:
# NAME USED IN THE CODE -> NAME USED IN THE MS

# PlotID            -> Plot
# MeanTotAbu        -> Arthropod Mean Total Abundance (MeanTotAbu in figures)
# Stab_Abu          -> Arthropod community stability (computed on arthropod abundance)
# eta_w             -> Synchrony index
# weighted_avg_var  -> Weighted average population variability
# LUI               -> Land-use intensity index (grasslands)
# Formi             -> Land-use intensity index (forests)
# rao               -> Multi-Trait Rao's Q index for arthropod functional diversity (FD in figures)
# PC1               -> First axis of the PCA (representing the dominant ecological strategy)
# PC2               -> Second axis of the PCA (representing the dominant ecological strategy)
# mean_richness     -> Mean species richness
# MeanTreeCov       -> Mean plot tree cover in forests
# PlantFD           -> Multi-Trait Rao's Q index for plant functional diversity 
# PlantPC1          -> First axis of the PCA computed on plant functional traits (representing the dominant ecological strategy in plants)
# PlantPC2          -> Second axis of the PCA computed on plant functional traits (representing the dominant ecological strategy in plants)


# 5.1 grasslands--------------
# Fitting hypothesised and reduced SEMs -------------------------------------
# Fitting hypothesised model
acroSEM_arthro_gr_hypo <- psem(
  lm(PlantPC1 ~ LUI, data = arthro_gr),
  lm(PlantPC2 ~ LUI, data = arthro_gr),
  lm(PlantFD ~ LUI, data = arthro_gr),
  lm(PC1 ~ LUI + PlantFD + PlantPC1 + PlantPC2, data = arthro_gr), 
  lm(PC2 ~ LUI + PlantFD + PlantPC1 + PlantPC2, data = arthro_gr),
  lm(mean_richness ~ LUI + PlantFD + PlantPC1 + PlantPC2 + PC1 + PC2, data = arthro_gr), 
  lm(rao ~ PC1 + PC2 + LUI + PlantFD + PlantPC1 + PlantPC2, data = arthro_gr),
  lm(eta_w ~ mean_richness + LUI + rao + PC1 + PC2 + PlantFD + PlantPC1 + PlantPC2, data = arthro_gr),  
  lm(MeanTotAbu ~ mean_richness + rao + LUI + eta_w + PlantFD + PlantPC1 + PlantPC2 + PC1 + PC2, data = arthro_gr), 
  lm(weighted_avg_var ~ LUI + mean_richness + MeanTotAbu + PC1 + PC2 + rao + PlantFD + PlantPC1 + PlantPC2, data = arthro_gr),
  lm(Stab_Abu ~ eta_w + mean_richness + MeanTotAbu + PC1 + PC2 + rao + LUI + weighted_avg_var + PlantFD + PlantPC1 + PlantPC2, data = arthro_gr),
  rao %~~% mean_richness,
  PlantPC1 %~~% PlantFD,
  PlantPC2 %~~% PlantFD,
  data = arthro_gr)

# Check model fit
summary(acroSEM_arthro_gr_hypo) 

# Model selection for each of the models in the SEM
lapply(list(
  PlantPC1 = lm(PlantPC1 ~ LUI, na.action = "na.fail", data = arthro_gr),
  PlantPC2 = lm(PlantPC2 ~ LUI, na.action = "na.fail", data = arthro_gr),
  PlantFD = lm(PlantFD ~ LUI, na.action = "na.fail", data = arthro_gr),
  PC1 = lm(PC1 ~ LUI + PlantFD + PlantPC1 + PlantPC2, na.action = "na.fail", data = arthro_gr), 
  PC2 = lm(PC2 ~ LUI + PlantFD + PlantPC1 + PlantPC2, na.action = "na.fail", data = arthro_gr), 
  RICH = lm(mean_richness ~ LUI + PlantFD + PlantPC1 + PlantPC2 + PC1 + PC2, na.action = "na.fail", data = arthro_gr),
  RAO = lm(rao ~ PC1 + PC2 + LUI + PlantFD + PlantPC1 + PlantPC2, na.action = "na.fail", data = arthro_gr),
  SYNC = lm(eta_w ~ mean_richness + LUI + rao + PC1 + PC2 + PlantFD + PlantPC1 + PlantPC2, na.action = "na.fail", data = arthro_gr), 
  MEANTOTABUND = lm(MeanTotAbu ~ mean_richness + rao + LUI + eta_w + PlantFD + PlantPC1 + PlantPC2 + PC1 + PC2, na.action = "na.fail", data = arthro_gr),  
  AVG_POP_VAR = lm(weighted_avg_var ~ LUI + mean_richness + MeanTotAbu + PC1 + PC2 + rao + PlantFD + PlantPC1 + PlantPC2, na.action = "na.fail", data = arthro_gr),
  STAB = lm(Stab_Abu ~ eta_w + mean_richness + MeanTotAbu + PC1 + PC2 + rao + LUI + weighted_avg_var + PlantFD + PlantPC1 + PlantPC2, na.action = "na.fail", data = arthro_gr)), Var_to_incl)

# Fitting reduced model
acroSEM_arthro_gr_rd <- psem(
  lm(PlantPC1 ~ LUI, data = arthro_gr),
  #lm(PlantPC2 ~ 1, data = arthro_gr),
  lm(PlantFD ~ LUI, data = arthro_gr),
  #lm(PC1 ~ 1, data = arthro_gr), 
  lm(PC2 ~ PlantPC1, data = arthro_gr),
  lm(mean_richness ~ PC1 + LUI + PlantPC1, data = arthro_gr), 
  lm(rao ~ PC1 + PlantPC1, data = arthro_gr),
  lm(eta_w ~ PC2 + LUI, data = arthro_gr),  
  lm(MeanTotAbu ~ mean_richness + eta_w + PC1 + rao, data = arthro_gr), 
  lm(weighted_avg_var ~ mean_richness + PC1 + rao + MeanTotAbu, data = arthro_gr),
  lm(Stab_Abu ~ eta_w + weighted_avg_var + MeanTotAbu, data = arthro_gr), 
  rao %~~% mean_richness,
  PlantPC1 %~~% PlantFD,
  PlantPC2 %~~% PlantFD,
  data = arthro_gr)

# Check model fit
summary(acroSEM_arthro_gr_rd) 

# comparing AICc values
AIC_psem(acroSEM_arthro_gr_hypo) 
AIC_psem(acroSEM_arthro_gr_rd) 

# Bootstrapping and calculating standardised effects------ 
# Bootstrapping model effects with 10,000 resamples for the reduced SEM
acroSEM_arthro_gr_boot <- bootEff(acroSEM_arthro_gr_rd, R = 10000, catch.err = FALSE, parallel = "snow", seed = 23)
# Calculating standardised effects (direct, indirect, total) for endogenous variables in the bootstrapped model
acroSEM_arthro_gr_boot_eff <- semEff(acroSEM_arthro_gr_boot)

# Summary of effects and confidence intervals for endogenous variables
summary(acroSEM_arthro_gr_boot_eff) 
# Same as above, for the stability model only
summary(acroSEM_arthro_gr_boot_eff, responses = "Stab.Abu") 
# Effect of LUI on stability and mediating effect of functional features 
summary(semEff(acroSEM_arthro_gr_boot, predictor = "LUI")) 


# 5.2 forests--------------
# Fitting hypothesised and reduced SEMs -------------------------------------
# Fitting hypothesised model
acroSEM_arthro_for_hypo <- psem(
  lm(PlantPC1 ~ Formi + MeanTreeCov, data = arthro_for),
  lm(PlantPC2 ~ Formi + MeanTreeCov, data = arthro_for),
  lm(PlantFD ~ Formi + MeanTreeCov, data = arthro_for),
  lm(PC1 ~ Formi + PlantFD + PlantPC1 + PlantPC2 + MeanTreeCov, data = arthro_for), 
  lm(PC2 ~ Formi + PlantFD + PlantPC1 + PlantPC2 + MeanTreeCov, data = arthro_for),
  lm(mean_richness ~ Formi + PC1 + PC2 + PlantFD + PlantPC1 + PlantPC2 + MeanTreeCov, data = arthro_for),  
  lm(rao ~ PC1 + PC2 + Formi + PlantFD + PlantPC1 + PlantPC2 + MeanTreeCov, data = arthro_for),  
  lm(eta_w ~  mean_richness + Formi + rao + PC1 + PC2 + PlantFD + PlantPC1 + PlantPC2 + MeanTreeCov, data = arthro_for), 
  lm(MeanTotAbu ~ mean_richness + rao + Formi + eta_w + PlantFD + PlantPC1 + PlantPC2 + PC1 + PC2 + MeanTreeCov, data = arthro_for), 
  lm(weighted_avg_var ~ Formi + mean_richness + MeanTotAbu + PC1 + PC2 + rao + PlantFD + PlantPC1 + PlantPC2 + MeanTreeCov, data = arthro_for), 
  lm(Stab_Abu ~ eta_w + mean_richness + MeanTotAbu + PC1 + PC2 + rao + Formi + weighted_avg_var + PlantFD + PlantPC1 + PlantPC2 + MeanTreeCov, data = arthro_for), 
  rao %~~% mean_richness,
  PlantPC1 %~~% PlantFD,
  PlantPC2 %~~% PlantFD,
  data = arthro_for)

# Check model fit
summary(acroSEM_arthro_for_hypo) 

# Model selection for each of the models in the SEM
lapply(list(
  PLANTPC1 = lm(PlantPC1 ~ Formi + MeanTreeCov, na.action = "na.fail", data = arthro_for),
  PLANTPC2 = lm(PlantPC2 ~ Formi + MeanTreeCov, na.action = "na.fail", data = arthro_for),
  PLANTFD = lm(PlantFD ~ Formi + MeanTreeCov, na.action = "na.fail", data = arthro_for),
  PC1 = lm(PC1 ~ Formi + PlantFD + PlantPC1 + PlantPC2 + MeanTreeCov, na.action = "na.fail", data = arthro_for),
  PC2 = lm(PC2 ~ Formi + PlantFD + PlantPC1 + PlantPC2 + MeanTreeCov, na.action = "na.fail", data = arthro_for),
  RICH = lm(mean_richness ~ Formi + PC1 + PC2 + PlantFD + PlantPC1 + PlantPC2 + MeanTreeCov, na.action = "na.fail", data = arthro_for),
  RAO = lm(rao ~ PC1 + PC2 + Formi + PlantFD + PlantPC1 + PlantPC2 + MeanTreeCov, na.action = "na.fail", data = arthro_for),
  SYNC = lm(eta_w ~  mean_richness + Formi + rao + PC1 + PC2 + PlantFD + PlantPC1 + PlantPC2 + MeanTreeCov, na.action = "na.fail", data = arthro_for),
  MEANTOTABU = lm(MeanTotAbu ~ mean_richness + rao + Formi + eta_w + PlantFD + PlantPC1 + PlantPC2 + PC1 + PC2 + MeanTreeCov, na.action = "na.fail", data = arthro_for),
  AVGPOPVAR = lm(weighted_avg_var ~ Formi + mean_richness + MeanTotAbu + PC1 + PC2 + rao + PlantFD + PlantPC1 + PlantPC2 + MeanTreeCov, na.action = "na.fail", data = arthro_for),
  STAB = lm(Stab_Abu ~ eta_w + mean_richness + MeanTotAbu + PC1 + PC2 + rao + Formi + weighted_avg_var + PlantFD + PlantPC1 + PlantPC2 + MeanTreeCov, na.action = "na.fail", data = arthro_for)),
  Var_to_incl)

# Fitting reduced model
acroSEM_arthro_for_rd <- psem(
  #lm(PlantPC1 ~ 1, data = arthro_for), # Uncomment to avoid fitting intercept-only model
  #lm(PlantPC2 ~ 1, data = arthro_for), # Uncomment to avoid fitting intercept-only model
  lm(PlantFD ~ Formi, data = arthro_for),
  lm(PC1 ~ Formi, data = arthro_for), 
  #lm(PC2 ~ Formi, data = arthro_for),  # Uncomment to avoid fitting intercept-only model
  lm(mean_richness ~ PC1 + Formi + PlantFD, data = arthro_for), # PlantFD included based on d-separation test
  lm(rao ~ PC2 + PC1, data = arthro_for),  
  lm(eta_w ~ rao + MeanTreeCov, data = arthro_for), # MeanTreeCov included based on d-separation test
  lm(MeanTotAbu ~ mean_richness + eta_w + PlantFD, data = arthro_for), # PlantFD included based on d-separation test
  lm(weighted_avg_var ~ rao + MeanTreeCov + PlantPC1 + Formi, data = arthro_for), 
  lm(Stab_Abu ~ weighted_avg_var + eta_w + mean_richness + MeanTotAbu + PlantFD, data = arthro_for), # PlantFD included based on d-separation test
  rao %~~% mean_richness,
  PlantPC1 %~~% PlantFD,
  PlantPC2 %~~% PlantFD,
  data = arthro_for)

# Check model fit
summary(acroSEM_arthro_for_rd) 

# comparing AICc values
AIC_psem(acroSEM_arthro_for_hypo) 
AIC_psem(acroSEM_arthro_for_rd) 

# Bootstrapping and calculating standardised effects------ 
# Bootstrapping model effects with 10,000 resamples for the reduced SEM
acroSEM_arthro_for_boot <- bootEff(acroSEM_arthro_for_rd, R = 10000, catch.err = FALSE, parallel = "snow", seed = 23)
# Calculating standardised effects (direct, indirect, total) for endogenous variables in the bootstrapped model
acroSEM_arthro_for_boot_eff <- semEff(acroSEM_arthro_for_boot)

# Summary of effects and confidence intervals for endogenous variables
summary(acroSEM_arthro_for_boot_eff) 
# Same as above, for the stability model only
summary(acroSEM_arthro_for_boot_eff, responses = "Stab.Abu") 
# Effect of Formi on stability and mediating effect of functional features 
summary(semEff(acroSEM_arthro_for_boot, predictor = "Formi")) 

# ---------------------------------------------------------------------------------
# 6. Herbivores-------
# ---------------------------------------------------------------------------------
# legend--------
# A legend linking variable names from the HERBIVORES dataframes (herbiv_gr; herbiv_for)
# to the terminology used in the manuscript:
# NAME USED IN THE CODE -> NAME USED IN THE MS

# PlotID            -> Plot
# MeanTotAbu_H      -> Herbivore Mean Total Abundance (MeanTotAbu in figures)
# Stab_Abu_H        -> Herbivore community stability (computed on arthropod abundance)
# eta_w             -> Synchrony index
# weighted_avg_var  -> Weighted average population variability
# LUI               -> Land-use intensity index (grasslands)
# Formi             -> Land-use intensity index (forests)
# rao               -> Multi-Trait Rao's Q index for herbivore functional diversity (FD in figures)
# PC1               -> First axis of the PCA (representing the dominant ecological strategy)
# PC2               -> Second axis of the PCA (representing the dominant ecological strategy)
# MeanRich_H        -> Herbivore mean species richness
# MeanTreeCov       -> Mean plot tree cover in forests
# PlantFD           -> Multi-Trait Rao's Q index for plant functional diversity 
# PlantPC1          -> First axis of the PCA computed on plant functional traits (representing the dominant ecological strategy in plants)
# PlantPC2          -> Second axis of the PCA computed on plant functional traits (representing the dominant ecological strategy in plants)

# 6.1 grasslands--------------
# Fitting hypothesised and reduced SEMs -------------------------------------
# Fitting hypothesised model
acroSEM_herb_gr_hypo <- psem(
  lm(PlantPC1 ~ LUI, data = herbiv_gr),
  lm(PlantPC2 ~ LUI, data = herbiv_gr),
  lm(PlantFD ~ LUI, data = herbiv_gr),
  lm(PC1 ~ LUI + PlantFD + PlantPC1 + PlantPC2, data = herbiv_gr), 
  lm(PC2 ~ LUI + PlantFD + PlantPC1 + PlantPC2, data = herbiv_gr),
  lm(MeanRich_H ~ LUI + PlantFD + PlantPC1 + PlantPC2 + PC1 + PC2, data = herbiv_gr),  
  lm(rao ~ PC1 + PC2 + LUI + PlantFD + PlantPC1 + PlantPC2, data = herbiv_gr),  
  lm(eta_w ~ MeanRich_H + LUI + rao + PC1 + PC2 + PlantFD + PlantPC1 + PlantPC2, data = herbiv_gr),  
  lm(MeanTotAbu_H ~ MeanRich_H + rao + LUI + eta_w + PlantFD + PlantPC1 + PlantPC2 + PC1 + PC2, data = herbiv_gr), 
  lm(weighted_avg_var ~ LUI + MeanRich_H + MeanTotAbu_H + PC1 + PC2 + rao + PlantFD + PlantPC1 + PlantPC2, data = herbiv_gr),
  lm(Stab_Abu_H ~ eta_w + MeanRich_H + MeanTotAbu_H + PC1 + PC2 + rao + LUI + weighted_avg_var + PlantFD + PlantPC1 + PlantPC2, data = herbiv_gr), 
  rao %~~% MeanRich_H,
  PlantPC1 %~~% PlantFD,
  PlantPC2 %~~% PlantFD,
  data = herbiv_gr)

# Check model fit
summary(acroSEM_herb_gr_hypo)

# Model selection for each of the models in the SEM
lapply(list(
  PLANTPC1 = lm(PlantPC1 ~ LUI, na.action = "na.fail", data = herbiv_gr),
  PLANTPC2 =lm(PlantPC2 ~ LUI, na.action = "na.fail", data = herbiv_gr),
  PLANTFD = lm(PlantFD ~ LUI, na.action = "na.fail", data = herbiv_gr),
  PC1 = lm(PC1 ~ LUI + PlantFD + PlantPC1 + PlantPC2, na.action = "na.fail", data = herbiv_gr),
  PC2 = lm(PC2 ~ LUI + PlantFD + PlantPC1 + PlantPC2, na.action = "na.fail", data = herbiv_gr),
  RICH = lm(MeanRich_H ~ LUI + PlantFD + PlantPC1 + PlantPC2 + PC1 + PC2, na.action = "na.fail", data = herbiv_gr),
  RAO = lm(rao ~ PC1 + PC2 + LUI + PlantFD + PlantPC1 + PlantPC2, na.action = "na.fail", data = herbiv_gr),
  SYNC = lm(eta_w ~ MeanRich_H + LUI + rao + PC1 + PC2 + PlantFD + PlantPC1 + PlantPC2, na.action = "na.fail", data = herbiv_gr),
  MEANTOTABU = lm(MeanTotAbu_H ~ MeanRich_H + rao + LUI + eta_w + PlantFD + PlantPC1 + PlantPC2 + PC1 + PC2, na.action = "na.fail", data = herbiv_gr),
  AVGPOPVAR = lm(weighted_avg_var ~ LUI + MeanRich_H + MeanTotAbu_H + PC1 + PC2 + rao + PlantFD + PlantPC1 + PlantPC2, na.action = "na.fail", data = herbiv_gr),
  STAB = lm(Stab_Abu_H ~ eta_w + MeanRich_H + MeanTotAbu_H + PC1 + PC2 + rao + LUI + weighted_avg_var + PlantFD + PlantPC1 + PlantPC2, na.action = "na.fail", data = herbiv_gr)),
  Var_to_incl)

# Fitting reduced model
acroSEM_herb_gr_rd <- psem(
  lm(PlantPC1 ~ LUI, data = herbiv_gr),
  #lm(PlantPC2 ~ 1, data = herbiv_gr), # Uncomment to avoid fitting intercept-only model
  lm(PlantFD ~ LUI, data = herbiv_gr),
  lm(PC1 ~ PlantPC2, data = herbiv_gr), 
  #lm(PC2 ~ 1, data = herbiv_gr), # Uncomment to avoid fitting intercept-only model
  lm(MeanRich_H ~ PC2 + LUI + PlantPC1 + PlantPC2, data = herbiv_gr),  
  lm(rao ~ PC2 + PC1 + PlantPC1 + PlantPC2, data = herbiv_gr),  
  lm(eta_w ~ rao + MeanRich_H, data = herbiv_gr),  
  lm(MeanTotAbu_H ~ MeanRich_H + eta_w + rao + LUI, data = herbiv_gr), 
  lm(weighted_avg_var ~ MeanRich_H + PC1 + MeanTotAbu_H + PC2 + rao + PlantPC1, data = herbiv_gr),
  lm(Stab_Abu_H ~ eta_w + weighted_avg_var + MeanTotAbu_H, data = herbiv_gr), 
  rao %~~% MeanRich_H,
  PlantPC1 %~~% PlantFD,
  PlantPC2 %~~% PlantFD,
  data = herbiv_gr)

# Check model fit
summary(acroSEM_herb_gr_rd) 

# comparing AICc values
AIC_psem(acroSEM_herb_gr_hypo) 
AIC_psem(acroSEM_herb_gr_rd) 

# Bootstrapping and calculating standardised effects------ 
# Bootstrapping model effects with 10,000 resamples for the reduced SEM
acroSEM_herb_gr_boot <- bootEff(acroSEM_herb_gr_rd, R = 10000, catch.err = FALSE, parallel = "snow", seed = 23)
# Calculating standardised effects (direct, indirect, total) for endogenous variables in the bootstrapped model
acroSEM_herb_gr_boot_eff <- semEff(acroSEM_herb_gr_boot)

# Summary of effects and confidence intervals for endogenous variables
summary(acroSEM_herb_gr_boot_eff)
# Same as above, for the stability model only
summary(acroSEM_herb_gr_boot_eff, responses = "Stab.Abu.H") 
# Effect of LUI on stability and mediating effect of functional features 
summary(semEff(acroSEM_herb_gr_boot, predictor = "LUI")) 

# 6.2 forests--------------
# Fitting hypothesised and reduced SEMs -------------------------------------
# Fitting hypothesised model
acroSEM_herb_for_hypo <- psem(
  lm(PlantPC1 ~ Formi + MeanTreeCov, data = herbiv_for),
  lm(PlantPC2 ~ Formi + MeanTreeCov, data = herbiv_for),
  lm(PlantFD ~ Formi + MeanTreeCov, data = herbiv_for),
  lm(PC1 ~ Formi + PlantFD + PlantPC1 + PlantPC2 + MeanTreeCov, data = herbiv_for), 
  lm(PC2 ~ Formi + PlantFD + PlantPC1 + PlantPC2 + MeanTreeCov, data = herbiv_for),
  lm(MeanRich_H ~ Formi + PC1 + PC2 + PlantFD + PlantPC1 + PlantPC2 + MeanTreeCov, data = herbiv_for), 
  lm(rao ~ PC1 + PC2 + Formi + PlantFD + PlantPC1 + PlantPC2 + MeanTreeCov, data = herbiv_for), 
  lm(eta_w ~  MeanRich_H + Formi + rao + PC1 + PC2 + PlantFD + PlantPC1 + PlantPC2 + MeanTreeCov, data = herbiv_for),
  lm(MeanTotAbu_H ~ MeanRich_H + rao + Formi + eta_w + PlantFD + PlantPC1 + PlantPC2 + PC1 + PC2 + MeanTreeCov, data = herbiv_for), 
  lm(weighted_avg_var ~ Formi + MeanRich_H + MeanTotAbu_H + PC1 + PC2 + rao + PlantFD + PlantPC1 + PlantPC2 + MeanTreeCov, data = herbiv_for), 
  lm(Stab_Abu_H ~ eta_w + MeanRich_H + MeanTotAbu_H + PC1 + PC2 + rao + Formi + weighted_avg_var + PlantFD + PlantPC1 + PlantPC2 + MeanTreeCov, data = herbiv_for), 
  rao %~~% MeanRich_H,
  PlantPC2 %~~% PlantFD, # included based on d-separation test
  data = herbiv_for)

# Check model fit
summary(acroSEM_herb_for_hypo) 

# Model selection for each of the models in the SEM
lapply(list(
  PlantPC1 = lm(PlantPC1 ~ Formi + MeanTreeCov, na.action = "na.fail", data = herbiv_for),
  PlantPC2 = lm(PlantPC2 ~ Formi + MeanTreeCov, na.action = "na.fail", data = herbiv_for),
  PlantFD = lm(PlantFD ~ Formi + MeanTreeCov, na.action = "na.fail", data = herbiv_for),
  PC1 = lm(PC1 ~ Formi + PlantFD + PlantPC1 + PlantPC2 + MeanTreeCov, na.action = "na.fail", data = herbiv_for), 
  PC2 = lm(PC2 ~ Formi + PlantFD + PlantPC1 + PlantPC2 + MeanTreeCov, na.action = "na.fail", data = herbiv_for), 
  RICH = lm(MeanRich_H ~ Formi + PC1 + PC2 + PlantFD + PlantPC1 + PlantPC2 + MeanTreeCov, na.action = "na.fail", data = herbiv_for), 
  RAO = lm(rao ~ PC1 + PC2 + Formi + PlantFD + PlantPC1 + PlantPC2 + MeanTreeCov, na.action = "na.fail", data = herbiv_for),
  SYNC = lm(eta_w ~  MeanRich_H + Formi + rao + PC1 + PC2 + PlantFD + PlantPC1 + PlantPC2 + MeanTreeCov, na.action = "na.fail", data = herbiv_for), 
  MEANTOTABUND = lm(MeanTotAbu_H ~ MeanRich_H + rao + Formi + eta_w + PlantFD + PlantPC1 + PlantPC2 + PC1 + PC2 + MeanTreeCov, na.action = "na.fail", data = herbiv_for), 
  AVG_POP_VAR = lm(weighted_avg_var ~ Formi + MeanRich_H + MeanTotAbu_H + PC1 + PC2 + rao + PlantFD + PlantPC1 + PlantPC2 + MeanTreeCov, na.action = "na.fail", data = herbiv_for),
  STAB = lm(Stab_Abu_H ~ eta_w + MeanRich_H + MeanTotAbu_H + PC1 + PC2 + rao + Formi + weighted_avg_var + PlantFD + PlantPC1 + PlantPC2 + MeanTreeCov, na.action = "na.fail", data = herbiv_for)),
  Var_to_incl)

# Fitting reduced model
acroSEM_herb_for_rd <- psem(
  # lm(PlantPC1 ~ 1, data = herbiv_for), # Uncomment to avoid fitting intercept-only model 
  # lm(PlantPC2 ~ 1, data = herbiv_for), # Uncomment to avoid fitting intercept-only model 
  lm(PlantFD ~ Formi, data = herbiv_for),
  #lm(PC1 ~ 1, data = herbiv_for), # Uncomment to avoid fitting intercept-only model 
  lm(PC2 ~ PlantFD, data = herbiv_for),
  lm(MeanRich_H ~ PC1, data = herbiv_for), # dead end after removing MeanTotAbu_H: removed from figure
  lm(rao ~ PC1, data = herbiv_for), # dead end after removing MeanTotAbu_H: removed from figure
  lm(eta_w ~ PlantPC2 + PC1, data = herbiv_for), 
  lm(MeanTotAbu_H ~ PC2 + MeanRich_H + rao, data = herbiv_for), # now dead end - removed from figure
  lm(weighted_avg_var ~ Formi, data = herbiv_for), 
  lm(Stab_Abu_H ~ weighted_avg_var + eta_w + PC2 + PC1, data = herbiv_for), 
  rao %~~% MeanRich_H,
  PlantPC2 %~~% PlantFD, # included based on d-separation test
  data = herbiv_for)

# Check model fit
summary(acroSEM_herb_for_rd)

# comparing AICc values
AIC_psem(acroSEM_herb_for_hypo) 
AIC_psem(acroSEM_herb_for_rd)

# Bootstrapping and calculating standardised effects------ 
# Bootstrapping model effects with 10,000 resamples for the reduced SEM
acroSEM_herb_for_boot <- bootEff(acroSEM_herb_for_rd, R = 10000, catch.err = FALSE, parallel = "snow", seed = 23)
# Calculating standardised effects (direct, indirect, total) for endogenous variables in the bootstrapped model
acroSEM_herb_for_boot_eff <- semEff(acroSEM_herb_for_boot)
# Summary of effects and confidence intervals for endogenous variables
summary(acroSEM_herb_for_boot_eff) 
# Same as above, for the stability model only
summary(acroSEM_herb_for_boot_eff, responses = "Stab.Abu.H") 
# Effect of Formi on stability and mediating effect of functional features 
summary(semEff(acroSEM_herb_for_boot, predictor = "Formi")) 

# ---------------------------------------------------------------------------------
# 7. Carnivores---------
# ---------------------------------------------------------------------------------
# legend--------
# A legend linking variable names from the CARNIVORES dataframes (carniv_gr; carniv_for)
# to the terminology used in the manuscript:
# NAME USED IN THE CODE -> NAME USED IN THE MS

# PlotID            -> Plot
# MeanTotAbu_Carn   -> Carnivore Mean total Abundance (MeanTotAbu in figures)
# Stab_Abu_Carn     -> Carnivore community stability (computed on carnivore abundance)
# eta_w             -> Synchrony index
# weighted_avg_var  -> Weighted average population variability
# LUI               -> Land-use intensity index (grasslands)
# Formi             -> Land-use intensity index (forests)
# rao               -> Multi-Trait Rao's Q index for carnivore functional diversity (FD in figures)
# PC1               -> First axis of the PCA (representing the dominant ecological strategy)
# PC2               -> Second axis of the PCA (representing the dominant ecological strategy)
# MeanRich_Carn     -> Carnivore mean species richness
# MeanTreeCov       -> Mean plot tree cover in forests
# HerbPC1           -> First axis of the PCA computed on herbivores functional traits (representing the dominant ecological strategy in herbivores)
# HerbPC2           -> Second axis of the PCA computed on herbivores functional traits (representing the dominant ecological strategy in herbivores)
# Herbrao           -> Multi-Trait Rao's Q index for herbivore functional diversity 
# PlantFD           -> Multi-Trait Rao's Q index for plant functional diversity 
# PlantPC1          -> First axis of the PCA computed on plant functional traits (representing the dominant ecological strategy in plants)
# PlantPC2          -> Second axis of the PCA computed on plant functional traits (representing the dominant ecological strategy in plants)

# 7.1 grasslands--------------
# Fitting hypothesised and reduced SEMs -------------------------------------
# Fitting hypothesised model
acroSEM_carn_gr_hypo <- psem(
  lm(PlantPC1 ~ LUI, data = carniv_gr),
  lm(PlantPC2 ~ LUI, data = carniv_gr),
  lm(PlantFD ~ LUI, data = carniv_gr),
  lm(HerbPC1 ~ PlantFD + PlantPC1 + PlantPC2 + LUI, data = carniv_gr),
  lm(HerbPC2 ~ PlantFD + PlantPC1 + PlantPC2 + LUI, data = carniv_gr),
  lm(Herbrao ~ HerbPC1 + HerbPC2 + LUI + PlantFD + PlantPC1 + PlantPC2, data = carniv_gr),
  lm(PC1 ~ HerbPC1 + HerbPC2 + Herbrao + PlantPC1 + PlantPC2 + LUI, data = carniv_gr),  # PlantPC1 and PlantPC2 included based on d-separation test
  lm(PC2 ~ HerbPC1 + HerbPC2 + Herbrao + LUI, data = carniv_gr),
  lm(MeanRich_Carn ~ LUI + PC1 + PC2 + HerbPC1 + HerbPC2 + Herbrao + PlantPC1 + PlantPC2, data = carniv_gr), # PlantPC1 and PlantPC2 included based on d-separation test
  lm(rao ~ PC1 + PC2 + LUI + HerbPC1 + HerbPC2 + Herbrao, data = carniv_gr),  
  lm(eta_w ~ MeanRich_Carn + LUI + rao + PC1 + PC2 + HerbPC1 + HerbPC2 + Herbrao, data = carniv_gr),  
  lm(MeanTotAbu_Carn ~ MeanRich_Carn + rao + LUI + eta_w + HerbPC1 + HerbPC2 + Herbrao + PC1 + PC2, data = carniv_gr), 
  lm(weighted_avg_var ~ LUI + MeanRich_Carn + MeanTotAbu_Carn + PC1 + PC2 + rao + HerbPC1 + HerbPC2 + Herbrao, data = carniv_gr),
  lm(Stab_Abu_Carn ~ eta_w + MeanRich_Carn + MeanTotAbu_Carn + PC1 + PC2 + rao + LUI + weighted_avg_var + HerbPC1 + HerbPC2 + Herbrao, data = carniv_gr),
  rao %~~% MeanRich_Carn,
  PlantPC1 %~~% PlantFD, # included based on d-separation test
  PlantFD %~~% PlantPC2, # included based on d-separation test
  data = carniv_gr)

# Check model fit
summary(acroSEM_carn_gr_hypo) 

# Model selection for each of the models in the SEM
lapply(list(
  PlantPC1 = lm(PlantPC1 ~ LUI, na.action = "na.fail", data = carniv_gr),
  PlantPC2 = lm(PlantPC2 ~ LUI, na.action = "na.fail", data = carniv_gr),
  PlantFD = lm(PlantFD ~ LUI, na.action = "na.fail", data = carniv_gr),
  HerbPC1 = lm(HerbPC1 ~ PlantFD + PlantPC1 + PlantPC2 + LUI, na.action = "na.fail", data = carniv_gr),
  HerbPC2 = lm(HerbPC2 ~ PlantFD + PlantPC1 + PlantPC2 + LUI, na.action = "na.fail", data = carniv_gr),
  Herbrao = lm(Herbrao ~ HerbPC1 + HerbPC2 + LUI + PlantFD + PlantPC1 + PlantPC2, na.action = "na.fail", data = carniv_gr),
  PC1 = lm(PC1 ~ HerbPC1 + HerbPC2 + Herbrao + PlantPC1 + PlantPC2 + LUI, na.action = "na.fail", data = carniv_gr),
  PC2 = lm(PC2 ~ HerbPC1 + HerbPC2 + Herbrao + LUI, na.action = "na.fail", data = carniv_gr),# 
  RICH = lm(MeanRich_Carn ~ LUI + PC1 + PC2 + HerbPC1 + HerbPC2 + Herbrao + PlantPC1 + PlantPC2, na.action = "na.fail", data = carniv_gr), 
  RAO = lm(rao ~ PC1 + PC2 + LUI + HerbPC1 + HerbPC2 + Herbrao, na.action = "na.fail", data = carniv_gr),
  SYNC = lm(eta_w ~ MeanRich_Carn + LUI + rao + PC1 + PC2 + HerbPC1 + HerbPC2 + Herbrao, na.action = "na.fail", data = carniv_gr), 
  MEANTOTABUND = lm(MeanTotAbu_Carn ~ MeanRich_Carn + rao + LUI + eta_w + HerbPC1 + HerbPC2 + Herbrao + PC1 + PC2, na.action = "na.fail", data = carniv_gr),
  AVG_POP_VAR = lm(weighted_avg_var ~ LUI + MeanRich_Carn + MeanTotAbu_Carn + PC1 + PC2 + rao + HerbPC1 + HerbPC2 + Herbrao, na.action = "na.fail", data = carniv_gr),
  STAB = lm(Stab_Abu_Carn ~ eta_w + MeanRich_Carn + MeanTotAbu_Carn + PC1 + PC2 + rao + LUI + weighted_avg_var + HerbPC1 + HerbPC2 + Herbrao, na.action = "na.fail", data = carniv_gr)),
  Var_to_incl)

# Fitting reduced model
acroSEM_carn_gr_rd <- psem(
  lm(PlantPC1 ~ LUI, data = carniv_gr),
  #lm(PlantPC2 ~ 1, data = carniv_gr), # Uncomment to avoid fitting intercept-only model
  lm(PlantFD ~ LUI, data = carniv_gr), # dead-end: removed from figure
  lm(HerbPC1 ~ PlantPC2, data = carniv_gr),
  #lm(HerbPC2 ~ 1, data = carniv_gr), # Uncomment to avoid fitting intercept-only model
  lm(Herbrao ~ HerbPC2 + HerbPC1 + PlantPC1 + PlantPC2, data = carniv_gr),
  lm(PC1 ~ PlantPC2 + PlantPC1 + Herbrao, data = carniv_gr),  
  lm(PC2 ~ HerbPC2 + LUI, data = carniv_gr),
  lm(MeanRich_Carn ~ PlantPC2 + LUI + PlantPC1 + PC2, data = carniv_gr),
  lm(rao ~ PC1, data = carniv_gr),  
  lm(eta_w ~ MeanRich_Carn, data = carniv_gr),  
  lm(MeanTotAbu_Carn ~ MeanRich_Carn + rao + PC1, data = carniv_gr), 
  lm(weighted_avg_var ~ MeanRich_Carn + LUI + rao + PC1, data = carniv_gr),
  lm(Stab_Abu_Carn ~ eta_w + MeanRich_Carn + MeanTotAbu_Carn + weighted_avg_var, data = carniv_gr), # eta_w included based on d-separation test
  rao %~~% MeanRich_Carn,
  PlantPC1 %~~% PlantFD, # included based on d-separation test
  PlantPC2 %~~% PlantFD, # included based on d-separation test
  data = carniv_gr)

# Check model fit
summary(acroSEM_carn_gr_rd) 

# comparing AICc values
AIC_psem(acroSEM_carn_gr_hypo) 
AIC_psem(acroSEM_carn_gr_rd) 

# Bootstrapping and calculating standardised effects------ 
# Bootstrapping model effects with 10,000 resamples for the reduced SEM--
acroSEM_carn_gr_boot <- bootEff(acroSEM_carn_gr_rd, R = 10000, catch.err = FALSE, parallel = "snow", seed = 23)
# Calculating standardised effects (direct, indirect, total) for endogenous variables in the bootstrapped model
acroSEM_carn_gr_boot_eff <- semEff(acroSEM_carn_gr_boot)

# Summary of effects and confidence intervals for endogenous variables
summary(acroSEM_carn_gr_boot_eff) 
# Same as above, for the stability model only
summary(acroSEM_carn_gr_boot_eff, responses = "Stab.Abu.Carn") 
# Effect of LUI on stability and mediating effect of functional features 
summary(semEff(acroSEM_carn_gr_boot, predictor = "LUI")) 


# 7.2 forests--------------
# Fitting hypothesised and reduced SEMs -------------------------------------
# Fitting hypothesised model
acroSEM_carn_for_hypo <- psem(
  lm(PlantPC1 ~ Formi + MeanTreeCov, data = carniv_for),
  lm(PlantPC2 ~ Formi + MeanTreeCov, data = carniv_for),
  lm(PlantFD ~ Formi + MeanTreeCov, data = carniv_for),
  lm(HerbPC1 ~ PlantFD + PlantPC1 + PlantPC2 + Formi + MeanTreeCov, data = carniv_for),
  lm(HerbPC2 ~ PlantFD + PlantPC1 + PlantPC2 + Formi + MeanTreeCov, data = carniv_for),
  lm(Herbrao ~ HerbPC1 + HerbPC2 + Formi + PlantFD + PlantPC1 + PlantPC2 + MeanTreeCov, data = carniv_for),
  lm(PC1 ~ HerbPC1 + HerbPC2 + Herbrao + Formi + MeanTreeCov, data = carniv_for),   
  lm(PC2 ~ HerbPC1 + HerbPC2 + Herbrao + Formi + MeanTreeCov, data = carniv_for),
  lm(MeanRich_Carn ~ Formi + PC1 + PC2 + HerbPC1 + HerbPC2 + Herbrao + MeanTreeCov, data = carniv_for), 
  lm(rao ~ PC1 + PC2 + Formi + HerbPC1 + HerbPC2 + Herbrao + MeanTreeCov, data = carniv_for),  
  lm(eta_w ~ MeanRich_Carn + Formi + rao + PC1 + PC2 + HerbPC1 + HerbPC2 + Herbrao + MeanTreeCov + PlantFD + PlantPC2, data = carniv_for), # PlantFD and PlantPC2 included based on d-separation test
  lm(MeanTotAbu_Carn ~ MeanRich_Carn + rao + Formi + eta_w + HerbPC1 + HerbPC2 + Herbrao + PC1 + PC2 + MeanTreeCov, data = carniv_for), 
  lm(weighted_avg_var ~ Formi + MeanRich_Carn + MeanTotAbu_Carn + PC1 + PC2 + rao + HerbPC1 + HerbPC2 + Herbrao + MeanTreeCov, data = carniv_for),
  lm(Stab_Abu_Carn ~ eta_w + MeanRich_Carn + MeanTotAbu_Carn + PC1 + PC2 + rao + Formi + weighted_avg_var + HerbPC1 + HerbPC2 + Herbrao + MeanTreeCov, data = carniv_for),
  rao %~~% MeanRich_Carn,
  PlantPC1 %~~% PlantFD, # included based on d-separation test
  PlantPC2 %~~% PlantFD, # included based on d-separation test
  data = carniv_for)

# Check model fit
summary(acroSEM_carn_for_hypo) 

# Model selection for each of the models in the SEM
lapply(list(
  PlantPC1 = lm(PlantPC1 ~ Formi + MeanTreeCov, na.action = "na.fail", data = carniv_for),
  PlantPC2 = lm(PlantPC2 ~ Formi + MeanTreeCov, na.action = "na.fail", data = carniv_for),
  PlantFD = lm(PlantFD ~ Formi + MeanTreeCov, na.action = "na.fail", data = carniv_for),
  HerbPC1 = lm(HerbPC1 ~ PlantFD + PlantPC1 + PlantPC2 + Formi + MeanTreeCov, na.action = "na.fail", data = carniv_for),
  HerbPC2 = lm(HerbPC2 ~ PlantFD + PlantPC1 + PlantPC2 + Formi + MeanTreeCov, na.action = "na.fail", data = carniv_for),
  Herbrao = lm(Herbrao ~ HerbPC1 + HerbPC2 + Formi + PlantFD + PlantPC1 + PlantPC2 + MeanTreeCov, na.action = "na.fail", data = carniv_for),
  PC1 = lm(PC1 ~ HerbPC1 + HerbPC2 + Herbrao + Formi + MeanTreeCov, na.action = "na.fail", data = carniv_for),
  PC2 = lm(PC2 ~ HerbPC1 + HerbPC2 + Herbrao + Formi + MeanTreeCov, na.action = "na.fail", data = carniv_for), 
  RICH = lm(MeanRich_Carn ~ Formi + PC1 + PC2 + HerbPC1 + HerbPC2 + Herbrao + MeanTreeCov, na.action = "na.fail", data = carniv_for), 
  RAO = lm(rao ~ PC1 + PC2 + Formi + HerbPC1 + HerbPC2 + Herbrao + MeanTreeCov, na.action = "na.fail", data = carniv_for),
  SYNC = lm(eta_w ~ MeanRich_Carn + Formi + rao + PC1 + PC2 + HerbPC1 + HerbPC2 + Herbrao + MeanTreeCov + PlantFD + PlantPC2, na.action = "na.fail", data = carniv_for), 
  MEANTOTABUND = lm(MeanTotAbu_Carn ~ MeanRich_Carn + rao + Formi + eta_w + HerbPC1 + HerbPC2 + Herbrao + PC1 + PC2 + MeanTreeCov, na.action = "na.fail", data = carniv_for), 
  AVG_POP_VAR = lm(weighted_avg_var ~ Formi + MeanRich_Carn + MeanTotAbu_Carn + PC1 + PC2 + rao + HerbPC1 + HerbPC2 + Herbrao + MeanTreeCov, na.action = "na.fail", data = carniv_for),
  STAB = lm(Stab_Abu_Carn ~ eta_w + MeanRich_Carn + MeanTotAbu_Carn + PC1 + PC2 + rao + Formi + weighted_avg_var + HerbPC1 + HerbPC2 + Herbrao + MeanTreeCov, na.action = "na.fail", data = carniv_for)),
  Var_to_incl)

# Fitting reduced model
acroSEM_carn_for_rd <- psem(
  #lm(PlantPC1 ~ 1, data = carniv_for), # Uncomment to avoid fitting intercept-only model
  #lm(PlantPC2 ~ 1, data = carniv_for), # Uncomment to avoid fitting intercept-only model
  lm(PlantFD ~ Formi, data = carniv_for), 
  #lm(HerbPC1 ~ 1, data = carniv_for), # Uncomment to avoid fitting intercept-only model
  lm(HerbPC2 ~ PlantFD, data = carniv_for),
  lm(Herbrao ~ HerbPC1, data = carniv_for),
  # lm(PC1 ~ 1, data = carniv_for), # Uncomment to avoid fitting intercept-only model
  lm(PC2 ~ HerbPC2, data = carniv_for),  
  lm(MeanRich_Carn ~ Herbrao, data = carniv_for), 
  #lm(rao ~ 1, data = carniv_for), # Uncomment to avoid fitting intercept-only model
  lm(eta_w ~ PlantPC2, data = carniv_for), 
  lm(MeanTotAbu_Carn ~ MeanRich_Carn + Formi, data = carniv_for), 
  lm(weighted_avg_var ~ HerbPC1 + MeanTotAbu_Carn, data = carniv_for),
  lm(Stab_Abu_Carn ~ eta_w + weighted_avg_var + MeanTotAbu_Carn + MeanRich_Carn, data = carniv_for), 
  #rao %~~% MeanRich_Carn, 
  data = carniv_for)

# Check model fit
summary(acroSEM_carn_for_rd) 

# cmparing AICc values
AIC_psem(acroSEM_carn_for_hypo) 
AIC_psem(acroSEM_carn_for_rd) 

# Bootstrapping and calculating standardised effects------ 
# Bootstrapping model effects with 10,000 resamples for the reduced SEM
acroSEM_carn_for_boot <- bootEff(acroSEM_carn_for_rd, R = 10000, catch.err = FALSE, parallel = "snow", seed = 23)
# Calculating standardised effects (direct, indirect, total) for endogenous variables in the bootstrapped model
acroSEM_carn_for_boot_eff <- semEff(acroSEM_carn_for_boot)

# Summary of effects and confidence intervals for endogenous variables
summary(acroSEM_carn_for_boot_eff) 
# Same as above, for the stability model only
summary(acroSEM_carn_for_boot_eff, responses = "Stab.Abu.Carn") 
# Effect of Formi on stability and mediating effect of functional features 
summary(semEff(acroSEM_carn_for_boot, predictor = "Formi")) # effect of Formi on stability + mediating eff of functional features

