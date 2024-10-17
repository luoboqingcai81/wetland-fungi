##### 加载R包
library(tidyverse)
library(sf)
library(terra)
library(purrr)
library(stars)
library(caret)
library(s2)
library(CAST)
library(fs)
library(parallel)
source("script/process_script/map_function.R")

covariable <-read_csv("asv_process/env/preds_names.csv") %>% pull(1)
### process model predict for richness
modelname = "asv_process/05.spatialcv_predict/05.spatialcv_richness_predict"
classProperty <- "richness"
svs_ffs_predict(
  modelname = modelname,
  classProperty = classProperty,
  GEE_point_env = "asv_process/env/GEE_point/richness_fc_agg_unique_NoTundra.csv",
  covariable
)

AOA_global_predict(
  modelname = modelname,
  classProperty = classProperty,
  variable_tiff = "~/Downloads/select_tiff2/all_predictors_renamed_masked.tiff"
)

### process model predict for shannon ##########
modelname2 = "asv_process/05.spatialcv_predict/06.spatialcv_shannon_predict"
classProperty <- "Shannon"
svs_ffs_predict(
  modelname = modelname2,
  classProperty = classProperty,
  GEE_point_env = "asv_process/env/GEE_point/Shannon_fc_agg_unique_NoTundra.csv",
  covariable
)

svs_model2 <-
  readRDS("asv_process/05.spatialcv_predict/06.spatialcv_shannon_predict/rfmodel.RDS")
svs_model2$selectedvars
svs_model2$results
coef_det(svs_model2$pred$obs, svs_model2$pred$pred)


### process model predict for Pathotrophs richness  ##########

modelname3 = "asv_process/05.spatialcv_predict/07.spatialcv_Pathotrophs_predict"
classProperty <- "richness"
svs_ffs_predict(
  modelname = modelname3,
  classProperty = classProperty,
  GEE_point_env = "asv_process/env/GEE_point/Pathotrophs_fc_agg_unique_NoTundra.csv",
  covariable
)

svs_model3 <-
  readRDS(
    "asv_process/05.spatialcv_predict/07.spatialcv_Pathotrophs_predict/rfmodel.RDS"
  )
svs_model3$selectedvars
svs_model3$results
coef_det(svs_model3$pred$obs, svs_model3$pred$pred)

### process model predict for Saprotrophs richness  ##########

modelname4 = "asv_process/05.spatialcv_predict/08.spatialcv_Saprotrophs_predict"
classProperty <- "richness"
svs_ffs_predict(
  modelname = modelname4,
  classProperty = classProperty,
  GEE_point_env = "asv_process/env/GEE_point/Saprotrophs_fc_agg_unique_NoTundra.csv",
  covariable
)

svs_model4 <-
  readRDS(
    "asv_process/05.spatialcv_predict/08.spatialcv_Saprotrophs_predict/rfmodel.RDS"
  )
svs_model4$selectedvars
svs_model4$results
coef_det(svs_model4$pred$obs, svs_model4$pred$pred)

### process model predict for Symbiotrophs richness  ##########

modelname5 = "asv_process/05.spatialcv_predict/09.spatialcv_Symbiotrophs_predict"
classProperty <- "richness"
svs_ffs_predict(
  modelname = modelname5,
  classProperty = classProperty,
  GEE_point_env = "asv_process/env/GEE_point/Symbiotrophs_fc_agg_unique_NoTundra.csv",
  covariable
)

svs_model5 <-
  readRDS(
    "asv_process/05.spatialcv_predict/09.spatialcv_Symbiotrophs_predict/rfmodel.RDS"
  )
svs_model5$selectedvars
svs_model5$results
coef_det(svs_model5$pred$obs, svs_model5$pred$pred)

############# all_predictors_bands_name ############
all_model_name <- c("asv_process/05.spatialcv_predict/05.spatialcv_richness_predict",
                    "asv_process/05.spatialcv_predict/06.spatialcv_shannon_predict",
                    "asv_process/05.spatialcv_predict/07.spatialcv_Pathotrophs_predict",
                    "asv_process/05.spatialcv_predict/08.spatialcv_Saprotrophs_predict",
                    "asv_process/05.spatialcv_predict/09.spatialcv_Symbiotrophs_predict",
                    "biomass/spatialcv_predict",
                    "biomass/FunB_spatialcv_predict")


all_model_name <- set_names(all_model_name) 

imap(all_model_name, var_importance) %>% 
  purrr::list_rbind(., names_to ="model") %>% 
  write_tsv(., file = "asv_process/results/all_model_results_importance.txt")
  
all_model_results <- purrr::map(all_model_name, model_tidy) %>% 
  purrr::list_rbind()

write_tsv(all_model_results, file = "asv_process/results/all_model_results.txt")

all_predictors_part1_BandName <- c(
  "EarthEnv_Texture_Variance_EVI",
  "SG_SOC_Density_0_5cm",
  "SBIO_0_5cm_Isothermality",
  "CHELSA_BIO_MeanSurfaceDownwellingShortwave",
  "CHELSA_BIO_Temperature_Annual_Range",
  "CHELSA_BIO_FrostChangeFrequency"
)

all_predictors_part2_BandName <-
  c(
    "EarthEnv_CloudCover_Mean",
    "CHELSA_BIO_Precipitation_of_Driest_Month",
    "CHELSA_BIO_MaximumMonthlyPotential_Evapotranspiration",
    "SBIO_0_5cm_Mean_Diurnal_Range",
    "CHELSA_BIO_Precipitation_of_Driest_Quarter"
  )

all_predictors_part3_BandName <-
  c(
    "CHELSA_BIO_MinimumMonthlyVaporPressureDeficit",
    "SG_CEC_0_5cm",
    "CHELSA_BIO_MeanMonthlyVaporPressureDeficit",
    "CHELSA_BIO_MeanMonthlyMoisture_Index",
    "CHELSA_BIO_Precipitation_of_Coldest_Quarter",
    "SG_Nitrogen_Content_0_5cm",
    "CHELSA_BIO_GrowingDegreeDays_gt10degrees"
  )

all_predictors_part4_BandName <-
  c(
    "CHELSA_BIO_Annual_Precipitation",
    "SG_SOC_Content_0_5cm",
    "SG_Soil_pH_0_5cm"
  )

