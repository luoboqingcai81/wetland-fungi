library(tidyverse)
library(lmerTest)
library(lme4)
library(glmm.hp)
library(effectsize)
### load R function

glmmHP_tidy <- function(data,property,output){
  
  envs <- c("CHELSA_BIO_Annual_Mean_Temperature","CHELSA_BIO_Temperature_Annual_Range","CHELSA_BIO_Temperature_Seasonality",
            "CHELSA_BIO_Annual_Precipitation","CHELSA_BIO_Precipitation_Seasonality","CHELSA_BIO_MeanMonthlyMoisture_Index",
            "CHELSA_BIO_AridityIndex","CHELSA_BIO_MeanMonthlyPotential_Evapotranspiration",
            "EarthEnv_Texture_Shannon_Index","CHELSA_BIO_NPP",
            "SG_Nitrogen_Content_0_5cm","SG_SOC_Content_0_5cm","SG_Soil_pH_0_5cm")
  
  envs_tibble <- tibble(Varibles = envs,
                        Type = rep(c("temperature","precipitation","aridityindex","plant","soil"), c(3,3,2,2,3)))
  
  Data <<-  read_csv(data) %>% 
    select(-`.geo`, - `system:index`) %>% 
    select(property, all_of(envs), Paper) %>% 
    drop_na() %>% 
    mutate(across(all_of(envs), ~ scale(.x)[,1]))
  
  formula_raw = as.formula(
    paste(property,"~",
          paste(c(envs,"(1|Paper)"),collapse = "+")
    ))
  
  lmer_raw_results <- lme4::lmer(formula = formula_raw,data = Data)
  property_factor <- effectsize(lmer_raw_results)[-1,] %>% as_tibble()
  
  lmer_gmmhp <<- lme4::lmer(formula = formula_raw,data = Data)
  lmer_gmmhp_results <- glmm.hp(lmer_gmmhp)
  a <- lmer_gmmhp_results$hierarchical.partitioning[,4] 
  A <- lmer_gmmhp_results$r.squaredGLMM
  
  aa <- a %>%
    base::as.data.frame() %>%
    tibble::rownames_to_column(var = "Varibles") %>%
    dplyr::rename(percent = 2)
  print(aa)
  all_results <- envs_tibble %>% left_join(property_factor, by =join_by(Varibles == Parameter)) %>%
    left_join(aa) %>%
    mutate(R2m = A[1,1], R2c = A[1,2])%>%
    write_csv(output)
}

all_gee_table_name_byPaper <-
  c("asv_process/env/GEE_point_byPaper/RawEnv_unique_site_richness_without_tundra_byPaper.csv",
    "asv_process/env/GEE_point_byPaper/RawEnv_unique_site_FunB_NoTundra_byPaper.csv")

all_classProperty <-
  c("richness","FunB")

all_output <-c("asv_process/results/5type_glmmhp_fungi_richness.csv",
    "asv_process/results/5type_glmmhp_FunB.csv")

pwalk(list(all_gee_table_name_byPaper, all_classProperty, all_output), glmmHP_tidy)



