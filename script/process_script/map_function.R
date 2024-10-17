
model_tidy <- function(model_name, output_file = NULL){
  model <- readRDS(paste0(model_name, "/rfmodel.RDS"))
  vars <- paste(model$selectedvars, collapse = ";")
  model_R2 <- model$results$Rsquared
  coef_R2 <- coef_det(model$pred$obs, model$pred$pred)
  Table <- tibble(name = model_name, factors = vars, model_R2 = model_R2, coef_R2 = coef_R2)
  Table
  if(!is.null(output_file)){
    readr::write_csv(Table, file = output_file)
  }
}

var_importance <- function(model_name, output_file = NULL){
  model <- readRDS(paste0(model_name, "/rfmodel.RDS"))
  feat_imp_df <- model$finalModel$variable.importance %>% 
    data.frame() %>% 
    rownames_to_column(var = "variable") %>%
    rename(importance = 2)
}

###Function borrowed and modified from ENMeval
get.block <- function(occs){
  rownames(occs) <- 1:nrow(occs)
  # SPLIT occs POINTS INTO FOUR SPATIAL GROUPS
  noccs <- nrow(occs)
  n1 <- ceiling(nrow(occs)/2)
  n2 <- floor(nrow(occs)/2)
  n3 <- ceiling(n1/2)
  n4 <- ceiling(n2/2)
  grpA <- occs[order(occs[, "Latitude"]),][1:n1,]
  grpB <- occs[rev(order(occs[, "Latitude"])),][1:n2,]
  grp1 <- grpA[order(grpA[, "Longitude"]),][1:(n3),]
  grp2 <- grpA[!rownames(grpA)%in%rownames(grp1),]
  grp3 <- grpB[order(grpB[, "Longitude"]),][1:(n4),]
  grp4 <- grpB[!rownames(grpB)%in%rownames(grp3),]
  r <- data.frame()
  if (nrow(grp1) > 0) grp1$grp <- 1; r <- rbind(r, grp1)
  if (nrow(grp2) > 0) grp2$grp <- 2; r <- rbind(r, grp2)
  if (nrow(grp3) > 0) grp3$grp <- 3; r <- rbind(r, grp3)
  if (nrow(grp4) > 0) grp4$grp <- 4; r <- rbind(r, grp4)
  occ.grp <- r[order(as.numeric(rownames(r))),]$grp
  return(occ.grp)
}

coef_det <- function(xtrue, xpred){
  return(1-sum((xtrue-xpred)^2)/sum((xtrue-mean(xtrue))^2))
}

Split_index <- function(split_num, split_Data, split_Var){
  split_Var <- sym(split_Var)
  index<-list()
  indexOut<-list()
  for (h in 1:split_num){
    indexOut[[h]]<-which(pull(split_Data, split_Var) == h)
    index[[h]]<-which(pull(split_Data, split_Var) != h)
  }
  indices<-list(index= index,indexOut= indexOut)
  indices
}

###### Creates index and indexOut for caret::trainControl from a fold vector
#' @param folds vector with fold labels
#' @return list, with index and indexOut
#' @author Marvin Ludwig

fold2index = function(fold){
  
  fold = data.frame(fold = fold)
  
  indOut = fold %>% dplyr::group_by(fold) %>%
    attr('groups') %>% dplyr::pull(.rows)
  
  ind = purrr::map(seq(length(indOut)), function(x){
    s = seq(nrow(fold))
    s = s[!s %in% indOut[[x]]]
    return(s)
  })
  return(
    list(
      index = ind,
      indexOut = indOut
    )
  )
  
}

spatial_variable_selection = function(modelname, training_samples, predictors, response, folds, hyperparameter, MinVar = 2){
  
  training_samples = training_samples %>% 
    dplyr::select(all_of(c(predictors, response))) %>% 
    st_drop_geometry()
  
  i = fold2index(folds)
  set.seed(4815)
  
  
  ffsModel = CAST::ffs(predictors = training_samples %>% dplyr::select(all_of(predictors)),
                       response = training_samples %>% pull(response),
                       method = "ranger",
                       minVar = MinVar,
                       tuneGrid = hyperparameter,
                       num.trees = 300,
                       trControl = caret::trainControl(method = "cv",number = 10,
                                                       index = i$index, indexOut = i$indexOut,
                                                       savePredictions = "final"),
                       importance = "impurity")
  
  saveRDS(ffsModel, paste0(modelname, "/rfmodel.RDS"))
  return(ffsModel)
  
}

create_folds_spatial = function(modelname, training_samples, n_folds, gridsize, seed){
  # create grid
  grid = world_grid(cellsize = gridsize)
  
  
  # spatially match points
  # remove grid cells without points
  grid = grid[lengths(st_intersects(grid, training_samples)) > 0,]
  
  # randomly create groups
  grid$fold = seq(nrow(grid)) %% n_folds
  set.seed(seed)
  grid$fold = sample(grid$fold)
  
  st_write(grid, paste0(modelname, "/spatial_folds_grid.gpkg"), append = FALSE)
  
  spatial_folds = st_join(training_samples, grid, left = TRUE) %>% pull(fold)
  
  saveRDS(spatial_folds, paste0(modelname, "/folds.RDS"))
  return(spatial_folds)
}

world_grid = function(extent = c(-180, -90, 180, 90), cellsize = 15, crs = 4326){
  world = sf::st_multipoint(x = matrix(extent, ncol = 2, byrow = TRUE)) %>%
    sf::st_sfc(crs = crs)
  world_grid = sf::st_make_grid(world, cellsize = cellsize)
  world_grid = sf::st_sf(fold = seq(length(world_grid)), world_grid)
  
  return(world_grid)
}

initModel = function(modelname){
  if(!dir.exists(modelname)){
    dir.create(modelname)
    dir.create(file.path(modelname, "results"))}}


# compute DI
pi_prediction = function(modelname, model, predictor_layers){
  # reduce predictor stack to model predictors
  mp = colnames(model$trainingData)[-length(model$trainingData)]
  predictor_layers = predictor_layers[[mp]]
  p = terra::predict(predictor_layers, model,na.rm =T)
  terra::writeRaster(p, file.path(modelname, "/prediction.tiff"), overwrite = TRUE)
  return(p)
}

pi_trainDI = function(modelname, model){
  tdi = CAST::trainDI(model)
  saveRDS(tdi, paste0(modelname, "/trainDI.RDS"))
  return(tdi)
}

pi_aoa = function(modelname, trainDI, model, predictor_layers){
  a = CAST::aoa(newdata = predictor_layers, model = model, trainDI = trainDI)
  saveRDS(a, paste0(modelname, "/aoa.RDS"))
  return(a)
}

mosaic_list <- function(list_rast){
  rsrc <- sprc(list_rast)
  m <- mosaic(rsrc)
  return(m)
}

# The function spatially aggregates the original raster
# it turns each aggregated cell into a polygon
# then the extent of each polygon is used to crop
# the original raster.
# The function returns a list with all the pieces
# in case you want to keep them in the memory. 
# it saves and plots each piece
# The arguments are:
# raster = raster to be chopped            (raster object)
# ppside = pieces per side                 (integer)
# save   = write raster                    (TRUE or FALSE)
# plot   = do you want to plot the output? (TRUE or FALSE)
SplitRas <- function(modelname,rast,ppside){
  h        <- ceiling(ncol(rast)/ppside)
  v        <- ceiling(nrow(rast)/ppside)
  agg      <- terra::aggregate(rast,fact=c(h,v))
  agg[]    <- 1:terra::ncell(agg)
  agg_poly <- terra::as.polygons(agg)
  names(agg_poly) <- "polis"
  r_list <- list()
  for(i in 1:terra::ncell(agg)){
    e1          <- terra::ext(agg_poly[agg_poly$polis==i,])
    cropped_raster <- terra::crop(rast,e1)
    terra::writeRaster(cropped_raster, paste0(modelname, "/results/tile_", i, ".tiff"))
  }
}

summary_tiff <- function(tiff_file, classProperty, modelname = NULL,degree = 1, Ext = c(-180,180,-60,90)){
  library(terra)
  library(dtplyr)
  library(dplyr)
  library(readr)
  if(is.null(modelname)){
    Path <- ""
  }else{
    Path <- paste0(modelname,"/result/")
  }
  preAOA <- terra::rast(tiff_file)
  df_preAOA <- terra::as.data.frame(preAOA, xy=T, na.rm = T)
  colnames(df_preAOA) <- c("Longitude","Latitude", "ClassProperty")
  write_csv(df_preAOA, file = paste0(Path,"preAOA_",classProperty,".csv"))
  
  dt_preAOA <- lazy_dt(df_preAOA)
  
  ### for Latitude data summarized by the degree
  breaks_Lat <- seq(Ext[3], Ext[4], degree)
  Label_breaks_Lat <- breaks_Lat[-length(breaks_Lat)]
  dt_preAOA %>%  
    mutate(Lat = cut(Latitude, breaks = breaks_Lat, labels = Label_breaks_Lat, 
                     include.lowest = TRUE) |> as.character() |> as.numeric()) %>% 
    group_by(Lat) %>% 
    dplyr::summarise(Mean = mean(ClassProperty),
                     Sd = sd(ClassProperty),
                     n =n()) %>% 
    mutate(Sd = replace_na(Sd,0)) %>% 
    mutate(Min = Mean - Sd,
           Max = Mean + Sd) %>% 
    as_tibble() %>% 
    write_csv(paste0(Path,"Latitude_by_", degree,"_",classProperty,".csv"))
  
  ### for Longitude data summarized by the degree
  breaks_Long <- seq(Ext[1], Ext[2], degree)
  Label_breaks_Long <- breaks_Long[-length(breaks_Long)]
  dt_preAOA %>%  
    mutate(Long = cut(Longitude, breaks = breaks_Long, labels = Label_breaks_Long, 
                      include.lowest = TRUE) |> as.character() |> as.numeric()) %>% 
    group_by(Long) %>% 
    dplyr::summarise(Mean = mean(ClassProperty),
                     Sd = sd(ClassProperty),
                     n = n()) %>% 
    mutate(Sd = replace_na(Sd,0)) %>% 
    mutate(Min = Mean - Sd,
           Max = Mean + Sd) %>% 
    as_tibble() %>% 
    write_csv(paste0(Path,"Longitude_by_", degree,"_",classProperty,".csv"))
  
  
  ### for Latitude data summarized by raw point
  dt_preAOA %>%
    group_by(Latitude) %>% 
    dplyr::summarise(Mean = mean(ClassProperty),
                     Sd = sd(ClassProperty),
                     n=n()) %>% 
    mutate(Sd = replace_na(Sd,0)) %>%  
    mutate(Min = Mean - Sd,
           Max = Mean + Sd) %>% 
    as_tibble() %>% 
    write_csv(paste0(Path,"Latitude_",classProperty,".csv"))
  
  ### for Longitude data summarized by raw point
  dt_preAOA %>%                                                               
    group_by(Longitude) %>% 
    dplyr::summarise(Mean = mean(ClassProperty),
                     Sd = sd(ClassProperty),
                     n=n()) %>% 
    mutate(Sd = replace_na(Sd,0)) %>%  
    mutate(Min = Mean - Sd,
           Max = Mean + Sd) %>% 
    as_tibble() %>% 
    write_csv(paste0(Path,"Longitude_",classProperty,".csv"))
} 

summary_tiff_by_LatLon <- function(tiff_file, classProperty, output_dir=NULL, degree = 1, Ext = c(-180,180,-60,90)){
  library(terra)
  library(dtplyr)
  library(dplyr)
  library(readr)
  
  preAOA <- terra::rast(tiff_file)
  df_preAOA <- terra::as.data.frame(preAOA, xy=T, na.rm = T)
  colnames(df_preAOA) <- c("Longitude","Latitude", "ClassProperty")
  write_csv(df_preAOA, file = paste0(output_dir,"preAOA_",classProperty,".csv"))
  
  # convert dataframe to dt
  dt_preAOA <- lazy_dt(df_preAOA)
  
  ### for Latitude data summarized by the degree
  breaks_Lat <- seq(Ext[3], Ext[4], degree)
  Label_breaks_Lat <- breaks_Lat[-length(breaks_Lat)]
  dt_preAOA %>%  
    mutate(Lat = cut(Latitude, breaks = breaks_Lat, labels = Label_breaks_Lat, 
                     include.lowest = TRUE) |> as.character() |> as.numeric()) %>% 
    group_by(Lat) %>% 
    dplyr::summarise(Mean = mean(ClassProperty),
                     Sd = sd(ClassProperty),
                     n =n()) %>% 
    mutate(Sd = replace_na(Sd,0)) %>% 
    mutate(Min = Mean - Sd,
           Max = Mean + Sd) %>% 
    as_tibble() %>% 
    write_csv(paste0(output_dir,"Latitude_by_", degree,"_",classProperty,".csv"))
  
  ### for Longitude data summarized by the degree
  breaks_Long <- seq(Ext[1], Ext[2], degree)
  Label_breaks_Long <- breaks_Long[-length(breaks_Long)]
  dt_preAOA %>%  
    mutate(Long = cut(Longitude, breaks = breaks_Long, labels = Label_breaks_Long, 
                      include.lowest = TRUE) |> as.character() |> as.numeric()) %>% 
    group_by(Long) %>% 
    dplyr::summarise(Mean = mean(ClassProperty),
                     Sd = sd(ClassProperty),
                     n = n()) %>% 
    mutate(Sd = replace_na(Sd,0)) %>% 
    mutate(Min = Mean - Sd,
           Max = Mean + Sd) %>% 
    as_tibble() %>% 
    write_csv(paste0(output_dir,"Longitude_by_", degree,"_",classProperty,".csv"))
  
  slice_mean <- slice_by_latitude(dt_preAOA)
  write_csv(slice_mean,paste0(output_dir,"Latitude_by_a1b1", degree,"_",classProperty,".csv"))
} 

#function for slice
slice_by_latitude <- function(Data, Ext=c(-60,90)){
  library(dplyr)
  library(dtplyr)
  library(purrr)
  library(tibble)
  library(readr)
  library(tidyr)
  ## create index
  breaks_Lat <- seq(Ext[1], Ext[2], 1)
  Range <- breaks_Lat[-length(breaks_Lat)]
  
  Slice_Index <- map(1:length(Range), ~{
    i <- .
    if (i == 1) {
      return(c(Range[i], Range[i+1]))
    } else if (i == length(Range)) {
      return(c(Range[i-1], Range[i]))
    } else {
      return(c(Range[i-1], Range[i], Range[i+1]))
    }})
  Slice_Index <- set_names(Slice_Index, Range)
  
  new_data <- Data %>% 
    mutate(Index = cut(Latitude, breaks = breaks_Lat, labels = Range, 
                       include.lowest = TRUE) |> as.character() |> as.numeric())
  
  slice_mean <- map(Slice_Index, \(x) new_data %>% filter(Index %in% x) %>% 
                      summarise(Mean = mean(ClassProperty), Sd = sd(ClassProperty), n =n()) %>% 
                      as_tibble) %>% 
    list_rbind(names_to = "Latitude") %>% 
    mutate(Sd = replace_na(Sd,0)) %>% 
    mutate(Min = Mean - Sd,
           Max = Mean + Sd)
  slice_mean
}


mess_map <- function(model_name, model, training_data, predicts_layers){
  var_names <- names(model$variable.importance)
  data <- as.matrix(training_data[,var_names])
  mess_map <- dismo::mess(predicts_layers[[var_names]], data, full = FALSE)
  writeRaster(mess_map,paste0(modelname,"/mess_ratio.tif"),overwrite=TRUE)
  rcl <- c(-10000,0,-1,0,0,0,0,10000,1)
  rcl <- matrix(rcl,ncol=3,byrow=T)
  mess_map_reclass <- raster::reclassify(mess_map,rcl)
  writeRaster(mess_map_reclass,paste0(modelname,"/mess_ratio_reclass.tif"),overwrite=TRUE)
}

tiles_aoa_f <- function(modelname,name, svs, tdi){
  tile_name <- paste0(modelname, "/results/tile_", name, ".tiff")
  tile <- terra::rast(tile_name)
  aoa_tile <- CAST::aoa(newdata = tile, model = svs, trainDI = tdi)
  terra::writeRaster(aoa_tile$AOA, paste0(modelname, "/results/AOA_tile", name ,".tiff"), overwrite = TRUE)
  terra::writeRaster(aoa_tile$DI, paste0(modelname, "/results/DI_tile", name ,".tiff"), overwrite = TRUE)
}

AOA_global_predict <- function(modelname, classProperty, variable_tiff, Ext = c(-180,180,-60,90), summarised=FALSE){
  library(terra)
  library(purrr)
  library(readr)
  library(CAST)
  library(fs)
  ## 1. read model part
  svs = readRDS(paste0(modelname,"/rfmodel.RDS"))
  tdi=readRDS(paste0(modelname,"/trainDI.RDS"))

  ## 2. process predict and save predicted tiff
  Variable_tiff <- terra::rast(variable_tiff, lyrs=svs$selectedvars)
  pre =terra::predict(Variable_tiff, svs, na.rm = T)
  terra::writeRaster(pre, paste0(modelname, "/results/prediction_", classProperty, ".tiff"), overwrite = TRUE)

  ## 3. process AOA analysis and save AOA tiff
  ### 3.0 split tiff to 16 pieces and save
  SplitRas(modelname = modelname,Variable_tiff,4)
  purrr::walk(1:16,\(x) tiles_aoa_f(name = x,modelname = modelname, svs = svs, tdi = tdi))

  ### 3.1 AOA
  file_all_aoa <- paste0(modelname,"/results/AOA_tile",1:16,".tiff")
  aoa_list <- map(file_all_aoa, \(x) rast(x))
  aoa_rsrc <- sprc(aoa_list)
  aoa_final <- mosaic(aoa_rsrc)

  ### 3.2 DI
  file_all_DI <- paste0(modelname,"/results/DI_tile",1:16,".tiff")
  di_list <- map(file_all_DI, \(x) rast(x))
  di_rsrc <- sprc(di_list)
  di_final <- mosaic(di_rsrc)

  ### 3.3 save AOA ,DI tiff and remove split AOA,ID tiff
  writeRaster(aoa_final, filename = paste0(modelname,"/results/AOA_", classProperty, ".tiff"), overwrite = TRUE)
  writeRaster(di_final, filename = paste0(modelname,"/results/DI_", classProperty, ".tiff"), overwrite = TRUE)

  fs::file_delete(paste0(modelname,"/results/AOA_tile",1:16,".tiff"))
  fs::file_delete(paste0(modelname,"/results/DI_tile",1:16,".tiff"))
  fs::file_delete(paste0(modelname,"/results/tile_",1:16,".tiff"))

  ## 4. crop the tiff by Ext, and mask the predited_tiff by AOA_tiff
  pre_croped <- terra::crop(pre, Ext)
  AOA_croped <- terra::crop(aoa_final, Ext)
  pre_AOA <- terra::mask(pre_croped,AOA_croped, maskvalues=0)

  writeRaster(pre_croped,paste0(modelname,"/results/prediction_croped_", classProperty, ".tiff"), overwrite = TRUE)
  writeRaster(AOA_croped,paste0(modelname,"/results/AOA_croped_", classProperty, ".tiff"), overwrite = TRUE)
  writeRaster(pre_AOA,paste0(modelname,"/results/prediction_AOA_croped_", classProperty, ".tiff"), overwrite = TRUE)
  
  pre_AOA <- terra::rast(paste0(modelname,"/results/prediction_AOA_croped_", classProperty, ".tiff"))
  
  ## 5. summarise tiff by Lat
  if(summarised){
    summary_tiff_by_LatLon(tiff_file = pre_AOA, classProperty = classProperty, output_dir=paste0(modelname,"/results/"), degree = 1, Ext = c(-180,180,-60,90))
  }
}

svs_ffs_predict <- function(modelname,classProperty, GEE_point_env, covariable, MinVar = 2){
  fs::dir_create(paste0(modelname,"/results"),recurse = TRUE)

  training_samples <-
    read_csv(GEE_point_env) %>%
    select(!!classProperty, `.geo`, all_of(covariable)) %>%
    mutate(locate = str_extract(`.geo`, "\\[.+\\]") |> str_sub(2,-2) , .before = 2, .keep = "unused") %>%
    separate_wider_delim(locate, delim = ",", names = c("Longitude", "Latitude"))
  
  training_samples_sf <- st_as_sf(training_samples, coords = c("Longitude", "Latitude"), crs = 4326)
  
  folds <-  create_folds_spatial(modelname, training_samples_sf, n_folds = 10, gridsize = 6, seed = 11)
  
  hyperparameter = expand.grid(mtry = 2,splitrule = "variance",min.node.size = 5)
  
  svs = spatial_variable_selection(modelname, training_samples,
                                   predictors = covariable, response = classProperty,
                                   folds = folds, hyperparameter = hyperparameter, MinVar = 2)
  
  tdi = pi_trainDI(modelname, svs)
  svs$results
  svs$selectedvars
  coef_det(svs$pred$obs,svs$pred$pred)
}

rfe_model_results <- function(modelname,classProperty, GEE_point_env, covariable, sizes =84){
  training_samples <-
    read_csv(GEE_point_env) %>%
    select(!!classProperty, `.geo`, all_of(covariable)) %>%
    mutate(locate = str_extract(`.geo`, "\\[.+\\]") |> str_sub(2,-2) , .before = 2, .keep = "unused") %>%
    separate_wider_delim(locate, delim = ",", names = c("Longitude", "Latitude"))
  x <- training_samples[, covariable] 
  y <- training_samples[,classProperty] %>% pull
  
  control <- rfeControl(functions = rfFuncs, method = "cv")
  set.seed(123)
  result <- rfe(x, y, sizes = sizes, rfeControl = control)# c(1:ncol(x))
  selected_features <- predictors(result)
  selected_features 
  formu <- paste0(classProperty, "~ .")
  model_caret <- ranger(formu , data = training_samples[,c(classProperty,selected_features)], mtry = 2, num.trees = 300, importance = "impurity")
  saveRDS(model_caret, paste0(modelname, "/model_caret.RDS"))
  
  coef_R2 <- coef_det(y, model_caret$predictions)
  
  all_results_caret <- model_caret$variable.importance %>% 
    as_tibble(rownames = "variable") %>% 
    rename(importance = value) %>% 
    mutate(r2 = model_caret$r.squared) %>% 
    mutate(coef = coef_R2)
  saveRDS(all_results_caret, paste0(modelname, "/results/all_results_caret.RDS"))
}

caret_model_tidy <- function(model_name, output_file = NULL){
  model <- readRDS(paste0(model_name, "/model_caret.RDS"))
  vars <- paste(model$selectedvars, collapse = ";")
  model_R2 <- model$results$Rsquared
  coef_R2 <- coef_det(model$pred$obs, model$pred$pred)
  Table <- tibble(name = model_name, factors = vars, model_R2 = model_R2, coef_R2 = coef_R2)
  Table
  if(!is.null(output_file)){
    readr::write_csv(Table, file = output_file)
  }
}

randomForest_model <- function(modelname,classProperty, GEE_point_env, covariable, sizes =84){
  library(rfPermute)
  training_samples <-
    read_csv(GEE_point_env) %>%
    dplyr::select(!!classProperty, all_of(covariable))
  print(head(training_samples))
  set.seed(123)
  formu <- as.formula(paste0(classProperty, "~ .")) 
  model_rfPermute <- rfPermute(formu , data = training_samples, ntrees = 300,num.cores =4)
  saveRDS(model_rfPermute, paste0(modelname, "/model_rfPermute.RDS"))
}

#####################
glmmHP_tidy <- function(data,property,output){
  envs <- c("CHELSA_BIO_Annual_Mean_Temperature","CHELSA_BIO_Temperature_Annual_Range","CHELSA_BIO_Temperature_Seasonality",
            "CHELSA_BIO_Annual_Precipitation","CHELSA_BIO_GrowingSeasonPrecipitationSum","CHELSA_BIO_Precipitation_Seasonality",
            "CHELSA_BIO_AridityIndex","CHELSA_BIO_MeanMonthlyMoisture_Index","CHELSA_BIO_MeanMonthlyPotential_Evapotranspiration","CHELSA_BIO_MeanMonthlyVaporPressureDeficit",
            "EarthEnv_Texture_Shannon_Index","EarthEnv_Texture_Simpson_Index","EarthEnv_Texture_Variance_EVI","CHELSA_BIO_NPP",
            "SG_Bulk_density_0_5cm","SG_CEC_0_5cm","SG_Nitrogen_Content_0_5cm","SG_SOC_Content_0_5cm","SG_Soil_pH_0_5cm")

  envs_tibble <- tibble(Varibles = envs,
                        Type = rep(c("temperature","precipitation","aridityindex","plant","soil"), c(3,3,4,4,5)))

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

glmmHP_tidy2 <- function(data,property,output){
  
  envs <- c("CHELSA_BIO_Annual_Mean_Temperature","CHELSA_BIO_Temperature_Annual_Range","CHELSA_BIO_Temperature_Seasonality",
            "CHELSA_BIO_Annual_Precipitation","CHELSA_BIO_GrowingSeasonPrecipitationSum","CHELSA_BIO_Precipitation_Seasonality",
            "CHELSA_BIO_AridityIndex","CHELSA_BIO_MeanMonthlyMoisture_Index","CHELSA_BIO_MeanMonthlyPotential_Evapotranspiration","CHELSA_BIO_MeanMonthlyVaporPressureDeficit",
            "EarthEnv_Texture_Shannon_Index","EarthEnv_Texture_Simpson_Index","EarthEnv_Texture_Variance_EVI","CHELSA_BIO_NPP",
            "SG_Bulk_density_0_5cm","SG_CEC_0_5cm","SG_Nitrogen_Content_0_5cm","SG_SOC_Content_0_5cm","SG_Soil_pH_0_5cm")
  
  envs_tibble <- tibble(Varibles = envs,
                        Type = rep(c("temperature","precipitation","aridityindex","plant","soil"), c(3,3,4,4,5)))
  Data <<-  read_csv(data) %>% 
    select(-`.geo`, - `system:index`) %>% 
    select(property, all_of(envs), Paper) %>% 
    drop_na() %>% 
    mutate(across(all_of(envs), ~ scale(.x)[,1]))
  
  formula_raw <<- as.formula(
    paste(property,"~",
          paste(envs,collapse = "+")
    ))
  print(formula_raw)
  lm_raw_results <- lm(formula = formula_raw,data = Data)
  property_factor <- effectsize(lm_raw_results)[-1,] %>% as_tibble()
  
  lm_gmmhp <<- lm(formula = formula_raw,data = Data)
  lm_gmmhp_results <- glmm.hp(lm_gmmhp)
  a <- lm_gmmhp_results$hierarchical.partitioning[,4] 
  A <- lm_gmmhp_results$r.squaredGLMM
  
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