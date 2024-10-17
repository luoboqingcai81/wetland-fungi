
library(argparser)
p <- arg_parser("summarise latitude mean value by degrees")
p <- add_argument(p, "--preTable",
                  help = "predicted csv")

p <- add_argument(p, "--property",
                  help = "the predicted classProperty, such as richness, biomass..")

p <- add_argument(p, "--output",
                  help = "the output dictionary")

argv <- parse_args(p)

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


summary_tiff_by_Lat_slice <- function(preTable, classProperty, output_dir=NULL, degree = 1, Ext = c(-180,180,-60,90)){
  library(terra)
  library(dtplyr)
  library(dplyr)
  library(readr)
  
  df_preAOA <- read_csv(preTable)
  colnames(df_preAOA )
  # convert dataframe to dt
  dt_preAOA <- lazy_dt(df_preAOA)
  
  slice_mean <- slice_by_latitude(dt_preAOA)
  write_csv(slice_mean,paste0(output_dir,"Latitude_by_a1b1_", degree,"_",classProperty,".csv"))
} 

summary_tiff_by_Lat_slice(preTable = argv$preTable, classProperty = argv$property)

