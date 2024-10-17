source("script/process_script/Rfunction.R")
library(tidyverse)
library(conflicted)
library(fs)
library(rlang)
library(patchwork)
library(vegan)
conflict_scout()
conflict_prefer("filter", "dplyr")
conflict_prefer("summarise", "dplyr")
conflict_prefer("map", "purrr")
conflict_prefer("rename", "dplyr")

############################ 1.0 generate table (unique_site_richness_without_tundra.csv)
site_table <-
  readxl::read_xlsx("asv_process/env/metadata_wetland_fungi.xlsx") %>%
  select(1,3,7,8,Wetland_type,Wetland_type_b,vegetation_Type,Plant_Species,Paper,depth) %>%
  dplyr::rename(sample = SRR_Run_num) %>%
  mutate(across(c(3, 4), ~ round(as.numeric(.x), 2))) %>%
  mutate(lat_lon = paste0(lon, "_", lat))  %>%
  mutate(across(c(Wetland_type, vegetation_Type), str_to_title)) %>%
  filter(Wetland_type != "Tundra")

####### get the sample_name for analysis
sample_name <- "asv_process/rawdata/" %>%
  dir_ls %>%
  path_file %>%
  path_filter(., glob = "unoise_otu_tab*") %>%
  str_remove("unoise_otu_tab_") %>%
  str_remove(".txt") %>%
  set_names(., .)

## 01 Merge the asv and tax tables of fungi to generate the asv_tax table, and remove ASVs that are not mapped to the fungal kingdom
dir_create("asv_process/tmp/01.ASV_Tax_tab")
walk(
  sample_name,
  ~ combineTaxAsvtable(
    paste0("asv_process/rawdata/unoise_otu_tab_", .x, ".txt"),
    paste0("asv_process/rawdata/zotus_classified_", .x , ".txt"),
    table_input = TRUE,
    is_fungi = TRUE,
    output_file = paste0("asv_process/tmp/01.ASV_Tax_tab/asv_tax_", .x, ".txt")
  )
)

### 1.2 Dilute the ASV_Tax table to generate a rarified table for 2000 sequences 
dir_create("asv_process/tmp/02.ASV_Tax_RR")

sample_filter_list <-map(sample_name, 
                         ~ filter_sample_By_readsNum(paste0("asv_process/tmp/01.ASV_Tax_tab/asv_tax_", .x, ".txt"),
                                                     MinNums = 2000,table_input = TRUE))

asv_tab_list <-imap(sample_filter_list, ~ raryfied_asv(.x, rarify_num = 2000,
         output_file = paste0("asv_process/tmp/02.ASV_Tax_RR/ASV_Tax_RR_", .y, ".txt")
       ))

## 1.3 combined all otu table into one table (unique_site_richness_without_tundra.csv)
### Global richness
global_fungi <- map(paste0("asv_process/tmp/02.ASV_Tax_RR/ASV_Tax_RR_", sample_name, ".txt"), read_tsv) %>% 
  map(., richness_div) %>% 
  bind_rows() %>% 
  left_join(site_table) %>% 
  filter(Wetland_type != "Tundra") %>% 
  write_csv(file = "asv_process/results/Global_fungi_richness_without_tundra.csv") %>% 
  mutate(fungi = "Fungi")

### Richness (unique site)
unique_site_richness <- global_fungi %>%
  group_by(lat_lon, Wetland_type) %>%
  summarise(richness = mean(richness),N = n()) %>%
  mutate(richness = round(richness,0)) %>%
  ungroup() %>%
  separate_wider_delim(lat_lon, delim = "_", names = c("Lon","Lat"), cols_remove = FALSE) %>%
  write_csv("asv_process/results/unique_site_richness_without_tundra.csv")

### Shannon index  (unique site)
unique_site_Shannon <- global_fungi %>%
  group_by(lat_lon, Wetland_type) %>%
  summarise(Shannon = mean(Shannon),N = n()) %>%
  mutate(Shannon = round(Shannon,3)) %>%
  ungroup() %>%
  separate_wider_delim(lat_lon, delim = "_", names = c("Lon","Lat"), cols_remove = FALSE) %>%
  write_csv("asv_process/results/unique_site_Shannon_without_tundra.csv")

#########################################################################################################
#2.0 Fungaltrait  ###################
dir_create("asv_process/tmp/03.ASV_fungaltrait")
fungal_trait <- read_tsv("asv_process/env/fungal_trait_v2.txt") %>% 
  mutate(ID = as.character(ID))

sample_name <- "asv_process/tmp/02.ASV_Tax_RR/" %>% dir_ls %>% path_file %>% path_filter(., glob = "ASV_Tax_RR_*") %>% 
  str_remove("ASV_Tax_RR_") %>% 
  str_remove(".txt") %>% 
  set_names(., .)

# combined fungal trait accroding to genus 
tax_list <- map(sample_name,\(x) read_tsv(paste0("asv_process/tmp/02.ASV_Tax_RR/ASV_Tax_RR_",x, ".txt")) %>% 
                  select(ASV_ID,taxonomy) %>% 
                  separate_wider_delim(taxonomy,delim = ";",names = c("kingdom","phylum","class","order","family","genus","species")) %>% 
                  select(ASV_ID,genus) %>% 
                  left_join(fungal_trait, by = "genus"))

# generate fungal trait table (numbers rather than relative abundance)
ASV_table <- map(sample_name,
                 \(x) read_tsv(paste0("asv_process/tmp/02.ASV_Tax_RR/ASV_Tax_RR_",x, ".txt")))
ASV_tax_num_fungaltrait <-
  map2(ASV_table,tax_list,\(x, y) x %>% 
         left_join(y, by = "ASV_ID") %>%mutate(across(c("trait", "Aquatic"), ~ replace_na(.x, "Unknown")))) %>%
  iwalk(\(x, idx) write_tsv(x,paste0("asv_process/tmp/03.ASV_fungaltrait/ASV_tax_num_",idx,"_fungaltrait.txt")))

# generate the relative abundance of fungal trait 
ASV_relative <- map(sample_name,\(x) read_tsv(paste0("asv_process/tmp/02.ASV_Tax_RR/ASV_Tax_RR_",x, ".txt")) %>% 
                      mutate(across(where(is.numeric), \(x) x/sum(x)))) 

ASV_tax_relative_fungaltrait <-map2(ASV_relative,tax_list, \(x, y) x %>% 
                                      left_join(y, by = "ASV_ID") %>%
                                      mutate(across(c("trait", "Aquatic"), ~ replace_na(.x, "Unknown")))) %>%
  iwalk(\(x, idx) write_csv(x,paste0("asv_process/tmp/03.ASV_fungaltrait/ASV_tax_relative_",idx,"_fungaltrait.txt")))

### combined all fungal trait table and generate richness or Shannon index for Symbiotrophs,Saprotrophs and Pathotrophs
tmp_metadata <- read_csv("asv_process/results/Global_fungi_richness_without_tundra.csv") %>% 
  select(sample,Wetland_type) 

fungal_trait_list <- map(paste0("asv_process/tmp/03.ASV_fungaltrait/ASV_tax_num_", sample_name, "_fungaltrait.txt"), read_tsv) %>%
  map(., ~ .x %>%
        select(-taxonomy,-ID,-Aquatic,-genus) %>% 
        group_nest(trait) %>%
        mutate(data = map(data, richness_div)) %>%
        unnest(data)) %>%
  bind_rows()%>%
  left_join(tmp_metadata) %>% 
  filter(trait %in% c("Symbiotrophs", "Saprotrophs","Pathotrophs")) %>% 
  filter(!is.na(Wetland_type)) %>% 
  mutate(fungi = "Fungi") %>% 
  write_tsv("asv_process/results/ASV_fungaltrait_richness_for_ggplot.txt")

#### Pathotrophs_richness
unique_Pathotrophs_richness <- read_tsv("asv_process/results/ASV_fungaltrait_richness_for_ggplot.txt") %>% 
  filter(trait == "Pathotrophs") %>% 
  left_join(site_table) %>% 
  summarise(richness = mean(richness), N = n(), .by = "lat_lon") %>% 
  mutate(richness = round(richness,0)) %>% 
  filter(!richness == 0) %>% 
  ungroup() %>%
  separate_wider_delim(lat_lon, delim = "_", names = c("Lon","Lat"), cols_remove = FALSE) %>% 
  write_csv("asv_process/results/unique_Pathotrophs_richness.csv")

#### Symbiotrophs_richness
unique_Symbiotrophs_richness <- read_tsv("asv_process/results/ASV_fungaltrait_richness_for_ggplot.txt") %>% 
  filter(trait == "Symbiotrophs") %>% 
  left_join(site_table) %>% 
  summarise(richness = mean(richness), N = n(), .by = "lat_lon") %>% 
  mutate(richness = round(richness,0)) %>% 
  filter(!richness == 0) %>% 
  ungroup() %>%
  separate_wider_delim(lat_lon, delim = "_", names = c("Lon","Lat"), cols_remove = FALSE) %>% 
  write_csv("asv_process/results/unique_Symbiotrophs_richness.csv")

#### Saprotrophs_richness
unique_Saprotrophs_richness <- read_tsv("asv_process/results/ASV_fungaltrait_richness_for_ggplot.txt") %>% 
  filter(trait == "Saprotrophs") %>% 
  left_join(site_table) %>% 
  summarise(richness = mean(richness), N = n(), .by = "lat_lon") %>% 
  mutate(richness = round(richness,0)) %>% 
  filter(!richness == 0) %>% 
  ungroup() %>%
  separate_wider_delim(lat_lon, delim = "_", names = c("Lon","Lat"), cols_remove = FALSE) %>% 
  write_csv("asv_process/results/unique_Saprotrophs_richness.csv")

### combined all fungal trait table and generate the relative abundance for all trait group
#### (Symbiotrophs,Saprotrophs and Pathotrophs) or ("aquatic","partly-aquatic","non-aquatic")
tmp_metadata <- read_csv("asv_process/results/Global_fungi_richness_without_tundra.csv") %>% 
  select(sample,Wetland_type) 

##### the relative abundance 
ASV_fungaltrait <- ASV_tax_relative_fungaltrait %>% 
  map(\(x) x %>% summarise(across(where(is.numeric), sum), .by = trait)) %>% 
  reduce(\(x,y) full_join(x,y,by = "trait")) %>% 
  pivot_rowtocol("trait","sample") %>% 
  mutate(across(where(is.numeric), \(x) replace_na(x,replace = 0))) %>% 
  right_join(tmp_metadata,by = "sample") %>%
  select(-Wetland_type) %>% 
  write_tsv("asv_process/results/ASV_fungaltrait.txt")

ASV_fungaltrait_global <- ASV_fungaltrait %>% 
  select(-sample) %>% 
  summarise(across(where(is.numeric), \(x) mean(x)*100)) %>% 
  mutate(Wetland_type = "Global")

ASV_fungaltrait_bytype <- ASV_fungaltrait %>% left_join(tmp_metadata,by = "sample") %>% 
  select(-sample) %>% 
  summarise(across(where(is.numeric), \(x) mean(x) * 100), .by = "Wetland_type") %>% 
  bind_rows(ASV_fungaltrait_global) %>% 
  pivot_rowtocol("Wetland_type","trait") %>% 
  mutate(trait = factor(trait,levels = c("Symbiotrophs","Saprotrophs","Pathotrophs","Other","unspecified","Unknown"))) %>% 
  pivot_longer(-trait, names_to = "Group") %>% 
  mutate(Group = factor(Group, levels = c("Coastal","Inland","Peatland", "Global"))) %>% 
  write_tsv("asv_process/results/ASV_fungaltrait_for_ggplot.txt")

##### Aquatic type for all fungal ASV
ASV_fungalAquatic <- ASV_tax_relative_fungaltrait %>% 
  map(\(x) x %>% summarise(across(where(is.numeric), sum), .by = "Aquatic")) %>% 
  reduce(\(x,y) full_join(x,y,by = "Aquatic")) %>% 
  pivot_rowtocol("Aquatic","sample") %>% 
  mutate(across(where(is.numeric), \(x) replace_na(x,replace = 0))) %>%
  #mutate(Unknown = Unknown + undefined, .keep = "unused") %>% 
  right_join(tmp_metadata,by = "sample") %>% 
  select(-Wetland_type) %>% 
  write_tsv("asv_process/results/ASV_fungalAquatic.txt")

ASV_fungalAquatic_global <- ASV_fungalAquatic %>% 
  select(-sample) %>% 
  summarise(across(where(is.numeric), \(x) mean(x)*100)) %>% 
  mutate(Wetland_type = "Global")

##### Aquatic type for fungal ASV by wetland type
ASV_fungalAquatic_bytype <- ASV_fungalAquatic %>% left_join(tmp_metadata,by = "sample") %>% 
  select(-sample) %>% 
  summarise(across(where(is.numeric), \(x) mean(x) * 100), .by = "Wetland_type") %>% 
  bind_rows(ASV_fungalAquatic_global) %>% 
  pivot_rowtocol("Wetland_type","Aquatic") %>% 
  mutate(Aquatic = factor(Aquatic,levels = c("aquatic","partly-aquatic","non-aquatic","unspecified","Unknown"))) %>% 
  arrange(Aquatic) %>% 
  write_tsv("asv_process/results/ASV_fungalAquatic_for_ggplot.txt")


##### add others
ASV_fungaltrait_num_for_otherTrait <- map(sample_name,\(x) read_tsv(
    paste0("asv_process/tmp/03.ASV_fungaltrait/ASV_tax_num_",x,"_fungaltrait.txt")) %>%
    filter(trait == "Other") %>%
    select(-ID) %>%
    summarise(across(where(is.numeric), sum), .by = primary_lifestyle)) %>%
  reduce(\(x, y) full_join(x, y, by = "primary_lifestyle")) %>% 
  mutate(across(where(is.numeric), \(x) replace_na(x, replace = 0))) %>% 
  select("primary_lifestyle", where( ~ is.numeric(.x) && sum(.x) > 0)) %>% 
  mutate(across(where(is.numeric), \(x) x/sum(x))) %>%
  pivot_rowtocol("primary_lifestyle", "sample") %>%
  right_join(tmp_metadata, by = "sample")

ASV_fungaltrait_num_for_otherTrait_global <- ASV_fungaltrait_num_for_otherTrait %>% 
  select(-sample,-Wetland_type) %>% 
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)*100)) %>% 
  mutate(Wetland_type = "Global")

ASV_fungaltrait_num_for_otherTrait_bytype <- ASV_fungaltrait_num_for_otherTrait %>% 
  select(-sample) %>% 
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)*100), .by = "Wetland_type") %>% 
  bind_rows(ASV_fungaltrait_num_for_otherTrait_global) %>% 
  pivot_rowtocol("Wetland_type","primary_lifestyle") %>% 
  arrange(desc(Global)) %>% 
  write_tsv("asv_process/results/ASV_fungalTrait_others_for_ggplot.txt")

asv_tab_list <- map(sample_name,\(x) read_tsv(paste0("asv_process/tmp/02.ASV_Tax_RR/ASV_Tax_RR_",x, ".txt")))

###### phylum_fungi_top12
phylum_fungi_top12 <- imap(asv_tab_list, ~ split_taxonomy_by_level(.x, level = "phylum" ,is_fungi = TRUE)) %>% 
  reduce(., full_join) %>% 
  mutate(across(-phylum, ~ replace_na(.x, 0))) %>% 
  ggtable_composition_fungi_phylum(., 12) 

phylum_fungi_top12_global <-  phylum_fungi_top12 %>% 
  pivot_longer(-phylum, names_to = "sample") %>%
  pivot_wider(names_from = phylum, values_from = value) %>% 
  left_join(site_table) %>% 
  select(-lat,-lon) %>% 
  filter(!is.na(Wetland_type)) %>% 
  summarise(across(where(is.numeric), ~ mean(.))) %>% 
  mutate(sample = "Global", .before = 1) %>% 
  pivot_longer(-sample, names_to = "Phylum") %>%
  pivot_wider(names_from = sample, values_from = value)

phylum_fungi_top12_byGroup <- phylum_fungi_top12 %>% 
  pivot_longer(-phylum, names_to = "sample") %>%
  pivot_wider(names_from = phylum, values_from = value) %>% 
  left_join(site_table) %>% 
  select(-lat,-lon) %>% 
  filter(!is.na(Wetland_type)) %>% 
  group_by(Wetland_type) %>%
  summarise(across(where(is.numeric), ~ mean(.))) %>% 
  pivot_longer(-Wetland_type, names_to = "Phylum") %>%
  pivot_wider(names_from = Wetland_type, values_from = value) %>% 
  left_join(phylum_fungi_top12_global) %>% 
  mutate(Phylum = factor(Phylum, levels = Phylum))

write_csv(phylum_fungi_top12_global, file = "asv_process/results/phylum_fungi_top12_global.csv")
write_csv(phylum_fungi_top12_byGroup, file = "asv_process/results/phylum_fungi_top12_byGroup.csv")

###### remove tmp file
dir_delete("asv_process/tmp/")



