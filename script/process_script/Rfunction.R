

####################   combined the tax table to asv table and generate a new table ######################
combineTaxAsvtable <- function(asv_table, rdp_table,table_input = FALSE, is_fungi =  FALSE, output_file = NULL){
  if(!is_fungi){
    if(table_input){
      ASV <- read_tsv(asv_table) %>% 
        rename(ASV_ID = `#OTU ID`)
      rdp_tax <- read_tsv(rdp_table) %>% 
        rename(ASV_ID = `OTU_ID`)
    }else{
      ASV <- asv_table %>% 
        rename(ASV_ID = `#OTU ID`)
      rdp_tax <- rdp_table %>% 
        rename(ASV_ID = `OTU_ID`)
    }
    
    colnames(rdp_tax) <- c("ASV_ID","kindom","phylum","class","order","family","genus")
    
    combined_ASV <- ASV %>% 
      left_join(rdp_tax) %>% 
      filter(class != "Chloroplast") %>% 
      filter(phylum != "Cyanobacteria/Chloroplast") %>% 
      filter(kindom == "Bacteria") %>%
      unite("taxonomy", c("kindom","phylum","class","order","family","genus"),sep = ";", remove = TRUE) %>%
      select(ASV_ID, where(is.numeric), taxonomy) 
    
    if(!is.null(output_file)){
      write_tsv(combined_ASV, file = output_file)
    }
  }else{
    if(table_input){
      ASV <- read_tsv(asv_table) %>% 
        rename(ASV_ID = `#OTU ID`)
      rdp_tax <- read_tsv(rdp_table, col_names = FALSE, show_col_types = FALSE)%>% 
        select(1,6,9,12,15,18,21,24) %>% 
        mutate(X24 = str_split_fixed(X24, "\\|", n=2)[,1]) 
    }else{
      ASV <- asv_table %>% 
        rename(ASV_ID = `#OTU ID`)
      rdp_tax <- rdp_table %>% 
        select(1,6,9,12,15,18,21,24) %>% 
        mutate(X24 = str_split_fixed(X24, "\\|", n=2)[,1]) 
    }
    
    colnames(rdp_tax) <- c("ASV_ID","kindom","phylum","class","order","family","genus","species")
    
    combined_ASV <- ASV %>% 
      left_join(rdp_tax) %>% 
      filter(!kindom == "NA") %>% 
      unite("taxonomy", c("kindom","phylum","class","order","family","genus","species"),sep = ";", remove = TRUE) %>%
      select(ASV_ID, where(is.numeric), taxonomy) 
    
    if(!is.null(output_file)){
      write_tsv(combined_ASV, file = output_file)
    }
    combined_ASV
  }
}

combineFunguildTax_asv <- function(asv_tax_rr, Funguild_tax ,table_input = FALSE, output_file = NULL){
  if(table_input){
    Asv_tax_rr <- read_tsv(asv_tax_rr)
    Funguild_tax <- read_tsv(Funguild_tax) %>% 
      rename(ASV_ID = OTU) %>% 
      select(ASV_ID, trophicMode, guild)
  }else{
    Asv_tax_rr <- asv_tax_rr
    Funguild_tax <- Funguild_tax
  }
  
  asv_tax_funguild <- Asv_tax_rr %>% 
    left_join(Funguild_tax)
  
  if(!is.null(output_file)){
    write_tsv(asv_tax_funguild, file = output_file)
  }
  asv_tax_funguild
}


####################   filter_ASV_ByTaxname
filter_ASV_ByTaxname <- function(ASV_file, tax_name, level, output_file){
  tax_name <- sym(tax_name)
  level <- sym(level)
  ASV <- read_tsv(ASV_file) %>% 
    tidyr::separate(taxonomy, sep = ";",into = c("kindom","phylum","class","order","family","genus","species")) %>%
    dplyr::filter(!(!!level == rlang::as_string(tax_name))) %>% 
    unite("taxonomy", c("kindom","phylum","class","order","family","genus","species"),sep = ";", remove = TRUE) %>%
    write_tsv(output_file)
}

####################   filter_sample_By_readsNum
filter_sample_By_readsNum <- function(asv_file, MinNums, table_input = FALSE, output_file = NULL){
  if(table_input){
    qualited_sample <- read_tsv(asv_file) %>% 
      select(-taxonomy) %>% 
      column_to_rownames(var = "ASV_ID") %>% 
      summarise(across(everything(), ~ sum(.x))) %>% 
      pivot_longer(everything(), names_to = "sample", values_to = "readsNum") %>% 
      filter(readsNum > MinNums) %>% 
      pull(sample)
    
    asv_filter <- read_tsv(asv_file) %>% 
      select(ASV_ID, all_of(qualited_sample), taxonomy) 
  }else{
    qualited_sample <- asv_file %>% 
      select(-taxonomy) %>% 
      column_to_rownames(var = "ASV_ID") %>% 
      summarise(across(everything(), ~ sum(.x))) %>% 
      pivot_longer(everything(), names_to = "sample", values_to = "readsNum") %>% 
      filter(readsNum > MinNums) %>% 
      pull(sample)
    
    asv_filter <- asv_file %>% 
      select(ASV_ID, all_of(qualited_sample), taxonomy)
  }

  if(!is.null(output_file)){
    write_tsv(asv_filter, file = output_file)
  }
  asv_filter
}

#################### raryfied_as
raryfied_asv <- function(asv_table, table_input = FALSE, rarify_num = NULL, output_file = NULL ){
  if(table_input){
    asv_file <- read_tsv(asv_table)
  }else{
    asv_file <- asv_table
  }
  
  tax_table <- asv_file %>% 
    select(ASV_ID, taxonomy)
  
  asv_table <- asv_file %>% 
    select(-taxonomy) %>%
    column_to_rownames(var = "ASV_ID")
  
  if(is.null(rarify_num)){
    Min_reads_Num <- min(colSums(asv_table))
    print(paste0("Min_rarified_number :" ,Min_reads_Num))
    
    asv_rarefied <- as.data.frame(t(vegan::rrarefy(t(asv_table), Min_reads_Num))) %>%
      filter(if_any(everything(), ~ .x != 0)) %>% 
      rownames_to_column(var = "ASV_ID") %>% 
      left_join(tax_table)
  }else{
    asv_rarefied <- as.data.frame(t(vegan::rrarefy(t(asv_table), rarify_num))) %>%
      filter(if_any(everything(), ~ .x != 0)) %>% 
      rownames_to_column(var = "ASV_ID") %>% 
      left_join(tax_table)
  }
  if(!is.null(output_file)){
    write_tsv(asv_rarefied, file = output_file)
  }
  asv_rarefied
}

####################
split_taxonomy_by_level <- function(asv_table, level = NULL, relative = FALSE, is_fungi =  FALSE, output_file = NULL){
  if(!is_fungi){
    if(is.null(level)){
      ASV_table <- asv_table %>%
        mutate(across(where(is.numeric), ~ .x / sum(.x, na.rm = TRUE)))
    }else{
      level <- sym(level)
      if(relative){
        ASV_table <- asv_table %>%
          mutate(across(where(is.numeric), ~ .x / sum(.x, na.rm = TRUE))) %>% 
          tidyr::separate(taxonomy, sep = ";",
                          into = c("kindom","phylum","class","order","family","genus")) %>%
          dplyr::select(!!level, where(is.numeric))%>%
          dplyr::group_by(!!level) %>%
          dplyr::summarise(dplyr::across(where(is.numeric),~ sum(.x, na.rm = TRUE)))
      }else{
        ASV_table <- asv_table %>%
          tidyr::separate(taxonomy, sep = ";",
                          into = c("kindom","phylum","class","order","family","genus")) %>%
          dplyr::select(!!level, where(is.numeric))%>%
          dplyr::group_by(!!level) %>%
          dplyr::summarise(dplyr::across(where(is.numeric),~ sum(.x, na.rm = TRUE)))
      }
    }
    
    if(!is.null(output_file)){
      write_tsv(ASV_table, file = output_file)
    }
    return(ASV_table)
  }else{
    if(is.null(level)){
      ASV_table <- asv_table %>%
        mutate(across(where(is.numeric), ~ .x / sum(.x, na.rm = TRUE)))
    }else{
      level <- sym(level)
      if(relative){
        ASV_table <- asv_table %>%
          mutate(across(where(is.numeric), ~ .x / sum(.x, na.rm = TRUE))) %>% 
          tidyr::separate(taxonomy, sep = ";",
                          into = c("kindom","phylum","class","order","family","genus","species")) %>%
          dplyr::select(!!level, where(is.numeric))%>%
          dplyr::group_by(!!level) %>%
          dplyr::summarise(dplyr::across(where(is.numeric),~ sum(.x, na.rm = TRUE)))
      }else{
        ASV_table <- asv_table %>%
          tidyr::separate(taxonomy, sep = ";",
                          into = c("kindom","phylum","class","order","family","genus","species")) %>%
          dplyr::select(!!level, where(is.numeric))%>%
          dplyr::group_by(!!level) %>%
          dplyr::summarise(dplyr::across(where(is.numeric),~ sum(.x, na.rm = TRUE)))
      }
    }
    
    if(!is.null(output_file)){
      write_tsv(ASV_table, file = output_file)
    }
    return(ASV_table)
  }
}

###################
reform_fungal_guild_table <- function(asv_tax_rr_funguild, input_table = FALSE, output_file = NULL){
  if(input_table){
    Asv_tax_rr_funguild <- read_tsv(asv_tax_rr_funguild)
  }else{
    Asv_tax_rr_funguild <- asv_tax_rr_funguild
  }
  Fungal_guild_reformed <- Asv_tax_rr_funguild %>% 
    filter(trophicMode != "na") %>% 
    mutate(Counts = str_count(trophicMode, "-") + 1, .before = guild) %>%
    mutate(across(where(is.numeric), ~ .x/Counts)) %>%
    separate_longer_delim(trophicMode,delim = "-") %>% 
    select(-Counts, -guild, -taxonomy)
  
  if(!is.null(output_file)){
    write_tsv(Fungal_guild_reformed, file = output_file)
  }
  Fungal_guild_reformed
}

######################
fungal_guild_abundance <- function(Fungal_guild_reformed, sample_reads,input_table = FALSE, output_file = NULL){
  if(input_table){
    fungal_guild_reformed <- read_tsv(Fungal_guild_reformed)
  }else{
    fungal_guild_reformed <- Fungal_guild_reformed
  }
  Fungal_guild_abundance <- fungal_guild_reformed %>% 
    mutate(across(where(is.numeric), ~ .x/sample_reads)) %>% 
    group_by(trophicMode) %>% 
    summarise(across(where(is.numeric), ~ sum(.x)))  %>% 
    pivot_longer(-trophicMode, names_to = "sample", values_to = "relative abundance") 
  
  if(!is.null(output_file)){
    write_tsv(Fungal_guild_abundance, file = output_file)
  }
  Fungal_guild_abundance
}

########################
alpha_diversity <- function(data_table, output_file = NULL) {
  new_table <- t(data_table)
  richness <- specnumber(new_table)
  Chao1 <- estimateR(new_table)[2, ]
  ACE <- estimateR(new_table)[4, ]
  Shannon <- diversity(new_table, index = 'shannon')
  Simpson <- diversity(new_table, index = 'simpson')
  Inverse_Simpson <- diversity(new_table, index = "inv")
  
  alpha_result <- data.frame(richness, ACE, Chao1, Shannon, Simpson, Inverse_Simpson) %>%
    rownames_to_column(var = "sample")
  
  if(!is.null(output_file)){
    write_tsv(alpha_result, file = output_file)
  }
  alpha_result 
}


richness_div <- function(data_table, output_file = NULL) {
  Data_table <- data_table %>% column_to_rownames(var = "ASV_ID") %>% 
    select(where(is.numeric)) 
  new_table <- t(Data_table)
  Shannon <- diversity(new_table, index = 'shannon')
  Simpson <- diversity(new_table, index = 'simpson')
  Inverse_Simpson <- diversity(new_table, index = "inv")
  
  alpha_result <- data.frame(Shannon, Simpson, Inverse_Simpson) %>%
    rownames_to_column(var = "sample")
  
  richness_result <- data_table %>% 
    summarise(across(where(is.numeric), ~ sum(.x > 0))) %>% 
    pivot_longer(everything(), names_to = "sample", values_to = "richness") %>% 
    left_join(alpha_result)
  
  if(!is.null(output_file)){
    write_tsv(richness_result, file = output_file)
  }
  richness_result 
}

ggtable_composition_fungi_phylum <- function(Combined_tax, Num, input_table = FALSE, output_file = NULL){
  if(input_table){
    combined_tax <- read_tsv(Combined_tax)
  }else{
    combined_tax <- Combined_tax
  }
  
  tax_name <- rlang::sym(colnames(combined_tax)[1])
  Tax_name <- colnames(combined_tax)[1]
  Sum <- NULL
  ggcomposition_tmp <- combined_tax %>%
    dplyr::mutate(Sum = purrr::pmap_dbl(.[,-1], ~ sum(c(...))))%>%
    dplyr::mutate(across(where(is.numeric), ~ .x/sum(.x)))
  
  top10_name <- ggcomposition_tmp %>% 
    filter(!phylum == "Fungi_phy_Incertae_sedis") %>% 
    dplyr::arrange(dplyr::desc(Sum)) %>% 
    dplyr::slice(1:10) %>%
    dplyr::pull(phylum)
  
  Col_other <- ggcomposition_tmp %>% 
    filter(!(phylum %in% c("Fungi_phy_Incertae_sedis",top10_name))) %>% 
    summarise(across(where(is.numeric), ~ sum(.x))) %>% 
    mutate(phylum = "Others", .before = 1)
  
  ggcomposition <- ggcomposition_tmp %>% 
    filter(phylum %in% c("Fungi_phy_Incertae_sedis",top10_name)) %>% 
    bind_rows(Col_other) %>% 
    mutate(across(where(is.numeric), ~ .x * 100)) %>% 
    mutate(phylum = factor(phylum, levels = c(top10_name,"Others","Fungi_phy_Incertae_sedis"))) %>% 
    arrange(phylum)
  
  if(!is.null(output_file)){
    write_tsv(ggcomposition, file = output_file)
  }
  
  return(ggcomposition)
}

pivot_rowtocol <- function(data, col, names_tofrom ="name"){
  data %>% pivot_longer(cols = -col, names_to = names_tofrom) %>% 
    pivot_wider(names_from = col)
}

#################################
ggtable_composition_num <- function(combined_tax, Num){
  tax_name <- rlang::sym(colnames(combined_tax)[1])
  Sum <- NULL
  ggcomposition <- combined_tax %>%
    dplyr::mutate(Sum = purrr::pmap_dbl(.[,-1], ~ sum(c(...)))) %>%
    dplyr::filter(!!tax_name != "other") %>%
    dplyr::arrange(dplyr::desc(Sum)) %>%
    dplyr::slice(1:Num) %>%
    dplyr::select(-Sum)
  Col_other <- combined_tax %>% 
    dplyr::mutate(Sum = purrr::pmap_dbl(.[,-1], ~ sum(c(...)))) %>%
    dplyr::filter(!(!!tax_name == "other")) %>%
    dplyr::arrange(dplyr::desc(Sum)) %>%
    dplyr::slice(Num+1: -1) %>% 
    summarise(across(where(is.numeric), sum)) %>% 
    mutate(!!tax_name := "others", .before = 1) %>% 
    dplyr::select(-Sum)
  
  ggcomposition2 <- ggcomposition %>%
    dplyr::rows_insert(Col_other)
  return(ggcomposition2)
}

##############
ggtable_composition <- function(Combined_tax, Num, input_table = FALSE, output_file = NULL){
  if(input_table){
    combined_tax <- read_tsv(Combined_tax)
  }else{
    combined_tax <- Combined_tax
  }
  
  tax_name <- rlang::sym(colnames(combined_tax)[1])
  Sum <- NULL
  ggcomposition <- combined_tax %>%
    dplyr::mutate(Sum = purrr::pmap(.[,-1], ~ sum(c(...)))) %>%
    dplyr::mutate(Sum = as.numeric(Sum)) %>%
    dplyr::filter(!(!!tax_name == "other")) %>%
    dplyr::mutate(across(where(is.numeric), ~ .x/sum(.x))) %>% 
    dplyr::arrange(dplyr::desc(Sum)) %>%
    dplyr::slice(1:Num) %>%
    dplyr::select(-Sum)
  Col_other <- dplyr::tibble(!!tax_name := "others") %>%
    dplyr::bind_cols(dplyr::summarise(ggcomposition,dplyr::across(where(is.numeric),~ 1-sum(.x))))
  ggcomposition2 <- ggcomposition %>%
    dplyr::rows_insert(Col_other) %>% 
    mutate(across(where(is.numeric), ~ .x * 100))
  
  if(!is.null(output_file)){
    write_tsv(ggcomposition2, file = output_file)
  }
  
  return(ggcomposition2)
}

###################
split_tax_from_asvtable <- function(asv_tax_table, input_table = FALSE, output_file = NULL){
  if(input_table){
    Asv_tax_table <- read_tsv(asv_tax_table)
  }else{
    Asv_tax_table <- asv_tax_table
  }
  tax <-  Asv_tax_table %>% 
    select(ASV_ID, taxonomy) %>% 
    tidyr::separate(taxonomy, sep = ";",
                    into = c("kindom","phylum","class","order","family","genus","species"))
  if(!is.null(output_file)){
    write_tsv(tax, file = output_file)
  }
  tax 
}

Se <- function(x) {
  c <- sd(x) / sqrt(length(x))
  return(c)
}

cat("Rfunction source process is running \n")
