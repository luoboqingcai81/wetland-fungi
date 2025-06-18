
# Wetland Fungal Diversity and Biomass Analysis

This repository contains data and code for investigating fungal diversity and biomass in anoxic wetlands, including analysis of spatial variation, environmental drivers, and latitudinal patterns.

## Repository Structure

### Data Directory
`/data/`
- Contains all raw data used for analyses and figure generation

### Main Visualization Script
`script/wetlandfungi_allfigures.R`
- Generates all main publication figures

### Processing Pipeline

**Processing Scripts:**
1. **ASV Table Processing**  
   `01.ASV_table_process.R`  
   - Cleans ASV tables  
   - Calculates fungal richness  
   - Generates community composition matrices  
   - *Requires:* `Rfunction.R`

2. **Spatial Variation Analysis**  
   `02.RF_predict.R`  
   - Models global spatial patterns using Random Forests  
   - Predicts fungal richness and biomass across anoxic wetlands  

3. **Environmental Drivers**  
   `03.stats.R`  
   - Analyzes effects of environmental factors on:  
     - Fungal richness  
     - Biomass carbon  

4. **Latitudinal Patterns**  
   `summarise_slice_lat.R`  
   - Calculates latitudinal means of fungal richness  
   - Computes standard deviations across gradients

**Supporting Functions:**
- `Rfunction.R`: Core functions for ASV processing
- `map_function.R`: Geographic visualization functions

