#=========================================Packages loading==========================
library(dplyr)
library(readr)
library(tidyverse)
library(terra)
library(tidyterra) 
library(ggplot2)
library(rgdal)
library(purrr)
library(RColorBrewer)
library(sf)
library(magrittr)
library(rcolors)
library(cowplot)
library(readr)
library(reshape2)
library(hrbrthemes)
library(viridis)
library(ggpubr)
library(ggdist)
library(gghalves) 
library(patchwork)
library(ggnewscale)
library(readxl)
library(ggtrendline)
library(fs)
library(ggthemes)
library(broom)
#==============================================0.Data loading==========================
setwd("C:\\Users\\11986\\OneDrive\\Desktop-office\\metafungi-LX\\WETLAND_FUNGI\\DATA\\")
map_richness<-read_xlsx("Wetland_fungi.xlsx",sheet="richness_site")
map_biomass<-read_xlsx("Wetland_fungi.xlsx",sheet="Fmass_site")
data_richness<-read_xlsx("Wetland_fungi.xlsx",sheet="richness_notundra_topsoil")
data_biomass<-read_xlsx("Wetland_fungi.xlsx",sheet="Fmass_notundra_topsoil")
lat_r<-read_xlsx("Wetland_fungi.xlsx",sheet="Latitude_by_a1b1_1_richness")
lat_b<-read_xlsx("Wetland_fungi.xlsx",sheet="Latitude_by_a1b1_1_FunB")
compare<-read_xlsx("Wetland_fungi.xlsx",sheet="Comparison")
area_richness<-read_xlsx("Wetland_fungi.xlsx",sheet="Latitude_by_a1b1_richness_area")
glmm_richness<-read_xlsx("Wetland_fungi.xlsx",sheet="glmmhp_lmer_fungi_richnes")
glmm_biomass<-read_xlsx("Wetland_fungi.xlsx",sheet="glmmhp_lmer_fungi_FunB")
phylum_fungi_top12_global <- read_csv("phylum_fungi_top12_global.csv")
phylum_fungi_top12_byGroup <- read_csv("phylum_fungi_top12_byGroup.csv")

ASV_fungaltrait_bytype_forPie <- read_tsv("ASV_fungaltrait_for_ggplot.txt")%>% 
  pivot_wider(names_from = Group, values_from = value) %>% 
  mutate(trait = factor(trait,levels = c("Symbiotrophs","Saprotrophs","Pathotrophs","Other","unspecified","Unknown"))) %>% 
  arrange(trait)

ASV_fungalAquatic_bytype <- read_tsv("ASV_fungalAquatic_for_ggplot.txt")%>% 
  mutate(Aquatic = factor(Aquatic,levels = c("aquatic","partly-aquatic","non-aquatic","unspecified","Unknown"))) %>% 
  arrange(Aquatic)


wetland_base <- terra::rast("GLWD-3\\glwd_3") 
wetland_one <- ifel(wetland_base >= 1 & wetland_base <= 12, 1, NA)
sf_fn <- st_read("borde/è¾¹ç•Œ.shx") 
crs_84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
st_crs(sf_fn) <- crs_84
#sf_fn_new <-  terra::project(sf_fn, crs_84)
sf_fn_new <- st_transform(sf_fn, crs = 4326)
world <- map_data("world")

#==============================================1.Plotting ==========================
##Figure 1A--richness map-------------------------------
p1.1<-ggplot()+
 geom_spatraster(data = wetland_one, show.legend = F) +
  scale_fill_gradientn(
    colours ="grey60" ,
    na.value = NA,
  ) +
  geom_spatvector(data = sf_fn, fill = NA, linewidth = 0.1) +
  geom_point(data=map_richness,aes(x=Long,y=Lat,size=N,color=richness),shape=16,alpha=0.6)+
  scale_color_gradient(low = "#aea536" , high = "#378000") +
  coord_sf() +
  ylab("Latitude (N)")+
  xlab("Longitude (E)")+
  theme_test()+
  theme(
    axis.text.x = element_text(color = "black", size = 12),  # ÏÔÊ¾xÖá×ø±ê¿Ì¶È
    axis.text.y = element_text(color = "black", size = 12),  # ÏÔÊ¾yÖá×ø±ê¿Ì¶È
    panel.grid.major = element_line(color = "gray", linetype = "dotted"),
    plot.background = element_rect(fill = "white"),
    panel.background = element_rect(fill = "white", colour = "black")
  )



##Figure 1B--Richness.total-------------------------------
p1.2_1<-ggplot(data=data_richness,aes(x=group,y=richness),color=group)+
  
  stat_halfeye(mapping = aes(fill=group),width = 0.4, justification = -0.7,.width = 0,point_colour = NA)+
  geom_jitter(aes(color=group),shape=21,position=position_jitter(0.2), size=2.5)+
  scale_color_manual(values=c("#AECFD4"))+
  scale_fill_manual(values=c("#AECFD4"))+
  geom_boxplot(aes(),position = position_nudge(),linetype="solid",width=0.2,size=1,fill=NA)+
  labs(x="", y="Fungal richness")+
  theme_classic()+
  theme(
    legend.position="none", 
    axis.title.x = element_text(face = "bold",color = "black",size =18),
    axis.title.y = element_text(face = "bold",color = "black",size = 18),
    axis.text.x = element_text(color = "black",size = 15),
    axis.text.y= element_text(color = "black",size = 15))

###Richness.compare
kruskal.test(richness ~ Wetland_type, data = data_richness)
pairwise.wilcox.test(data_richness$richness, data_richness$Wetland_type,
                     p.adjust.method = "BH")

kruskal.test(richness ~ salted, data = data_richness)
pairwise.wilcox.test(data_richness$richness, data_richness$salted,
                     p.adjust.method = "BH")

data_richness$Wetland_type<-factor(data_richness$Wetland_type,levels = c("Coastal","Inland","Peatland"), ordered = TRUE)

p1.2_2<-ggplot(data_richness,mapping = aes(x=Wetland_type,y=richness,color=Wetland_type)) +
  stat_halfeye(mapping = aes(fill=Wetland_type),width = 0.4, justification = -0.7,.width = 0,point_colour = NA)+
  geom_jitter(shape=21,position=position_jitter(0.2), size=2.5)+
  scale_color_manual(values=c("#C5dff4","#aeb2d1","#d9b9d4"))+
  scale_fill_manual(values=c("#C5dff4","#aeb2d1","#d9b9d4"))+
  labs(x="", y="Shannon diversity")+theme_bw(base_line_size = 1.05,base_rect_size =1.05)+
  geom_boxplot(aes(),position = position_nudge(),linetype="solid",width=0.2,size=1,color="black",fill=NA)+
  theme_bw()+ 
  theme(
    legend.position="none",
    axis.title.x = element_text(face = "bold",color = "black",size =18),
    axis.title.y = element_text(face = "bold",color = "black",size = 18),
    axis.text.x = element_text(color = "black",size = 15),
    axis.text.y= element_text(color = "black",size = 15),
    panel.background = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    plot.background = element_blank())+
    coord_flip() 

data_richness$salted<-factor(data_richness$salted,levels = c("yes","no"), ordered = TRUE)

p1.2_3<-ggplot(data_richness,mapping = aes(x=salted,y=richness,color=salted)) +
  stat_halfeye(mapping = aes(fill=salted),width = 0.4, justification = -0.7,.width = 0,point_colour = NA)+
  geom_jitter(shape=21,position=position_jitter(0.2), size=2.5)+
  geom_boxplot(aes(),position = position_nudge(),linetype="solid",width=0.2,size=1,color="black",fill=NA)+
  scale_color_manual(values=c("#daa87c","#f4eeac"))+
  scale_fill_manual(values=c("#daa87c","#f4eeac"))+
  labs(x="", y="richness")+theme_bw(base_line_size = 1.05,base_rect_size =1.05)+
  coord_flip() +
  theme_classic() 

##Figure 1C--The relationship between fungal richness and wetland area from study sites-------------------------------
area_richness_site<-read.csv("raw_richness_area_mean.csv",sep=",")
area_richness_site<-area_richness_site[-c(1,2,4),]
y<-area_richness_site$richness
x<-area_richness_site$area

p1.3<-ggtrendline(x, y, model = "log2P",linecolor = "#e99e00",linetype = 1,linewidth = 1.3,CI.fill = "#237cee",CI.alpha = 0.3,CI.lty = 0)  +
  xlab("Wetland area (km2)")+
  ylab("Fungal richness")+
  theme_bw()+
  theme(
    legend.position="none", 
    axis.title.x = element_text(face = "bold",color = "black",size =16),
    axis.title.y = element_text(face = "bold",color = "black",size = 16),
    axis.text.x = element_text(color = "black",size = 12),
    axis.text.y= element_text(color = "black",size = 12),
    panel.background = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    plot.background = element_blank())

##Figure 1D--Fungal composition-------------------------------
### Loading Pie Chart Custom Functions
p_pie <- function(df, group, Colors_set) {
  group = ensym(group)
  string_group = as_string(group)
  Tax_name <- colnames(df)[1]
  tax_name <- sym(Tax_name)
  tax_name
  group_f <- df %>%
    select(!!tax_name, !!group) %>%
    mutate(prop = round(!!group, digits = 1)) %>%
    mutate(lab_ypos = c(rep(100, length(prop)) -
                          c(0, cumsum(prop)[-length(prop)]) - 0.5 * prop)) %>% 
    mutate(!!tax_name := factor(!!tax_name, levels = !!tax_name))
  group_f
  ggplot(group_f, aes(x = "", y = prop, fill = !!tax_name)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0, direction = 1) +
    theme_void() +
    geom_text(aes(y = lab_ypos, x = sum(prop) / 75, label = prop),
              color = "white", size = 3) +
    scale_fill_manual(values = Colors_set) +
    theme(panel.border = element_rect(fill = NA, colour = NA)) +
    ggtitle(string_group) +
    theme(plot.title = element_text(hjust=0.5, face = "bold", size = 18),
          legend.title = element_text(size = 18, face = "bold"))
}


### Creating Forms
dir_create("asv_process/picture")

### Fungal composition,phylum level

Colors_set <- c("#3399cc", "#99cccc","#6699cc" , "#99ccff" , 
                "#9999cc", "#ccccff", "#cc99cc", "#ff99cc", "#ffcccc", 
                "#669966", "#66cc99", "#9E9E9E")

Colors_set <- set_names(Colors_set, nm = phylum_fungi_top12_byGroup$Phylum )

as_string <- function(x) {
  as.character(x)
}

gg_global_fungi <- p_pie(phylum_fungi_top12_byGroup, Global,Colors_set)
gg_peatland_fungi <- p_pie(phylum_fungi_top12_byGroup, Peatland,Colors_set)
gg_inland_fungi <- p_pie(phylum_fungi_top12_byGroup, Inland,Colors_set)
gg_coastal_fungi <- p_pie(phylum_fungi_top12_byGroup, Coastal,Colors_set)

layout_pie <- c(
  area(1,1,2,3),
  area(3,1),
  area(3,2),
  area(3,3)
)

p1.4 <- gg_global_fungi + gg_peatland_fungi + gg_inland_fungi + gg_coastal_fungi  + plot_layout(guides = 'collect', design = layout_pie)
p1.4
ggsave(p1.4, filename = "01.Phylum_fungi_composition.pdf", height = 10, width = 16)


##Figure 1E--Fungal habitat trait ---------------------------------------
Colors_set <- c("#66cc99","#3399cc", "#99cccc","#6699cc" , "#9E9E9E")
Colors_set <- set_names(Colors_set, nm = c("aquatic","partly-aquatic","non-aquatic","unspecified","Unknown"))

gg_global_fungi <- p_pie(ASV_fungalAquatic_bytype, Global,Colors_set)
gg_peatland_fungi <- p_pie(ASV_fungalAquatic_bytype, Peatland,Colors_set)
gg_inland_fungi <- p_pie(ASV_fungalAquatic_bytype, Inland,Colors_set)
gg_coastal_fungi <- p_pie(ASV_fungalAquatic_bytype, Coastal,Colors_set)

layout_pie <- c(
  area(1,1,2,3),
  area(3,1),
  area(3,2),
  area(3,3)
)

p1.5 <- gg_global_fungi + gg_peatland_fungi + gg_inland_fungi + gg_coastal_fungi  + plot_layout(guides = 'collect', design = layout_pie)
p1.5

## Figure 1F--Fungal trait composition --------------------------------
Colors_set <- c("#66cc99", "#3399cc","#ffcccc" , "#99ccff", "#9999cc","#9E9E9E")
Colors_set <- set_names(Colors_set, nm = c("Symbiotrophs","Saprotrophs","Pathotrophs","Other","unspecified","Unknown"))

gg_global_fungitrait <- p_pie(ASV_fungaltrait_bytype_forPie, Global,Colors_set)
gg_peatland_fungitrait <- p_pie(ASV_fungaltrait_bytype_forPie, Peatland,Colors_set)
gg_inland_fungitrait <- p_pie(ASV_fungaltrait_bytype_forPie, Inland,Colors_set)
gg_coastal_fungitrait <- p_pie(ASV_fungaltrait_bytype_forPie, Coastal,Colors_set)

layout_pie <- c(
  area(1,1,2,3),
  area(3,1),
  area(3,2),
  area(3,3)
)

p1.6 <- gg_global_fungitrait + gg_peatland_fungitrait + gg_inland_fungitrait + gg_coastal_fungitrait  + plot_layout(guides = 'collect', design = layout_pie)
p1.6

##Figure 2A--Biomass map-------------------------------
p2.1<-ggplot()+
  geom_spatraster(data = wetland_one, show.legend = F) +
  scale_fill_gradientn(
    colours ="grey60" ,
    na.value = NA,
  ) +
  geom_spatvector(data = sf_fn, fill = NA, linewidth = 0.1) +
  geom_point(data=map_biomass,aes(x=Long,y=Lat,size=N,color=Fmass),shape=16,alpha=0.6)+
  scale_color_gradient(low = "#aea536" , high = "#378000") +
  ylab("Latitude (N)")+
  xlab("Longitude (E)")+
  coord_sf() +
  theme_test()+
  theme(
    axis.text.x = element_text(color = "black", size = 12),  # ÏÔÊ¾xÖá×ø±ê¿Ì¶È
    axis.text.y = element_text(color = "black", size = 12),  # ÏÔÊ¾yÖá×ø±ê¿Ì¶È
    panel.grid.major = element_line(color = "gray", linetype = "dotted"),
    plot.background = element_rect(fill = "white"),
    panel.background = element_rect(fill = "white", colour = "black")
  )


##Figure 2B--Biomass total ----------------------------------------------
p2.2_1<-ggplot(data=data_biomass,aes(x=group,y=Fmass),color=group)+
  stat_halfeye(mapping = aes(fill=group),width = 0.4, justification = -0.7,.width = 0,point_colour = NA)+
  geom_jitter(aes(color=group),shape=21,position=position_jitter(0.2), size=2.5)+
  scale_color_manual(values=c("#AECFD4"))+
  scale_fill_manual(values=c("#AECFD4"))+
  ylim(0,2000)+
  ylab(expression(Fungal~biomass~C~(mg~kg^-1~soil)))+
  xlab("")+
  geom_boxplot(aes(),position = position_nudge(),linetype="solid",width=0.2,size=1,fill=NA)+
  theme_classic()+
  theme(
    legend.position="none", 
    axis.title.x = element_text(face = "bold",color = "black",size =18),
    axis.title.y = element_text(face = "bold",color = "black",size = 18),
    axis.text.x = element_text(color = "black",size = 15),
    axis.text.y= element_text(color = "black",size = 15))

### Biomass compare
data_biomass2<-subset(data_biomass,data_biomass$Fmass<2000)

kruskal.test(Fmass ~ Type, data = data_biomass2)
pairwise.wilcox.test(data_biomass2$Fmass, data_biomass2$Type,
                     p.adjust.method = "BH")

kruskal.test(Fmass ~ Salinity, data = data_biomass2)
pairwise.wilcox.test(data_biomass2$Fmass, data_biomass2$Salinity,
                     p.adjust.method = "BH")

data_biomass2$Type<-factor(data_biomass2$Type,levels = c("Coastal","Inland","Peatland"), ordered = TRUE)

p2.2_2<-ggplot(data_biomass2,mapping = aes(x=Type,y=Fmass,color=Type)) +
  stat_halfeye(mapping = aes(fill=Type),width = 0.4, justification = -0.7,.width = 0,point_colour = NA)+
  geom_jitter(shape=21,position=position_jitter(0.2), size=2.5)+
  scale_color_manual(values=c("
                              ","#aeb2d1","#d9b9d4"))+
  scale_fill_manual(values=c("#C5dff4","#aeb2d1","#d9b9d4"))+
  geom_boxplot(aes(),position = position_nudge(),linetype="solid",width=0.2,size=1,color="black",fill=NA)+
  
  labs(x="", y="Log (F:B biomass ratio)")+theme_bw(base_line_size = 1.05,base_rect_size =1.05)+
  theme_test()+ 
  coord_flip()+
  theme(
    legend.position="none",
    axis.title.x = element_text(face = "bold",color = "black",size =18),
    axis.title.y = element_text(face = "bold",color = "black",size = 18),
    axis.text.x = element_text(color = "black",size = 15),
    axis.text.y= element_text(color = "black",size = 15),
    panel.background = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    plot.background = element_blank())

data_biomass2$Salinity<-factor(data_biomass2$Salinity,levels = c("yes","no"), ordered = TRUE)

p2.2_3<-ggplot(data_biomass2,mapping = aes(x=Salinity,y=Fmass,color=Salinity)) +
  stat_halfeye(mapping = aes(fill=Salinity),width = 0.4, justification = -0.7,.width = 0,point_colour = NA)+
  geom_jitter(shape=21,position=position_jitter(0.2), size=2.5)+
  geom_boxplot(aes(),position = position_nudge(),linetype="solid",width=0.2,size=1,color="black",,fill=NA)+
  scale_color_manual(values=c("#daa87c","#f4eeac"))+
  scale_fill_manual(values=c("#daa87c","#f4eeac"))+
  labs(x="", y="Fungal biomass carbon")+theme_bw(base_line_size = 1.05,base_rect_size =1.05)+
  coord_flip() +
  theme_classic()

##Figure 2C--Relationships betweem SOC and fungal biomass carbon ----------------------------------------------
setwd("C:\\Users\\11986\\OneDrive\\Desktop-office\\metafungi-LX\\WETLAND_FUNGI\\DATA\\")
biomass<-read.csv("soc_bio.csv")

fit<-lm(SOC~FunBC,biomass)
model_summary<-summary(fit)
r_squared <- model_summary$r.squared
p_value <- coef(model_summary)[2, 4]
r_squared_text <- paste("R^2 = ", round(r_squared, 3))
p_value_text <- paste("P = ", format.pval(p_value, digits = 3))

label_data <- data.frame(
  x = Inf,
  y = Inf,
  label = paste(r_squared_text, "\n", p_value_text)
)

p2.3<-ggplot(data=biomass,aes(x=FunBC, y=SOC))+
  geom_point(data=biomass,aes(x=FunBC,y=SOC),shape=21,color="#e99e00",stroke=2)+
  geom_smooth(method="lm",se=T,linewidth=1.5,color="black")+
  ylab("Soil organic carbon (g-1 kg)")+
  xlab("Fungal biomass carbon (mg-1 kg)")+
  geom_text(data = label_data, aes(x = x, y = y, label = label),
            hjust = 1, vjust = 1, size = 5, color = "black") +
  theme_bw()+
  theme(
    panel.background = element_blank(), 
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(), 
    plot.background = element_blank())+
  theme(axis.text.x = element_text(size=16,color="black"),
        axis.text.y = element_text(size=16,color="black"),
        axis.title.y = element_text(size=18,color="black"),
        axis.title.x = element_text(size=18,color="black"),
        legend.title=element_blank(),
        legend.text = element_text(size=14))

##Figure 3A--Effect of factors on richness-------------------------------
glmm_richness$Varibles<-factor(glmm_richness$Varibles,levels = c("Annual_Mean_Temperature",
                                                                 "Temperature_Annual_Range",
                                                                 "Temperature_Seasonality",
                                                                 "Annual_Precipitation",
                                                                 "Precipitation_Seasonality",
                                                                 "MeanMonthlyMoisture_Index",
                                                                 "AridityIndex",
                                                                 "Potential_Evapotranspiration",
                                                                 "Shannon_Index",
                                                                 "NPP",
                                                                 "Nitrogen_Content",
                                                                 "SOC_Content",
                                                                 "Soil_pH"), ordered = TRUE)
glmm_richness$Type<-factor(glmm_richness$Type,levels = c("soil","plant","aridityindex","precipitation","temperature"), ordered = TRUE)

p3.1_1<-ggplot(data=glmm_richness, aes(x=Std_Coefficient, y=Varibles,color=Type))+ 
  geom_vline(xintercept=0, color = "grey",linetype="dashed", lwd=1)+
  geom_point(size=6)+
  geom_errorbar(aes(xmin=CI_low, xmax=CI_high), 
                width=.2,lwd=1,color="grey40") +
  scale_color_manual(values =c("grey","#C2DCBF","#E7D8FA","#BAD8FB","#FBE5BE"))+ 
  theme_classic()+ 
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(size = 14, color = "black"))+
  theme(axis.text.y = element_text(size = 14, color = "black"))+
  theme(title=element_text(size=16))+
  labs(title = "Fungal richness",x = "Standardized coefficients",y="")+
  scale_y_discrete(position = "right")+
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank())

p3.1_2<-ggplot(glmm_richness, aes(x="",y=percent,fill=Type))+
  geom_bar(stat="identity")+ 
  theme(axis.text.x = element_blank())+labs(x = "",y = "Relative effect (%)")+
  theme_classic()+
  scale_fill_manual(values =c("grey","#C2DCBF","#E7D8FA","#BAD8FB","#FBE5BE"))+
  theme(legend.position = "none",
        axis.text.x=element_blank(),
        axis.line.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(axis.text.x = element_text(size = 14, color = "black"))+
  theme(axis.text.y = element_text(size = 14, color = "black"))+
  theme(title=element_text(size=16))

p3.1_2+p3.1_1+plot_layout(widths = c(2,3))


##Figure 3B--Effect of factors on biomass -------------------------------
glmm_biomass$Varibles<-factor(glmm_biomass$Varibles,levels = c("Annual_Mean_Temperature",
                                                               "Temperature_Annual_Range",
                                                               "Temperature_Seasonality",
                                                               "Annual_Precipitation",
                                                               "Precipitation_Seasonality",
                                                               "MeanMonthlyMoisture_Index",
                                                               "AridityIndex",
                                                               "Potential_Evapotranspiration",
                                                               "Shannon_Index",
                                                               "NPP",
                                                               "Nitrogen_Content",
                                                               "SOC_Content",
                                                               "Soil_pH"), ordered = TRUE)

glmm_biomass$Type<-factor(glmm_biomass$Type,levels = c("soil","plant","aridityindex","precipitation","temperature"), ordered = TRUE)

p3.2_1<-ggplot(data=glmm_biomass, aes(x=Std_Coefficient, y=Varibles,color=Type))+ 
  geom_vline(xintercept=0, color = "grey",linetype="dashed", lwd=1)+
  geom_point(size=6)+
  geom_errorbar(aes(xmin=CI_low, xmax=CI_high), 
                width=.2,lwd=1,color="grey40") +
  scale_color_manual(values =c("grey","#C2DCBF","#E7D8FA","#BAD8FB","#FBE5BE"))+ 
  theme_classic()+ 
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(size = 14, color = "black"))+
  theme(axis.text.y = element_text(size = 14, color = "black"))+
  theme(title=element_text(size=16))+
  labs(title = "Fungal biomass",x = "Standardized coefficients",y="")+
  scale_y_discrete(position = "right")+
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank())

p3.2_2<-ggplot(glmm_biomass, aes(x="",y=percent,fill=Type))+
  geom_bar(stat="identity")+ 
  theme(axis.text.x = element_blank())+labs(x = "",y = "Relative effect (%)")+
  theme_classic()+
  scale_fill_manual(values =c("grey","#C2DCBF","#E7D8FA","#BAD8FB","#FBE5BE"))+
  theme(legend.position = "none",
        axis.text.x=element_blank(),
        axis.line.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(axis.text.x = element_text(size = 14, color = "black"))+
  theme(axis.text.y = element_text(size = 14, color = "black"))+
  theme(title=element_text(size=16))

p3.2_2+p3.2_1+plot_layout(widths = c(2,3))

##Figure 4A--Changes in latitude of richness-------------------------------
p4.1 <- ggplot(lat_r,aes(x= Latitude, y= Mean)) +
  geom_ribbon(aes(ymin= Min, ymax= Max), fill="mediumpurple2",alpha=0.7) +
  geom_line(col="blue", lwd=1) +
  coord_flip()+
  theme(
    panel.border = element_rect(color="black",linewidth = 1,fill=NA),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(face = "bold", color = "black", size = 15),
    axis.text.y = element_text(face = "bold", color = "black", size = 15),
    axis.line.x = element_line(color = "black", size = 1),
    axis.line.y = element_line(color = "black", size = 1),
    axis.ticks.length.x = unit(0.20, "cm"),
    axis.ticks.length.y = unit(0.20, "cm"),
    panel.background = element_blank(),
    legend.position = "bottom"
  ) +
  scale_x_continuous(name=expression("Latitude (" * degree * ")"),
                     breaks=c( -60,-45, -30, -15, 0,15, 30, 45, 60,75),
                     labels=c(  -60,-45, -30, -15, 0,15, 30, 45, 60,75),
                     limits=c(-60, 80))


##Figure 4B--Changes in latitude of biomass -----------------------------
p4.2 <- ggplot(lat_b,aes(x= Latitude, y= Mean)) +
  geom_ribbon(aes(ymin= Min, ymax= Max), fill="mediumpurple2",alpha=0.7) +
  geom_line(col="blue", lwd=1) +
  coord_flip()+
  theme(
    panel.border = element_rect(color="black",linewidth = 1,fill=NA),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(face = "bold", color = "black", size = 15),
    axis.text.y = element_text(face = "bold", color = "black", size = 15),
    axis.line.x = element_line(color = "black", size = 1),
    axis.line.y = element_line(color = "black", size = 1),
    
    axis.ticks.length.x = unit(0.20, "cm"),
    axis.ticks.length.y = unit(0.20, "cm"),
    panel.background = element_blank(),
    legend.position = "bottom"
  ) +
  scale_x_continuous(name=expression("Latitude (" * degree * ")"),
                     breaks=c( -60,-45, -30, -15, 0,15, 30, 45, 60,75),
                     labels=c(  -60,-45, -30, -15, 0,15, 30, 45, 60,75),
                     limits=c(-60, 80))

##Figure S4--Comparing richness and biomass carbon-------------------------------
Sp4.1 = ggplot(compare, aes(x=factor(biome,levels = c("Boreal forest","Temperate forest","Tropical/subtropical forest","Shrub","Grassland","Wetland_Woody","Wetland_Shrub/Moss","Wetland_Herbage","Wetland_Bare","Wetland_Total"),order=TRUE),  # è½¬åŒ–ä¸ºå› å­ï¼Œç›®çš„æ˜¯æ˜¾ç¤ºé¡ºåºä¸Žæ–‡ä»¶é¡ºåºç›¸åŒï¼Œå¦åˆ™æŒ‰ç…§å­—æ¯é¡ºåºæŽ’åº?
                            y=-richness), 
)+
  labs(
    x="",   
    y="Fungal richness",) +
  geom_bar(data = compare,aes(fill=biome),
           stat="identity",
           position=position_dodge(),width=0.7)+ 
  scale_fill_manual(values = c("grey80","grey80","grey80","grey80","grey80","#009583","#009583","#009583","#009583","#009583"))+
  geom_hline(yintercept = c(-186.74),linetype="dashed",colour="grey40",size=0.6)+
  coord_flip()+
  theme_pander()+ 
  theme(
    legend.position="none", 
    axis.title.x = element_text(face = "bold",color = "black",size =18),
    axis.title.y = element_text(face = "bold",color = "black",size = 18),
    axis.text.x = element_text(color = "black",size = 12),
    axis.text.y= element_text(color = "black",size = 12),
    panel.background = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    plot.background = element_blank())

###Biomass C
sp4.2 = ggplot(compare, aes(x=factor(biome,levels = c("Boreal forest","Temperate forest","Tropical/subtropical forest","Shrub","Grassland","Wetland_Woody","Wetland_Shrub/Moss","Wetland_Herbage","Wetland_Bare","Wetland_Total"),order=TRUE),  # è½¬åŒ–ä¸ºå› å­ï¼Œç›®çš„æ˜¯æ˜¾ç¤ºé¡ºåºä¸Žæ–‡ä»¶é¡ºåºç›¸åŒï¼Œå¦åˆ™æŒ‰ç…§å­—æ¯é¡ºåºæŽ’åº?
                            y=FC), 
)+
  labs(
    x="",   
    y="Fungal biomass C (mg kg -1 soil)",) +
  geom_bar(data = compare,aes(fill=biome),
           stat="identity",
           position=position_dodge(),width=0.7)+ 
  scale_fill_manual(values = c("grey80","grey80","grey80","grey80","grey80","#a4b3e8","#a4b3e8","#a4b3e8","#a4b3e8","#a4b3e8"))+
  geom_hline(yintercept = c(425.34),linetype="dashed",colour="grey40",size=0.6)+
  coord_flip()+
  theme_pander()+ 
  theme(
    legend.position="none", 
    axis.title.x = element_text(face = "bold",color = "black",size =18),
    axis.title.y = element_text(face = "bold",color = "black",size = 18),
    axis.text.x = element_text(color = "black",size = 12),
    axis.text.y= element_text(color = "black",size = 12),
    panel.background = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    plot.background = element_blank())








