library(sf)
library(stars)
library(tidyverse)
library(ggplot2)
library(sp)
library(tmap)
library(spdep)
library(raster)
library(fasterize)
library(snow)
library(RColorBrewer)
library(rasterVis)    
library(colorspace)
library(readr)
library(gridExtra)
library(readxl)
library(MMWRweek)
library(maps)
library(plotrix)
library(terra)
library(readxl)
library(tidygeocoder)
library(writexl)
library(phylin)
library(spatstat)
library(tidycensus)
library(prism)
library(reshape)
library(stars)
library(rsample)
library(caret)
library(h2o)
library(ranger)
library(recipes)
library(Metrics)
library(vip)
library(vivid) 

setwd('/Users/nathanialodell/Library/CloudStorage/OneDrive-TexasTechUniversity/JRM_Tarrant County Aedes/NODell Tarrant/Modeling (again)')
tx <- read_sf("DFW No Residuals.shp")

### extracting variables and renaming them for the sake of plotting
tx <- data.frame(aeg = tx$AEGYPTI, alb = tx$ALBOPIC, 
                 lat = tx$lat, lon = tx$lon, wk = tx$Epi_wk,
                 fMonth = tx$fMonth)

assign_season <- function(month) {
  if (month %in% c("December", "January", "Feburary")) return("Winter")
  if (month %in% c("March", "April", "May")) return("Spring")
  if (month %in% c("June", "July", "August")) return("Summer")
  if (month %in% c("September", "October", "November")) return("Fall")
}

### add a season column to the data
tx <- tx %>%
  mutate(Season = sapply(fMonth, assign_season))

### Split data into subsets by season
process_season_data <- function(season_data) {
  season_data %>%
    group_by(lat,lon) %>%
    summarise(
      AEGYPTI = sum(aeg),
      ALBOPIC = sum(alb),
      Effort = n()
    ) %>%
    mutate(
      AEGYPTI = AEGYPTI / Effort,
      ALBOPIC = ALBOPIC / Effort,
      Dominant_Species = case_when(
        AEGYPTI > ALBOPIC ~ "AEGYPTI",
        ALBOPIC > AEGYPTI ~ "ALBOPICTUS",
        TRUE ~ "EQUAL"
      )
      )%>%
    data.frame()
}

### split data by season and process
tx_spring <- process_season_data(tx %>% filter(Season == "Spring"))%>%
  st_as_sf(coords = c("lon", "lat"))
tx_summer <- process_season_data(tx %>% filter(Season == "Summer"))%>%
  st_as_sf(coords = c("lon", "lat"))
tx_fall <- process_season_data(tx %>% filter(Season == "Fall"))%>%
  st_as_sf(coords = c("lon", "lat"))

tx_spring$Season = "Spring"
tx_fall$Season = "Fall"
tx_summer$Season = "Summer"
ov <- rbind(tx_spring, tx_summer, tx_fall)

### overall species dominance counts 
table(ov$Season, ov$Dominant_Species)

plot_season <- function(data, season_name) {
  ggplot(data) +
    geom_sf(aes(color = Dominant_Species), size = 5) +
    scale_color_manual(values = c("AEGYPTI" = "pink", "ALBOPICTUS" = "dodgerblue",
                                  "EQUAL" = "grey")) +
    theme_minimal() +
    labs(title = paste(season_name),
         size = "Count",
         color = "Dominant Species") +
    theme(legend.position = "right")
}

### plot for each season
p_spring <- plot_season(tx_spring, "Predicted Spring")
p_summer <- plot_season(tx_summer, "Predicted Summer")
p_fall <- plot_season(tx_fall, "Predicted Fall")

### display plots
print(p_spring)
print(p_summer)
print(p_fall)
print(p_winter)
library(ggpubr)
ggarrange(p_spring, p_summer, p_fall)