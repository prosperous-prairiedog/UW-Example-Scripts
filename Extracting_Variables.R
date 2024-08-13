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

setwd('/Users/nathanialodell/Library/CloudStorage/OneDrive-TexasTechUniversity/JRM_Tarrant County Aedes/NODell Tarrant/Tarrant pseudo-ab')
tx1 <- read_sf("Aeg Point Resid.shp") 
tx2 <- read_sf("Alb Point Resid.shp")
tx <- bind_rows(tx1, tx2)
state_map <- map('county', fill = TRUE, plot = FALSE) %>% st_as_sf() 
######################################
### landcover
setwd("/Users/nathanialodell/Library/CloudStorage/OneDrive-TexasTechUniversity/JRM_Tarrant County Aedes")
wrd.pp.05 <- raster("gpw_v4_population_density_adjusted_to_2015_unwpp_country_totals_rev11_2005_30_sec.tif")
wrd.pp.10 <- raster("gpw_v4_population_density_adjusted_to_2015_unwpp_country_totals_rev11_2010_30_sec.tif")
wrd.pp.15 <- raster("gpw_v4_population_density_adjusted_to_2015_unwpp_country_totals_rev11_2015_30_sec.tif")
wrd.pp.20 <- raster("gpw_v4_population_density_adjusted_to_2015_unwpp_country_totals_rev11_2020_30_sec.tif")

tx_sf1 <- st_transform(tx, crs(wrd.pp.05))
state_map2 <- st_transform(state_map, crs(wrd.pp.05))
e <- extent(tx_sf1)

tx.pop.05 <- crop(wrd.pp.05, e)
tx.pop.10 <- crop(wrd.pp.10, e)
tx.pop.15 <- crop(wrd.pp.15, e)
tx.pop.20 <- crop(wrd.pp.20, e)

ns.line = st_sfc(st_linestring(rbind(c(-97.53708, 32.49527), 
                                     c(-97.53708, 32.99781))), crs = crs(tx))
st_length(ns.line)/1000

ew.line = st_sfc(st_linestring(rbind(c(-97.5311, 32.49527), 
                                     c(-96.61046, 32.49527))), crs = crs(tx))
st_length(ew.line)/1000  

r <- raster(e, nrow=86, ncol=56, vals=1,crs=crs(tx))

### load available variables if needed to look through things
v.21 <- load_variables(2021, "acs1", cache = TRUE)

dens.fun <- function(data, YR) {
  pop.mean.15 <- raster::extract(x=data, y=tx%>%filter(Year == YR))
  pop.mean.15.2a <- unlist(lapply(pop.mean.15, mean))
  
  tx_sf1.x <- tx %>% filter(Year==YR)
  tx_sf1.x$Dens <- pop.mean.15.2a
  return(tx_sf1.x)
}

dens.15 <- dens.fun(tx.pop.15, 2015)
dens.16 <- dens.fun(tx.pop.15, 2016)
dens.17 <- dens.fun(tx.pop.15, 2017)
dens.18 <- dens.fun(tx.pop.15, 2018)
dens.19 <- dens.fun(tx.pop.15, 2019)
dens.20 <- dens.fun(tx.pop.15, 2020)
dens.21 <- dens.fun(tx.pop.15, 2021)
dens.22 <- dens.fun(tx.pop.15, 2022)
dens <- rbind(dens.15, dens.16, dens.17, dens.18, dens.19, dens.20, dens.21,
              dens.22)


##### need function for impervious surface  
imp.fun <- function(data, YR) {  
  imp <- raster::extract(x=data, y=tx%>%filter(Year == YR))
  imp.2 <- unlist(lapply(imp, mean))
  
  tx_sf1.x <- tx %>% filter(Year==YR)
  tx_sf1.x$Imp <- imp.2
  return(tx_sf1.x)
}

imp.15 <- imp.fun(data=r.imp.13.3, YR=2015)
imp.16 <- imp.fun(data=r.imp.16.3, YR=2016)
imp.17 <- imp.fun(data=r.imp.16.3, YR=2017)
imp.18 <- imp.fun(data=r.imp.16.3, YR=2018)
imp.19 <- imp.fun(data=r.imp.19.3, YR=2019)
imp.20 <- imp.fun(data=r.imp.19.3, YR=2020)
imp.21 <- imp.fun(data=r.imp.21.3, YR=2021)
imp.22 <- imp.fun(data=r.imp.21.3, YR=2022)

imp <- rbind(imp.15, imp.16, imp.17, imp.18, imp.19, imp.20, imp.21,
             imp.22)

#### function to capture acs data based on yr
acs.fun <- function(variable,state,year) {
  
  ### grabbing data  
  acs.dat <- get_acs(
    geography = "tract", 
    variables = variable,
    state = state, 
    year = year,
    geometry = TRUE
  )
  
  ## transform crs and clip to extent of Tarrant date
  acs.dat2 <- st_transform(acs.dat, crs(tx))
  acs.dat2.clip <- st_crop(acs.dat2, xmin= -97.5311, ymin= 32.49527, xmax= -96.61046, ymax= 32.99781)
  acs.cast <- st_cast(acs.dat2.clip, to="MULTIPOLYGON")
  
  ### rasterize and extract value to points codes of CB mp
  acs.r    <- rasterize(acs.cast, r, field="estimate", fun="mean")
  acs.var  <- raster::extract(x=acs.r, y=tx %>% filter(Year==year))
  acs.var2 <- unlist(lapply(acs.var, mean))
  
  tx_sf1.x <- tx %>% filter(Year==year)
  tx_sf1.x$Md.inc <- acs.var2
  return(tx_sf1.x)
  
}

#### median income
Inc.15 <- acs.fun(variable="B19013_001", state="TX", year=2015)
Inc.16 <- acs.fun(variable="B19013_001", state="TX", year=2016)
Inc.17 <- acs.fun(variable="B19013_001", state="TX", year=2017)
Inc.18 <- acs.fun(variable="B19013_001", state="TX", year=2018)
Inc.19 <- acs.fun(variable="B19013_001", state="TX", year=2019)
Inc.20 <- acs.fun(variable="B19013_001", state="TX", year=2020)
Inc.21 <- acs.fun(variable="B19013_001", state="TX", year=2021)
Inc.22 <- acs.fun(variable="B19013_001", state="TX", year=2022)
Inc.all <- bind_rows(Inc.15,Inc.16,Inc.17,Inc.18,
                          Inc.19,Inc.20,Inc.21,Inc.22)



#### median population age
age.15 <- acs.fun(variable="B01002_001", state="TX", year=2015)
age.16 <- acs.fun(variable="B01002_001", state="TX", year=2016)
age.17 <- acs.fun(variable="B01002_001", state="TX", year=2017)
age.18 <- acs.fun(variable="B01002_001", state="TX", year=2018)
age.19 <- acs.fun(variable="B01002_001", state="TX", year=2019)
age.20 <- acs.fun(variable="B01002_001", state="TX", year=2020)
age.21 <- acs.fun(variable="B01002_001", state="TX", year=2021)
age.22 <- acs.fun(variable="B01002_001", state="TX", year=2022)
age.all <- bind_rows(age.15,
                          age.16,age.17,age.18,age.19,age.20,
                          age.21,age.22)

#### median housing age
hage.15 <- acs.fun(variable="B25035_001", state="TX", year=2015)
hage.16 <- acs.fun(variable="B25035_001", state="TX", year=2016)
hage.17 <- acs.fun(variable="B25035_001", state="TX", year=2017)
hage.18 <- acs.fun(variable="B25035_001", state="TX", year=2018)
hage.19 <- acs.fun(variable="B25035_001", state="TX", year=2019)
hage.20 <- acs.fun(variable="B25035_001", state="TX", year=2020)
hage.21 <- acs.fun(variable="B25035_001", state="TX", year=2021)
hage.22 <- acs.fun(variable="B25035_001", state="TX", year=2022)
hage.all <- bind_rows(hage.15,
                           hage.16,
                           hage.17,
                           hage.18,
                           hage.19,
                           hage.20,
                           hage.21,
                           hage.22)

#### median housing value
hInc.15 <- acs.fun(variable="B25077_001", state="TX", year=2015)
hInc.16 <- acs.fun(variable="B25077_001", state="TX", year=2016)
hInc.17 <- acs.fun(variable="B25077_001", state="TX", year=2017)
hInc.18 <- acs.fun(variable="B25077_001", state="TX", year=2018)
hInc.19 <- acs.fun(variable="B25077_001", state="TX", year=2019)
hInc.20 <- acs.fun(variable="B25077_001", state="TX", year=2020)
hInc.21 <- acs.fun(variable="B25077_001", state="TX", year=2021)
hInc.22 <- acs.fun(variable="B25077_001", state="TX", year=2022)
hInc.all <- bind_rows(hInc.15,hInc.16,hInc.17,
                           hInc.18,hInc.19,
                           hInc.20,hInc.21,hInc.22)

### pop below poverty line
pov.15 <- acs.fun(variable="B17001_002", state="TX", year=2015)
pov.16 <- acs.fun(variable="B17001_002", state="TX", year=2016)
pov.17 <- acs.fun(variable="B17001_002", state="TX", year=2017)
pov.18 <- acs.fun(variable="B17001_002", state="TX", year=2018)
pov.19 <- acs.fun(variable="B17001_002", state="TX", year=2019)
pov.20 <- acs.fun(variable="B17001_002", state="TX", year=2020)
pov.21 <- acs.fun(variable="B17001_002", state="TX", year=2021)
pov.22 <- acs.fun(variable="B17001_002", state="TX", year=2022)
pov.all <- bind_rows(pov.15, pov.16, pov.17, pov.18,
                          pov.19, pov.20,
                          pov.21, pov.22)

### pop size
pop.15 <- acs.fun(variable="B01003_001", state="TX", year=2015)
pop.16 <- acs.fun(variable="B01003_001", state="TX", year=2016)
pop.17 <- acs.fun(variable="B01003_001", state="TX", year=2017)
pop.18 <- acs.fun(variable="B01003_001", state="TX", year=2018)
pop.19 <- acs.fun(variable="B01003_001", state="TX", year=2019)
pop.20 <- acs.fun(variable="B01003_001", state="TX", year=2020)
pop.21 <- acs.fun(variable="B01003_001", state="TX", year=2021)
pop.22 <- acs.fun(variable="B01003_001", state="TX", year=2022)
pop.all <- bind_rows(pop.15, pop.16, pop.17, pop.18, 
                          pop.19, pop.20, pop.21, pop.22)


### diploma having
ed.15 <- acs.fun(variable="B15003_017", state="TX", year=2015)
ed.16 <- acs.fun(variable="B15003_017", state="TX", year=2016)
ed.17 <- acs.fun(variable="B15003_017", state="TX", year=2017)
ed.18 <- acs.fun(variable="B15003_017", state="TX", year=2018)
ed.19 <- acs.fun(variable="B15003_017", state="TX", year=2019)
ed.20 <- acs.fun(variable="B15003_017", state="TX", year=2020)
ed.21 <- acs.fun(variable="B15003_017", state="TX", year=2021)
ed.22 <- acs.fun(variable="B15003_017", state="TX", year=2022)
ed.all <- bind_rows(ed.15, ed.16, ed.17, ed.18, 
                     ed.19, ed.20, ed.21, ed.22)
## insurance
in.15 <- acs.fun(variable="B27001_001", state="TX", year=2015)
in.16 <- acs.fun(variable="B27001_001", state="TX", year=2016)
in.17 <- acs.fun(variable="B27001_001", state="TX", year=2017)
in.18 <- acs.fun(variable="B27001_001", state="TX", year=2018)
in.19 <- acs.fun(variable="B27001_001", state="TX", year=2019)
in.20 <- acs.fun(variable="B27001_001", state="TX", year=2020)
in.21 <- acs.fun(variable="B27001_001", state="TX", year=2021)
in.22 <- acs.fun(variable="B27001_001", state="TX", year=2022)
in.all <- bind_rows(in.15, in.16, in.17, in.18, 
                    in.19, in.20, in.21, in.22)

### percent below poverty line
blw.pov = pov.all$Md.inc/pop.all$Md.inc
House.age = 2024 - hage.all$Md.inc
diploma = ed.all$Md.inc/pop.all$Md.inc

### for joining
tx <- tx %>%
  arrange(Year)

tx$Mn.Inc = Inc.all$Md.inc
tx$Pop.Age = age.all$Md.inc
tx$Hous_Age = hage.all$Md.inc
tx$Mn.House.Val = hInc.all$Md.inc
tx$Pop.Sz = pop.all$Md.inc
tx$Pov.Pct = blw.pov
tx$Pop.Dns = dens$Dens
tx$Imp.Sf = imp$Imp
tx$diploma = diploma



### save! 
setwd('/Users/nathanialodell/Library/CloudStorage/OneDrive-TexasTechUniversity/JRM_Tarrant County Aedes/NODell Tarrant/Modeling (again)')
write_sf(tx, "No Buffers with Var.shp")

tx <- read_sf("No Buffers with Var.shp")
#########################################################
#Zhao, X., Sun, Y., Xu, J. et al. Urban landscape pattern mainly drives the 
#early epidemic distribution of dengue fever in Hangzhou, China. Landsc Ecol 39, 116 (2024)
install.packages('osmdata')
library(osmdata)
### define bb
dfw <- c(-97.53708, 32.49527, -96.61046, 32.99781)
dfwbb <- matrix(dfw, nrow = 2, ncol = 2)
rownames(dfwbb) <- c("x", "y")
colnames(dfwbb) <- c("min", "max")
### list of feature we can grab
available_features()
### tags for features

### testing things.. 
dfw_hospitals <-
  opq(bbox = dfwbb) %>%
  add_osm_feature(key = "amenity", value = "hospital") %>%
  osmdata_sf()

ggplot() +
  geom_sf(data = dfw_hospitals$osm_polygons) #polygons of hospitals... what if we found distance of points to hospitals? 
#something like distance to closest hospital... 
library(ggmap)
plot(st_geometry(tx), col = "red")
plot(st_geometry(dfw_hospitals$osm_points), add = T) #interesting... 

distances <- st_distance(tx, dfw_hospitals$osm_points)
# Extract the minimum distance for each point in 'tx'
min_distances <- apply(distances, 1, min)

# Add these minimum distances as a new column in 'tx'
tx$nearest_hospital_distance <- min_distances/1000 #hmmm... keep going with this!

#roads and density... 
available_tags("highway")
#extract road data from OpenStreetMaps
rd_dfw <- opq(bbox = dfwbb) %>%
  add_osm_feature(key = "highway", value = c("tertiary", "unclassified", 
                                             "residential"))%>% #non-highways
  osmdata_sf()

# only want line data for plotting and density calculations
rd_lines <- rd_dfw$osm_lines
#save!
write_sf(rd_lines%>%dplyr::select(geometry), "DFW Roads.shp")

#load roads line geometry sf
roads <- read_sf("DFW Roads.shp")
# Function to calculate road density for different buffer distances

calculate_road_density <- function(tx_data, roads, dist) {
  
  # Generate centroids for each unique address
  tx_unique <- tx_data %>%
    group_by(Address) %>%
    slice(1) %>%
    ungroup()
  
  centroids <- st_centroid(tx_unique)
  
  # Create buffer around each centroid
  buffer <- st_buffer(centroids, dist = dist)
  
  # Find the roads that intersect the buffer
  roads_intersect <- st_intersection(roads, buffer)
  
  # Calculate the length of the intersecting roads
  roads_intersect$length <- st_length(roads_intersect)
  
  # Aggregate the total road length by address
  roads_aggregate <- aggregate(roads_intersect$length, by = list(roads_intersect$Address), FUN = sum)
  colnames(roads_aggregate) <- c('Address', 'length')
  
  # Calculate the area of the buffer
  buffer$area <- st_area(buffer) / 10000 # Convert to km^2
  
  # Merge the buffer data with the road lengths
  buffers_merge <- merge(buffer, roads_aggregate, by = 'Address', all = TRUE)
  
  # Fill NAs with 0 for lengths
  buffers_merge$length[is.na(buffers_merge$length)] <- 0
  
  # Convert length to kilometers
  buffers_merge$length <- buffers_merge$length / 1000
  
  # Calculate road density (km of road per km^2)
  buffers_merge$road_density <- buffers_merge$length / buffers_merge$area
  
  # Select columns from buffers_merge
  Address  <- buffers_merge$Address
  area <- buffers_merge$area
  length <- buffers_merge$length
  road_density <- buffers_merge$road_density
  
  # Create a new data frame that has measurements
  new_df <- data.frame(Address, area, length, road_density) 
  
  #now we append!
  tx.rd <- merge(tx, new_df[, c("Address", "road_density")], by = "Address", all.x = TRUE)
  
  return(tx.rd)
}
tx.rd.bufs.10 <- calculate_road_density(tx, roads, 100)
tx.rd.bufs.25 <- calculate_road_density(tx, roads, 250)
tx.rd.bufs.50 <- calculate_road_density(tx, roads, 500)

#append
tx.rd.bufs.10$road_density_250 <- tx.rd.bufs.25$road_density 
tx.rd.bufs.10$road_density_500 <- tx.rd.bufs.50$road_density 
tx.rd.bufs.10 <- tx.rd.bufs.10 %>%
  arrange(Year)

setwd('/Users/nathanialodell/Library/CloudStorage/OneDrive-TexasTechUniversity/JRM_Tarrant County Aedes/NODell Tarrant/Modeling (again)')
write_sf(tx.rd.bufs.10, "DFW point with LC CLIM.shp")
#for parks we probably want to use a polygon for intersections. let's make the function more robust as well
tx <- read_sf("DFW point with LC CLIM.shp")%>%
  arrange(Year)
#define bb
dfw <- c(-97.53708, 32.49527, -96.61046, 32.99781)
dfwbb <- matrix(dfw, nrow = 2, ncol = 2)
rownames(dfwbb) <- c("x", "y")
colnames(dfwbb) <- c("min", "max")

calculate_density <- function(tx_data, dist, osmdat) {
  
  var.shape <- osmdat$osm_polygons #store polygon shapefile
  
  # Generate centroids for each unique address
  tx_unique <- tx_data %>%
    group_by(Address) %>%
    slice(1) %>%
    ungroup()
  
  centroids <- st_centroid(tx_unique)
  
  # Create buffer around each centroid
  buffer <- st_buffer(centroids, dist = dist)
  
  # Find the shapes that intersect the buffer
  intersect <- st_intersection(var.shape, buffer)
  
  # Calculate the length of the intersecting shapes
  intersect$shp_area <- st_area(intersect)
  
  # Aggregate the total area by address
  agg <- aggregate(intersect$shp_area, by = list(intersect$Address), FUN = sum)
  colnames(agg) <- c('Address', 'shp_area')
  
  # Calculate the area of the buffer
  buffer$pt_area <- st_area(buffer) / 1000000 # Convert to km^2
  
  # Merge the buffer data with the road lengths
  buffers_merge <- merge(buffer, agg, by = 'Address', all = TRUE)
  
  # Fill NAs with 0 for lengths
  buffers_merge$shp_area[is.na(buffers_merge$shp_area)] <- 0
  
  # Convert length to kilometers
  buffers_merge$shp_area <- buffers_merge$shp_area / 1000000
  
  # Calculate percent density
  buffers_merge$density <- buffers_merge$shp_area / buffers_merge$pt_area 
  
  # Select columns from buffers_merge
  Address  <- buffers_merge$Address
  pt_area <- buffers_merge$pt_area
  shp_area <- buffers_merge$shp_area
  density <- buffers_merge$density
  
  # Create a new data frame that has measurements
  new_df <- data.frame(Address, pt_area, shp_area, density)
  
  #now we append!
  tx.ap <- merge(tx_data, new_df[, c("Address", "density")], by = "Address", all.x = TRUE)%>%
    arrange(Year) #this is how tx is arranged
  
  return(tx.ap)
}
#retrieve features 
parks <-
  opq(bbox = dfwbb) %>%
  add_osm_feature(key = "leisure", value = c("dog_park", "park", "playground")) %>%
  osmdata_sf()

tx.pk.bufs.10 <- calculate_density(tx, 100, parks)
tx.pk.bufs.25 <- calculate_density(tx, 250, parks)
tx.pk.bufs.10 <- calculate_density(tx, 500, parks) #note that these are percentages
#add
colnames(tx.pk.bufs.10)[colnames(tx.pk.bufs.10) == 'density'] <- 'pk_den_100'
tx.pk.bufs.10$pk_den_250 <- tx.pk.bufs.25$density
tx.pk.bufs.10$pk_den_500 <- tx.pk.bufs.50$density

#shopping centers
shops <-
  opq(bbox = dfwbb) %>%
  add_osm_feature(key = "shop") %>%
  osmdata_sf()

tx.shop.10 <- calculate_density(tx, 100, shops)
tx.shop.25 <- calculate_density(tx, 250, shops)
tx.shop.50 <- calculate_density(tx, 500, shops)

tx.pk.bufs.10$shop_den_100 <- tx.shop.10$density
tx.pk.bufs.10$shop_den_250 <- tx.shop.25$density
tx.pk.bufs.10$shop_den_500 <- tx.shop.50$density

#schools
pub.ser <- 
  opq(bbox = dfwbb) %>%
  add_osm_feature(key = "amenity", value = c("school", "library" ))%>%
  osmdata_sf()
tx.pub.10 <- calculate_density(tx, 100, pub.ser)
tx.pub.25 <- calculate_density(tx, 250, pub.ser)
tx.pub.50 <- calculate_density(tx, 500, pub.ser)
tx$pub_den_100 <- tx.pub.10$density
tx$pub_den_250 <- tx.pub.25$density
tx$pub_den_500 <- tx.pub.50$density

# entertainment
en <- 
  opq(bbox = dfwbb) %>%
  add_osm_feature(key = "amenity", value = c("arts_centre", "cinema", "community_centre",
                                             "events_venue", "music_venue",
                                             "social_centre", "theatre"))%>%
  osmdata_sf()
tx.en.10 <- calculate_density(tx, 100, en)
tx.en.25 <- calculate_density(tx, 250, en)
tx.en.50 <- calculate_density(tx, 500, en)
tx$en_den_100 <- tx.en.10$density
tx$en_den_250 <- tx.en.25$density
tx$en_den_500 <- tx.en.50$density
write_sf(tx, "DFW point with LC CLIM.shp")

#landscape
tx <- read_sf("No Buffers with Var.shp")
setwd("/Users/nathanialodell/Library/CloudStorage/OneDrive-TexasTechUniversity/JRM_Tarrant County Aedes")
wrd.pp <- raster("gpw_v4_population_density_adjusted_to_2015_unwpp_country_totals_rev11_2020_30_sec.tif")

txa <- st_transform(tx, crs(wrd.pp))


########## you don't need to run any of this#################
#### creating warp raster for reprojection of MMCD lnd.cvr
e <- extent(txa)

#### creating raster - I think this code is about 350m X 350m  grids
r <- raster(e, ncol=50,nrow=93, crs=crs(wrd.pp), vals=1)
e2 <- as(extent(r), 'SpatialPolygons')

### land cover data for DFW
setwd("/Users/nathanialodell/Library/CloudStorage/OneDrive-TexasTechUniversity/JRM_Tarrant County Aedes/DFW landcover images")
r.lnd.11 <- rast("NLCD_2021_Land_Cover_L48_20230630_Y52lVFP1jJFtNVyhTYKr.tiff")

### changing utm to lat lon...
##exporting
r.lnd.11.2 <- project(r.lnd.11, crs(txa))

#### counting pixels ######

rasterdf <- function(x, aggregate = 1) {
  resampleFactor <- aggregate
  inputRaster <- x
  inCols <- ncol(inputRaster)
  inRows <- nrow(inputRaster)
  # Compute numbers of columns and rows in the new raster for mapping
  resampledRaster <- raster(ncol=(inCols / resampleFactor),
                            nrow=(inRows / resampleFactor))
  # Match to the extent of the original raster
  extent(resampledRaster) <- extent(inputRaster)
  # Resample data on the new raster
  y <- resample(inputRaster,resampledRaster,method='ngb')

  # Extract cell coordinates
  coords <- xyFromCell(y, seq_len(ncell(y)))
  dat <- stack(as.data.frame(getValues(y)))
  # Add names - 'value' for data, 'variable' to indicate different raster layers
  # in a stack
  names(dat) <- c('value', 'variable')
  dat <- cbind(coords, dat)
  dat
}

r.lnd2x <- raster(r.lnd.11.2)
r.lnd.df <- rasterdf(r.lnd2x)
r.lnd.tb <- as_tibble(r.lnd.df)


### pulling out individual values
wtr <- r.lnd.tb %>% filter(value==11)
d1 <- r.lnd.tb %>% filter(value==21)
d2 <- r.lnd.tb %>% filter(value==22)
d3 <- r.lnd.tb %>% filter(value==23)
d4 <- r.lnd.tb %>% filter(value==24)
bar <- r.lnd.tb %>% filter(value==31)
decid <- r.lnd.tb %>% filter(value==41)
conf <- r.lnd.tb %>% filter(value==42)
mix <- r.lnd.tb %>% filter(value==43)
scrub <- r.lnd.tb %>% filter(value==52)
grass <- r.lnd.tb %>% filter(value==71)
pastr <- r.lnd.tb %>% filter(value==81)
crop <- r.lnd.tb %>% filter(value==82)
wetland_f <- r.lnd.tb %>% filter(value==90)
wetland_e <- r.lnd.tb %>% filter(value==95)

### function from stack overflow

pointcount = function(r, pts){
  # make a raster of zeroes like the input
  r2 = r
  r2[] = 0
  # get the cell index for each point and make a table:
  counts = table(cellFromXY(r,pts))
  # fill in the raster with the counts from the cell index:
  r2[as.numeric(names(counts))] = counts
  return(r2)
}

wtr.cnt <- pointcount(r=r, pts=as.data.frame(wtr[,1:2]))
d1.cnt <- pointcount(r=r, pts=as.data.frame(d1[,1:2]))
d2.cnt <- pointcount(r=r, pts=as.data.frame(d2[,1:2]))
d3.cnt <- pointcount(r=r, pts=as.data.frame(d3[,1:2]))
d4.cnt <- pointcount(r=r, pts=as.data.frame(d4[,1:2]))
bar.cnt <- pointcount(r=r, pts=as.data.frame(bar[,1:2]))
decid.cnt <- pointcount(r=r, pts=as.data.frame(decid[,1:2]))
conf.cnt <- pointcount(r=r, pts=as.data.frame(conf[,1:2]))
mix.cnt <- pointcount(r=r, pts=as.data.frame(mix[,1:2]))
scrub.cnt <- pointcount(r=r, pts=as.data.frame(scrub[,1:2]))
grass.cnt <- pointcount(r=r, pts=as.data.frame(grass[,1:2]))
pastr.cnt <- pointcount(r=r, pts=as.data.frame(pastr[,1:2]))
crop.cnt <- pointcount(r=r, pts=as.data.frame(crop[,1:2]))
wet_f.cnt <- pointcount(r=r, pts=as.data.frame(wetland_f[,1:2]))
wet_e.cnt <- pointcount(r=r, pts=as.data.frame(wetland_e[,1:2]))

###### calculating percent

wtr.mat <- as.matrix(wtr.cnt)
d1.mat <- as.matrix(d1.cnt)
d2.mat <- as.matrix(d2.cnt)
d3.mat <- as.matrix(d3.cnt)
d4.mat <- as.matrix(d4.cnt)
bar.mat <- as.matrix(bar.cnt)
decid.mat <- as.matrix(decid.cnt)
conf.mat <- as.matrix(conf.cnt)
mix.mat <- as.matrix(mix.cnt)
scrub.mat <- as.matrix(scrub.cnt)
grass.mat <- as.matrix(grass.cnt)
pastr.mat <- as.matrix(pastr.cnt)
crop.mat <- as.matrix(crop.cnt)
wet_f.mat <- as.matrix(wet_f.cnt)
wet_e.mat <- as.matrix(wet_e.cnt)

tot.mat <- wtr.mat + d1.mat + d2.mat + d3.mat + d4.mat +
  bar.mat + decid.mat + conf.mat + mix.mat + scrub.mat +
  grass.mat + pastr.mat + crop.mat + wet_f.mat + wet_e.mat

### convert to %
wtr.pct <-wtr.mat/tot.mat
d1.pct <-d1.mat/tot.mat
d2.pct <-d2.mat/tot.mat
d3.pct <-d3.mat/tot.mat
d4.pct <-d4.mat/tot.mat
bar.pct <-bar.mat/tot.mat
decid.pct <-decid.mat/tot.mat
conf.pct <-conf.mat/tot.mat
mix.pct <-mix.mat/tot.mat
scrub.pct <-scrub.mat/tot.mat
grass.pct <-grass.mat/tot.mat
pastr.pct <-pastr.mat/tot.mat
crop.pct <- crop.mat/tot.mat
wet_f.pct <- wet_f.mat/tot.mat
wet_e.pct <- wet_e.mat/tot.mat


#### create raster stack of land cover
r.wtr <- raster(e, ncol=50,nrow=93, crs=crs(r.lnd2x), vals=1)
values(r.wtr) <- wtr.pct
r.d1 <- raster(e, ncol=50,nrow=93, crs=crs(r.lnd2x), vals=1)
values(r.d1) <- d1.pct
r.d2 <- raster(e, ncol=50,nrow=93, crs=crs(r.lnd2x), vals=1)
values(r.d2) <- d2.pct
r.d3 <- raster(e, ncol=50,nrow=93, crs=crs(r.lnd2x), vals=1)
values(r.d3) <- d3.pct
r.d4 <- raster(e, ncol=50,nrow=93, crs=crs(r.lnd2x), vals=1)
values(r.d4) <- d4.pct
r.bar <- raster(e, ncol=50,nrow=93, crs=crs(r.lnd2x), vals=1)
values(r.bar) <- bar.pct
r.decid <- raster(e, ncol=50,nrow=93, crs=crs(r.lnd2x), vals=1)
values(r.decid) <- decid.pct
r.conf <- raster(e, ncol=50,nrow=93, crs=crs(r.lnd2x), vals=1)
values(r.conf) <- conf.pct
r.mix <- raster(e, ncol=50,nrow=93, crs=crs(r.lnd2x), vals=1)
values(r.mix) <- mix.pct
r.scrub <- raster(e, ncol=50,nrow=93, crs=crs(r.lnd2x), vals=1)
values(r.scrub) <- scrub.pct
r.grass <- raster(e, ncol=50,nrow=93, crs=crs(r.lnd2x), vals=1)
values(r.grass) <- grass.pct
r.pastr <- raster(e, ncol=50,nrow=93, crs=crs(r.lnd2x), vals=1)
values(r.pastr) <- pastr.pct
r.crop <- raster(e, ncol=50,nrow=93, crs=crs(r.lnd2x), vals=1)
values(r.crop) <- crop.pct
r.wet_f <- raster(e, ncol=50,nrow=93, crs=crs(r.lnd2x), vals=1)
values(r.wet_f) <- wet_f.pct
r.wet_e <- raster(e, ncol=50,nrow=93, crs=crs(r.lnd2x), vals=1)
values(r.wet_e) <- wet_e.pct

dfw.lnd.stck.11 <- stack(r.wtr, r.d1, r.d2, r.d3, r.d4,
                         r.bar, r.decid, r.conf, r.mix,
                         r.scrub, r.grass, r.pastr, r.crop,
                         r.wet_f, r.wet_e)
names(dfw.lnd.stck.11) <- c("Water","Develop_1","Develop_2","Develop_3","Develop_4",
                            "Barren","Deciduous","Coniferous","Mixed","Scrub",
                            "Grass","Pasture","Crop","Wetland_F","Wetland_E")
e2 <- as(extent(r), 'SpatialPolygons')
dfw.lnd2 <- crop(dfw.lnd.stck.11, e2)
dfw.lnd3 <- resample(dfw.lnd2, r, "bilinear")
plot(dfw.lnd3)

setwd("/Users/nathanialodell/Library/CloudStorage/OneDrive-TexasTechUniversity/JRM_Tarrant County Aedes/DFW landcover images")
writeRaster(dfw.lnd3, filename="DFW 2021 percent lnd cover.tif", options="INTERLEAVE=BAND", overwrite=TRUE)
############################################################################################

#### add values to txa file - this works, just need to make sure we get the years aligned with the images
setwd("/Users/nathanialodell/Library/CloudStorage/OneDrive-TexasTechUniversity/JRM_Tarrant County Aedes")
wrd.pp <- raster("gpw_v4_population_density_adjusted_to_2015_unwpp_country_totals_rev11_2020_30_sec.tif")

txa <- st_read("No Buffers with Var.shp")
txa <- st_transform(txa, crs(wrd.pp))

setwd("/Users/nathanialodell/Library/CloudStorage/OneDrive-TexasTechUniversity/JRM_Tarrant County Aedes/DFW landcover images")
lnd.13.15 <- rast("DFW 2013 percent lnd cover.tif")
txa.13.15 <- txa %>% filter(Year<2016)

lnd.cvr.13.15 <- matrix(nrow=nrow(txa.13.15), ncol=15)
for (i in 1:15) {
  dfw <- terra::extract(x=lnd.13.15[[i]], y=txa.13.15, fun="mean", method="bilinear")
  lnd.cvr.13.15[,i] = dfw[,2]
}

txa.13.15 <- txa.13.15 %>%
  mutate(Water = lnd.cvr.13.15[,1],
         Develop_1 = lnd.cvr.13.15[,2],
         Develop_2 = lnd.cvr.13.15[,3],
         Develop_3 = lnd.cvr.13.15[,4],
         Develop_4 = lnd.cvr.13.15[,5],
         Barren = lnd.cvr.13.15[,6],
         Deciduous = lnd.cvr.13.15[,7],
         Coniferous = lnd.cvr.13.15[,8],
         Mixed = lnd.cvr.13.15[,9],
         Scrub = lnd.cvr.13.15[,10],
         Grass = lnd.cvr.13.15[,11],
         Pasture = lnd.cvr.13.15[,12],
         Crop = lnd.cvr.13.15[,13],
         Wetland_F = lnd.cvr.13.15[,14],
         Wetland_E = lnd.cvr.13.15[,15])

lnd.16.18 <- rast("DFW 2016 percent lnd cover.tif")
txa.16.18 <- txa %>% filter(Year > 2015 & Year < 2019)

lnd.cvr.16.18 <- matrix(nrow=nrow(txa.16.18), ncol=15)
for (i in 1:15) {
  dfw <- terra::extract(x=lnd.16.18[[i]], y=txa.16.18, fun="mean", method="bilinear")
  lnd.cvr.16.18[,i] = dfw[,2]
}

txa.16.18 <- txa.16.18 %>%
  mutate(Water = lnd.cvr.16.18[,1],
         Develop_1 = lnd.cvr.16.18[,2],
         Develop_2 = lnd.cvr.16.18[,3],
         Develop_3 = lnd.cvr.16.18[,4],
         Develop_4 = lnd.cvr.16.18[,5],
         Barren = lnd.cvr.16.18[,6],
         Deciduous = lnd.cvr.16.18[,7],
         Coniferous = lnd.cvr.16.18[,8],
         Mixed = lnd.cvr.16.18[,9],
         Scrub = lnd.cvr.16.18[,10],
         Grass = lnd.cvr.16.18[,11],
         Pasture = lnd.cvr.16.18[,12],
         Crop = lnd.cvr.16.18[,13],
         Wetland_F = lnd.cvr.16.18[,14],
         Wetland_E = lnd.cvr.16.18[,15])

lnd.19.20 <- rast("DFW 2019 percent lnd cover.tif")
txa.19.20 <- txa %>% filter(Year > 2018 & Year < 2021)

lnd.cvr.19.20 <- matrix(nrow=nrow(txa.19.20), ncol=15)
for (i in 1:15) {
  dfw <- terra::extract(x=lnd.19.20[[i]], y=txa.19.20, fun="mean", method="bilinear")
  lnd.cvr.19.20[,i] = dfw[,2]
}

txa.19.20 <- txa.19.20 %>%
  mutate(Water = lnd.cvr.19.20[,1],
         Develop_1 = lnd.cvr.19.20[,2],
         Develop_2 = lnd.cvr.19.20[,3],
         Develop_3 = lnd.cvr.19.20[,4],
         Develop_4 = lnd.cvr.19.20[,5],
         Barren = lnd.cvr.19.20[,6],
         Deciduous = lnd.cvr.19.20[,7],
         Coniferous = lnd.cvr.19.20[,8],
         Mixed = lnd.cvr.19.20[,9],
         Scrub = lnd.cvr.19.20[,10],
         Grass = lnd.cvr.19.20[,11],
         Pasture = lnd.cvr.19.20[,12],
         Crop = lnd.cvr.19.20[,13],
         Wetland_F = lnd.cvr.19.20[,14],
         Wetland_E = lnd.cvr.19.20[,15])

lnd.21 <- rast("DFW 2021 percent lnd cover.tif")
txa.21 <- txa %>% filter(Year > 2020)

lnd.cvr.21 <- matrix(nrow=nrow(txa.21), ncol=15)
for (i in 1:15) {
  dfw <- terra::extract(x=lnd.21[[i]], y=txa.21, fun="mean", method="bilinear")
  lnd.cvr.21[,i] = dfw[,2]
}

txa.21 <- txa.21 %>%
  mutate(Water = lnd.cvr.21[,1],
         Develop_1 = lnd.cvr.21[,2],
         Develop_2 = lnd.cvr.21[,3],
         Develop_3 = lnd.cvr.21[,4],
         Develop_4 = lnd.cvr.21[,5],
         Barren = lnd.cvr.21[,6],
         Deciduous = lnd.cvr.21[,7],
         Coniferous = lnd.cvr.21[,8],
         Mixed = lnd.cvr.21[,9],
         Scrub = lnd.cvr.21[,10],
         Grass = lnd.cvr.21[,11],
         Pasture = lnd.cvr.21[,12],
         Crop = lnd.cvr.21[,13],
         Wetland_F = lnd.cvr.21[,14],
         Wetland_E = lnd.cvr.21[,15])

txa.V2 <- rbind(txa.13.15, txa.16.18, txa.19.20, txa.21) %>%
  na.omit() #removes 12 addresses 
# save 
setwd('/Users/nathanialodell/Library/CloudStorage/OneDrive-TexasTechUniversity/JRM_Tarrant County Aedes/NODell Tarrant/Modeling (again)')
write_sf(txa.V2, "DFW point with LC.shp")

tx <- read_sf("DFW point with LC.shp")
#let's go ahead and pivot now
tx <- tx %>%
  pivot_wider(names_from = "SpcsTyp", values_from = "Total")
tx$ALBOPICTUS <- tx$ALBOPICTUS %>%
  replace_na(0)
tx$AEGYPTI <- tx$AEGYPTI %>%
  replace_na(0)
write_sf(tx, "DFW point with LC.shp")


#CLIMATE
#################### CLIMATE ########## 
tx <- read_sf("DFW point with LC.shp")
library(reshape)
library(prism)
library(stars)
#### function of above for weekly averages 
pavg.fun <- function(variable, dir) {
  prism_set_dl_dir(dir)
  get_prism_dailys(
    type = variable,
    minDate = "2015-03-03",
    maxDate = "2022-11-02",
    keepZip = FALSE
  )
  # Stack and subset the data
  stackp <- pd_stack(
    prism_archive_subset(
      type = variable,
      temp_period = "daily",
      minDate = "2015-03-03",
      maxDate = "2022-11-02"
    )
  )
  
  # Define the bounding box and crop the stack
  e <- st_bbox(tx)
  stack_crop <- crop(stackp, e)
  
  # Extract sites
  sites <- unique(tx[, c("Address", "geometry")])
  sites$ID <- 1:nrow(sites) # Numerical labels
  
  # Extract data for the specified sites
  tx_mn <- terra::extract(
    stack_crop, sites
  )
  
  # Convert to long format
  tx_mn_long <- as.data.frame(tx_mn)
  tx_mn_long$site <- rownames(tx_mn_long)
  tx_mn_long <- pivot_longer(tx_mn_long, -site, names_to = "date", values_to = "mn.var")
  tx_mn_long$site <- as.numeric(tx_mn_long$site)
  
  # Merge Address from sites into tx_mn_long
  tx_mn_long <- tx_mn_long %>%
    left_join(sites %>% dplyr::select(ID, Address), by = c("site" = "ID"))
  
  # Get the date and week
  tx_mn_long$date <- pd_get_date(tx_mn_long$date)
  stck.dates <- MMWRweek(tx_mn_long$date)
  tx_mn_long$Epi_wk <- stck.dates$MMWRweek
  tx_mn_long$Year <- stck.dates$MMWRyear
  
  # Remove NA values
  tx_mn_long <- tx_mn_long %>% drop_na(mn.var)
  
  # Calculate the weekly average for each site
  weekly_avg <- tx_mn_long %>%
    group_by(Address, Year, Epi_wk) %>%
    summarize(mn_value = mean(mn.var, na.rm = TRUE))
  
  # Append to the main dataset
  tx_updated <- left_join(tx, weekly_avg %>% dplyr::select(Address, Year, Epi_wk, mn_value), 
                        by = c("Address", "Year", "Epi_wk"))
  
  return(tx_updated)
}
tx_dew <- pavg.fun("tdmean", dir = "/Users/nathanialodell/Library/CloudStorage/OneDrive-TexasTechUniversity/JRM_Tarrant County Aedes/NODell Tarrant/Modeling (again)/Climate Extraction/Dew")
tx_tmean <- pavg.fun("tmean", "/Users/nathanialodell/Library/CloudStorage/OneDrive-TexasTechUniversity/JRM_Tarrant County Aedes/NODell Tarrant/Modeling (again)/Climate Extraction")
tx$mn.dew <- tx_dew$mn_value
tx$tmean <- tx_tmean$mn_value

#####function for weekly mins 
pdailymin.fun <- function(variable, dir) {
  prism_set_dl_dir(dir)
  # Stack and subset the data
  stackp <- pd_stack(
    prism_archive_subset(
      type = variable,
      temp_period = "daily",
      minDate = "2015-03-03",
      maxDate = "2022-11-02"
    )
  )
  
  # Define the bounding box and crop the stack
  e <- st_bbox(tx)
  stack_crop <- crop(stackp, e)
  
  # Extract sites
  sites <- unique(tx[, c("Address", "geometry")])
  sites$ID <- 1:nrow(sites) # Numerical labels
  
  # Extract data for the specified sites
  tx_mn <- terra::extract(
    stack_crop, sites
  )
  
  # Convert to long format
  tx_mn_long <- as.data.frame(tx_mn)
  tx_mn_long$site <- rownames(tx_mn_long)
  tx_mn_long <- pivot_longer(tx_mn_long, -site, names_to = "date", values_to = "mn.var")
  tx_mn_long$site <- as.numeric(tx_mn_long$site)
  
  # Merge Address from sites into tx_mn_long
  tx_mn_long <- tx_mn_long %>%
    left_join(sites %>% dplyr::select(ID, Address), by = c("site" = "ID"))
  
  # Get the date and week
  tx_mn_long$date <- pd_get_date(tx_mn_long$date)
  stck.dates <- MMWRweek(tx_mn_long$date)
  tx_mn_long$Epi_wk <- stck.dates$MMWRweek
  tx_mn_long$Year <- stck.dates$MMWRyear
  
  # Remove NA values
  tx_mn_long <- tx_mn_long %>% drop_na(mn.var)
  
  # Calculate the weekly average for each site
  weekly_avg <- tx_mn_long %>%
    group_by(Address, Year, Epi_wk) %>%
    summarize(mn_value = min(mn.var)) #min or max
  
  # Append to the main dataset
  tx_updated <- left_join(tx, weekly_avg %>% dplyr::select(Address, Year, Epi_wk, mn_value), 
                          by = c("Address", "Year", "Epi_wk"))
  
  return(tx_updated)
}
#####function for weekly max
pdailymax.fun <- function(variable, dir) {
  prism_set_dl_dir(dir)
  # Stack and subset the data
  stackp <- pd_stack(
    prism_archive_subset(
      type = variable,
      temp_period = "daily",
      minDate = "2015-03-03",
      maxDate = "2022-11-02"
    )
  )
  
  # Define the bounding box and crop the stack
  e <- st_bbox(tx)
  stack_crop <- crop(stackp, e)
  
  # Extract sites
  sites <- unique(tx[, c("Address", "geometry")])
  sites$ID <- 1:nrow(sites) # Numerical labels
  
  # Extract data for the specified sites
  tx_mn <- terra::extract(
    stack_crop, sites
  )
  
  # Convert to long format
  tx_mn_long <- as.data.frame(tx_mn)
  tx_mn_long$site <- rownames(tx_mn_long)
  tx_mn_long <- pivot_longer(tx_mn_long, -site, names_to = "date", values_to = "mn.var")
  tx_mn_long$site <- as.numeric(tx_mn_long$site)
  
  # Merge Address from sites into tx_mn_long
  tx_mn_long <- tx_mn_long %>%
    left_join(sites %>% dplyr::select(ID, Address), by = c("site" = "ID"))
  
  # Get the date and week
  tx_mn_long$date <- pd_get_date(tx_mn_long$date)
  stck.dates <- MMWRweek(tx_mn_long$date)
  tx_mn_long$Epi_wk <- stck.dates$MMWRweek
  tx_mn_long$Year <- stck.dates$MMWRyear
  
  # Remove NA values
  tx_mn_long <- tx_mn_long %>% drop_na(mn.var)
  
  # Calculate the weekly average for each site
  weekly_avg <- tx_mn_long %>%
    group_by(Address, Year, Epi_wk) %>%
    summarize(mn_value = max(mn.var)) #min or max
  
  # Append to the main dataset
  tx_updated <- left_join(tx, weekly_avg %>% dplyr::select(Address, Year, Epi_wk, mn_value), 
                          by = c("Address", "Year", "Epi_wk"))
  
  return(tx_updated)
}
tx_tmin <- pdailymin.fun("tmin", "/Users/nathanialodell/Library/CloudStorage/OneDrive-TexasTechUniversity/JRM_Tarrant County Aedes/NODell Tarrant/Modeling (again)/Climate Extraction/tminmax")
tx_vpdmin <- pdailymin.fun("vpdmin", "/Users/nathanialodell/Library/CloudStorage/OneDrive-TexasTechUniversity/JRM_Tarrant County Aedes/NODell Tarrant/Modeling (again)/Climate Extraction/VPD")
tx_tmax <- pdailymax.fun("tmax", "/Users/nathanialodell/Library/CloudStorage/OneDrive-TexasTechUniversity/JRM_Tarrant County Aedes/NODell Tarrant/Modeling (again)/Climate Extraction/tminmax")
tx_vpdmax <- pdailymax.fun("vpdmax", "/Users/nathanialodell/Library/CloudStorage/OneDrive-TexasTechUniversity/JRM_Tarrant County Aedes/NODell Tarrant/Modeling (again)/Climate Extraction/VPD")
tx$tmin <- tx_tmin$mn_value
tx$vpdmin <- tx_vpdmin$mn_value
tx$tmax <- tx_tmax$mn_value
tx$vpdmax <- tx_vpdmax$mn_value


####function for lags... incorporate it with the first function at some point... 
plag.fun <- function(variable, dir, n) {
  prism_set_dl_dir(dir)
  # Stack and subset the data
  stackp <- pd_stack(
    prism_archive_subset(
      type = variable,
      temp_period = "daily",
      minDate = "2015-03-03",
      maxDate = "2022-11-02"
    )
  )
  
  # Define the bounding box and crop the stack
  e <- st_bbox(tx)
  stack_crop <- crop(stackp, e)
  
  # Extract sites
  sites <- unique(tx[, c("Address", "geometry")])
  sites$ID <- 1:nrow(sites) # Numerical labels
  
  # Extract data for the specified sites
  tx_mn <- terra::extract(
    stack_crop, sites
  )
  
  # Convert to long format
  tx_mn_long <- as.data.frame(tx_mn)
  tx_mn_long$site <- rownames(tx_mn_long)
  tx_mn_long <- pivot_longer(tx_mn_long, -site, names_to = "date", values_to = "mn.var")
  tx_mn_long$site <- as.numeric(tx_mn_long$site)
  
  # Merge Address from sites into tx_mn_long
  tx_mn_long <- tx_mn_long %>%
    left_join(sites %>% dplyr::select(ID, Address), by = c("site" = "ID"))
  
  # Get the date and week
  tx_mn_long$date <- pd_get_date(tx_mn_long$date)
  stck.dates <- MMWRweek(tx_mn_long$date)
  tx_mn_long$Epi_wk <- stck.dates$MMWRweek
  tx_mn_long$Year <- stck.dates$MMWRyear
  
  # Remove NA values
  tx_mn_long <- tx_mn_long %>% drop_na(mn.var)
  
  # Calculate the weekly average for each site
  weekly_avg <- tx_mn_long %>%
    group_by(Address, Year, Epi_wk) %>%
    summarize(mn_value = mean(mn.var, na.rm = TRUE))
  
  weekly_avg_lag <- weekly_avg %>%
    group_by(Address)%>%
    mutate(mn_value_lag = lag(mn_value, n))
  # Append to the main dataset
  tx_updated <- left_join(tx, weekly_avg_lag %>% dplyr::select(Address, Year, Epi_wk, mn_value_lag), 
                          by = c("Address", "Year", "Epi_wk"))
  
  return(tx_updated)
}
minlag.fun <- function(variable, dir, n) {
  prism_set_dl_dir(dir)
  # Stack and subset the data
  stackp <- pd_stack(
    prism_archive_subset(
      type = variable,
      temp_period = "daily",
      minDate = "2015-03-03",
      maxDate = "2022-11-02"
    )
  )
  
  # Define the bounding box and crop the stack
  e <- st_bbox(tx)
  stack_crop <- crop(stackp, e)
  
  # Extract sites
  sites <- unique(tx[, c("Address", "geometry")])
  sites$ID <- 1:nrow(sites) # Numerical labels
  
  # Extract data for the specified sites
  tx_mn <- terra::extract(
    stack_crop, sites
  )
  
  # Convert to long format
  tx_mn_long <- as.data.frame(tx_mn)
  tx_mn_long$site <- rownames(tx_mn_long)
  tx_mn_long <- pivot_longer(tx_mn_long, -site, names_to = "date", values_to = "mn.var")
  tx_mn_long$site <- as.numeric(tx_mn_long$site)
  
  # Merge Address from sites into tx_mn_long
  tx_mn_long <- tx_mn_long %>%
    left_join(sites %>% dplyr::select(ID, Address), by = c("site" = "ID"))
  
  # Get the date and week
  tx_mn_long$date <- pd_get_date(tx_mn_long$date)
  stck.dates <- MMWRweek(tx_mn_long$date)
  tx_mn_long$Epi_wk <- stck.dates$MMWRweek
  tx_mn_long$Year <- stck.dates$MMWRyear
  
  # Remove NA values
  tx_mn_long <- tx_mn_long %>% drop_na(mn.var)
  
  # Calculate the weekly average for each site
  weekly_avg <- tx_mn_long %>%
    group_by(Address, Year, Epi_wk) %>%
    summarize(mn_value = min(mn.var))
  
  weekly_avg_lag <- weekly_avg %>%
    group_by(Address)%>%
    mutate(mn_value_lag = lag(mn_value, n))
  # Append to the main dataset
  tx_updated <- left_join(tx, weekly_avg_lag %>% dplyr::select(Address, Year, Epi_wk, mn_value_lag), 
                          by = c("Address", "Year", "Epi_wk"))
  
  return(tx_updated)
}
maxlag.fun <- function(variable, dir, n) {
  prism_set_dl_dir(dir)
  # Stack and subset the data
  stackp <- pd_stack(
    prism_archive_subset(
      type = variable,
      temp_period = "daily",
      minDate = "2015-03-03",
      maxDate = "2022-11-02"
    )
  )
  
  # Define the bounding box and crop the stack
  e <- st_bbox(tx)
  stack_crop <- crop(stackp, e)
  
  # Extract sites
  sites <- unique(tx[, c("Address", "geometry")])
  sites$ID <- 1:nrow(sites) # Numerical labels
  
  # Extract data for the specified sites
  tx_mn <- terra::extract(
    stack_crop, sites
  )
  
  # Convert to long format
  tx_mn_long <- as.data.frame(tx_mn)
  tx_mn_long$site <- rownames(tx_mn_long)
  tx_mn_long <- pivot_longer(tx_mn_long, -site, names_to = "date", values_to = "mn.var")
  tx_mn_long$site <- as.numeric(tx_mn_long$site)
  
  # Merge Address from sites into tx_mn_long
  tx_mn_long <- tx_mn_long %>%
    left_join(sites %>% dplyr::select(ID, Address), by = c("site" = "ID"))
  
  # Get the date and week
  tx_mn_long$date <- pd_get_date(tx_mn_long$date)
  stck.dates <- MMWRweek(tx_mn_long$date)
  tx_mn_long$Epi_wk <- stck.dates$MMWRweek
  tx_mn_long$Year <- stck.dates$MMWRyear
  
  # Remove NA values
  tx_mn_long <- tx_mn_long %>% drop_na(mn.var)
  
  # Calculate the weekly average for each site
  weekly_avg <- tx_mn_long %>%
    group_by(Address, Year, Epi_wk) %>%
    summarize(mn_value = max(mn.var))
  
  weekly_avg_lag <- weekly_avg %>%
    group_by(Address)%>%
    mutate(mn_value_lag = lag(mn_value, n))
  # Append to the main dataset
  tx_updated <- left_join(tx, weekly_avg_lag %>% dplyr::select(Address, Year, Epi_wk, mn_value_lag), 
                          by = c("Address", "Year", "Epi_wk"))
  
  return(tx_updated)
}
tx_updated <- plag.fun("ppt", dir = "/Users/nathanialodell/Library/CloudStorage/OneDrive-TexasTechUniversity/JRM_Tarrant County Aedes/NODell Tarrant/Modeling (again)/Climate Extraction/ppt", 1)
tx_pp2 <- plag.fun("ppt", dir = "/Users/nathanialodell/Library/CloudStorage/OneDrive-TexasTechUniversity/JRM_Tarrant County Aedes/NODell Tarrant/Modeling (again)/Climate Extraction/ppt", 2)
txd2 <- plag.fun("tdmean", dir = "/Users/nathanialodell/Library/CloudStorage/OneDrive-TexasTechUniversity/JRM_Tarrant County Aedes/NODell Tarrant/Modeling (again)/Climate Extraction/Dew", 2)
vpdmi1 <- minlag.fun("vpdmin", "/Users/nathanialodell/Library/CloudStorage/OneDrive-TexasTechUniversity/JRM_Tarrant County Aedes/NODell Tarrant/Modeling (again)/Climate Extraction/VPD", 1)
vpdmi2 <- minlag.fun("vpdmin", "/Users/nathanialodell/Library/CloudStorage/OneDrive-TexasTechUniversity/JRM_Tarrant County Aedes/NODell Tarrant/Modeling (again)/Climate Extraction/VPD", 2)
vpdma1 <- maxlag.fun("vpdmax", "/Users/nathanialodell/Library/CloudStorage/OneDrive-TexasTechUniversity/JRM_Tarrant County Aedes/NODell Tarrant/Modeling (again)/Climate Extraction/VPD", 1)
vpdma2 <- maxlag.fun("vpdmax", "/Users/nathanialodell/Library/CloudStorage/OneDrive-TexasTechUniversity/JRM_Tarrant County Aedes/NODell Tarrant/Modeling (again)/Climate Extraction/VPD", 2)
tmi1 <- minlag.fun("tmin", "/Users/nathanialodell/Library/CloudStorage/OneDrive-TexasTechUniversity/JRM_Tarrant County Aedes/NODell Tarrant/Modeling (again)/Climate Extraction/tminmax", 1)
tmi2 <- minlag.fun("tmin", "/Users/nathanialodell/Library/CloudStorage/OneDrive-TexasTechUniversity/JRM_Tarrant County Aedes/NODell Tarrant/Modeling (again)/Climate Extraction/tminmax", 2)
tma1 <- maxlag.fun("tmax", "/Users/nathanialodell/Library/CloudStorage/OneDrive-TexasTechUniversity/JRM_Tarrant County Aedes/NODell Tarrant/Modeling (again)/Climate Extraction/tminmax", 1)
tma2 <- maxlag.fun("tmax", "/Users/nathanialodell/Library/CloudStorage/OneDrive-TexasTechUniversity/JRM_Tarrant County Aedes/NODell Tarrant/Modeling (again)/Climate Extraction/tminmax", 2)
tx$pptL1 <- tx_updated$mn_value_lag
tx$mn_dewL1 <- tx_updated$mn_value_lag
tx$pptL2 <- tx_pp2$mn_value_lag
tx$mn_dewL2 <- txd2$mn_value_lag
tx$vpdminL1 <- vpdmi1$mn_value_lag
tx$vpdminL2 <- vpdmi2$mn_value_lag
tx$vpdmaxL1 <- vpdma1$mn_value_lag
tx$vpdmaxL2 <- vpdma2$mn_value_lag
tx$tminL1 <- tmi1$mn_value_lag
tx$tminL2 <- tmi2$mn_value_lag
tx$tmaxL1 <- tma1$mn_value_lag
tx$tmaxL2 <- tma2$mn_value_lag

getwd()
write_sf(tx, "DFW point with LC CLIM.shp")

### Meteorological data ###
######################## Meteorological data ####
library(openmeteo)
setwd('/Users/nathanialodell/Library/CloudStorage/OneDrive-TexasTechUniversity/JRM_Tarrant County Aedes/NODell Tarrant/Modeling (again)')
tx1 <- read_sf("DFW point with LC CLIM.shp") %>%
  st_drop_geometry()

## function for retrieving variables
meteo_daily_fun <- function(x, var, start, end) {
  
  # Extract unique combinations of lat and lon
  unique_locations <- unique(x[c("lat", "lon")])
  
  # Initialize empty dataframe to store results
  all_data <- data.frame()
  
  # Loop through each unique location
  for (i in 1:nrow(unique_locations)) {
    # Retrieve location coordinates
    location <- as.numeric(unique_locations[i, ])
    
    # Retrieve daily meteorological data for the location and specified variables
    dum <- weather_history(
      location = location,
      start = start,
      end = end,
      daily = var
    )
    colnames(dum)[2] <- "value" # value is stored in second column; change the name, 
    # otherwise it's dynamic to variable being retrieved
    
    # Add lat and lon columns to 'dum'
    dum$lat <- location[1]
    dum$lon <- location[2]
    
    # Append 'dum' to 'all_data'
    all_data <- rbind(all_data, dum)
  }
  
  # Calculate MMWR week and year for each date
  dates <- MMWRweek(all_data$date)
  all_data$Epi_wk <- dates$MMWRweek
  all_data$Year <- dates$MMWRyear
  
  # Calculate the weekly average for each site
  weekly_avg <- all_data %>%
    group_by(lon, lat, Year, Epi_wk) %>%
    summarize(mn_value = mean(value, na.rm = TRUE), .groups = 'drop') %>%
    ungroup() %>%
    group_by(lon, lat) %>%
    mutate(mn_value_L1 = lag(mn_value, 1),
           mn_value_L2 = lag(mn_value, 2)) %>% # for certain variables (windspeed) the lag should actually be MORE important than present values
    ungroup()
  
  return(weekly_avg)
}
meteo_hour_fun <- function(x, var, start, end) {
  
  # Extract unique combinations of lat and lon
  unique_locations <- unique(x[c("lat", "lon")])
  
  # Initialize empty dataframe to store results
  all_data <- data.frame()
  
  # Loop through each unique location
  for (i in 1:nrow(unique_locations)) {
    # Retrieve location coordinates
    location <- as.numeric(unique_locations[i, ])
    
    # Retrieve daily meteorological data for the location and specified variables
    dum <- weather_history(
      location = location,
      start = start,
      end = end,
      hourly = var # instead of daily
    ) %>%
      na.omit()
    colnames(dum)[2] <- "value" # value is stored in second column; change the name, 
    # otherwise it's dynamic to variable being retrieved
    
    ## the alteration to the daily fun is this piece right here
    dum$datetime <- as.Date(dum$datetime)
    
    # Add lat and lon columns to 'dum'
    dum$lat <- location[1]
    dum$lon <- location[2]
    
    # Append 'dum' to 'all_data'
    all_data <- rbind(all_data, dum)
  }
  
  # Calculate MMWR week and year for each date
  dates <- MMWRweek(all_data$datetime)
  all_data$Epi_wk <- dates$MMWRweek
  all_data$Year <- dates$MMWRyear
  
  # Calculate the weekly average for each site
  weekly_avg <- all_data %>%
    group_by(lon, lat, Year, Epi_wk) %>%
    summarize(mn_value = mean(value, na.rm = TRUE), .groups = 'drop') %>%
    ungroup() %>%
    group_by(lon, lat) %>%
    mutate(mn_value_L1 = lag(mn_value, 1),
           mn_value_L2 = lag(mn_value, 2)) %>%
    ungroup()
  
  return(weekly_avg)
}

test <- meteo_daily_fun(tx, "windspeed_10m_max", "2015-09-01", "2015-09-30") # looks like it works... 
test <- meteo_hour_fun(tx, "cloudcover", "2022-03-03", "2022-03-04") # works as well!

# so we want daily values of windspeed, and then hourly of every thing else 
# NOTE: we have to do yearly, otherwise it overloads the API
wnd.tx.15 <- meteo_daily_fun(tx, "windspeed_10m_max", "2015-03-02", "2015-12-25")
wnd.tx.16 <- meteo_daily_fun(tx, "windspeed_10m_max", "2016-03-02", "2016-12-25")
wnd.tx.17 <- meteo_daily_fun(tx, "windspeed_10m_max", "2017-03-02", "2017-12-25")
wnd.tx.18 <- meteo_daily_fun(tx, "windspeed_10m_max", "2018-03-02", "2018-12-25")
wnd.tx.19 <- meteo_daily_fun(tx, "windspeed_10m_max", "2019-03-02", "2019-12-25")
wnd.tx.20 <- meteo_daily_fun(tx, "windspeed_10m_max", "2020-03-02", "2020-12-25")
wnd.tx.21 <- meteo_daily_fun(tx, "windspeed_10m_max", "2021-03-02", "2021-12-25")
wnd.tx.22 <- meteo_daily_fun(tx, "windspeed_10m_max", "2022-03-02", "2022-12-25")
weekly_avg <- rbind(wnd.tx.15,wnd.tx.16,wnd.tx.17,wnd.tx.18,wnd.tx.19,wnd.tx.20,
                    wnd.tx.21,wnd.tx.22)
tx_updated <- left_join(tx, weekly_avg %>% dplyr::select(lat, lon, Year, Epi_wk, mn_value, mn_value_L1,
                                                         mn_value_L2), 
                        by = c("lat", "lon", "Year", "Epi_wk"))

# Cloudcover
cc.15 <- meteo_hour_fun(tx, "cloudcover", "2015-03-02", "2015-12-25")
cc.16 <- meteo_hour_fun(tx, "cloudcover", "2016-03-02", "2016-12-25")
cc.17 <- meteo_hour_fun(tx, "cloudcover", "2017-03-02", "2017-12-25")
cc.18 <- meteo_hour_fun(tx, "cloudcover", "2018-03-02", "2018-12-25")
cc.19 <- meteo_hour_fun(tx, "cloudcover", "2019-03-02", "2019-12-25")
cc.20 <- meteo_hour_fun(tx, "cloudcover", "2020-03-02", "2020-12-25")
cc.21 <- meteo_hour_fun(tx, "cloudcover", "2021-03-02", "2021-12-25")
cc.22 <- meteo_hour_fun(tx, "cloudcover", "2022-03-02", "2022-12-25")
weekly_avg.cc <- rbind(cc.15,cc.16,cc.17,cc.18,cc.19,cc.20,
                    cc.21,cc.22)
tx_updated.cc <- left_join(tx, weekly_avg.cc %>% dplyr::select(lat, lon, Year, Epi_wk, mn_value, mn_value_L1,
                                                               mn_value_L2), 
                        by = c("lat", "lon", "Year", "Epi_wk"))

# atmospheric pressure
ap.15 <- meteo_hour_fun(tx, "pressure_msl", "2015-03-02", "2015-12-25")
ap.16 <- meteo_hour_fun(tx, "pressure_msl", "2016-03-02", "2016-12-25")
ap.17 <- meteo_hour_fun(tx, "pressure_msl", "2017-03-02", "2017-12-25")
ap.18 <- meteo_hour_fun(tx, "pressure_msl", "2018-03-02", "2018-12-25")
ap.19 <- meteo_hour_fun(tx, "pressure_msl", "2019-03-02", "2019-12-25")
ap.20 <- meteo_hour_fun(tx, "pressure_msl", "2020-03-02", "2020-12-25")
ap.21 <- meteo_hour_fun(tx, "pressure_msl", "2021-03-02", "2021-12-25")
ap.22 <- meteo_hour_fun(tx, "pressure_msl", "2022-03-02", "2022-12-25")
weekly_avg.ap <- rbind(ap.15,ap.16,ap.17,ap.18,ap.19,ap.20,
                       ap.21,ap.22)
tx_updated.ap <- left_join(tx, weekly_avg.ap %>% dplyr::select(lat, lon, Year, Epi_wk, mn_value, mn_value_L1,
                                                               mn_value_L2), 
                           by = c("lat", "lon", "Year", "Epi_wk"))

# surface soil temp
so.15 <- meteo_hour_fun(tx, "soil_temperature_0_to_7cm", "2015-03-02", "2015-12-25")
so.16 <- meteo_hour_fun(tx, "soil_temperature_0_to_7cm", "2016-03-02", "2016-12-25")
so.17 <- meteo_hour_fun(tx, "soil_temperature_0_to_7cm", "2017-03-02", "2017-12-25")
so.18 <- meteo_hour_fun(tx, "soil_temperature_0_to_7cm", "2018-03-02", "2018-12-25")
so.19 <- meteo_hour_fun(tx, "soil_temperature_0_to_7cm", "2019-03-02", "2019-12-25")
so.20 <- meteo_hour_fun(tx, "soil_temperature_0_to_7cm", "2020-03-02", "2020-12-25")
so.21 <- meteo_hour_fun(tx, "soil_temperature_0_to_7cm", "2021-03-02", "2021-12-25")
so.22 <- meteo_hour_fun(tx, "soil_temperature_0_to_7cm", "2022-03-02", "2022-12-25")
weekly_avg.so <- rbind(so.15,so.16,so.17,so.18,so.19,so.20,
                       so.21,so.22)
tx_updated.so <- left_join(tx, weekly_avg.so %>% dplyr::select(lat, lon, Year, Epi_wk, mn_value, mn_value_L1,
                                                               mn_value_L2), 
                           by = c("lat", "lon", "Year", "Epi_wk"))

# surface soil moisture
som.15 <- meteo_hour_fun(tx, "soil_moisture_0_to_7cm", "2015-03-02", "2015-12-25")
som.16 <- meteo_hour_fun(tx, "soil_moisture_0_to_7cm", "2016-03-02", "2016-12-25")
som.17 <- meteo_hour_fun(tx, "soil_moisture_0_to_7cm", "2017-03-02", "2017-12-25")
som.18 <- meteo_hour_fun(tx, "soil_moisture_0_to_7cm", "2018-03-02", "2018-12-25")
som.19 <- meteo_hour_fun(tx, "soil_moisture_0_to_7cm", "2019-03-02", "2019-12-25")
som.20 <- meteo_hour_fun(tx, "soil_moisture_0_to_7cm", "2020-03-02", "2020-12-25")
som.21 <- meteo_hour_fun(tx, "soil_moisture_0_to_7cm", "2021-03-02", "2021-12-25")
som.22 <- meteo_hour_fun(tx, "soil_moisture_0_to_7cm", "2022-03-02", "2022-12-25")
weekly_avg.som <- rbind(som.15,som.16,som.17,som.18,som.19,som.20,
                       som.21,som.22)
tx_updated.som <- left_join(tx, weekly_avg.som %>% dplyr::select(lat, lon, Year, Epi_wk, mn_value, mn_value_L1,
                                                                 mn_value_L2), 
                           by = c("lat", "lon", "Year", "Epi_wk"))

# relative humidity
hum.15 <- meteo_hour_fun(tx, "relative_humidity_2m", "2015-03-02", "2015-12-25")
hum.16 <- meteo_hour_fun(tx, "relative_humidity_2m", "2016-03-02", "2016-12-25")
hum.17 <- meteo_hour_fun(tx, "relative_humidity_2m", "2017-03-02", "2017-12-25")
hum.18 <- meteo_hour_fun(tx, "relative_humidity_2m", "2018-03-02", "2018-12-25")
hum.19 <- meteo_hour_fun(tx, "relative_humidity_2m", "2019-03-02", "2019-12-25")
hum.20 <- meteo_hour_fun(tx, "relative_humidity_2m", "2020-03-02", "2020-12-25")
hum.21 <- meteo_hour_fun(tx, "relative_humidity_2m", "2021-03-02", "2021-12-25")
hum.22 <- meteo_hour_fun(tx, "relative_humidity_2m", "2022-03-02", "2022-12-25")
weekly_avg.hum <- rbind(hum.15,hum.16,hum.17,hum.18,hum.19,hum.20,
                        hum.21,hum.22)
tx_updated.hum <- left_join(tx, weekly_avg.hum %>% dplyr::select(lat, lon, Year, Epi_wk, mn_value, mn_value_L1,
                                                                 mn_value_L2), 
                            by = c("lat", "lon", "Year", "Epi_wk"))

# irradiance (Brownsville paper)
irr.15 <- meteo_hour_fun(tx, "direct_normal_irradiance", "2015-03-02", "2015-12-25")
irr.16 <- meteo_hour_fun(tx, "direct_normal_irradiance", "2016-03-02", "2016-12-25")
irr.17 <- meteo_hour_fun(tx, "direct_normal_irradiance", "2017-03-02", "2017-12-25")
irr.18 <- meteo_hour_fun(tx, "direct_normal_irradiance", "2018-03-02", "2018-12-25")
irr.19 <- meteo_hour_fun(tx, "direct_normal_irradiance", "2019-03-02", "2019-12-25")
irr.20 <- meteo_hour_fun(tx, "direct_normal_irradiance", "2020-03-02", "2020-12-25")
irr.21 <- meteo_hour_fun(tx, "direct_normal_irradiance", "2021-03-02", "2021-12-25")
irr.22 <- meteo_hour_fun(tx, "direct_normal_irradiance", "2022-03-02", "2022-12-25")
weekly_avg.irr <- rbind(irr.15,irr.16,irr.17,irr.18,irr.19,irr.20,
                        irr.21,irr.22)
tx_updated.irr <- left_join(tx, weekly_avg.irr %>% dplyr::select(lat, lon, Year, Epi_wk, mn_value, mn_value_L1,
                                                                 mn_value_L2), 
                            by = c("lat", "lon", "Year", "Epi_wk"))

# bind all these guys
tx$mn.wnd.spd <- tx_updated$mn_value
tx$mn.wnd.spd.L1 <- tx_updated$mn_value_L1
tx$mn.wnd.spd.L2 <- tx_updated$mn_value_L2
tx$cld.cvr <- tx_updated.cc$mn_value
tx$cld.cvr.L1 <- tx_updated.cc$mn_value_L1
tx$cld.cvr.L2 <- tx_updated.cc$mn_value_L2
tx$atmos.ps <- tx_updated.ap$mn_value
tx$atmo.ps.L1 <- tx_updated.ap$mn_value_L1
tx$atmo.ps.L2 <- tx_updated.ap$mn_value_L2
tx$soil.tmp <- tx_updated.so$mn_value
tx$soil.tmp.L1 <- tx_updated.so$mn_value_L1
tx$soil.tmp.L2 <- tx_updated.so$mn_value_L2
tx$soil.moi <- tx_updated.som$mn_value
tx$soil.moi.L1 <- tx_updated.som$mn_value_L1
tx$soil.moi.L2 <- tx_updated.som$mn_value_L2
tx$rel.hum <- tx_updated.hum$mn_value
tx$rel.hum.L1 <- tx_updated.hum$mn_value_L1
tx$rel.hum.L2 <- tx_updated.hum$mn_value_L2
tx$rad <- tx_updated.irr$mn_value
tx$rad.L1 <- tx_updated.irr$mn_value_L1
tx$rad.L2 <- tx_updated.irr$mn_value_L2

tx <- tx %>%
  st_as_sf(coords = c("lon", "lat"))

write_sf(tx, "DFW point with LC CLIM METEO.shp")
