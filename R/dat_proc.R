library(tidyverse)
library(tbeptools)
library(sf)
library(raster)

source('R/funcs.R')

# dem data ----------------------------------------------------------------

# utm <- '+proj=utm +zone=17 +datum=NAD83 +units=m +no_defs'
# 
# # https://www.ngdc.noaa.gov/mgg/bathymetry/estuarine/
# dem <- raster('~/Desktop/TBEP/tampa_bay_G070_2017.nc')
# dem <- readAll(dem)
# save(dem, file = 'data/dem.RData', compress = 'xz')

# transect data -----------------------------------------------------------

# import entire transect dataset as JSON
transect <- read_transect(training = FALSE) %>% 
  dplyr::select(-Crew, -MonitoringAgency) %>% 
  arrange(desc(Date))

# get transect species occurrence summaries
transectocc <- anlz_transectocc(transect) %>% 
  arrange(desc(Date))

save(transect, file = 'data/transect.RData', compress = 'xz')
save(transectocc, file = 'data/transectocc.RData', compress = 'xz')

# dem depth at each point -------------------------------------------------

load(file = 'data/transect.RData')
load(file = 'data/dem.RData')

utm <- '+proj=utm +zone=17 +datum=NAD83 +units=m +no_defs'
spp <- c('Halodule', 'Syringodium', 'Thalassia', 'Ruppia', 'Halophila')

# subset transct, get pa for five major seagrass speies
trn <- transect %>% 
  dplyr::filter(var %in% 'Abundance') %>% 
  dplyr::filter(Savspecies %in% !!c(spp, 'No Cover')) %>% 
  group_by(Transect, Site, Date, Depth_obs = Depth) %>% 
  summarise(
    pa = sum(aveval, na.rm = T), 
    meanabu = mean(aveval, na.rm = T),
    .groups = 'drop'
    ) %>% 
  mutate(
    pa = ifelse(pa > 0, 1, 0),
    patxt = factor(pa, levels = c('1', '0'), labels = c('present', 'absent'))
    )

# get transect angle
# location of starting point in UTM
lns <- trnlns %>% 
  st_transform(crs = utm) %>% 
  group_by(Site) %>% 
  nest() %>% 
  mutate(
    LONG_M = purrr::map(data, function(x) st_coordinates(x)[1, 1]),
    LAT_M = purrr::map(data, function(x) st_coordinates(x)[1, 2]), 
    bearing = purrr::map(data, function(x) x$bearing)
  ) %>% 
  dplyr::select(-data) %>% 
  unnest(c('LONG_M', 'LAT_M', 'bearing')) %>% 
  ungroup() %>% 
  dplyr::rename(Transect = Site) 

# get location of transect points
# extract depth from dem using locations
# get location of points by angle and distance
transectdem <- trn %>% 
  inner_join(., lns, by = 'Transect') %>% 
  dplyr::select(Transect, Date, Site, Depth_obs, meanabu, pa, patxt, LAT_Mstr = LAT_M, LONG_Mstr = LONG_M, bearing) %>% 
  mutate(
    Site = as.numeric(Site),
    LONG_M = Site * sin(bearing * pi / 180),
    LAT_M = LONG_M / tan(bearing * pi / 180), 
    LONG_M = LONG_M + LONG_Mstr,
    LAT_M = LAT_M + LAT_Mstr
  ) %>% 
  st_as_sf(coords = c('LONG_M', 'LAT_M'), crs = utm) %>% 
  st_transform(crs = 4326) %>% 
  mutate(
    Depth_dem = raster::extract(dem, .), 
    Depth_dem = 100 * Depth_dem, 
    Date = lubridate::floor_date(Date, unit = 'month')
  )

save(transectdem, file = 'data/transectdem.RData', compress = 'xz')
