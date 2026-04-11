library(sf)
library(ggplot2)
library(ggspatial)
library(tidyverse)
library(data.table)
library(patchwork)



load_catchments = function() {
  # load ----
  ara_shp = sf::read_sf("data/spatial/230214_ARA_BAG/230214_ARA_BAG.shp") %>% 
    dplyr::select(ara_id=ara_id,geometry) %>% 
    dplyr::filter(ara_id!="100000", ara_id!="296300")
  
  
  
  return(ara_shp)
  
  
}

get_plz_areas = function(file_path = 'data/spatial/plz/PLZO_SHP_LV95/PLZO_PLZ.shp', crs_required=NULL) {
  # read file
  plz_area = st_read(file_path)
  
  # set to required crs
  if(!is.null(crs_required)){
    plz_area = st_transform(plz_area, crs = crs_required)
  }
  
  #merge split polygons
  plz_area = rmapshaper::ms_dissolve(plz_area,'PLZ')
  
  
  #return shapes in sf object
  plz_area
  
}

load_plz_centroids = function(file_path = 'data/spatial/plz/PLZO_SHP_LV95/PLZO_PLZ.shp', crs_required=NULL){
  
  #load plz areas
  plz_area = get_plz_areas(file_path, crs_required)
  
  #calculate centroids
  plz_pos = st_centroid(plz_area)
  
  #return centroids in sf object
  plz_pos
}



