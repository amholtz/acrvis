#' function to create leaflet map from mcc tree

#' @param mcc - mcc tree from BEAST
#' @param date - date of last displayed parent 
#' @param gps - GPS coordinates for loctions
#' @param burn - burn-in
#' @param mean - Poisson prior mean and offset
#' @param offset - Poisson prior mean and offset
#' @param mrst_par - most recent sampled tip


maps_mcc <- function(mcc, date, gps, burn, mean, offset, mrst_par) {
  
  library(tidyverse, quietly = T, verbose = F, warn.conflicts = F)
  library(sf, quietly = T, verbose = F)
  library(tidytree, warn.conflicts = F, quietly = T)
  library(leaflet)
  library(dplyr)
  
  
  
  # configuration -----------------------------------------------------------
  
  # Additional functions
  hpd = function(x, conf, b) {
    
    L = length(x)
    S = sort(x)
    
    # Number of values of x representing conf length of x  
    extent = floor(L * conf)
    
    # If the length of x is too small, the 
    # hpd interval corresponds to x
    if (L==0) {
      return(NaN)
      
    } else if (L == 1) {
      return(S)
      
    } else {
      hpdIndex = 1
      hpdRange = 0
      
      # Test the length of all the intervals 
      # starting from the first value of S
      for (i in 1:(L - extent)) {
        lower = S[i]
        upper = S[i + extent]
        currRange = upper - lower
        
        # Update the hpd interval and 
        # the index of the left-bound value
        if (hpdRange == 0) hpdRange = currRange
        if (currRange < hpdRange) {
          hpdRange = currRange
          hpdIndex = i
        }
      }
      
      if (b == "lower") return(S[hpdIndex])
      if (b == "upper") return(S[hpdIndex + extent])
    }
  }
  
  # 1. Specify parameters------------------------------
  # Burnin
  #b = 0.1
  b = burn
  
  # Poisson prior mean and offset
  #p_mean = log(2)
  #p_offset = 5 
  p_mean = log(mean)
  p_offset = offset 
  
  
  # Most recent sampled tip
  #mrst = 2005.5
  mrst = mrst_par
  
  
  # 2. Load files--------------------------------------
  #log = read.table("/Volumes/NGS_Viroscreen/workshop/PHINDAccess/batRABV.state.rates.log", header = T, sep = "\t")

  #gps = read.table("/Volumes/NGS_Viroscreen/workshop/PHINDAccess/locationStates.txt", sep = "\t", header = F)
  colnames(gps) = c("location", "lat", "lon")
  locations = gps$location
  n = length(locations)
  
  #map = read_sf("/Volumes/NGS_Viroscreen/workshop/PHINDAccess/gadm40_USA_1.shp")
  #map = map[!map$NAME_1 %in% c("Alaska", "Hawaii", "Rhode Island"),]
  # write_sf(map, "/Volumes/NGS_Viroscreen/workshop/PHINDAccess/usa_states.shp")
  map = read_sf("nc.shp")
  
  # 4. MCC tree---------------------------------------
  # Create dataframe with info for all pairs of 
  # parent and child in the mcc tree
  mcc@data$calendar_date = mrst - as.numeric(mcc@data$height)
  edge_df = data.frame(parent = mcc@phylo$edge[,1], child = mcc@phylo$edge[,2])
  
  edge_df$parent_loc = sapply(edge_df$parent, function(x) mcc@data$state[mcc@data$node == x])
  edge_df$parent_lat = sapply(edge_df$parent_loc, function(x) gps$lat[gps$location == x])
  edge_df$parent_lon = sapply(edge_df$parent_loc, function(x) gps$lon[gps$location == x])
  edge_df$parent_date = sapply(edge_df$parent, function(x) mcc@data$calendar_date[mcc@data$node == x])
  
  edge_df$child_loc = sapply(edge_df$child, function(x) mcc@data$state[mcc@data$node == x])
  edge_df$child_lat = sapply(edge_df$child_loc, function(x) gps$lat[gps$location == x])
  edge_df$child_lon = sapply(edge_df$child_loc, function(x) gps$lon[gps$location == x])
  edge_df$child_date = sapply(edge_df$parent_loc, function(x) mcc@data$calendar_date[mcc@data$node == x])
  
  # Remove edge where no migration happened
  edge_df = edge_df %>%
    dplyr::filter(parent_lat != child_lat & parent_lon != child_lon)
  row.names(edge_df) = NULL
  
  # Create multilinestring object
  parent <- edge_df %>%
    select(parent, parent_lon, parent_lat) %>%
    st_as_sf(coords = c('parent_lon', 'parent_lat'), crs = 4326)
  
  child <- edge_df %>%
    select(child, child_lon, child_lat) %>%
    st_as_sf(coords = c('child_lon', 'child_lat'), crs = 4326)
  
  point_geoms = cbind(parent, child)
  mcc_tree = st_sfc(mapply(function(a,b){st_cast(st_union(a,b),"LINESTRING")},
                           point_geoms$geometry, point_geoms$geometry.1, SIMPLIFY=FALSE))
  mcc_tree = st_sf(parent = edge_df$parent, geom = mcc_tree)
  mcc_tree$child = edge_df$child
  mcc_tree$child_date = edge_df$child_date
  mcc_tree$parent_date = edge_df$parent_date
  st_crs(mcc_tree) = 4326
  
  mcc_tree <- mcc_tree %>% filter(parent_date < date)

  
  # Plot MCC tree
  ggplot() +
    geom_sf(data = map) +
    geom_sf(data = mcc_tree) +
    geom_point(data = gps, aes(x = lon, y = lat), size = 2) +
    geom_label(data = gps, aes(x = lon, y = lat-1, label = location), fontface = "bold")
  
  #Popup_text for leaflet map
  date_text <- paste(
    "parent Date: ", round(mcc_tree$parent_date,2),"<br/>", 
    sep="") %>%
    lapply(htmltools::HTML)
  
  state_text <- paste(
    "State: ", map$name,"<br/>", 
    sep="") %>%
    lapply(htmltools::HTML)
  
  bins <- c(1807, 1850, 1960, 1970, 1980, 1990, 2000, 2005)
  pal <- colorBin("YlGnBu", domain = mcc_tree$parent_date, bins = bins)
  
  
  #leaflet map
  
  # migration rates
  map_mcc<- leaflet() %>%
    addProviderTiles(providers$CartoDB.Positron) %>% 
    addPolygons(data = map,
                highlightOptions = highlightOptions(stroke = 2, weight = 3),
                stroke=TRUE, 
                color="grey", 
                weight=1,
                label = state_text,
                labelOptions = labelOptions( 
                  style = list("font-weight" = "normal", padding = "3px 8px"), 
                  textsize = "13px", 
                  direction = "auto"
                )) %>% 
    addPolylines(data=mcc_tree, color = ~pal(parent_date), weight = 3, opacity =1,
                 label = date_text) %>% 
    addLegend(data = mcc_tree, pal = pal, values = ~parent_date, opacity = 0.7,
              title = "Parent Date", position = "bottomleft", na.label = 'no data',
              labFormat = labelFormat(big.mark="")) %>% 
    addCircles(lng = gps$lon, lat = gps$lat, weight = 20,  popup = gps$location)
    
  
  
  return(map_mcc)
}