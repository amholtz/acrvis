#' function to create leaflet map for incidence and death

#' @param language - language for figure
#' @param maps_type - notification rate or change in cases
#' @param maps_format - WHO or ECDC


maps_migrationrates <- function(log,
                                gps) {
  
  library(tidyverse, quietly = T, verbose = F, warn.conflicts = F)
  library(sf, quietly = T, verbose = F)
  library(tidytree, warn.conflicts = F, quietly = T)
  library(leaflet)

  
  
  # configuration -----------------------------------------------------------
  
  #install.packages("USAboundariesData", repos = "https://ropensci.r-universe.dev", type = "source")
  
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
  b = 0.1
  
  # Poisson prior mean and offset
  p_mean = log(2)
  p_offset = 5 
  
  # Most recent sampled tip
  mrst = 2005.5
  
  
  # 2. Load files--------------------------------------
  #log = read.table("/Volumes/NGS_Viroscreen/workshop/PHINDAccess/batRABV.state.rates.log", header = T, sep = "\t")
  log = log[ceiling(nrow(log)*b):nrow(log),]
  
  #gps = read.table("/Volumes/NGS_Viroscreen/workshop/PHINDAccess/locationStates.txt", sep = "\t", header = F)
  colnames(gps) = c("location", "lat", "lon")
  locations = gps$location
  n = length(locations)
  
  #map = read_sf("/Volumes/NGS_Viroscreen/workshop/PHINDAccess/gadm40_USA_1.shp")
  #map = map[!map$NAME_1 %in% c("Alaska", "Hawaii", "Rhode Island"),]
  # write_sf(map, "/Volumes/NGS_Viroscreen/workshop/PHINDAccess/usa_states.shp")
  map = read_sf("nc.shp")
  
  # 3. Migration rates---------------------------------
  # Number of rates estimated in the analysis
  if (sum(grepl("indicators", colnames(log))) == n * (n-1) /2) {
    nRates = n * (n-1) /2
  } else {
    nRates = n * (n-1)
  }
  
  out = data.frame()
  indcols = colnames(log)[grepl("indicators", colnames(log))]
  
  for (ri in indcols) {
    
    # Posterior distribution of the migration rate
    r = gsub("indicators", "rates", ri)
    posterior = log[[r]][log[[ri]] == 1]
    
    # BF
    qk = (p_mean + p_offset) / nRates
    prior_or = qk / (1-qk)
    
    pk = mean(log[[ri]])
    if (pk == 1) pk = 1 - 1/nrow(log)
    post_or = pk / ( 1 - pk )
    bf = post_or / prior_or
    
    # Compile all rates
    temp_out = data.frame(
      rate = r,
      BF = bf,
      min = min(posterior),
      q025 = quantile(posterior, 0.025),
      hpd025 = hpd(posterior, 0.95, "lower"),
      median = median(posterior),
      mean = mean(posterior),
      q975 = quantile(posterior, 0.975),
      hpd975 = hpd(posterior, 0.95, "upper"),
      max = max(posterior)
    )
    out = bind_rows(out, temp_out)
  }
  row.names(out) = NULL
  
  # Merge GPS data with summary statistics 
  out$source = combn(locations, 2)[1,]
  out$destination = combn(locations, 2)[2,]
  out = out %>%
    left_join(., gps, by = c("source" = "location")) %>%
    dplyr::rename(lat_source = lat, lon_source = lon) %>%
    left_join(., gps, by = c("destination" = "location")) %>%
    dplyr::rename(lat_dest = lat, lon_dest = lon)
  
  # Create linestring object
  sources <- out %>%
    filter(BF >= 3) %>%
    select(source, lon_source, lat_source) %>%
    st_as_sf(coords = c('lon_source', 'lat_source'), crs = 4326)
  
  destinations <- out %>%
    filter(BF >= 3) %>%
    select(destination, lon_dest, lat_dest) %>%
    st_as_sf(coords = c('lon_dest', 'lat_dest'), crs = 4326)
  
  point_geoms = cbind(sources, destinations)
  mig_rates = st_sfc(mapply(function(a,b){st_cast(st_union(a,b),"LINESTRING")},
                            point_geoms$geometry, point_geoms$geometry.1, SIMPLIFY=FALSE))
  mig_rates = st_sf(NAME = out$rate[out$BF >= 3], BF = out$BF[out$BF >= 3], geom = mig_rates)
  st_crs(mig_rates) = 4326
  
  # Plot migration rates with BF >= 3
  x <- ggplot() +
    geom_sf(data = map) +
    geom_sf(data = mig_rates, aes(col = BF)) +
    geom_point(data = gps, aes(x = lon, y = lat), size = 2) +
    geom_label(data = gps, aes(x = lon, y = lat-1, label = location), fontface = "bold")
  
  
  #Popup_text for leaflet map
  rate_text <- paste(
    "Bayes Factor: ", round(mig_rates$BF,2),"<br/>", 
    sep="") %>%
    lapply(htmltools::HTML)
  
  state_text <- paste(
    "State: ", map$name,"<br/>", 
    sep="") %>%
    lapply(htmltools::HTML)
  
  bins <- c(3, 10, 30, 90, 200, 500, 1000, 5000, Inf)
  pal <- colorBin("YlGnBu", domain = mig_rates$BF, bins = bins)
  
  
  #leaflet map
  
  # migration rates
  map_migration <- leaflet() %>%
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
    addPolylines(data=mig_rates, color = ~pal(BF), weight = 3, opacity =1,
                 label = rate_text) %>% 
    addLegend(data = mig_rates, pal = pal, values = ~BF, opacity = 0.7,
              title = "Bayes Factors > 3", position = "bottomleft", na.label = 'no data') %>% 
    addCircles(lng = gps$lon, lat = gps$lat, weight = 20,  popup = gps$location)
  
               
return(map_migration)
}




