

# Setup things ----

library(tidyverse)
library(ggmap)
register_google(key = "AIzaSyBVNJV0Myoy3M50IBBDIjXZhUyU8Mnagvs")
scope_uname <- "scope"
scope_pw <- "SCOPE2014"

stations <- c(62, 64, 77, 80)

# Grab data and format nicely ----

dataScraper <- function(filename){
  url <- paste0("ftp://ftp.soest.hawaii.edu/dkarl/scope/ctd/fk180310/", filename)
  read_table(url, skip=3, 
             col_names = c("Depth", 
                           "Temp", 
                           "Salinity", 
                           "O2", 
                           "Chl", 
                           "PAR", 
                           "N1", "N2", "N3"),
             col_types = cols(
               Depth = col_double(),
               Temp = col_double(),
               Salinity = col_double(),
               O2 = col_double(),
               Chl = col_double(),
               PAR = col_double(),
               N1 = col_double(),
               N2 = col_double(),
               N3 = col_double()
             )) %>%
    select(c("Depth", "Temp", "Salinity", "O2", "Chl", "PAR"))
}

dataMaker <- function(station){
  raw_data_dn <- dataScraper(paste0("fk180310s", station, "c1dn.ctd")) %>% 
    mutate("Dir"="Down")
  raw_data_up <- dataScraper(paste0("fk180310s", station, "c1up.ctd")) %>%
    mutate("Dir"="Up")
  rbind(raw_data_up, raw_data_dn) %>% mutate("Station"=station)
}

gpxing <- function(coordinate){
  v <- strsplit(coordinate, " ")
  return(as.numeric(v[[1]][1])+as.numeric(v[[1]][length(v[[1]])])/60)
}



# Visualize cruise track ----

url <- paste0("http://", scope_uname, ":", scope_pw,
              "@scope.soest.hawaii.edu/collaborators/datainventory/Data/SCOPEcore/FK_CTDsummary_current.txt")
station_locations <- read_table(url, skip = 2, col_types = cols(
  Stn = col_double(), Cast = col_double(),  Month = col_character(),
  Date = col_character(),  Time = col_time(format = ""),
  Latitude = col_character(),  NorthSouth = col_character(),
  Longitude = col_character(),  EastWest = col_character(),
  Depth = col_double(),  Bottles = col_double()),
  col_names = c("Stn", "Cast", "Month", "Date", "Time", "Latitude", "NorthSouth",
                "Longitude", "EastWest", "Depth", "Bottles")) %>%
  slice(-1) %>%
  mutate(DecimalLat=sapply(Latitude, gpxing, USE.NAMES = F)) %>%
  mutate(DecimalLon=sapply(Longitude, gpxing, USE.NAMES = F)*-1) %>%
  mutate(StnName=paste("Station", Stn))




gyre1_map <- get_map(location = c(-156.5, 22.5), zoom = 9, maptype = "satellite", source = "google")
gyre1 <- ggmap(gyre1_map) +
  geom_point(data = station_locations, aes(x=DecimalLon, y=DecimalLat), col=rgb(1,0,0), size=2)+
  geom_path(data = station_locations, aes(x=DecimalLon, y=DecimalLat), col=rgb(1,0,0), lwd=1) +
  geom_point(data = filter(station_locations, Stn%in%stations), 
             aes(x=DecimalLon, y=DecimalLat), col=rgb(0,1,0), size=3, pch=17)
gyre2_map <- get_map(location = c(-160.7, 24.5), zoom = 10, maptype = "satellite", source = "google")
gyre2 <- ggmap(gyre2_map) +
  geom_point(data = station_locations, aes(x=DecimalLon, y=DecimalLat), col=rgb(1,0,0), size=2)+
  geom_path(data = station_locations, aes(x=DecimalLon, y=DecimalLat), col=rgb(1,0,0), lwd=1) +
  geom_point(data = filter(station_locations, Stn%in%stations), 
             aes(x=DecimalLon, y=DecimalLat), col=rgb(0,1,0), size=3, pch=17)


bbox_center <- c(mean(station_locations$DecimalLon), 
                 mean(station_locations$DecimalLat)) + c(-2, 2)
#bbox_center <- unlist(geocode("Hawaii"))+c(-3,5)
hawaii_map <- get_map(location = bbox_center, zoom = 6, maptype = "satellite", source = "google")
ggmap(hawaii_map) + 
  geom_point(data = station_locations, aes(x=DecimalLon, y=DecimalLat), col=rgb(1,0,0), size=2)+
  geom_path(data = station_locations, aes(x=DecimalLon, y=DecimalLat), col=rgb(1,0,0), lwd=1) +
  geom_point(data = filter(station_locations, Stn%in%stations), 
             aes(x=DecimalLon, y=DecimalLat), col=rgb(0,1,0), size=3, pch=17) +
  inset(ggplotGrob(gyre2), xmin = -166, xmax = -160, ymin = 25, ymax = 31) +
  inset(ggplotGrob(gyre1), xmin = -159, xmax = -152.5, ymin = 24, ymax = 31)

# Interactive plots ----
g <- c(scope = 'world', showframe = T, 
       showland = T, landcolor = rgb(1,1,1),
       showcoastlines = T, coastlinecolor = rgb(0,0,1),
       showocean=T, oceancolor=rgb(0.8,0.8,1),
       resolution = 50, projection = list(type = 'Mercator'),
       list(lonaxis = list(range = c(-162, -154))),
       list(lataxis = list(range = c(19, 25))),
       list(domain = list(x = c(0, 1), y = c(0, 1))))

plot_geo(station_locations, mode="markers",
         lat = ~DecimalLat, lon = ~DecimalLon) %>%
  add_paths(x=~DecimalLat, y=~DecimalLon, color=I("grey50"),
            showlegend=F) %>%
  add_markers(x=~DecimalLat, y=~DecimalLon,
              color=~Stn, hoverinfo="text", text=~StnName,
              showlegend=F) %>%
  layout(geo = g)

ocean_surface <- matrix(500, ncol = 2, nrow = 2)
colnames(ocean_surface) <- c(min(station_locations$DecimalLat), max(station_locations$DecimalLat))
rownames(ocean_surface) <- c(min(station_locations$DecimalLon), max(station_locations$DecimalLon))
plot_ly() %>%
  add_trace(data = station_locations, 
            x=~DecimalLat, y=~DecimalLon, z=~Depth*-1,
            type = "scatter3d", mode="markers",
            showlegend=F,
            hoverinfo="text",
            text=~StnName,
            color=~Stn) %>%
  add_surface(z = ocean_surface,
              x = c(min(station_locations$DecimalLat), max(station_locations$DecimalLat)),
              y = c(min(station_locations$DecimalLon), max(station_locations$DecimalLon)),
              color=0,
              colorscale=list(c(0,1), c("blue", "blue")),
              opacity = 0.5,
              showscale=F,
              hoverinfo="none") %>%
  layout(scene = list(xaxis = list(title = "Latitude", showspikes=FALSE),
                      yaxis = list(title = "Longitude", showspikes=FALSE),
                      zaxis = list(title = "Depth", showspikes=FALSE)))



# Visualize my stations ----

my_stations <- station_locations %>%
  filter(Stn%in%stations) %>%
  mutate(TrueDate=format(as.Date(paste(Month, Date), format = "%b %d %Y"), "%Y/%m/%d")) %>%
  mutate(URL=paste0("https://earth.nullschool.net/#",
                    TrueDate,
                    "/0000Z/ocean/surface/currents/overlay=significant_wave_height/orthographic=",
                    DecimalLon, ",",
                    DecimalLat, "/loc=",
                    DecimalLon, ",",
                    DecimalLat))
#sapply(my_stations$URL, browseURL)
#zoom in all the way, Screen2Gif ~50 frames & save as S[station number].gif"




# Visualize my stations CTD data ----

station_data <- do.call(rbind, lapply(stations, dataMaker))
station_data$Chl <- ifelse(c(0, diff(station_data$Chl))<sd(station_data$Chl), station_data$Chl, NA)

temp_plot <- ggplot(station_data) + 
  geom_line(aes(x=Temp, y=Depth, group=Dir, color=Temp), lwd=1) + 
  facet_grid(rows = ~Station) +
  scale_y_reverse() +
  scale_color_continuous(low = rgb(0,0,1), high = rgb(1,0,0)) +
  theme_bw()

chl_plot <- ggplot(station_data) + 
  geom_path(aes(x=Chl, y=Depth, group=Dir, color=Chl), lwd=1) + 
  facet_grid(rows = ~Station) +
  scale_y_reverse() +
  scale_color_continuous(low = rgb(0,0,0.3), high = rgb(0,1,0)) +
  theme_bw()

sal_plot <- ggplot(station_data) + 
  geom_path(aes(x=Salinity, y=Depth, group=Dir, color=Salinity), lwd=1) + 
  facet_grid(rows = ~Station) +
  scale_y_reverse() +
  scale_color_continuous(low = rgb(0.1,0.1,0.1), high = rgb(0.9,0.9,0.9)) +
  theme_bw()

o2_plot <- ggplot(station_data) + 
  geom_path(aes(x=O2, y=Depth, group=Dir, color=O2), lwd=1) + 
  facet_grid(rows = ~Station) +
  scale_y_reverse() +
  scale_color_continuous(low = rgb(0,0,1), high = rgb(1,0,0)) +
  theme_bw()


grid.arrange(temp_plot, sal_plot, o2_plot, chl_plot, nrow = 4)
