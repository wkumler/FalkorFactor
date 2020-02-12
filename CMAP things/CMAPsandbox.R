
# Setup things ----
library(cmap4r)
library(dplyr)
library(plotly)

#set_authorization()

catalog <- get_catalog() %>%
  select(Variable, Table_Name, Unit, Sensor, Unit)

# Modeled chlorophyll ----
# From Pisces model
get_head("tblPisces_NRT")
as.data.frame(get_var_coverage(tableName = "tblPisces_NRT", varName = "CHL"))
as.data.frame(get_var_resolution(tableName = "tblPisces_NRT", varName = "CHL"))
dat <- exec_manualquery("SELECT * FROM tblPisces_NRT
         WHERE
         [time] BETWEEN '2017-06-01' AND '2017-07-01' AND
         lat BETWEEN 0 AND 70 AND
         lon BETWEEN -180 AND -80 AND
         depth BETWEEN 0 AND 0.5
         ORDER BY [time], lat, lon, depth")
plot_ly(data=dat, x=~lon, y=~lat, z=~log10(CHL), type = "heatmap")


# Satellite chlorophyll! ----
# From MODIS Aqua (?)
as.data.frame(get_var_coverage(tableName = "tblCHL_REP", varName = "chl"))
as.data.frame(get_var_resolution(tableName = "tblCHL_REP", varName = "chl"))

dat <- exec_manualquery("SELECT * FROM tblCHL_REP
                        WHERE 
                        [time] BETWEEN '2017-06-18' AND '2017-07-18' AND
                        lat BETWEEN 0 AND 70 AND
                        lon BETWEEN -180 AND -80
                        ORDER BY [time], lat, lon")
plot_ly(data=dat, x=~lon, y=~lat, z=~log10(chl), type = "heatmap")


# ARGO chlorophyll? ----
unname(unlist(get_columns('tblArgoMerge_REP')))
as.data.frame(get_var_coverage(tableName = "tblArgoMerge_REP", 
                               varName = "argo_merge_chl"))
dat <- exec_manualquery("SELECT lat, lon, time, depth,  argo_merge_chl FROM tblArgoMerge_REP
                        WHERE 
                        [time] BETWEEN '2017-06-18' AND '2017-07-18' AND
                        lat BETWEEN 0 AND 70 AND
                        lon BETWEEN -180 AND -80
                        ORDER BY [time], lat, lon")
plot_ly(data=dat, x=~lon, y=~lat, color=~argo_merge_chl, type = "scatter", mode="markers")
