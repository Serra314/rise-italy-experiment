
# read the parameters file in the input folder
par.raw <- readLines(file('parameters'))
# retrieve information from the parameters file
for(i in 1:length(par.raw)){
  row.i <- par.raw[[i]]
  if(grepl('forecastStart', row.i)){
    start.fore.date <- as.POSIXct(substr(row.i, nchar(row.i) - 18, nchar(row.i)))
  }
  if(grepl('forecastEnd', row.i)){
    end.fore.date <- as.POSIXct(substr(row.i, nchar(row.i) - 18, nchar(row.i)))
  }
  if(grepl('PathCatalogData', row.i)){
    path.catalogue <- substr(row.i, 22, nchar(row.i))
  }
  if(grepl('PathOutputData', row.i)){
    path.fore <- substr(row.i, 21, nchar(row.i))
  }
}


# load functions
source('source/source_inlabru_utilities.R')

# load catalogue
iside.sel.cat <- read.table(file = path.catalogue,
                            sep = ',', header = FALSE, skip = 1) 

names(iside.sel.cat) <- strsplit(readLines(file('input/ISIDE_catalog_selected_PyCSEP'), 
                                           n = 1), ',')[[1]]
# filter catalogue for magnitude and create time_string
iside.sel.cat <- iside.sel.cat %>%
  mutate(time_date = as.POSIXct(gsub('T', ' ', time_string))) %>%
  filter(iside.sel.cat$M > 2.5)
# sort events by time
iside.sel.cat <- iside.sel.cat[order(iside.sel.cat$time_date),]

# take the last 200 points for model fitting and build polygon representing study region
df.500 <- tail(iside.sel.cat, 201)
x.lims <- c(min(iside.sel.cat$lon), max(iside.sel.cat$lon)) + c(-0.01, 0.01)
y.lims <- c(min(iside.sel.cat$lat), max(iside.sel.cat$lat)) + c(-0.01, 0.01)
x_coords <- c(x.lims[1],x.lims[2],x.lims[2],x.lims[1],x.lims[1])
y_coords <- c(y.lims[1],y.lims[1],y.lims[2],y.lims[2],y.lims[1])
poly1 <- sp::Polygon(cbind(x_coords,y_coords))
bdy <- sp::Polygons(list(poly1), ID = "A")
bdy <- sp::SpatialPolygons(list(bdy))
bdy <- as(bdy, 'SpatialPolygonsDataFrame')

# prepare data.frame for inla.bru
start.date = min(df.500$time_date)
# convert it in days from starting event
T1.fore <- as.numeric(difftime(start.fore.date, start.date, units = 'days'))
T2.fore <- as.numeric(difftime(end.fore.date, start.date, units = 'days'))


df.bru <- data.frame(x = df.500$lon, y = df.500$lat, 
                     ts = as.numeric(difftime(df.500$time_date, 
                                              start.date,
                                              units = 'days')),
                     mags = df.500$M)

df.bru <- df.bru[-1,]

# set time limits
T1 = 0
T2 = max(df.bru$ts) + 0.001

# to get the initial values
# bru_options_set(bru_method=list(max_step= 0.5))
# fit_8 <- ETAS.fit(sample.s = df.bru,
#                   N.breaks.min = 30, 
#                   max.length = (T2 - T1)/250,
#                   M0 = 2.5, T1 = T1, T2 = T2, bdy = bdy,
#                   prior.mean = c(0,0,0,0,0,-3),
#                   prior.prec = c(5,3,3,3,3,15),
#                   bru.opt = list(bru_verbose = 3,
#                                  bru_max_iter = 40)) 

# initial values
init_v <- c(-0.3734126, 0.5474037, -0.2237425, -2.0321369, 0.1422050, -5.3822198)
# fit the model
fit_ <- ETAS.fit(sample.s = df.bru,
                  N.breaks.min = 30, 
                  max.length = (T2 - T1)/250,
                  M0 = 2.5, T1 = T1, T2 = T2, bdy = bdy,
                  prior.mean = c(0,0,0,0,0,-3),
                  prior.prec = c(5,3,3,3,3,15),
                  bru.opt = list(bru_verbose = 3,
                                 bru_max_iter = 20,
                                 bru_initial = list(th1 = init_v[1],
                                                    th2 = init_v[2],
                                                    th3 = init_v[3],
                                                    th4 = init_v[4],
                                                    th5 = init_v[5],
                                                    th6 = init_v[6])))


# sample 100000 from the posterior of the parameters
sample.param <- foreach(i = 1:100, .combine = cbind) %do% {
  generate(fit_, data.frame(x = 0, y = 0, ts = 0, mags = 0), 
                         ~ c(th1, th2, th3, th4, th5, th6), 1000)}

# simulate a catalogue for each posterior sample 
dailyforecast <- foreach(i = 1:ncol(sample.param), .combine = rbind) %do% {
  Sigma.p <- matrix(c(exp(sample.param[6,i]), 0, 0, exp(sample.param[6,i])), 
                    byrow = TRUE, ncol = 2)
  
  ss <- sample.ETAS(th.p = sample.param[1:5, i], 
                    beta.p = 2.3, M0 = 2.5, T1 = T1.fore, T2 = T2.fore, 
                    bdy = bdy, Sigma = Sigma.p, Ht = df.bru)
  ss <- bind_rows(ss)
  ss <- ss[order(ss$ts),]
  if(nrow(ss) == 0){
    data.frame(ts = NA, x = NA, y = NA, mags = NA, gen = NA, catalog_id = i, event_id = 1)
  }
  else{
    ss$catalog_id = i - 1
    ss$event_id = 1:nrow(ss)
    ss
  }
}


## cut for magnitude in [4, 9]
dailyforecast <- dailyforecast[dailyforecast$mags >= 4 & dailyforecast$mags <= 9, ]

# set up collection area as sf object
coll.area <- read.table(file = 'input/collection.area.dat')
coll.area.poly <- create_poly(coll.area)
coll.area.sf <- st_as_sf(coll.area.poly)
st_crs(coll.area.sf) <- "+proj=utm +zone=33 +datum=WGS84"
coll.area.sf <- st_transform(coll.area.sf, "+proj=utm +zone=33 +datum=WGS84")

# retain only simulated catalogues with at least one event
dailyforecast.points <- dailyforecast[!is.na(dailyforecast$x),]

# transform point in sf object
coordinates(dailyforecast.points) <- c('x', 'y')
dailyforecast.sf <- st_as_sf(dailyforecast.points)
st_crs(dailyforecast.sf) <- "+proj=utm +zone=33 +datum=WGS84"
dailyforecast.sf <- st_transform(dailyforecast.sf, "+proj=utm +zone=33 +datum=WGS84")
# retain points inside collection area
dailyforecast.sf.filt <- st_intersection(dailyforecast.sf, coll.area.sf)

# set up time_string for the output
output.time.string <- start.date + (60*60*24*dailyforecast.sf.filt$ts)
output.time.string.char <- format(output.time.string, 
                                  format = '%Y-%m-%d-%dT%H:%M:%S',
                                  tz = 'GMT')

start.date.char <- format(start.date, 
                          format = '%Y-%m-%d-%dT%H:%M:%S',
                          tz = 'GMT')

# set up output
output <- data.frame(longitude = st_coordinates(dailyforecast.sf.filt)[,1],
                     latitude = st_coordinates(dailyforecast.sf.filt)[,2],
                     magnitude = dailyforecast.sf.filt$mags,
                     time_string = output.time.string.char,
                     depth = NA, 
                     catalog_id = dailyforecast.sf.filt$catalog_id,
                     event_id = dailyforecast.sf.filt$event_id)

# write output
write.csv(output, paste0(path.fore,'/ETAS_Inlabru.', start.date.char, '.csv'),
          row.names = FALSE)

