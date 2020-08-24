# install necessary packages
if(!require(raster)){install.packages('raster'); library(raster)}
if(!require(mgcv)){install.packages('mgcv'); library(mgcv)}

# function to predict across the entire raster 
predfun <- function(model, data) {
  v <- predict(model, data, se.fit=TRUE)
  cbind(p=as.vector(v$fit), se=as.vector(v$se.fit))
}

# specify paths in and out 
path.in = '../../data/eswatini_inference'
path.out = '../../data/figs/fig_4/'

# load the rasters 
load(file = paste(path.in, '/rasters.RData', sep = ''))

# load the nodes and extract coordinates
nodes = read.csv(list.files(path = path.in, pattern = 'nodes.csv', full.names = T))
coords = data.frame(lon = nodes$Lon, lat = nodes$Lat)

lon = nodes$Lon
lat = nodes$Lat

# extract list of mean Rc values for each node under each data type assumption
files <- list.files(path = path.in, pattern = 'networks_processed.RData', recursive = T, full.names = T)
files <- files[!grepl(pattern = 'replicate', files)]

Rc = list()
prob.imported = list()
for(ii in 1:length(files))
{
  load(files[ii])
  Rc[[ii]] = rowSums(network.data[[1]]$edgeProbs[-1,-1])
  prob.imported[[ii]] = network.data[[1]]$edgeProbs[1,-1]
}

# generate data frame with each row denoting a node and each column denoting a covariate or output
r.points = rasterToPoints(pop, spatial=TRUE)
r.lon = coordinates(r.points)[,1]
r.lat = coordinates(r.points)[,2]
r.df = data.frame(r.lon, r.lat)
names(r.df) <- c('lon', 'lat')

Rc.1 = log(Rc[[1]] + 1e-5)
Rc.2 = log(Rc[[2]] + 1e-5)
Rc.3 = log(Rc[[3]] + 1e-5)
Rc.4 = log(Rc[[4]] + 1e-5)
Rc.5 = log(Rc[[5]] + 1e-5)

# transform data to ensure that 1's and 0's are pushed inward slightly. This allows us to regress on the importation probability.
# taken from Smithson and Verkuilen (2006)
imported.1 = (prob.imported[[1]] * (length(prob.imported[[1]]) - 1) + 0.5)  / length(prob.imported[[1]])
imported.2 = (prob.imported[[2]] * (length(prob.imported[[2]]) - 1) + 0.5)  / length(prob.imported[[2]])
imported.3 = (prob.imported[[3]] * (length(prob.imported[[3]]) - 1) + 0.5)  / length(prob.imported[[3]])
imported.4 = (prob.imported[[4]] * (length(prob.imported[[4]]) - 1) + 0.5)  / length(prob.imported[[4]])
imported.5 = (prob.imported[[5]] * (length(prob.imported[[5]]) - 1) + 0.5)  / length(prob.imported[[5]])

# fit GP to each of the 5 outputs corresponding to the different data type assumptions
gam.Rc.1 = gam(Rc.1 ~ s(lon, lat, bs = 'gp', m = 1), family = 'gaussian')
gam.Rc.2 = gam(Rc.2 ~ s(lon, lat, bs = 'gp', m = 1), family = 'gaussian')
gam.Rc.3 = gam(Rc.3 ~ s(lon, lat, bs = 'gp', m = 1), family = 'gaussian')
gam.Rc.4 = gam(Rc.4 ~ s(lon, lat, bs = 'gp', m = 1), family = 'gaussian')
gam.Rc.5 = gam(Rc.5 ~ s(lon, lat, bs = 'gp', m = 1), family = 'gaussian')

gam.imported.1 = gam(imported.1 ~ s(lon, lat, bs = 'gp', m = 1), family = betar(link = 'logit'))
gam.imported.2 = gam(imported.2 ~ s(lon, lat, bs = 'gp', m = 1), family = betar(link = 'logit'))
gam.imported.3 = gam(imported.3 ~ s(lon, lat, bs = 'gp', m = 1), family = betar(link = 'logit'))
gam.imported.4 = gam(imported.4 ~ s(lon, lat, bs = 'gp', m = 1), family = betar(link = 'logit'))
gam.imported.5 = gam(imported.5 ~ s(lon, lat, bs = 'gp', m = 1), family = betar(link = 'logit'))

# generate predictions of Rc 
Rc.coord.pred.1 <- predict(gam.Rc.1, r.df, se.fit = T, type = 'response')
Rc.coord.pred.1 <- exp(Rc.coord.pred.1[[1]] + ((Rc.coord.pred.1[[2]])^2)/2) - 1e-5
Rc.coord.pred.1 = data.frame(r.lon, r.lat, Rc.coord.pred.1)
Rc.coord.pred.1 <- rasterFromXYZ(Rc.coord.pred.1)

Rc.coord.pred.2 <- predict(gam.Rc.2, r.df, se.fit = T)
Rc.coord.pred.2 <- exp(Rc.coord.pred.2[[1]] + ((Rc.coord.pred.2[[2]])^2)/2) - 1e-5
Rc.coord.pred.2 = data.frame(r.lon, r.lat, Rc.coord.pred.2)
Rc.coord.pred.2 <- rasterFromXYZ(Rc.coord.pred.2)

Rc.coord.pred.3 <- predict(gam.Rc.3, r.df, se.fit = T)
Rc.coord.pred.3 <- exp(Rc.coord.pred.3[[1]] + ((Rc.coord.pred.3[[2]])^2)/2) - 1e-5
Rc.coord.pred.3 = data.frame(r.lon, r.lat, Rc.coord.pred.3)
Rc.coord.pred.3 <- rasterFromXYZ(Rc.coord.pred.3)

Rc.coord.pred.4 <- predict(gam.Rc.4, r.df, se.fit = T)
Rc.coord.pred.4 <- exp(Rc.coord.pred.4[[1]] + ((Rc.coord.pred.4[[2]])^2)/2) - 1e-5
Rc.coord.pred.4 = data.frame(r.lon, r.lat, Rc.coord.pred.4)
Rc.coord.pred.4 <- rasterFromXYZ(Rc.coord.pred.4)

Rc.coord.pred.5 <- predict(gam.Rc.5, r.df, se.fit = T)
Rc.coord.pred.5 <- exp(Rc.coord.pred.5[[1]] + ((Rc.coord.pred.5[[2]])^2)/2) - 1e-5
Rc.coord.pred.5 = data.frame(r.lon, r.lat, Rc.coord.pred.5)
Rc.coord.pred.5 <- rasterFromXYZ(Rc.coord.pred.5)

imported.coord.pred.1 <- predict(gam.imported.1, r.df)
imported.coord.pred.1 <- exp(imported.coord.pred.1) / (1+ exp(imported.coord.pred.1))
imported.coord.pred.1 <- (length(prob.imported[[1]]) * imported.coord.pred.1 - 0.5) / length(prob.imported[[1]])
imported.coord.pred.1 <- data.frame(r.lon, r.lat, imported.coord.pred.1)
imported.coord.pred.1 <- rasterFromXYZ(imported.coord.pred.1)

imported.coord.pred.2 <- predict(gam.imported.2, r.df)
imported.coord.pred.2 <- exp(imported.coord.pred.2) / (1+ exp(imported.coord.pred.2))
imported.coord.pred.2 <- (length(prob.imported[[2]]) * imported.coord.pred.2 - 0.5) / length(prob.imported[[2]])
imported.coord.pred.2 <- data.frame(r.lon, r.lat, imported.coord.pred.2)
imported.coord.pred.2 <- rasterFromXYZ(imported.coord.pred.2)

imported.coord.pred.3 <- predict(gam.imported.3, r.df)
imported.coord.pred.3 <- exp(imported.coord.pred.3) / (1+ exp(imported.coord.pred.3))
imported.coord.pred.3 <- (length(prob.imported[[3]]) * imported.coord.pred.3 - 0.5) / length(prob.imported[[3]])
imported.coord.pred.3 <- data.frame(r.lon, r.lat, imported.coord.pred.3)
imported.coord.pred.3 <- rasterFromXYZ(imported.coord.pred.3)

imported.coord.pred.4 <- predict(gam.imported.4, r.df)
imported.coord.pred.4 <- exp(imported.coord.pred.4) / (1+ exp(imported.coord.pred.4))
imported.coord.pred.4 <- (length(prob.imported[[4]]) * imported.coord.pred.4 - 0.5) / length(prob.imported[[4]])
imported.coord.pred.4 <- data.frame(r.lon, r.lat, imported.coord.pred.4)
imported.coord.pred.4 <- rasterFromXYZ(imported.coord.pred.4)

imported.coord.pred.5 <- predict(gam.imported.5, r.df)
imported.coord.pred.5 <- exp(imported.coord.pred.5) / (1+ exp(imported.coord.pred.5))
imported.coord.pred.5 <- (length(prob.imported[[5]]) * imported.coord.pred.5 - 0.5) / length(prob.imported[[5]])
imported.coord.pred.5 <- data.frame(r.lon, r.lat, imported.coord.pred.5)
imported.coord.pred.5 <- rasterFromXYZ(imported.coord.pred.5)

# save predicted rasters
save(Rc.coord.pred.1, imported.coord.pred.1,
     Rc.coord.pred.2, imported.coord.pred.2,
     Rc.coord.pred.3, imported.coord.pred.3,
     Rc.coord.pred.4, imported.coord.pred.4,
     Rc.coord.pred.5, imported.coord.pred.5,
     file = paste(path.out, 'fig4.RData', sep = ''))
