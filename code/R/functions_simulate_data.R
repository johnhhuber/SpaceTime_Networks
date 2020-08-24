# load necessary packages
if(!require(igraph)){install.packages('igraph'); library(igraph)}
if(!require(maptools)){install.packages('maptools'); library(maptools)}
if(!require(pomp)){install.packages('pomp'); library(pomp)}
if(!require(spatstat)){install.packages('spatstat');library(spatstat)}

# load necessary data
load('../../data/network_sim/STS.RData')
load('../../data/network_sim/Huber_MalariaJournal/Fig2.RData')
load('../../data/network_sim/WorldPopSwazi.RData')

map.swazi = readRDS('../../data/network_sim/SWZ_adm0.rds')
polygon.swazi = SpatialPolygons(map.swazi@polygons, proj4string = map.swazi@proj4string)

# specify the maximum length of the IDP's
max.IDP.S = 44
max.IDP.A = 365

# load the spatial variances that are used in the inference framework 
var.STS = read.csv('../src/variancesSTS.csv', header = F)
var.STA = read.csv('../src/variancesSTA.csv', header = F)
var.ATS = read.csv('../src/variancesATS.csv', header = F)
var.ATA = read.csv('../src/variancesATA.csv', header = F)
var.AUS = read.csv('../src/variancesAUS.csv', header = F)
var.AUA = read.csv('../src/variancesAUA.csv', header = F)

# function to compute spatial variance for STS
variance.STS = function(SI)
{
  if((SI + length(IDP.S)) > dim(Pr_GI.SI)[2] | (SI + length(IDP.S)) < 0)
  {
    return(0)
  }
  else
  {
    sum((GI.vec) * Pr_GI.SI[GI.vec, SI + length(IDP.S)])
  }
}

variance.spatial = function(SI, type)
{
  if(type == 'STS')
  {
    index = SI + max.IDP.S 
    if(index < 1 | index > ncol(var.STS))
    {
      variance = 0
    }
    else 
    {
      variance  = var.STS[1, index]
    }
  }
  if(type == 'STA')
  {
    index = SI + max.IDP.S
    if(index < 1 | index > ncol(var.STA))
    {
      variance = 0
    }
    else 
    {
      variance = var.STA[1, index]
    }
  }
  if(type == 'ATS')
  {
    index = SI + max.IDP.A
    if(index < 1 | index > ncol(var.ATS))
    {
      variance = 0
    }
    else 
    {
      variance = var.ATS[1, index]
    }
  }
  if(type == 'ATA')
  {
    index = SI + max.IDP.A
    if(index < 1 | index > ncol(var.ATA))
    {
      variance = 0
    }
    else 
    {
      variance = var.ATA[1, index]
    }
  }
  if(type == 'AUS')
  {
    index = SI + max.IDP.A
    if(index < 1 | index > ncol(var.AUS))
    {
      variance = 0
    }
    else 
    {
      variance = var.AUS[1, index]
    }
  }
  if(type == 'AUA')
  {
    index = SI + max.IDP.A
    if(index < 1 | index > ncol(var.AUA))
    {
      variance = 0
    }
    else 
    {
      variance = var.AUA[1, index]
    }
  }
  return(variance)
}

loglikSpace = function(node1, node2)
{
  SI = dat$Time[node2] - dat$Time[node1]
  sd = sqrt(10 * variance.STS(SI))
  
  dLat = (dat$Lat[node1] - dat$Lat[node2]) * pi / 180
  dLon = (dat$Lon[node1] - dat$Lon[node2]) * pi / 180
  
  dx = dLat * 6371
  dy = dLon * 6371 * cos(pi * dat$Lat[node1] / 180)
  print(c(dx, dy))
  
  loglik = dnorm(dx, sd = sd, log = TRUE) + dnorm(dy, sd = sd, log = TRUE) 
  return(loglik)
}



# function to generate a single point on a disc with 90% the radius of Swaziland
point.poisson = function(r)
{
  rate = 1 / (pi * r ^ 2)
  ps = rpoispp(rate, win = disc(radius = r))
  while(ps$n != 1)
  {
   #ps = rpoispp(rate, win = disc(radius = 0.9 * sqrt(17364.0 / pi)))
    ps = rpoispp(rate, nsim = n, win = disc(radius = r))
  }
  return(c(ps$x,ps$y))
}

# function to generate a point process of according to a Matern clustering algorithm 
points.clustered <- function(n, frac.cluster, clusteredness = 0.25, r = sqrt(pi * 17360 * 0.258))
{
  # find points on unit disc 
  n.clusters <- n * frac.cluster
  mean.distance <- clusteredness / sqrt(n.clusters)
  mean.points <- n / n.clusters
  r.unit <- sqrt(1 / pi)
  
  ps = rMatClust(n.clusters, mean.distance, mean.points, win = disc(radius = r.unit))
  while(ps$n != n)
  {
    ps = rMatClust(n.clusters, mean.distance, mean.points, win = disc(radius = r.unit))
  } # end while
  
  # transform points to disc of specified radius 
  r.original <- sqrt(ps$x ^ 2 + ps$y ^ 2)
  theta <- base::atan2(ps$y, ps$x)
  ps$x <- r.original * (r / r.unit) * cos(theta)
  ps$y <- r.original * (r / r.unit) * sin(theta)
  
  return(return(rbind(ps$x, ps$y)))
}

# function to generate a single point by sampling from the population density of swaziland
point.swazi = function()
{
  num.coordinates = dim(coordinates.swazi)[1]
  index.sampled = sample(1:num.coordinates, size = 1, prob = popdens.swazi, replace = TRUE)
  return(as.numeric(c(coordinates.swazi[index.sampled,1], coordinates.swazi[index.sampled,2])))
}

# function to convert founder points to lat lon coordinates
convert.coordinates = function(coordinate, ref.lat = 0, ref.lon = 0)
{
  EARTH.RADIUS = 6371.0
  dLat = coordinate[1] / EARTH.RADIUS
  dLon = coordinate[2] / (EARTH.RADIUS * cos(pi * ref.lat / 180))
  
  lat = ref.lat - (180.0 / pi) * dLat
  lon = ref.lon - (180.0 / pi) * dLon
  
  return(c(lat,lon))
}

# function to simulate spatial coordinates of offspring
locs.offspring = function(loc.parent,number.offspring,sigma.offspring){
  x.offspring = rnorm(number.offspring, loc.parent[1], sigma.offspring)
  y.offspring = rnorm(number.offspring, loc.parent[2], sigma.offspring)
  return(cbind(x.offspring,y.offspring))
}



# function to write a set of transmission networks according to a parameter sweep
write.sweep = function(
  reps = 10,
  num.mcmcs = 1,
  diffusion.coef = c(-0.5679933,-0.5679933),
  founder.prob = c(0.574654, 0.574654),
  max_founder_date = c(1,365.0*10.0),
  radius = c(sqrt(1/pi), 10 * sqrt(17364.0 / pi)),
  clustered = TRUE,
  frac.cluster = c(0,1),
  prob.symp = c(1.0,1.0),
  data.dir)
{
  scalars = sobolDesign(lower = c(diffusion.coef = diffusion.coef[1],
                                  founder.prob = founder.prob[1],
                                  max_founder_date = max_founder_date[1],
                                  radius = radius[1], 
                                  prob.symp = prob.symp[1],
                                  frac.cluster = frac.cluster[1]),
                        upper = c(diffusion.coef = diffusion.coef[2],
                                  founder.prob = founder.prob[2],
                                  max_founder_date = max_founder_date[2],
                                  radius = radius[2], 
                                  prob.symp = prob.symp[2],
                                  frac.cluster = frac.cluster[2]),
                        nseq = reps)
  
  for(ii in 1:reps)
  {
    print(ii)
    write.sim.data(
      file = ii - 1,
      num.mcmcs = 1,
      diffusion.coef = 10 ^ scalars$diffusion.coef[ii],
      tau.s = 0.8,
      tau.l = 0.2,
      founder.prob = scalars$founder.prob[ii],
      max_founder_date = scalars$max_founder_date[ii],
      radius = scalars$radius[ii],
      cluster = TRUE,
      frac.cluster = scalars$frac.cluster[ii],
      imported.prob = 1.0,
      max_cases = 199,
      max_tries = 200,
      prob.symp = scalars$prob.symp[ii],
      prob.trtd.symp = 1.0,
      prob.trtd.asym = 0.0,
      data.dir = data.dir)
  }
}



# function to simulate a single transmission network and write data to file
write.sim.data = function(
  file = 1,
  num.mcmcs = 1,
  diffusion.coef = 1.0,
  tau.s = 1.0,
  tau.l = 0.0,
  founder.prob = 0.574654,
  max_founder_date = 500.0,
  radius = sqrt(1/pi),
  clustered = FALSE,
  frac.cluster = 0.5,
  imported.prob = 1.0,
  max_cases = 1524,
  max_tries = 1524,
  prob.symp = 1.0,
  prob.trtd.symp = 1.0,
  prob.trtd.asym = 0.0,
  data.dir)
{
  # determine number of founders
  max_infections = max_cases
  if(prob.symp < 1.0){
    max_infections = max_cases / (prob.trtd.symp * prob.symp + prob.trtd.asym * (1 - prob.symp))
  }
  founders = max(1,round(founder.prob * max_infections))
  print(founders)
  # simulate network
  repeat{
    # set up some variables for storage
    edges = matrix(1, founders, 2)
    edges[,2] = (1:founders)+1
    dates.infection = runif(founders) * max_founder_date
    
    # simulate properties of the founders
    symp = rbinom(founders, 1, prob.symp)
    trtd = rbinom(founders, 1, ifelse(symp, prob.trtd.symp, prob.trtd.asym))
    imprtd = rbinom(founders, 1, imported.prob)
    hist = rbinom(founders, 1, ifelse(imprtd, tau.s, tau.l))
    dates.detection = dates.infection +
      ifelse(
        symp,
        sample(length(IDP.S),length(dates.infection),prob=IDP.S,replace=T),
        sample(length(IDP.A),length(dates.infection),prob=IDP.A,replace=T))
    
    if(clustered)
    {
      coordinates.founders = points.clustered(n = founders, frac.cluster = frac.cluster)
      latlon.founders = (apply(coordinates.founders, 2, convert.coordinates))
    }else{
      coordinates.founders = replicate(founders, c(0,0))
      latlon.founders = replicate(founders, point.swazi())
    }
    x.infection = coordinates.founders[1,]
    y.infection = coordinates.founders[2,]
    lats.infection = latlon.founders[1,]
    lons.infection = latlon.founders[2,]
    
    # set up parameters for simulations
    tovisit = (1:founders)+1
    R0.local = (max_infections - founders) / founders / (1 + (max_infections - founders) / founders)
    
    # repeat process of simulating offspring for each infection
    while(length(tovisit) > 0){
      # set who the current parent is and how many offspring they have
      parent = tovisit[1]
      number.offspring = rpois(1, R0.local)

      # simulate various statuses of offspring
      symp.offspring = rbinom(number.offspring, 1, prob.symp)
      symp = c(symp, symp.offspring)
      trtd = c(trtd, rbinom(number.offspring, 1, ifelse(tail(symp, number.offspring), prob.trtd.symp, prob.trtd.asym)))
      imprtd = c(imprtd,0)
      hist = c(hist, rbinom(number.offspring, 1, tau.l))
      
      # simulate any offspring that there might be
      if(number.offspring > 0)
      {
        # simulate timing of offspring
        dates.infection.offspring = rep(dates.infection[parent-1], number.offspring) +
          ifelse(
            rep(trtd[parent-1],number.offspring),
            sample(length(GT.T.data), number.offspring, prob = c(GT.T.data), replace = T),
            sample(length(GT.UT.data), number.offspring, prob = c(GT.UT.data), replace = T))
        dates.detection.offspring = dates.infection.offspring +
          ifelse(
            symp.offspring,
            sample(length(IDP.S), length(dates.infection.offspring), prob = IDP.S, replace = T),
            sample(length(IDP.A), length(dates.infection.offspring), prob = IDP.A, replace = T))
        dates.infection = c(dates.infection, dates.infection.offspring)
        dates.detection = c(dates.detection, dates.detection.offspring)
        
        # simulate whereabouts of offspring
        type = paste(paste(ifelse(symp[parent-1], 'S','A'), ifelse(trtd[parent-1], 'T', 'U'), sep = ''),
                     ifelse(tail(symp, n = number.offspring), 'S', 'A'), sep = '')
        sigma.space = sqrt(diffusion.coef * sapply(1:number.offspring, function(x){return(variance.spatial(dates.detection.offspring[x] - dates.detection[parent-1], type[x]))}))
        repeat
        {
          locs.new = locs.offspring(c(0,0), number.offspring, sigma.space)
          latlon.new = t(apply(locs.new, 1, convert.coordinates, ref.lat = lats.infection[parent-1], ref.lon = lons.infection[parent-1]))
          lats = latlon.new[,1]
          lons = latlon.new[,2]
          coords.df = data.frame(lats, lons)
          coords.offspring = SpatialPointsDataFrame(cbind(coords.df$lats, coords.df$lons), data = coords.df)
          proj4string(coords.offspring) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ")
          break
        }
        x.infection = c(x.infection, locs.new[,1])
        y.infection = c(y.infection, locs.new[,2])
        lats.infection = c(lats.infection, latlon.new[,1])
        lons.infection = c(lons.infection, latlon.new[,2])
        rm(type)
        rm(locs.new)
        rm(latlon.new)
        rm(coords.df)
        rm(lats)
        rm(lons)
        
        # add offspring to records
        tovisit = c(tovisit, (max(edges)+1):(max(edges)+number.offspring))
        edges = rbind(
          edges,
          cbind(
            rep(parent, number.offspring),
            (max(edges)+1):(max(edges)+number.offspring)))
      }

      # scratch parent off the list and check if we're done
      tovisit = tovisit[-1]
      if(sum(trtd) >= max_cases + 1){break}
    }

    # convert edge list to network and exit if possible    
    if(sum(trtd) == max_cases){
      network = NULL
      network = matrix(0, nrow = max(edges), ncol = max(edges))
      network[edges] = 1
      break
    }
  }
  
  # write true parameter values to file
  write.table(
    paste(
      'diffusion.coef,', diffusion.coef, '\n',
      'founder.prob,', founder.prob, '\n',
      'max_founder_date,', max_founder_date, '\n',
      'prob.symp,', prob.symp, '\n',
      'R0.local,', R0.local, '\n'),
    file = paste(data.dir, 'scalars_true_', file, '.csv', sep=''),
    quote = FALSE, sep = '', row.names = FALSE, col.names = FALSE)

  # write data about nodes to file
  write(
    paste(
      'Time', 'Type', 'Lat', 'Lon', 'History','L1',
      sep = ','),
    paste(data.dir, 'nodes_', file, '_FULL.csv', sep = ''))
  for(nn in 1:(nrow(network)-1)){
    write(
      paste(
        floor(dates.detection[nn]),
        c('AU','SU','AT','ST')[1 + symp[nn] + 2 * trtd[nn]],
        as.numeric(lats.infection)[nn],
        as.numeric(lons.infection)[nn],
        as.numeric(hist)[nn],
        "00",sep=','),
      paste(data.dir, 'nodes_', file, '_FULL.csv',sep=''),
      append=TRUE)
  }
  
  # write simulated network to file
  string.to.write = ''
  network.ind = which(network > 0, arr.ind = T)
  for(ii in 1:length(which(network > 0))){
    if(network.ind[ii,1] == 1){
      string.to.write = paste(
        string.to.write,
        's-',
        as.character(network[network.ind[ii,1],network.ind[ii,2]]), '-',
        as.character(network.ind[ii,2]-1), ';', sep = '')
    }
    if(network.ind[ii,1] > 1){
      string.to.write = paste(
        string.to.write,
        as.character(network.ind[ii,1]-1), '-',
        as.character(network[network.ind[ii,1],network.ind[ii,2]]), '-',
        as.character(network.ind[ii,2]-1), ';', sep = '')
    }
  }
  for(jj in 1:num.mcmcs){
    write.table(
      string.to.write,
      file=paste(data.dir, 'network_true_', file,'_FULL.csv',sep=''),
      quote=FALSE,sep='',row.names=FALSE,col.names=FALSE)
  }
  
  # subsample network
  network.detected = network
  if(sum(!trtd) > 0){
    undetected = which(!trtd)
    undetected.keep = undetected
    num.undetected = length(undetected)
    for(ii in 1:num.undetected){
      parent.undetected = which(network.detected[,undetected[ii]+1]>0)
      offspring.undetected = which(network.detected[undetected[ii]+1,]>0)
      network.detected[parent.undetected,offspring.undetected] =
        network.detected[parent.undetected,undetected[ii]+1] + 1
      network.detected = network.detected[-(undetected[ii]+1),]
      network.detected = network.detected[,-(undetected[ii]+1)]
      undetected[which(undetected>undetected[ii])] = undetected[which(undetected>undetected[ii])] - 1
    }
  }
  
  # write data about detected nodes to file
  write(
    paste(
      'Time','Type','Lat','Lon','History', 'L1',sep=','),
    paste(data.dir, 'nodes_', file, '_DETECTED.csv',sep=''))
  for(nn in 1:(nrow(network.detected)-1)){
    write(
      paste(
        floor(dates.detection[trtd>0][nn]),
        c('AU','SU','AT','ST')[1 + symp[trtd>0][nn] + 2 * trtd[trtd>0][nn]],
        as.numeric(lats.infection[trtd>0])[nn],
        as.numeric(lons.infection[trtd>0])[nn],
        as.numeric(hist[trtd>0])[nn],
        "00",sep=','),
      paste(data.dir, 'nodes_', file, '_DETECTED.csv', sep=''),
      append=TRUE)
  }

  # write simulated network to file based on detected cases
  string.to.write = ''
  network = network.detected
  network.ind = which(network > 0, arr.ind = T)
  for(ii in 1:length(which(network > 0))){
    if(network.ind[ii,1] == 1){
      string.to.write = paste(
        string.to.write,
        's-',
        as.character(network[network.ind[ii,1],network.ind[ii,2]]), '-',
        as.character(network.ind[ii,2]-1), ';', sep = '')
    }
    if(network.ind[ii,1] > 1){
      string.to.write = paste(
        string.to.write,
        as.character(network.ind[ii,1]-1), '-',
        as.character(network[network.ind[ii,1],network.ind[ii,2]]), '-',
        as.character(network.ind[ii,2]-1), ';', sep = '')
    }
  }
  for(jj in 1:num.mcmcs){
    write.table(
      string.to.write,
      file=paste(data.dir, 'network_true_', file, '_DETECTED.csv',sep=''),
      quote=FALSE,sep='',row.names=FALSE,col.names=FALSE)
  }
}


distance = function(point.1, point.2)
{
  sqrt((point.1[1]- point.2[1])^2 + (point.1[2] - point.2[2])^2)
}

