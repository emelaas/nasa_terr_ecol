##This script will generate long-term mean and annual 32km maps of spring and autumn phenology
##for a given Landast overlap scene

system.time({
  library(ncdf4)
  require(rgdal)
  library(raster)
  library(foreach)
  library(iterators)
  library(doParallel)
  
  args = commandArgs(trailingOnly=T)
  stack_loc = args[1]  
  
  #stack_loc <- "sup_memphis"
  
  #Register the parallel backend
  registerDoParallel(16)
  
  deg2rad <- function(deg) {(deg * pi) / (180)}
  
  setwd('/projectnb/modislc/users/emelaas/scratch15/NASA_TE/NARR/')
  lon <- raster('air.2m.1979.nc',varname='lon')
  lat <- raster('air.2m.1979.nc',varname='lat')
  
  #Native NARR projection
  narr_crs <- CRS("+proj=lcc +lat_1=50 +lat_0=50 +lon_0=-107 +k_0=1 +x_0=5632642.22547 +y_0=4612545.65137 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
  
  #Generate lat/lon points from rasters and reproject to native NARR projection
  plat <- rasterToPoints(lat)
  plon <- rasterToPoints(lon)
  lonlat <- cbind(plon[,3], plat[,3])
  lonlat <- SpatialPoints(lonlat, proj4string = CRS("+proj=longlat +datum=WGS84"))
  
  plonlat <- spTransform(lonlat, CRSobj = narr_crs)
  
  projection(lat) <- narr_crs
  extent(lat) <- extent(plonlat)
  projection(lon) <- narr_crs
  extent(lon) <- extent(plonlat)
  
  #Import overlap shapefile
  setwd(paste('/projectnb/modislc/projects/te_phenology/landsat_stacks/',
    stack_loc,'/SHP/',sep=""))
  o <- readOGR(paste('/projectnb/modislc/projects/te_phenology/landsat_stacks/',
    stack_loc,'/SHP/',sep=""),stack_loc)
  
  #Reproject overlap shapefile
  o_reproj <- spTransform(o,CRS("+proj=lcc +lat_1=50 +lat_0=50 +lon_0=-107 +k_0=1 +x_0=5632642.22547 +y_0=4612545.65137 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
  
  lat.crop <- crop(lat,extent(o_reproj))
  lat_vals = getValues(lat.crop)
  lon.crop <- crop(lon,extent(o_reproj))
  lon_vals = getValues(lon.crop)
  
  setwd(paste('/projectnb/modislc/projects/te_phenology/landsat_stacks/',stack_loc,'/Daymet/',sep=""))
  load('narr2landsat_pcam_s50_blumel')
  pred.blumel <- pred.SPR
  load('narr2landsat_pcam_s50_sprwarm')
  pred.sprwarm <- pred.SPR
  
  load('narr2landsat_obs_s50_dbf')  
  obs.SPR <- round(t(obs.SPR[,-c(1:2)]))
  
  
  
  
  repmat = function(X,m,n){
    mx = dim(X)[1]
    nx = dim(X)[2]
    matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)} 
  
  chillDays <- matrix(NA,32,length(lon_vals))
  for (yr in 1983:2013){
    print(yr)
    
    #Load in NetCDF file data
    setwd(paste('/projectnb/modislc/projects/te_phenology/landsat_stacks/',stack_loc,'/Daymet/',sep=""))
    load(paste('tmean_',yr,sep=""))
    tmean2 <- tmean
    load(paste('tmean_',yr-1,sep=""))
    
    tmean <- rbind(tmean,tmean2)
    
    m <- matrix(1:nrow(tmean))
    time <- repmat(m,1,ncol(tmean))
    
    tmean[tmean<0] <- 0
    tmean[tmean>5] <- 0
    tmean[tmean!=0] <- 1
    
    tmean[1:244,] <- 0
    cum_tmean <- cumsum(tmean)
    
    tLF <- round(obs.SPR[(yr-1981),])
    
    for (i in ncol(tmean)){
      if (is.na(tLF[i]) == 0){
        if ((yr-1)%%4 == 0){
          chillDays[(yr-1981),] <- cum_tmean[(tLF[i]+366),i]
        } else {
          chillDays[(yr-1981),] <- cum_tmean[(tLF[i]+365),i]
        }
      }
    }
  }
  
  #Convert overlap raster to polygon in order to extract 1km spatial mean phenology dates
  tmp <- seq(1,length(lon.crop),1)
  tmp.sub <- setValues(lon.crop,tmp)
  
  setwd(paste('/projectnb/modislc/projects/te_phenology/landsat_stacks/',stack_loc,'/Daymet/',sep=""))
  save(tmp.sub,obs.SPR,chillDays,pred.sprwarm,pred.blumel,file = "narr2landsat_chillDays")
  
  blumel.diff <- abs(obs.SPR-pred.blumel)
  sprwarm.diff <- abs(obs.SPR-pred.sprwarm)
  blumel.adv <- sprwarm.diff-blumel.diff
  
  ###FIGURES###
  
  #Compare observed and predicted anomalies
  blumel.adv <- rowMeans(blumel.adv)
  chillDays <- rowMeans(chillDays)
  
  lm.pred <- lm(blumel.adv~chillDays)
  lm.pred_coef <- round(coef(lm.pred),3)
  r2 <- summary(lm.pred)$adj.r.squared
  #x11(h=6,w=6)
  pdf(h=6,w=6,paste('/projectnb/modislc/projects/te_phenology/figures/',stack_loc,'_32km_blumel_adv_vs_chillDays.pdf',sep=""))
  plot(chillDays,blumel.adv,ylab='Photoperiod Advantage (days)',
    xlab='No. Chill Days')
  abline(lm.pred)
  rp = vector('expression',2)
  rp[1] = substitute(expression(italic(R)^2 == MYVALUE),
    list(MYVALUE = format(r2,dig=2)))[2]
  rp[2] = substitute(expression(y == MYVALUE2*x+MYVALUE3), 
    list(MYVALUE2 = format(lm.pred_coef[2], digits = 4),
      MYVALUE3 = format(lm.pred_coef[1], digits = 4)))[2]
  legend('bottomright', legend = rp, bty = 'n')
  dev.off()
})