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
  
  #stack_loc <- "ma_tallahassee"
  
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
  
  for (yr in 1982:2013) {
    print(yr)
    
    #Load in NetCDF file data
    setwd('/projectnb/modislc/users/emelaas/scratch15/NASA_TE/NARR/')
    
    tair <- stack(paste('air.2m.',yr,'.nc',sep=""), varname="air")
    projection(tair) <- narr_crs
    extent(tair) <- extent(plonlat)
    tair.crop <- crop(tair,extent(o_reproj))
    tair_vals <- getValues(tair.crop)
    if (yr%%4==0) day <- rep(1:366,each=8)
    if (yr%%4!=0) day <- rep(1:365,each=8)
    tmean <- aggregate(t(tair_vals),by=list(day),FUN=mean)
    tmean <- tmean[,-1]
    tmean <- tmean-273.15
    
    apcp <- stack(paste('apcp.',yr,'.nc',sep=""), varname="apcp")
    dswrf <- stack(paste('dswrf.',yr,'.nc',sep=""), varname="dswrf")
    
    projection(apcp) <- narr_crs
    extent(apcp) <- extent(plonlat)
    projection(dswrf) <- narr_crs
    extent(dswrf) <- extent(plonlat)
    
    #Crop NARR files using overlap extent
    apcp.crop <- crop(apcp,extent(o_reproj))
    dswrf.crop <- crop(dswrf,extent(o_reproj))
  
    apcp <- getValues(apcp.crop)
    dswrf <- getValues(dswrf.crop)
      
    setwd(paste('/projectnb/modislc/projects/te_phenology/landsat_stacks/',stack_loc,'/Daymet/',sep=""))
    #save(tmean,file=paste('tmean_',yr,sep=''))
    save(apcp,file=paste('apcp_',yr,sep=''))
    save(dswrf,file=paste('dswrf_',yr,sep=''))
  }
  
  #Convert overlap raster to polygon in order to extract 1km spatial mean phenology dates
  tmp <- seq(1,length(lon.crop),1)
  tmp.sub <- setValues(lon.crop,tmp)
  tmp_poly <- rasterToPolygons(tmp.sub)  
  
  ##### Calculate zonal mean of phenology across Daymet polygon map #####
  
  #Import annual phenology dates
  setwd(paste('/projectnb/modislc/projects/te_phenology/landsat_stacks/',stack_loc,sep=""))
  if (file.exists('EOSD')==1){
    lc <- raster(paste('/projectnb/modislc/projects/te_phenology/landsat_stacks/',
      stack_loc,'/EOSD/eosd_clip.bip',sep=""))
    #w <- c(0,11,12,20,31,32,33,40,51,52,82,83,100)
    w <- c(0,11,12,20,31,32,33,40,51,52,81,82,83,100,211,212,213)
  } else if (file.exists('NLCD')==1) {
    lc <- raster(paste('/projectnb/modislc/projects/te_phenology/landsat_stacks/',
      stack_loc,'/NLCD/nlcd_clip.bip',sep=""))
    #w <- c(11,12,21,24,31,51,52,71,72,73,74,81,82,95)
    w <- c(11,12,21,22,23,24,31,42,51,52,71,72,73,74,81,82,90,95)
  }
  
  tmp_poly_proj <- spTransform(tmp_poly,projection(lc))
  
  load(paste('/projectnb/modislc/projects/te_phenology/landsat_stacks/',stack_loc,
    '/phenology/pheno_nodist_mat',sep=""))
  pheno_mat[pheno_mat==0]<-NA
  pheno_mat[which(pheno_mat[,2] %in% w),8:73] <- NA
  pheno_mat[which(is.na(pheno_mat[,2])==1),8:73] <- NA
  
  obs.SPR <- matrix(NA,length(tmp_poly),34) #ltmSPR,ltmAUT,32yrsSPR,32yrsAUT
  for (col in 8:41){
    pheno <- matrix(NA,ncell(lc),1)
    pheno[pheno_mat[,1]] <- pheno_mat[,col]
    pheno_hdr <- setValues(lc,pheno)
    
    #Zonal stats for LTM spring/autumn and each inidiv. year
    poly_means <- foreach(i = 1:length(tmp_poly), .combine = rbind) %dopar% {
      print(c(i,col))
      tmp_poly_subset <- tmp_poly_proj[tmp_poly_proj$longitude.coordinate==i, ]
      poly <- extract(pheno_hdr,tmp_poly_subset)
      poly_mean <- unlist(lapply(poly, 
        function(x) if (!is.null(x) && length(which(!is.na(unlist(x))))>500) 
          mean(x, na.rm=TRUE) else NA))
      ## length(which(!is.na(unlist(x))))/length(unlist(x))>0.05
    }
    
    obs.SPR[,col-7]<-poly_means
  }
  
  setwd(paste('/projectnb/modislc/projects/te_phenology/landsat_stacks/',stack_loc,'/Daymet/',sep=""))
  save(tmp.sub,obs.SPR,file = "narr2landsat_obs_s50_dbf")
  
})