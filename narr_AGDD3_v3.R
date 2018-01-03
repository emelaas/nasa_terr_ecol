##This script will generate long-term mean and annual 32km maps of spring and autumn phenology
##for a given Landast overlap scene

system.time({
  library(ncdf4)
  require(rgdal)
  library(raster)
  library(foreach)
  library(iterators)
  library(doParallel)
  library(zyp)
  library(Kendall)
  
  site_cor3 <- matrix(NA,75,24)
  site_ST3 <- matrix(NA,75,24)
  site_z <- matrix(NA,75,24)
  
  #Load in list of sidelap scene names
  setwd('/usr3/graduate/emelaas/Code/R/landsat_phenology/')
  scenes <- read.table('overlap_scenes.txt',header=FALSE)
  scenes <- as.character(scenes[,1])
  scenes <- scenes[-c(1,2,7,46,50,59,79,80,82,83,84)]
  
  #Loop through each site (i) and each time range (j) in increments of 15 between 15-120
  for (i in 1:75){
    for (j in 1:24){
      print(c(i,j))
      
      stack_loc <- scenes[i]
      
      #Load in grid-cell averaged observed Landsat phenology
      setwd(paste('/projectnb/modislc/projects/te_phenology/landsat_stacks/',stack_loc,'/Daymet/',sep=""))
      load('narr2landsat_obs_s50')
      obs.SPR <- round(obs.SPR)
      
      #For each year, calculate the mean preseason temperature prior to average leafout date (column 1)
      pred.SPR <- matrix(NA,32,nrow(obs.SPR))
      for (yr in 1984:2013){
        
        #Load in NetCDF file data
        setwd(paste('/projectnb/modislc/projects/te_phenology/landsat_stacks/',stack_loc,'/Daymet/',sep=""))
        load(paste('tmean_',yr,sep=""))
        tmean2 <- tmean
        load(paste('tmean_',yr-1,sep=""))
        tmean1 <- tmean
        tmean <- rbind(tmean1,tmean2)
        
        if ((yr-1)%%4 == 0){
          obs <- obs.SPR+366
        } else {
          obs <- obs.SPR+365
        }
        
        #Calculate Pre-season mean temperature N days before SOS
        N <- j*5
        
        #Calculate GDDs
        tmean[tmean<0] <- 0
        
        preTmean <- matrix(NA,1,ncol(tmean))
        for (c in 1:ncol(tmean)){
          if (is.na(obs[c,1])==0){
            tmean_sub <- tmean[(obs[c,1]-N):obs[c,1],c]
            preTmean[1,c] <- sum(tmean_sub)
          }
        }
        
        pred.SPR[(yr-1982+1),] <- preTmean
      }
      
      obs.SPR <- obs.SPR[,-c(1:2)]
      obs.SPR <- obs.SPR[,-c(33:64)]
      
      pred.SPR <- t(pred.SPR)
      
      pred.SPR.mean <- apply(pred.SPR,2,mean,na.rm=TRUE)
      obs.SPR.mean <- apply(obs.SPR,2,mean,na.rm=TRUE)
      
      #Calculate correlation coefficient between preseason temperature and SOS for each grid cell       
      if (length(which(is.na(obs.SPR.mean)==0))>5){
        lm.pred <- lm(obs.SPR.mean~pred.SPR.mean)
        cell_ST3 <- round(coef(lm.pred),3)[2]
        cell_r3 <- cor(obs.SPR.mean,pred.SPR.mean,use='pairwise.complete.obs')
        
        mk <- MannKendall(pred.SPR.mean)
        cell_mk <- as.numeric(mk$sl)
        
        y <- as.numeric(pred.SPR.mean)
        x <- as.integer(seq(1,32,1))
        z <- zyp.sen(y~x)
        cell_z <- z$coefficients[2] 
      }
      
      #Average across grid cells to obtain site-period specific correlation and sensitivity coefficient (days/degC)
      site_ST3[i,j] <- cell_ST3
      site_cor3[i,j] <- cell_r3
      
      #cell_z[cell_mk>0.05] <- NA
      site_z[i,j] <- cell_z
    }   
  }
  
  site_cor3[is.na(site_cor3)==1] <- 9999
  w <- apply(site_cor3,1,which.min)
  ST3 <- site_ST3[cbind(seq_along(w), w)]
  COR3 <- site_cor3[cbind(seq_along(w), w)]
  COR3[COR3==9999] <- NA
  
  Z <- site_z[cbind(seq_along(w), w)]
  
  setwd('/projectnb/modislc/users/emelaas/scratch15/NASA_TE/R_files/')
  save(site_ST3,site_cor3,ST3,COR3,scenes,site_z,Z,w,file = "preseason_temp_cor_gdd2")
  
})
