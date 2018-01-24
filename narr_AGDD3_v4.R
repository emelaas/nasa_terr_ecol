##This script will generate long-term mean and annual 32km maps of spring and autumn phenology
##for a given Landast overlap scene

system.time({
  library(ncdf4)
  require(rgdal)
  library(raster)
  library(zyp)
  library(Kendall)
  library(ppcor)
  
  site_cor <- matrix(NA,75,24)
  site_pcor <- matrix(NA,75,24)
  site_ST <- matrix(NA,75,24)
  site_z <- matrix(NA,75,24)
  
  #Load in list of sidelap scene names
  setwd('/usr3/graduate/emelaas/Code/R/landsat_phenology/')
  scenes <- read.table('overlap_scenes.txt',header=FALSE)
  scenes <- as.character(scenes[,1])
  
  setwd('/projectnb/modislc/projects/te_phenology/landsat_stacks')
  
  #Loop through each site (i) and each time range (j) in increments of 5 days between 5 and 120
  for (i in 1:75){
    for (j in 1:24){
      print(c(i,j))
      
      stack_loc <- scenes[i]
      
      #Load in grid-cell averaged observed Landsat phenology
      setwd(paste('/projectnb/modislc/projects/te_phenology/landsat_stacks/',stack_loc,'/Daymet/',sep=""))
      load('narr2landsat_obs_s50_dbf')
      obs.SPR <- round(obs.SPR)
      
      #For each year, calculate the mean preseason temperature prior to average leafout date (column 1)
      preseasonT <- matrix(NA,32,nrow(obs.SPR))
      preseasonP <- matrix(NA,32,nrow(obs.SPR))
      preseasonR <- matrix(NA,32,nrow(obs.SPR))
      for (yr in 1984:2013){
                
        # Load in NetCDF file data
        # Temperature
        setwd(paste('/projectnb/modislc/projects/te_phenology/landsat_stacks/',stack_loc,'/Daymet/',sep=""))
        load(paste('tmean_',yr,sep=""))
        tmean2 <- tmean
        load(paste('tmean_',yr-1,sep=""))
        tmean1 <- tmean
        tmean <- rbind(tmean1,tmean2)
        
        # Precipitation
        load(paste('apcp_',yr,sep=""))
        apcp2 <- apcp
        load(paste('apcp_',yr-1,sep=""))
        apcp1 <- apcp
        apcp <- rbind(t(apcp1),t(apcp2))
        
        # Radiation
        load(paste('dswrf_',yr,sep=""))
        dswrf2 <- dswrf
        load(paste('dswrf_',yr-1,sep=""))
        dswrf1 <- dswrf
        dswrf <- rbind(t(dswrf1),t(dswrf2))
        
        
        obs <- obs.SPR[,1]+nrow(tmean1)
        
        #Calculate Pre-season mean temperature N days before SOS
        N <- j*5
        
        #Calculate GDDs
        tmean[tmean<0] <- 0
        
        preTsum <- matrix(NA,1,ncol(tmean)) # Temperature
        prePsum <- matrix(NA,1,ncol(tmean)) # Precipitation
        preRsum <- matrix(NA,1,ncol(tmean)) # Radiation
        for (c in 1:ncol(tmean)){
          if (is.na(obs[c])==0){
            tmean_sub <- tmean[(obs[c]-N):obs[c],c]
            preTsum[1,c] <- sum(tmean_sub)
            apcp_sub <- apcp[(obs[c]-N):obs[c],c]
            prePsum[1,c] <- sum(apcp_sub)
            dswrf_sub <- dswrf[(obs[c]-N):obs[c],c]
            preRsum[1,c] <- sum(dswrf_sub)
          }
        }
        
        preseasonT[(yr-1982+1),] <- preTsum
        preseasonP[(yr-1982+1),] <- prePsum
        preseasonR[(yr-1982+1),] <- preRsum
      }
      
      obs.SPR <- obs.SPR[,-c(1:2)]
            
      #Calculate correlation coefficient between preseason temperature and SOS for each grid cell       
      cell_r <- matrix(NA,nrow(obs.SPR),1) # correlation coefficient
      cell_pr <- matrix(NA,nrow(obs.SPR),1) # partial correlation coefficient
      cell_ST <- matrix(NA,nrow(obs.SPR),1) # temperature sensitivity
      cell_mk <- matrix(NA,nrow(obs.SPR),1) # Mann-Kendall statistic
      cell_z <- matrix(NA,nrow(obs.SPR),1) # trend in AGDD using Theil-Sen
      for (k in 1:nrow(obs.SPR)){
        if (length(which(is.na(obs.SPR[k,1:32])==0))>5){
          lm.pred <- lm(obs.SPR[k,1:32]~preseasonT[1:32,k])
          cell_ST[k] <- round(coef(lm.pred),3)[2]
          cell_r[k] <- cor(obs.SPR[k,1:32],preseasonT[1:32,k],use='pairwise.complete.obs')
          
          # Partial correlation analysis
          df <- data.frame(obs.SPR[k,1:32],preseasonT[1:32,k],preseasonP[1:32,k],preseasonR[1:32,k])
          colnames(df) <- c('SOS','preT','preP','preR')
          w <- which(is.na(df$SOS)==0 & is.na(df$preT)==0 & is.na(df$preP)==0 & is.na(df$preR)==0)
          df <- df[w,]
          cell_pcor <- pcor(df,method='pearson')
          cell_pr[k] <- cell_pcor$estimate[1,2]
          
          mk <- MannKendall(preseasonT[1:32,k])
          cell_mk[k] <- as.numeric(mk$sl)
          
          y <- as.numeric(preseasonT[1:32,k])
          x <- as.integer(seq(1,32,1))
          z <- zyp.sen(y~x)
          cell_z[k] <- z$coefficients[2] 
        }
      }
      
      #Average across grid cells to obtain site-period specific correlation and sensitivity coefficient (days/degC)
      site_ST[i,j] <- mean(cell_ST,na.rm=TRUE)
      site_cor[i,j] <- mean(cell_r,na.rm=TRUE)
      site_pcor[i,j] <- mean(cell_pr,na.rm=TRUE)
      
      #cell_z[cell_mk>0.05] <- NA
      site_z[i,j] <- mean(cell_z,na.rm=TRUE)
    }   
  }
  
  site_cor[is.na(site_cor)==1] <- 9999
  w <- apply(site_cor,1,which.min)
  ST <- site_ST[cbind(seq_along(w), w)]
  COR <- site_cor[cbind(seq_along(w), w)]
  COR[COR==9999] <- NA
  
  site_pcor[is.na(site_pcor)==1] <- 9999
  w2 <- apply(site_pcor,1,which.min)
  ST2 <- site_ST[cbind(seq_along(w2), w2)]
  PCOR <- site_pcor[cbind(seq_along(w2), w2)]
  PCOR[PCOR==9999] <- NA
  
  Z <- site_z[cbind(seq_along(w), w)]
  Z2 <- site_z[cbind(seq_along(w2), w2)]
  
  setwd('/projectnb/modislc/users/emelaas/scratch15/NASA_TE/R_files/')
  save(site_ST,site_cor,site_pcor,ST,ST2,COR,PCOR,scenes,site_z,Z,Z2,w,w2,
    file = "preseason_temp_cor_gdd_v4")
  
})
