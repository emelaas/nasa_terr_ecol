#This version includes 0.02 subtraction from ETM+ data according to DSM's findings

system.time({
  
  #args = commandArgs(trailingOnly=T)
  #stack_loc = args[1]  
  
  stack_loc <- "ah_hubbard"
  
  library('rgdal')
  library("raster", lib.loc="/project/earth/packages/R-3.1.0/lib64/R/library")
  library("foreach", lib.loc="/project/earth/packages/R-3.1.0/lib64/R/library")
  library("iterators", lib.loc="/project/earth/packages/R-3.1.0/lib64/R/library")
  library("doParallel", lib.loc="/project/earth/packages/R-3.1.0/lib64/R/library")
  
  .libPaths('/usr3/graduate/emelaas/R/x86_64-unknown-linux-gnu-library/3.1')
  library("Kendall", lib.loc='/usr3/graduate/emelaas/R/x86_64-unknown-linux-gnu-library/3.1')
  library("plm", lib.loc='/usr3/graduate/emelaas/R/x86_64-unknown-linux-gnu-library/3.1')
  library("lmtest", lib.loc='/usr3/graduate/emelaas/R/x86_64-unknown-linux-gnu-library/3.1')
  library("sandwich", lib.loc='/usr3/graduate/emelaas/R/x86_64-unknown-linux-gnu-library/3.1')
  library("zyp", lib.loc='/usr3/graduate/emelaas/R/x86_64-unknown-linux-gnu-library/3.1')
  
  #Register the parallel backend
  registerDoParallel(16)
    
  ##########################################################################################
  #LANDSAT PHENOLOGY ALGORITHM RUN ON EACH PIXEL WITH SUFFICIENT NO. OF OBSERVATIONS
  
  #Phenology algorithm parameters:
  #evi_row - EVI time series
  #doy - Julian day of year
  #year - Year
  #sensor - Landat sensor (4,5(TM),7(ETM+))
  #prd - three-year period 
  #dYR - disturbance year 
  #pID - pixel ID
  #pLC - pixel land cover
  #pCDL - pixel cropland data layer
  #mID - MODIS panel ID
  #aID - AVHRR panel ID
  Landsat_Phenology <- function(evi_row,doy,yr,sensor,prd,dYR,pID,pLC,pCDL,mID,aID){
    
    #normalized EVI threshold for long-term mean phenology date used to calculate annual phenology
    thresh <- 0.5
    
    #Initialize Phenology matrices
    phenoSPR <- matrix(0,1,32) #annual spring phenology (1982-2013)
    phenoAUT <- matrix(0,1,32) #annual autumn phenology
    ltmSPR <- matrix(0,1,1) #long-term mean spring phenology date
    ltmAUT <- matrix(0,1,1) #long-term mean autumn phenology date
    bline <- matrix(0,1,1) #long-term mean winter background EVI
    hline <- matrix(0,1,1) #long-term mean summer maximum EVI
    rsmooth <- matrix(0,1,1) #correlation between observed and smoothed EVI
    nobs <- matrix(0,1,1) #number of available EVI observations used
    
    SPR_thresh = matrix(0,1,1)
    AUT_thresh = matrix(0,1,1)
    SPRsmooth_max = matrix(0,1,1)
    SPRsmooth_min = matrix(0,1,1)
    AUTsmooth_max = matrix(0,1,1)
    AUTsmooth_min = matrix(0,1,1)
    
    info_spr <- NA
    info_aut <- NA
    
    #Process phenology dates for each pixel in EVI matrix
    
    #Apply calibration adjustment for TM vs. ETM+ (subtract 0.02 from ETM+ observations)
    #Damien's analysis
    w <- which(sensor==7)
    evi_row[w] <- (evi_row[w]-0.019)/1.038
    
    #Set erroneous EVI values to NA
    EVI <- evi_row
    EVI[EVI==-Inf | EVI==Inf]<-NA
    EVI[EVI<0 | EVI>1]<-NA
    
    #Remove NAs and, if necessary, observations after disturbance year (dYR)
    if (dYR!=0) pos <- which(is.na(EVI)==0 & yr<dYR)
    if (dYR==0) pos <- which(is.na(EVI)==0)
    EVI <- EVI[pos]
    DOY <- doy[pos]
    YR <- yr[pos]
    PRD <- prd[pos]
    SEN <- sensor[pos]
    nobs <- length(which(is.na(EVI)==0)) #count number of available observations
    
    #Normalize EVI time series using 10% and 90% quantiles for EVI during three-year windows
    DF <- data.frame(PRD,EVI)
    colnames(DF) <- c('yr','evi')
    q1 <- with(DF,tapply(EVI,PRD,quantile,probs=0.10))
    q2 <- with(DF,tapply(EVI,PRD,quantile,probs=0.90))
    
    quant1 <- matrix(NA,10,1)
    quant2 <- matrix(NA,10,1)
    quant1[as.numeric(names(q1))] <- q1
    quant2[as.numeric(names(q2))] <- q2
    
    EVImax <- quant2[PRD]
    EVImin <- quant1[PRD]
    EVInorm <- (EVI-EVImin)/(EVImax-EVImin)
    
    EVInorm[EVInorm==Inf | EVInorm==-Inf]<-NA
    
    pos <- which(is.na(EVInorm)==0 & DOY<=365)
    EVInorm <- EVInorm[pos]
    DOY <- DOY[pos]
    YR <- YR[pos]
    SEN <- SEN[pos]
    
    if (nobs > 100){
      #Sort DOY/EVI by DOY and compute winter background EVI for improved smoothing spline fit
      x <- cbind(DOY,EVInorm) #Spring background
      x <- x[order(x[,1]),]
      start <- seq(1,x[1,1],1) 
      end <- seq(x[nrow(x),1],365,1)
      SENextra <- c(SEN,rep(3,each=length(start)),rep(3,each=length(end)))
      YRextra <- c(YR,rep(2012,each=length(start)),rep(2012,each=length(end)))
      DOYextra <- c(DOY,start,end)
      EVIextra <- c(EVInorm,matrix(0,length(start)+length(end),1))
      
      #Fit smoothing spline through EVI and DOY data
      fit <- smooth.spline(DOYextra,EVIextra,spar=0.55)
      EVIsmooth <- data.frame(predict(fit,x=1:365))
      rsmooth <- cor(EVIsmooth[DOY,2],EVInorm)
      
      if (rsmooth>0.85){
        #Separate spline into spring and autumn segments using annual maximum      
        pkval <- which.max(EVIsmooth[,2])
        SPRsmooth <- EVIsmooth[1:pkval,2]
        AUTsmooth <- EVIsmooth[(pkval+1):365,2];
        
        #Compute half-maximum of spring logistic for "ruling in" image dates (points) for
        #anamoly calculation
        SPR_thresh <- which.min(abs(SPRsmooth-thresh))
        SPR_halfmax <- which.min(abs(SPRsmooth-0.5))
        
        #Find anomalies inside of designated box
        SPRsmooth_max <- 1
        SPRsmooth_min <- 0
        box_max <- SPRsmooth_max-0.2*(SPRsmooth_max-SPRsmooth_min)
        box_min <- SPRsmooth_min+0.2*(SPRsmooth_max-SPRsmooth_min)
        
        #Generate a matrix with candidate spring phenology observations
        if (is.na(box_max)==0 && is.na(box_min)==0){
          info <- cbind(YR,DOY,SEN,EVInorm,matrix(0,length(DOY),1))
          info <- info[order(info[,1],info[,2]),]
          info <- rbind(matrix(0,1,5),info)
          pos <- which(is.na(info[,4])==0 & info[,4]>box_min & info[,4]<box_max &
              info[,2]<SPR_halfmax+20 & info[,2]>SPR_halfmax-20)
          if (length(pos)>4){
            k <- 1
            for (p in 1:length(pos)){
              smooth_ratio <- abs(SPRsmooth-info[pos[p],4])
              info[pos[p],5] <- which.min(smooth_ratio)
              
              if (k==1)
                info_spr <- info[pos[p],]
              else
                info_spr <- rbind(info_spr,info[pos[p],])
              
              k <- k+1
            }
            w <- which((info[pos,4]<info[(pos-1),4] & info[pos,1]==info[(pos-1),1])==1)
            if (length(w)>0) info_spr <- info_spr[-w,]
          }
        }
        
        #Compute half-maximum of spring logistic for "ruling in" image dates (points) for
        #anamoly calculation
        AUT_thresh <- which.min(abs(AUTsmooth-thresh))+pkval
        AUT_halfmax <- which.min(abs(AUTsmooth-0.5))+pkval
        
        #Find anomalies inside of designated box
        AUTsmooth_max <- 1
        AUTsmooth_min <- 0
        box_max <- AUTsmooth_max-0.3*(AUTsmooth_max-AUTsmooth_min)
        box_min <- AUTsmooth_min+0.2*(AUTsmooth_max-AUTsmooth_min)
        
        #Generate matrix with candidate autumn phenology observations
        if (is.na(box_max)==0 && is.na(box_min)==0){
          info <- cbind(YR,DOY,SEN,EVInorm,matrix(0,length(DOY),1))
          info <- info[order(info[,1],info[,2]),]
          info <- rbind(info,matrix(0,1,5))
          pos <- which(is.na(info[,4])==0 & info[,4]>box_min & info[,4]<box_max &
              info[,2]<AUT_halfmax+20 & info[,2]>AUT_halfmax-20)
          if (length(pos)>4){
            k <- 1
            for (p in 1:length(pos)){
              smooth_ratio <- abs(AUTsmooth-info[pos[p],4])
              info[pos[p],5] <- which.min(smooth_ratio)+pkval
              
              if (k==1)
                info_aut <- info[pos[p],]
              else
                info_aut <- rbind(info_aut,info[pos[p],])
              
              k <- k+1
            }
            w <- which((info[pos,4]<info[(pos+1),4] & info[pos,1]==info[(pos+1),1])==1)
            if (length(w)>0) info_aut <- info_aut[-w,]
          }
        }
        
        #Calculate interannual phenology dates by taking the distance
        #between each candidate observation and where the same magnitude 
        #of EVI occurs on the spline
        if (exists('info_spr') == 1 && exists('info_aut') == 1){
          if (ncell(info_spr) > 20 & ncell(info_aut) > 20) {
            info_spr <- cbind(info_spr,SPR_thresh+(info_spr[,2]-info_spr[,5]))
            info_aut <- cbind(info_aut,AUT_thresh+(info_aut[,2]-info_aut[,5]))
            
            for (y in 1982:2013){
              pos1 <- which(info_spr[,1] == y)
              pos2 <- which(info_aut[,1] == y)
              
              if (length(pos1) > 0 && ncol(info_spr)==6){
                w <- which.min(abs(info_spr[pos1,4]-0.5))
                phenoSPR[1,y-1981] <- ceiling(mean(info_spr[pos1[w],6]))
              }
              if (length(pos2) > 0 && ncol(info_aut)==6){
                w <- which.min(abs(info_aut[pos2,4]-0.5))
                phenoAUT[1,y-1981] <- ceiling(mean(info_aut[pos2[w],6]))
              }
            }
          }
          
          ltmSPR <- SPR_thresh
          ltmAUT <- AUT_thresh
          #bline <- median(x[1:15,2])
          #hline <- SPRsmooth_max
          
          remove(info_spr,info_aut) 
        }
      }
    }
    
    pheno_matrix <- cbind(pID,pLC,pCDL,mID,aID,nobs,rsmooth,ltmSPR,ltmAUT,phenoSPR,phenoAUT)
    
    if (is.character(pheno_matrix) == 1){
      pheno_matrix <- matrix(NA,1,73)
      print(paste("Pixel ID",pID,"had an error!",sep=''))
    } 
    
    return(pheno_matrix)
  }
  
  #Collect information re: Rdata files in scene folder and formulate table
  ## Damien's method:
  #in_dirs <- list.files(path=data_loc,pattern=glob2rx("L*"),full.names=T,include.dirs=T)
  
  setwd(paste('/projectnb/modislc/projects/te_phenology/landsat_stacks/',stack_loc,'/extract_R_dat',sep=""))
  dir1 <- dir()
  g <- grep(paste("^",stack_loc,'_evi_c',sep=""),dir1)
  g2 <- grep(paste("^",stack_loc,'_dist_r',sep=""),dir1)
  
  dir1_data <- dir1[g]
  dir1_dist <- dir1[g2]
  
  nums <- as.numeric(unlist(strsplit(unlist(dir1_data),"[^0-9]+")))
  dim(nums)<-c(length(nums)/length(g),length(g))
  nums <- data.frame(g,t(nums[2:3,]))
  nums <- nums[order(nums[,3]),]
  colnames(nums) <- c('g','c','r')
  
  nums2 <- unique(na.omit(as.numeric(unlist(strsplit(unlist(dir1_dist),"[^0-9]+")))))
  dim(nums2)<-c(length(nums2)/length(g2),length(g2))
  nums2 <- data.frame(g2,t(nums2))
  nums2 <- nums2[order(nums2[,2]),]
  colnames(nums2) <- c('g2','r')
  
  #LOAD IN EVI TIME SERIES CHUNKS
  for (chunk in 1:max(nums$r)){
    
    setwd(paste('/projectnb/modislc/projects/te_phenology/landsat_stacks/',stack_loc,
      '/extract_R_dat',sep=""))
    
    #Load disturbance information
    load(dir1[nums2$g2[chunk]])
    assign('dist_mat',out) 
    
    #Load in chunks of EVI time series for overlap scene
    nums_chunk <- nums[nums$r==chunk,]
    
    for (i in 1:nrow(nums_chunk)){
      print(c(chunk,i))  
      load(dir1[nums_chunk$g[i]])
      if (i == 1){
        nam <- "header_all"
        assign(nam, header)  
        nam <- "evi_mat"
        assign(nam, out)
      } else {
        header_all <- cbind(header_all,header[,6:ncol(header)])
        evi_mat <- cbind(evi_mat,out[,6:ncol(out)])
      }
      rm(out,header)
    }
    
    #Only keep pixels with 0 or 1 disturbances and the disturbance occurring 
    #after 1999
    w <- which(dist_mat[,5]==0 | dist_mat[,5]==1 & dist_mat[,3]>1999)
    dist_mat <- dist_mat[w,]
    evi_mat <- evi_mat[w,]
    
    #Remove AG pixels
    setwd(paste('/projectnb/modislc/projects/te_phenology/landsat_stacks/',stack_loc,sep=""))
    if (file.exists('EOSD')==1){
      w <- which(evi_mat[,2]!=100)
    } else if (file.exists('NLCD')==1) {
      w <- which(evi_mat[,2]!=81 & evi_mat[,2]!=82)  
    }    
    dist_mat <- dist_mat[w,]
    evi_mat <- evi_mat[w,]  

    dist_yr <- dist_mat[,3] #Year of disturbance
    pixID <- evi_mat[,1] #Pixel ID (cell location)
    pixLC <- evi_mat[,2] #Pixel LC (NLCD)
    pixCDL <- evi_mat[,3] #Pixel Cropland Datatype
    modisID <- evi_mat[,4] #MODIS panel ID
    avhrrID <- evi_mat[,5] #AVHRR panel ID
    evi_mat <- evi_mat[,6:ncol(evi_mat)]
    
    if (nrow(evi_mat)<10000){
      block_width <- nrow(evi_mat)-1
    } else {
      block_width <- 10000
    }
    nblocks <- nrow(evi_mat)%/%block_width
    bs_start <- seq(1,nblocks*block_width+1,block_width)
    bs_end <- c(seq(block_width,nblocks*block_width,block_width),nrow(evi_mat))
    
    #Compile Year, DOY, Sensor and Path information for each image
    yr <- as.numeric(header_all[1,6:ncol(header_all)])
    doy <- as.numeric(header_all[2,6:ncol(header_all)])
    info1 <- data.frame(yr,doy,seq(1,length(doy)))
    colnames(info1) <- c('yr','doy','ind')
    
    sceneIDs <- read.table(paste("/projectnb/modislc/projects/te_phenology/landsat_stacks/",stack_loc,"/sceneIDs/sceneIDs.txt",sep=""))
    sceneIDs <- as.character(as.matrix(sceneIDs))
    info_sceneIDs <- data.frame(as.numeric(substr(sceneIDs,10,13)),as.numeric(substr(sceneIDs,14,16)),
      as.numeric(substr(sceneIDs,3,3)),as.numeric(substr(sceneIDs,4,6)),as.numeric(substr(sceneIDs,7,9)))
    colnames(info_sceneIDs) <- c('yr','doy','sensor','path','row')
    
    if (stack_loc=='ah_hubbard'){
      w <- which(info_sceneIDs$path==12 & info_sceneIDs$row==30 | info_sceneIDs$path==13 & info_sceneIDs$row==29)
      info_sceneIDs <- info_sceneIDs[w,]
    }
    if (stack_loc=='mwp_harvard'){
      w <- which(info_sceneIDs$path==12 & info_sceneIDs$row==31 | info_sceneIDs$path==13 & info_sceneIDs$row==30)
      info_sceneIDs <- info_sceneIDs[w,]
    }
    
    info <- merge(info1,info_sceneIDs,by=c('yr','doy'),sort=FALSE)
    yr <- info$yr; doy <- info$doy; sensor <- info$sensor; path <- info$path
    
    if (ncol(evi_mat) > length(doy)){
      evi_mat <- evi_mat[,info$ind]
    }
    
    #Create "period" index for normalizing EVI time series 
    prd <- sensor
    prd[yr>=1982 & yr<=1986]<-1
    prd[yr>=1987 & yr<=1989]<-2
    prd[yr>=1990 & yr<=1992]<-3
    prd[yr>=1993 & yr<=1995]<-4
    prd[yr>=1996 & yr<=1998]<-5
    prd[yr>=1999 & yr<=2001]<-6
    prd[yr>=2002 & yr<=2004]<-7
    prd[yr>=2005 & yr<=2007]<-8
    prd[yr>=2008 & yr<=2010]<-9
    prd[yr>=2011 & yr<=2014]<-10
    
    pheno_mat_all <- matrix(NA,nrow(evi_mat),73)
    for (j in 1:length(bs_start)){
      #Use foreach to apply function(s) to r, each node will process a block until there are no more
      pheno_mat <- foreach(i = bs_start[j]:bs_end[j], .combine = rbind) %dopar% {
        if (i%%10000==0) print(c(chunk,i))
        
        dYR <- dist_yr[i] #year of first disturbance
        pID <- pixID[i]
        pLC <- pixLC[i]
        pCDL <- pixCDL[i]
        mID <- modisID[i]
        aID <- avhrrID[i]
        evi_row <- evi_mat[i,]
        pheno_matrix <- try(Landsat_Phenology(evi_row,doy,yr,sensor,prd,dYR,pID,pLC,pCDL,mID,aID)) 
        
#         if(!inherits(pheno_matrix,"try-error")){ 
#           #print(paste("Row",i,"was compiled successfully."))
#         } else{
#           print(paste("Row",i,"had an error."))
#           pheno_matrix <- matrix(NA,1,73)
#         }
        
      }
      
      pheno_mat_all[bs_start[j]:bs_end[j],] <- pheno_mat
    }
    
    w <- which(pheno_mat_all[,7]>0.85)
    pheno_mat_all <-pheno_mat_all[w,]
    evi_mat <- evi_mat[w,]
    
    print(paste("Saving pheno_matrix for",stack_loc))
    setwd(paste('/projectnb/modislc/projects/te_phenology/landsat_stacks/',stack_loc,'/phenology',sep=""))
    save(pheno_mat_all,file = paste("pheno_nodist_mat",chunk,sep=""))
    #save(evi_mat,doy,yr,sensor,prd,file = paste("evi_mat",chunk,sep=""))
    
    nam <- paste("pheno_mat_all",chunk,sep="")
    assign(nam, pheno_mat_all)
    
    rm(evi_mat,pheno_mat_all)
  }
  
  pheno_mat <- rbind(pheno_mat_all1,pheno_mat_all2,pheno_mat_all3,pheno_mat_all4,pheno_mat_all5,
    pheno_mat_all6,pheno_mat_all7,pheno_mat_all8,pheno_mat_all9,pheno_mat_all10)
  
  print(paste("Saving pheno_matrix for",stack_loc))
  setwd(paste('/projectnb/modislc/projects/te_phenology/landsat_stacks/',stack_loc,'/phenology',sep=""))
  save(pheno_mat,doy,yr,sensor,prd,file = "pheno_nodist_mat")
  
  ##########################################################################################
  #PANEL ANALYSIS BY PANEL ID (MODIS,AVHRR,ETC.)
  Panel_Analysis <- function(spr_panel,aut_panel,panel_num){
    
    plm1 <- data.frame(matrix(NA,1,1)); colnames(plm1) <- "coefficients"    
    plm2 <- data.frame(matrix(NA,1,1)); colnames(plm2) <- "coefficients"    
    plm1_pvalue <- matrix(NA,1,1)
    plm2_pvalue <- matrix(NA,1,1)
    
    #If there are at least ten deciduous pixels in a panel, make a panel
    if (length(spr_panel)>28*32){
      yr <- rep(1982:2013,times=ncol(spr_panel))
      pixID <- rep(1:ncol(spr_panel),each=32)
      panelID <- rep(panel_num,each=32*ncol(spr_panel))
      
      #dum1984 <- yr; dum1984[yr!=1984]<-0; dum1984[yr==1984]<-1
      #dum2012 <- yr; dum2012[yr!=2012]<-0; dum2012[yr==2012]<-1
      
      dim(spr_panel) <- c(length(spr_panel),1)
      
      panel1 <- cbind(panelID,pixID,yr,spr_panel)
      colnames(panel1) <- c("panelID","pixID","yr","springdoy")
      panel1 <- data.frame(panel1)
      
      #Fixed effects regression (panel)
      #names(panel1)[4] <- "springdoy"
      panel1$t <- panel1$yr - 1982
      
      plm1 <- plm(springdoy ~ t, index=c('pixID', 'yr'), model="within", data=panel1)
      plm1_pvalue <- coeftest(plm1, vcov.=vcovHC(plm1, method="arellano", cluster=c("group", "time")))[1, 4]      
    }
    
    #If there are at least ten deciduous pixels in a panel, make a panel
    if (length(aut_panel)>28*32){
      yr <- rep(1982:2013,times=ncol(aut_panel))
      pixID <- rep(1:ncol(aut_panel),each=32)
      panelID <- rep(panel_num,each=32*ncol(aut_panel))
      
      #dum1984 <- yr; dum1984[yr!=1984]<-0; dum1984[yr==1984]<-1
      #dum2012 <- yr; dum2012[yr!=2012]<-0; dum2012[yr==2012]<-1
      
      dim(aut_panel) <- c(length(aut_panel),1)
      
      panel2 <- cbind(panelID,pixID,yr,aut_panel)
      colnames(panel2) <- c("panelID","pixID","yr","DOY")
      panel2 <- data.frame(panel2)
      
      #Fixed effects regression (panel)
      names(panel2)[4] <- "autumndoy"
      panel2$t <- panel2$yr - 1982
      
      plm2 <- plm(autumndoy ~ t, index=c('pixID', 'yr'), model="within", data=panel2)
      plm2_pvalue <- coeftest(plm2, vcov.=vcovHC(plm2, method="arellano", cluster=c("group", "time")))[1, 4]      
    } 
    
    p <- t(matrix(c(panel_num,as.numeric(plm1$coefficients),plm1_pvalue,
      as.numeric(plm2$coefficients),plm2_pvalue)))
    return(p)
  }
  
  #Include only deciduous pixels with successful retrievals and non-zero panel IDs
  load(paste('/projectnb/modislc/projects/te_phenology/landsat_stacks/',stack_loc,'/phenology/pheno_nodist_mat',sep=""))
  pheno_mat <- pheno_mat[pheno_mat[,7]>=0.85,] 
  pheno_mat[pheno_mat==0]<-NA
  pheno_mat <- pheno_mat[is.na(pheno_mat[,8])==0 & is.na(pheno_mat[,9])==0,]
  
  #Extract all unique panel IDs from pheno matrix (col4 = modis; col5 = AVHRR)
  u <- unique(pheno_mat[,4])
  
  #Loop through each panel ID and perform panel analysis
  panel_stats <- foreach(i = 1:length(u), .combine = rbind) %dopar% {
    print(i)
    
    #Which Landsat pixels with successful retrievals have this panel ID?
    panel_num <- u[i]
    w <- which(pheno_mat[,4]==panel_num)
    
    if (length(w)>28){
      panel <- pheno_mat[w,]
      
      w <- which(is.na(panel[,5])==1)
      panel[w,5] <- 1
      avhrr_pix <- as.numeric(names(table(panel[,5]))[which.max(table(panel[,5]))])
      
      spr_panel <- t(panel[,10:41])
      aut_panel <- t(panel[,42:73])
      
      spr_mean <- apply(spr_panel,1,mean,na.rm=TRUE)
      aut_mean <- apply(aut_panel,1,mean,na.rm=TRUE)
      
      p <- Panel_Analysis(spr_panel,aut_panel,panel_num)
      if(!inherits(p,"try-error")){ 
        print(paste("Panel",i,"was compiled successfully."))
      } else{
        print(paste("Panel",i,"had an error."))
        p <- matrix(NA,1,5)
      }
      
      p <- cbind(p,t(spr_mean),t(aut_mean),avhrr_pix)
    } else {
      p <- cbind(panel_num,matrix(0,1,68),1)
    }
  }
  
  #Save panel stats matrix
  setwd(paste('/projectnb/modislc/projects/te_phenology/landsat_stacks/',stack_loc,'/phenology/',sep=""))
  save(panel_stats,file="panel_nodist_stats")
  
  #Extract all unique panel IDs from pheno matrix (AVHRR)
  u <- unique(panel_stats[,70])
  
  #Loop through each panel ID and perform panel analysis
  panel_stats_avhrr <- foreach(i = 1:length(u), .combine = rbind) %dopar% {
    print(i)
    
    #Which Landsat pixels with successful retrievals have this panel ID?
    panel_num <- u[i]
    w <- which(panel_stats[,70]==panel_num)
    
    if (length(w)>28){
      panel <- panel_stats[w,]
            
      spr_panel <- t(panel[,6:37])
      aut_panel <- t(panel[,38:69])
      
      spr_mean <- apply(spr_panel,1,mean,na.rm=TRUE)
      aut_mean <- apply(aut_panel,1,mean,na.rm=TRUE)
      
      p <- Panel_Analysis(spr_panel,aut_panel,panel_num)
      if(!inherits(p,"try-error")){ 
        print(paste("Panel",i,"was compiled successfully."))
      } else{
        print(paste("Panel",i,"had an error."))
        p <- matrix(NA,1,5)
      }
      
      p <- cbind(p,t(spr_mean),t(aut_mean))
    } else {
      p <- cbind(panel_num,matrix(0,1,68))
    }
  }
  
  #Save panel stats matrix
  setwd(paste('/projectnb/modislc/projects/te_phenology/landsat_stacks/',stack_loc,'/phenology/',sep=""))
  save(panel_stats_avhrr,file="panel_nodist_stats_avhrr")
  
  
  #Generate maps of panel stats
  setwd(paste('/projectnb/modislc/projects/te_phenology/landsat_stacks/',stack_loc,'/PANEL/',sep=""))
  modis_hdr <- raster('mod_id.bip')
  modis_vals <- getValues(modis_hdr)
  
  setwd(paste('/projectnb/modislc/projects/te_phenology/landsat_stacks/',stack_loc,'/PANEL/',sep=""))
  avhrr_hdr <- raster('avhrr_id.bip')
  avhrr_vals <- getValues(avhrr_hdr)
  
  #Match panel IDs with panel stat matrix
  m <- match(modis_vals,panel_stats[,1])
  w <- which(is.na(m)==0) 
  
  trend_spr <- panel_stats[m[w],2]
  pvalue_spr <- panel_stats[m[w],3]
  trend_aut <- panel_stats[m[w],4]
  pvalue_aut <- panel_stats[m[w],5]
  
  trend_spr[pvalue_spr>0.05 | pvalue_spr==0]<-NA
  trend_spr[trend_spr==0]<-NA
  trend_aut[pvalue_aut>0.05 | pvalue_aut==0]<-NA
  trend_aut[trend_aut==0]<-NA
  
  pixID <- matrix(NA,ncell(modis_hdr),1)
  pixID[w] <- trend_spr
  pixID_hdr <- setValues(modis_hdr,pixID)
  
  pixID2 <- matrix(NA,ncell(modis_hdr),1)
  pixID2[w] <- trend_aut
  pixID2_hdr <- setValues(modis_hdr,pixID2)
  
  my.colors = colorRampPalette(c("green4","greenyellow","yellow","gray24","cadetblue1","cadetblue3","navyblue"))
  pdf(h=8,w=6,paste('/projectnb/modislc/projects/te_phenology/figures/',stack_loc,'_nodist_spr_trend.pdf',sep=""))
  plot(pixID_hdr,colNA='black',col=my.colors(150),zlim=c(-0.5,0.5))
  dev.off()
  
  pdf(h=8,w=6,paste('/projectnb/modislc/projects/te_phenology/figures/',stack_loc,'_nodist_aut_trend.pdf',sep=""))
  plot(pixID2_hdr,colNA='black',col=my.colors(150),zlim=c(-0.5,0.5))
  dev.off()
  
  #Match panel IDs with panel stat matrix
  a <- match(avhrr_vals,panel_stats_avhrr[,1])
  w <- which(is.na(a)==0) 
  
  trend_spr <- panel_stats_avhrr[a[w],2]
  pvalue_spr <- panel_stats_avhrr[a[w],3]
  trend_aut <- panel_stats_avhrr[a[w],4]
  pvalue_aut <- panel_stats_avhrr[a[w],5]
  
  trend_spr[pvalue_spr>0.05 | pvalue_spr==0]<-NA
  trend_spr[trend_spr==0]<-NA
  trend_aut[pvalue_aut>0.05 | pvalue_aut==0]<-NA
  trend_aut[trend_aut==0]<-NA
  
  pixID <- matrix(NA,ncell(avhrr_hdr),1)
  pixID[w] <- trend_spr
  pixID_hdr <- setValues(avhrr_hdr,pixID)
  
  pixID2 <- matrix(NA,ncell(avhrr_hdr),1)
  pixID2[w] <- trend_aut
  pixID2_hdr <- setValues(avhrr_hdr,pixID2)
  
  my.colors = colorRampPalette(c("green4","greenyellow","yellow","gray24","cadetblue1","cadetblue3","navyblue"))
  pdf(h=8,w=6,paste('/projectnb/modislc/projects/te_phenology/figures/',stack_loc,'_nodist_spr_trend_avhrr.pdf',sep=""))
  plot(pixID_hdr,colNA='black',col=my.colors(150),zlim=c(-0.5,0.5))
  dev.off()
  
  pdf(h=8,w=6,paste('/projectnb/modislc/projects/te_phenology/figures/',stack_loc,'_nodist_aut_trend_avhrr.pdf',sep=""))
  plot(pixID2_hdr,colNA='black',col=my.colors(150),zlim=c(-0.5,0.5))
  dev.off()
  
  ##########################################################################################
  #MAKE STUDY AREA MAPS
  
  #Include only deciduous pixels with successful retrievals and non-zero panel IDs
  load(paste('/projectnb/modislc/projects/te_phenology/landsat_stacks/',stack_loc,'/phenology/pheno_nodist_mat',sep=""))
  pheno_mat <- pheno_mat[pheno_mat[,7]>=0.9,] 
  pheno_mat[pheno_mat==0]<-NA
  pheno_mat <- pheno_mat[is.na(pheno_mat[,4])==0,]
  
  setwd(paste('/projectnb/modislc/projects/te_phenology/landsat_stacks/',stack_loc,sep=""))
  if (file.exists('EOSD')==1){
    lc <- raster(paste('/projectnb/modislc/projects/te_phenology/landsat_stacks/',stack_loc,'/EOSD/eosd_clip.bip',sep=""))
  } else if (file.exists('NLCD')==1) {
    lc <- raster(paste('/projectnb/modislc/projects/te_phenology/landsat_stacks/',stack_loc,'/NLCD/nlcd_clip.bip',sep=""))
  }
  
  temp <- matrix(NA,ncell(lc),1)
  temp[pheno_mat[,1]] <- pheno_mat[,1]
  pixID <- setValues(lc,temp)
  
  writeRaster(pixID,filename=paste('/projectnb/modislc/projects/te_phenology/geotiffs/',
    stack_loc,'_pixID.tif',sep=""),format='GTiff',overwrite=TRUE)
  
  temp <- matrix(NA,ncell(lc),1)
  temp[pheno_mat[,1]] <- pheno_mat[,8]
  pSPR <- setValues(lc,temp)
  
  writeRaster(pSPR,filename=paste('/projectnb/modislc/projects/te_phenology/geotiffs/',
    stack_loc,'_pSPR_dbf.tif',sep=""),format='GTiff',overwrite=TRUE)
  
  my.colors3 = colorRampPalette(c("red","orange","yellow","green","darkgreen","blue"))
  
  pdf(h=8,w=6,paste('/projectnb/modislc/projects/te_phenology/figures/',stack_loc,'_nodist_spring_LTM.pdf',sep=""))
  plot(pSPR,colNA='black',col=my.colors3(150),
    zlim=c(quantile(temp,0.025,na.rm=TRUE),quantile(temp,0.975,na.rm=TRUE)))
  dev.off()
  
  temp <- matrix(NA,ncell(lc),1)
  temp[pheno_mat[,1]] <- pheno_mat[,9]
  pAUT <- setValues(lc,temp)
  
  writeRaster(pAUT,filename=paste('/projectnb/modislc/projects/te_phenology/geotiffs/',
    stack_loc,'_pAUT.tif',sep=""),format='GTiff',overwrite=TRUE)
  
  my.colors3 = colorRampPalette(c("blue","darkgreen","green","yellow","orange","red"))
  
  pdf(h=8,w=6,paste('/projectnb/modislc/projects/te_phenology/figures/',stack_loc,'_nodist_autumn_LTM.pdf',sep=""))
  plot(pAUT,colNA='black',col=my.colors3(150),
    zlim=c(quantile(temp,0.025,na.rm=TRUE),quantile(temp,0.975,na.rm=TRUE)))
  dev.off()
  
})