#### Perform single panel regression for each sidelap region using 500 m panels

require(rgdal)
library(raster)
library(foreach)
library(iterators)
library(doParallel)
library(zyp)
library(Kendall)
library(plm)
library(lmtest)
library(sandwich)

setwd('/usr3/graduate/emelaas/Code/R/landsat_phenology/')
scenes <- read.table('overlap_scenes.txt',header=FALSE)
scenes <- as.character(scenes[,1])
scenes <- scenes[-c(1,2,7,46,50,59,79,80,82,83,84)]

#Register the parallel backend
registerDoParallel(16)

##########################################################################################
#PANEL ANALYSIS BY PANEL ID (MODIS,AVHRR,ETC.)
Panel_Analysis <- function(spr_panel,aut_panel,panel_num,num_pix){
  
  plm1 <- data.frame(matrix(NA,1,1)); colnames(plm1) <- "coefficients"    
  plm2 <- data.frame(matrix(NA,1,1)); colnames(plm2) <- "coefficients"    
  plm1_pvalue <- matrix(NA,1,1)
  plm2_pvalue <- matrix(NA,1,1)
  
  #If there are at least ten deciduous pixels in a panel, make a panel
  if (ncol(spr_panel)>1){
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
  if (ncol(aut_panel)>1){
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


#Loop through each panel ID and perform panel analysis
panel_stats <- foreach(i = 1:length(scenes), .combine = rbind) %dopar% {
  
  print(i)
  
  stack_loc <- scenes[i]
  
  setwd(paste('/projectnb/modislc/projects/te_phenology/landsat_stacks/',stack_loc,'/phenology/',sep=""))
  load(file="panel_nodist_stats_MK-TS_50_dbf")
  
  #Which Landsat pixels with successful retrievals have this panel ID 
  panel_num <- 1
  w <- which(is.na(panel_stats[,2])==0)
  num_pix <- length(w)
  
  #If more than 1, count number of SOS/EOS observations for each year with
  #intent of excluding years with fewer than 50% of available observations
  if (length(w)>5){
    tmp <- panel_stats[w,]
    tmp[tmp==0] <- NA
    tmp[is.na(tmp)==0] <- 1
    tmp_colsum <- colSums(tmp,na.rm=TRUE)
    
    panel <- panel_stats[w,]
    panel[,which(tmp_colsum < num_pix/10)] <- NA
    
#     w <- which(is.na(panel[,5])==1)
#     panel[w,5] <- 1
#     avhrr_pix <- as.numeric(names(table(panel[,5]))[which.max(table(panel[,5]))])
    
    spr_panel <- t(panel[,9:40])
    aut_panel <- t(panel[,41:72])
    
    spr_mean <- apply(spr_panel,1,mean,na.rm=TRUE)
    aut_mean <- apply(aut_panel,1,mean,na.rm=TRUE)
    
    num_spr <- length(which(is.na(spr_mean)==0))
    num_aut <- length(which(is.na(aut_mean)==0))
    
    if (num_spr > 10 & num_aut > 10) {
      
      y <- as.numeric(spr_mean); x <- seq(1,32,1)
      mk <- MannKendall(y)
      mk_spr <- as.numeric(mk$sl)
      z <- zyp.sen(y~x)
      z_spr <- z$coefficients[2] 
      
      y <- as.numeric(aut_mean); x <- seq(1,32,1)
      mk <- MannKendall(y)
      mk_aut <- as.numeric(mk$sl)
      z <- zyp.sen(y~x)
      z_aut <- z$coefficients[2] 
      
      p <- Panel_Analysis(spr_panel,aut_panel,panel_num,num_pix)
      if(!inherits(p,"try-error")){ 
        print(paste("Panel",i,"was compiled successfully."))
      } else{
        print(paste("Panel",i,"had an error."))
        p <- matrix(NA,1,5)
      }
      
      avhrr_pix <- 1
      
      p <- cbind(p,num_pix,num_spr,num_aut,t(spr_mean),t(aut_mean),avhrr_pix,z_spr,mk_spr,z_aut,mk_aut)
    } else {
      p <- cbind(panel_num,NA,NA,NA,NA,num_pix,num_spr,num_aut,matrix(NA,1,64),1,NA,NA,NA,NA)
    }
  } else {
    num_spr <- NA
    num_aut <- NA
    p <- cbind(panel_num,NA,NA,NA,NA,num_pix,num_spr,num_aut,matrix(NA,1,64),1,NA,NA,NA,NA)
  }
}

#Save panel stats matrix
setwd('/projectnb/modislc/users/emelaas/scratch15/NASA_TE/R_files')
save(panel_stats,file="panel_sidelap_stats_MK-TS_50_dbf")