#### Calculate sidelap-specific statistics for mapping 

setwd('/usr3/graduate/emelaas/Code/R/landsat_phenology/GRL')
scenes <- read.table('overlap_scenes.txt',header=FALSE)
scenes <- as.character(scenes[,1])

scenes2 <- scenes
scenes2[22] <- 'mws_green_bay'
scenes2[68] <- 'oa_springfield'
scenes2[74] <- 'bp_kakwa'
scenes2[75] <- 'bp_williston'
L1 <- substring(scenes2,1,3)
L1_col <- L1
L1_col[which(L1 %in% c('ah_','ma_','mwp','oa_','sup'))] <- 'ETF'
L1_col[which(L1 %in% c('bp_','hp_','mws','ss_'))] <- 'NF'
L1_col[which(L1 %in% c('tc_','tp_','ts_'))] <- 'T'

loc <- NA
lat <- NA
lon <- NA
npan_spr <- NA
npan_aut <- NA
npan_tot <- NA
npix_tot <- NA
nyrs_spr <- NA
nyrs_aut <- NA
perc_spr_sig_adv <- NA 
perc_spr_sig_del <- NA
perc_aut_sig_adv <- NA
perc_aut_sig_del <- NA

k <- 1

for (i in 1:length(scenes)){
  print(i)
  
  setwd(paste('/projectnb/modislc/projects/te_phenology/landsat_stacks/',scenes[i],'/phenology',sep=""))
  if (file.exists('pheno_nodist_mat')==1 & file.exists('panel_nodist_stats_MK-TS_50_dbf')==1){
    load('pheno_nodist_mat')
    
    if (nrow(pheno_mat)>1000){
      
      setwd(paste('/projectnb/modislc/projects/te_phenology/landsat_stacks/',scenes[i],sep=""))
      if (file.exists('EOSD')==1){
        lc <- raster(paste('/projectnb/modislc/projects/te_phenology/landsat_stacks/',
          scenes[i],'/EOSD/eosd_clip.bip',sep=""))
        w <- c(221,222,223,231,232,233)
      } 
      if (file.exists('NLCD')==1) {
        lc <- raster(paste('/projectnb/modislc/projects/te_phenology/landsat_stacks/',
          scenes[i],'/NLCD/nlcd_clip.bip',sep=""))
        w <- c(41,43)
      }
      pheno_mat2 <- pheno_mat[which(pheno_mat[,2] %in% w),]
      npix_tot[k] <- nrow(pheno_mat2)
      
      assign(paste('pSPR',k,sep=""),pheno_mat2[,8])
      assign(paste('pAUT',k,sep=""),pheno_mat2[,9])
      assign(paste('pNOBS',k,sep=""),pheno_mat2[,6])
      remove(pheno_mat2,pheno_mat)
      
      setwd(paste('/projectnb/modislc/projects/te_phenology/landsat_stacks/',scenes[i],'/phenology',sep=""))
      load('panel_nodist_stats_MK-TS_50_dbf')
      w1 <- which(panel_stats[,75] <= 1 & panel_stats[,75] != 0)
      w1a <- which(panel_stats[,75] <= 0.05 & panel_stats[,75] != 0)
      panel_stats1 <- panel_stats[w1a,]
      assign(paste('tSPR',k,sep=""),panel_stats1[,74])
      
      w2 <- which(panel_stats[,77] <= 1 & panel_stats[,77] != 0)
      w2a <- which(panel_stats[,77] <= 0.05 & panel_stats[,77] != 0)
      panel_stats2 <- panel_stats[w2a,]
      if (length(w2a) == 1){
        panel_stats2 <- t(matrix(panel_stats2))
      }
      
      assign(paste('tAUT',k,sep=""),panel_stats2[,76])
      
      npan_tot[k] <- nrow(panel_stats)
      remove(panel_stats)
      
      #Obtain lat/lon coordinates of overlap region's centroid
      o <- readOGR(paste('/projectnb/modislc/projects/te_phenology/landsat_stacks/',
        scenes[i],'/SHP/',sep=""),scenes[i])
      o_proj <- spTransform(o,CRS("+proj=longlat +datum=WGS84"))
      lat[k] <- coordinates(o_proj)[2]
      lon[k] <- coordinates(o_proj)[1]
      
      loc[k] <- scenes[i]
      
      npan_spr[k] <- length(w1)
      npan_aut[k] <- length(w2)
      
      nyrs_spr[k] <- mean(panel_stats1[,7])
      nyrs_aut[k] <- mean(panel_stats1[,8])
      
      perc_spr_sig_adv[k] <- length(which(panel_stats1[,74]<0 & panel_stats1[,75]<0.05))/npan_spr[k]
      perc_spr_sig_del[k] <- length(which(panel_stats1[,74]>0 & panel_stats1[,75]<0.05))/npan_spr[k]
      perc_aut_sig_adv[k] <- length(which(panel_stats1[,76]<0 & panel_stats1[,77]<0.05))/npan_aut[k]
      perc_aut_sig_del[k] <- length(which(panel_stats1[,76]>0 & panel_stats1[,77]<0.05))/npan_aut[k]
      
      k <- k+1
    }
  }
}

spr_list <- lapply(ls(pattern='pSPR'),get)
aut_list <- lapply(ls(pattern='pAUT'),get)
spr2_list <- lapply(ls(pattern='tSPR'),get)
aut2_list <- lapply(ls(pattern='tAUT'),get)
nobs_list <- lapply(ls(pattern='pNOBS'),get)

setwd('/projectnb/modislc/users/emelaas/scratch15/NASA_TE/R_files/')
save('spr_list','aut_list','spr2_list','aut2_list','loc','npan_spr','npan_aut','npan_tot','npix_tot','nyrs_spr','nyrs_aut',
  file='scene_panel_mats_all')

spr_median <- unlist(lapply(spr_list,median,na.rm=TRUE))
aut_median <- unlist(lapply(aut_list,median,na.rm=TRUE))
spr2_median <- unlist(lapply(spr2_list,median,na.rm=TRUE))
aut2_median <- unlist(lapply(aut2_list,median,na.rm=TRUE))
nobs_median <- unlist(lapply(nobs_list,median,na.rm=TRUE))

index <- unlist(lapply(ls(pattern='pSPR'), function(x) as.numeric(substr(x,5,6))))
lat <- lat[index]
lon <- lon[index]
npan_spr <- npan_spr[index]
npan_aut <- npan_aut[index]
npan_tot <- npan_tot[index]
npix_tot <- npix_tot[index]
nyrs_spr <- nyrs_spr[index]
nyrs_aut <- nyrs_aut[index]
perc_spr_sig_adv <- perc_spr_sig_adv[index]
perc_spr_sig_del <- perc_spr_sig_del[index]
perc_aut_sig_adv <- perc_aut_sig_adv[index]
perc_aut_sig_del <- perc_aut_sig_del[index]

index <- data.frame(cbind(index,spr_median,aut_median,spr2_median,aut2_median,nobs_median,lat,lon,npan_spr,npan_aut,npan_tot,npix_tot,
  perc_spr_sig_adv,perc_spr_sig_del,perc_aut_sig_adv,perc_aut_sig_del,nyrs_spr,nyrs_aut))
colnames(index) <- c('ind','sLTM','aLTM','strend','atrend','nobs','lat','lon','npan_spr','npan_aut','npan_tot','npix_tot',
  'perc_spr_sig_adv','perc_spr_sig_del','perc_aut_sig_adv','perc_aut_sig_del','nyrs_spr','nyrs_aut')
index <- index[order(index[,1]),]

setwd('/projectnb/modislc/users/emelaas/scratch15/NASA_TE/R_files/')
save('loc','index',file='meta_analysis_nodist_all7')
