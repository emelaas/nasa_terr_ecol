require(rgdal)
require(raster)
require(RColorBrewer)
require(classInt)
require(vioplot)
require(ncdf4)
require(lmodel2)
require(zyp)
require(Kendall)

sidelaps <- readOGR('/projectnb/modislc/users/emelaas/scratch15/NASA_TE/SHP/','NASA_TE_sidelaps2')
usa_shp <- readOGR('/projectnb/modislc/users/emelaas/scratch32/DAYMET/USA_Boundary/states_21basic/states.dbf','states')
can_shp <- readOGR('/projectnb/modislc/users/emelaas/scratch32/DAYMET/Canada/','Canada')
mex_shp <- readOGR('/projectnb/modislc/users/emelaas/scratch32/DAYMET/Mexico/mexstates/','mexstates')
ecoreg_shp <- readOGR('/projectnb/modislc/projects/te_phenology/ecoregions/only_l1_eco/','only_l1_eco')
usa_shp_proj <- spTransform(usa_shp, CRS("+proj=laea +lat_0=50 +lon_0=-100 +x_0=0 
  +y_0=0 +ellps=WGS84 +units=m +no_defs"))
can_shp_proj <- spTransform(can_shp, CRS("+proj=laea +lat_0=50 +lon_0=-100 +x_0=0 
  +y_0=0 +ellps=WGS84 +units=m +no_defs"))
mex_shp_proj <- spTransform(mex_shp, CRS("+proj=laea +lat_0=50 +lon_0=-100 +x_0=0 
  +y_0=0 +ellps=WGS84 +units=m +no_defs"))
ecoreg_shp_proj <- spTransform(ecoreg_shp, CRS("+proj=laea +lat_0=50 +lon_0=-100 +x_0=0 
  +y_0=0 +ellps=WGS84 +units=m +no_defs"))

setwd('/usr3/graduate/emelaas/Code/R/landsat_phenology/GRL')
scenes <- read.table('overlap_scenes.txt',header=FALSE)
scenes <- as.character(scenes[,1])
scenes <- scenes[-c(1,2,7,46,50,59,79,80,82,83,84)]

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


#### Figure 1 ####

setwd('/projectnb/modislc/users/emelaas/scratch15/NASA_TE/R_files/')
load(file='meta_analysis_nodist_all7')
load(file='panel_sidelap_stats_MK-TS_50_dbf')
load('preseason_temp_cor_gdd')

#Generate map showing percentage of 500 m patches with significant spr/aut trends in each sidelap
cuts_s <- classIntervals(panel_stats[,2],dataPrecision=3,style='fixed',
  fixedBreaks=c(-0.8,-0.4,-0.2,-0.1,0.1,0.2,0.4,0.8))
plotclr <- rev(brewer.pal(9,'RdBu'))
colcode_s <- findColours(cuts_s,plotclr)

cuts_m <- classIntervals(index$sLTM,dataPrecision=3,style='fixed',
  fixedBreaks=c(90,100,110,120,130,140,150,160,170,180))
plotclr <- brewer.pal(9,'RdBu')
colcode_m <- findColours(cuts_m,plotclr)

index$ST3 <- ST3*100

#Plot sensitivity between mean preseason temperature and SOS (ST; days/degC)
cuts_ST3 <- classIntervals(index$ST3,dataPrecision=3,style='fixed',
  fixedBreaks=c(-12,-9,-6,-3))
plotclr <- rev(brewer.pal(5,'OrRd'))
colcode_ST3 <- findColours(cuts_ST3,plotclr)

cuts_z <- classIntervals(Z,dataPrecision=3,style='fixed',
  fixedBreaks=c(-8,-4,-2,-1,1,2,4,8))
plotclr <- rev(brewer.pal(7,'RdBu'))
colcode_z <- findColours(cuts_z,plotclr)

x11(h=6,w=40/3)
#pdf(h=6,w=40/3,'/projectnb/modislc/projects/te_phenology/figures/Storyboard_4/Figure1.pdf')
par(mfrow=c(1,2),mar=c(2.6,3,2,1))
plot(sidelaps,pch=21,cex=0.01,col=colcode_m,bg=colcode_m,
  main='SOS Mean (DOY)',cex.main=1)
lim <- par("usr")
rect(lim[1],lim[3],lim[2],lim[4],border="lightblue",col="lightblue")
plot(usa_shp_proj,axes=FALSE,add=TRUE,col='gray70',lwd=0.1,border='gray70')
plot(can_shp_proj,axes=FALSE,add=TRUE,col='gray70',lwd=0.1,border='gray70')
plot(mex_shp_proj,axes=FALSE,add=TRUE,col='gray70',lwd=0.1,border='gray70')
plot(ecoreg_shp,axes=FALSE,add=TRUE,col='white',lwd=0.2,border='gray40')
plot(sidelaps,pch=21,cex=1.6,col=colcode_m,bg=colcode_m,add=TRUE,lwd=0.5)
mtext('a',side=3,line=0.5,cex=1,adj=0,font=2)
legend("bottomright",title='',cex=0.9,legend=names(attr(colcode_m,'table')),
  fill=attr(colcode_m,'palette'),bty='n')
rect(grconvertX(0.025, from='npc'), grconvertY(0.025, from='npc'),
  grconvertX(0.325, from='npc'), grconvertY(0.55, from='npc'), 
  col="gray90", border=NA)

plot(sidelaps,pch=21,cex=0.01,col=colcode_s,bg=colcode_s,
  main=substitute(bold(paste('SOS Trend - 1984-2013 (days ', yr ^ -1, ")"), list(yr = ""))),cex.main=1)
lim <- par("usr")
rect(lim[1],lim[3],lim[2],lim[4],border="lightblue",col="lightblue")
plot(usa_shp_proj,axes=FALSE,add=TRUE,col='gray70',lwd=0.1,border='gray70')
plot(can_shp_proj,axes=FALSE,add=TRUE,col='gray70',lwd=0.1,border='gray70')
plot(mex_shp_proj,axes=FALSE,add=TRUE,col='gray70',lwd=0.1,border='gray70')
plot(ecoreg_shp,axes=FALSE,add=TRUE,col='white',lwd=0.2,border='gray40')
plot(sidelaps,pch=21,cex=1.6,col=colcode_s,bg=colcode_s,add=TRUE,lwd=0.5)
mtext('b',side=3,line=0.5,cex=1,adj=0,font=2)
legend("bottomright",title='',cex=0.9,legend=names(attr(colcode_s,'table')),
  fill=attr(colcode_s,'palette'),bty='n')
rect(grconvertX(0.025, from='npc'), grconvertY(0.025, from='npc'),
  grconvertX(0.325, from='npc'), grconvertY(0.55, from='npc'), 
  col="gray90", border=NA)

plotdim <- par("plt")
par(fig = c(plotdim[1]/2.3,plotdim[2]/4.7,plotdim[3]*1.2,plotdim[4]/1.65), new=TRUE, bty='n')
boxplot(index$sLTM~L1_col,col='white', main='Level 1 Ecoregions',cex.main=0.8)
abline(h=0,lty=2)

par(fig = c(plotdim[1]/2.3+(1/2),plotdim[2]/4.7+(1/2),plotdim[3]*1.2,plotdim[4]/1.65), new=TRUE, bty='n')
boxplot(panel_stats[,2]~L1_col,col='white')
abline(h=0,lty=2)


dev.off()  


#### Figure 2 ####

setwd('/usr3/graduate/emelaas/Code/R/landsat_phenology/GRL')
scenes <- read.table('overlap_scenes.txt',header=FALSE)
scenes <- as.character(scenes[,1])
scenes <- scenes[-c(1,2,7,46,50,59,79,80,82,83,84)]

scenes2 <- scenes
scenes2[22] <- 'mws_green_bay'
scenes2[68] <- 'oa_springfield'
scenes2[74] <- 'bp_kakwa'
scenes2[75] <- 'bp_williston'
L2 <- substring(scenes2,1,3)
#L2_tab <- names(table(L2))
L2_tab <- c('tc_','tp_','ts_','hp_','bp_','ss_','mws','ah_','oa_','mwp','ma_','sup')
L2_col <- L2_tab
L2_col[which(L2_tab %in% c('ah_','ma_','mwp','oa_','sup'))] <- 'blue'
L2_col[which(L2_tab %in% c('bp_','hp_','mws','ss_'))] <- 'orange'
L2_col[which(L2_tab %in% c('tc_','tp_','ts_'))] <- 'purple'
L2_nam <- c('Taiga Cordillera','Taiga Plain','Taiga Shield','Hudson Plain','Boreal Plain','Softwood Shield',
  'Mixed Wood Shield','Atlantic Highlands','Ozark/Ouachita Appalachian Forest','Mixed Wood Plains',
  'Mid Atlantic Coastal Plain','Southeastern USA Plains')
panel <- c('a','b','c','d','e','f','g','h','i','j','k','l')

x11(h=8,w=11)
#pdf(h=8,w=11,'/projectnb/modislc/projects/te_phenology/figures/Storyboard_3/Figure2.pdf')
par(mar=c(3.3,2,1,2))
par(mfrow=c(3,4))
for (i in 1:length(L2_tab)){
  w <- which(L2_tab[i]==L2)
  
  o.SPR_anom_all <- matrix(NA,1,30)
  for (j in 1:length(w)){
    setwd(paste('/projectnb/modislc/projects/te_phenology/landsat_stacks/',scenes[w[j]],'/phenology',sep=""))
    load('panel_nodist_stats_MK-TS_50_dbf')
    o.SPR_pixmean <- rowMeans(panel_stats[,11:40],na.rm=TRUE)
    o.SPR_anom <- panel_stats[,11:40]-replicate(30,o.SPR_pixmean)
    o.SPR_anom_all <- rbind(o.SPR_anom_all,o.SPR_anom)
  }
  
  o.SPR_mean <- colMeans(o.SPR_anom_all,na.rm=TRUE)
  o.SPR_std <- apply(o.SPR_anom_all,2,sd,na.rm=TRUE)
  o.SPR_lx <- o.SPR_mean-o.SPR_std
  o.SPR_ux <- o.SPR_mean+o.SPR_std
  
  x <- c(seq(1984,2013,1),seq(2013,1984,-1))
  y <- c(o.SPR_ux,rev(o.SPR_lx))
  plot(seq(1984,2013,1),o.SPR_mean,main="",
    ylim=c(-18,18),xlab="",ylab="",cex=0,bty='n',yaxt='n')
  if (L2_col[i]=='blue'){
    polygon(x,y,col=rgb(0,0,1,0.5),border=rgb(0,0,1,1),lty=0,lwd=2)
  } else if (L2_col[i]=='orange') {
    polygon(x,y,col=rgb(1,(165/255),0,0.5),border=rgb(1,(165/255),0,1),lty=0,lwd=2)
  } else {
    polygon(x,y,col=rgb((128/255),0,(128/255),0.35),border=rgb((128/255),0,(128/255),1),lty=0,lwd=2)
  }
  
  lines(seq(1984,2013,1),o.SPR_mean,lwd=5,col=L2_col[i])
  
  mtext(panel[i],side=2,las=1,adj=1,padj=-8,font=2)
  grid(nx=NA,ny=NULL,col='gray')
  text(par("usr")[1]-2,-11,srt=90,adj = 0,labels = "SOS Anomaly (days)",cex=1,xpd = TRUE)
  text(1984,seq(-15,15,5),seq(-15,15,5),cex=1)
  text(1996.5,par("usr")[3]-7,srt=0,adj = 0,labels = "Year",cex=1,xpd = TRUE)
  text(1990,par("usr")[4],srt=0,adj = c(0.25,NA),labels = L2_nam[i],cex=1,xpd = TRUE)
  
  x <- seq(1984,2013,1)
  y <- o.SPR_mean
  mk <- MannKendall(y)
  mk_spr <- as.numeric(mk$sl)
  x <- x[which(is.na(y)==0)]
  y <- y[which(is.na(y)==0)]
  fit <- zyp.sen(y~x)
  segments(x0=1984,y0=fit$coefficients[2]*1984+fit$coefficients[1],x1=2013,
    y1=fit$coefficients[2]*2013+fit$coefficients[1])
  text(1999,-17,labels=paste('trend = ',signif(fit$coefficients[2],2),', p = ',signif(mk_spr,2),sep=''))
}

dev.off()


#### Figure 3 ####

setwd('/projectnb/modislc/users/emelaas/scratch15/NASA_TE/R_files/')
load(file='meta_analysis_nodist_all7')
load(file='panel_sidelap_stats_MK-TS_50_dbf')
load('preseason_temp_cor_gdd')

index$ST3 <- ST3*100

#Plot sensitivity between mean preseason temperature and SOS (ST; days/degC)
cuts_ST3 <- classIntervals(index$ST3,dataPrecision=3,style='fixed',
  fixedBreaks=c(-12,-9,-6,-3))
plotclr <- rev(brewer.pal(5,'OrRd'))
colcode_ST3 <- findColours(cuts_ST3,plotclr)

cuts_z <- classIntervals(Z,dataPrecision=3,style='fixed',
  fixedBreaks=c(-8,-4,-2,-1,1,2,4,8))
plotclr <- rev(brewer.pal(7,'RdBu'))
colcode_z <- findColours(cuts_z,plotclr)

#x11(h=6,w=40/3)
pdf(h=6,w=40/3,'/projectnb/modislc/projects/te_phenology/figures/Storyboard_4/Figure3.pdf')
par(mfrow=c(1,2),mar=c(2.6,3,2,1))
plot(sidelaps,pch=21,cex=0.01,col=colcode_ST3,bg=colcode_ST3,
  main=expression(bold(paste(S[T],' (days/100 AGDD)',sep=""))),cex.main=1)
lim <- par("usr")
rect(lim[1],lim[3],lim[2],lim[4],border="lightblue",col="lightblue")
plot(usa_shp_proj,axes=FALSE,add=TRUE,col='gray70',lwd=0.1,border='gray70')
plot(can_shp_proj,axes=FALSE,add=TRUE,col='gray70',lwd=0.1,border='gray70')
plot(mex_shp_proj,axes=FALSE,add=TRUE,col='gray70',lwd=0.1,border='gray70')
plot(ecoreg_shp,axes=FALSE,add=TRUE,col='white',lwd=0.2,border='gray40')
plot(sidelaps,pch=21,cex=1.6,col=colcode_ST3,bg=colcode_ST3,add=TRUE,lwd=0.5)
mtext('a',side=3,line=0.5,cex=1,adj=0,font=2)
legend("bottomright",title='',cex=1.1,legend=names(attr(colcode_ST3,'table')),
  fill=attr(colcode_ST3,'palette'),bty='n')
rect(grconvertX(0.025, from='npc'), grconvertY(0.025, from='npc'),
  grconvertX(0.325, from='npc'), grconvertY(0.575, from='npc'), 
  col="gray90", border=NA)

plot(sidelaps,pch=21,cex=0.01,col=colcode_z,bg=colcode_z,
  main='Optimal Preseason Temperature Trend (AGDD/yr)',cex.main=1)
lim <- par("usr")
rect(lim[1],lim[3],lim[2],lim[4],border="lightblue",col="lightblue")
plot(usa_shp_proj,axes=FALSE,add=TRUE,col='gray70',lwd=0.1,border='gray70')
plot(can_shp_proj,axes=FALSE,add=TRUE,col='gray70',lwd=0.1,border='gray70')
plot(mex_shp_proj,axes=FALSE,add=TRUE,col='gray70',lwd=0.1,border='gray70')
plot(ecoreg_shp,axes=FALSE,add=TRUE,col='white',lwd=0.2,border='gray40')
plot(sidelaps,pch=21,cex=1.6,col=colcode_z,bg=colcode_z,add=TRUE,lwd=0.5)
mtext('b',side=3,line=0.5,cex=1.1,adj=0,font=2)
legend("bottomright",title='',cex=1.1,legend=names(attr(colcode_z,'table')),
  fill=attr(colcode_z,'palette'),bty='n')
rect(grconvertX(0.025, from='npc'), grconvertY(0.025, from='npc'),
  grconvertX(0.325, from='npc'), grconvertY(0.55, from='npc'), 
  col="gray90", border=NA)

plotdim <- par("plt")
par(fig = c(plotdim[1]/2.3,plotdim[2]/4.7,plotdim[3]*1.15,plotdim[4]/1.565), new=TRUE, bty='n')
boxplot(index$ST3~L1_col,col='white', main='Level 1 Ecoregions',cex.main=0.8)

par(fig = c(plotdim[1]/2.3+(1/2),plotdim[2]/4.7+(1/2),plotdim[3]*1.2,plotdim[4]/1.65), new=TRUE, bty='n')
boxplot(Z~L1_col,col='white')
abline(h=0,lty=2)

dev.off()


#### Figure 4 ####

setwd('/projectnb/modislc/users/emelaas/scratch15/NASA_TE/R_files/')
load('preseason_temp_cor_gdd')

setwd('/usr3/graduate/emelaas/Code/R/landsat_phenology/GRL')
scenes <- read.table('overlap_scenes.txt',header=FALSE)
scenes <- as.character(scenes[,1])
scenes <- scenes[-c(1,2,7,46,50,59,79,80,82,83,84)]
scenes2 <- scenes
scenes2[68] <- 'oa_springfield'
scenes2[74] <- 'bp_kakwa'
scenes2[75] <- 'bp_williston'
L2 <- substring(scenes2,1,3)
#L2_tab <- names(table(L2))
L2_tab <- c('tc_','tp_','ts_','hp_','bp_','ss_','mws','ah_','oa_','mwp','ma_','sup')
L2_col <- L2_tab
L2_col[which(L2_tab %in% c('ah_','ma_','mwp','oa_','sup'))] <- 'blue'
L2_col[which(L2_tab %in% c('bp_','hp_','mws','ss_'))] <- 'orange'
L2_col[which(L2_tab %in% c('tc_','tp_','ts_'))] <- 'purple'
L2_nam <- c('Taiga Cordillera','Taiga Plain','Taiga Shield','Hudson Plain','Boreal Plain','Softwood Shield',
  'Mixed Wood Shield','Atlantic Highlands','Ozark/Ouachita Appalachian Forest','Mixed Wood Plains',
  'Mid Atlantic Coastal Plain','Southeastern USA Plains')
panel <- c('a','b','c','d','e','f','g','h','i','j','k','l')

all.obs.SPR.anom.mean <- matrix(NA,12,32)
all.obs.SPR.anom.sd <- matrix(NA,12,32)
all.pred.SPR.anom.mean <- matrix(NA,12,32)
all.pred.SPR.anom.sd <- matrix(NA,12,32)
for (i in 1:length(L2_tab)){
  wL2 <- which(L2 %in% L2_tab[i])
  
  L2.obs.SPR.anom <- matrix(NA,1,32)
  L2.pred.SPR.anom <- matrix(NA,1,32)
  for (j in 1:length(wL2)){    
    print(c(i,j))
    
    stack_loc <- scenes[wL2[j]]
    
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
      N <- w[wL2[j]]*5
      
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
    
    obs.SPR.mean <- replicate(ncol(obs.SPR),apply(obs.SPR,1,mean,na.rm=TRUE))
    obs.SPR.anom <- obs.SPR-obs.SPR.mean
    L2.obs.SPR.anom <- rbind(L2.obs.SPR.anom,obs.SPR.anom)
    
    pred.SPR.mean <- replicate(ncol(pred.SPR),apply(pred.SPR,1,mean,na.rm=TRUE))
    pred.SPR.anom <- pred.SPR-pred.SPR.mean
    L2.pred.SPR.anom <- rbind(L2.pred.SPR.anom,pred.SPR.anom)
  }
  
  all.obs.SPR.anom.mean[i,] <- apply(L2.obs.SPR.anom,2,mean,na.rm=TRUE)
  all.obs.SPR.anom.sd[i,] <- apply(L2.obs.SPR.anom,2,sd,na.rm=TRUE)
  all.pred.SPR.anom.mean[i,] <- apply(L2.pred.SPR.anom,2,mean,na.rm=TRUE)
  all.pred.SPR.anom.sd[i,] <- apply(L2.pred.SPR.anom,2,sd,na.rm=TRUE)
}

x11(h=8,w=11)
#pdf(h=8,w=11,'/projectnb/modislc/projects/te_phenology/figures/Storyboard_3/tmp11.pdf')
par(mar=c(3.3,3,2,3))
par(mfrow=c(3,4))
for (i in 1:12){
  x <- all.pred.SPR.anom.mean[i,]
  xsd <- all.pred.SPR.anom.sd[i,]
  y <- all.obs.SPR.anom.mean[i,]
  ysd <- all.obs.SPR.anom.sd[i,]
  plot(x, y, xlim = c(-200,200), ylim = c(-20,20), pch = 21, col = 'black', bg = 'green', cex = 0 ,xlab = "",ylab = "",bty = 'n')
  arrows(x - xsd, y, x + xsd, y, length = 0.01, angl = 90, code = 3, col = 'gray50')
  arrows(x, y - ysd, x, y + ysd, length = 0.01, angl = 90, code = 3, col = 'gray50')
  points(x, y, xlim = c(-200,200), ylim = c(-20,20), pch = 21, col = 'black', bg = L2_col[i], cex = 1.45,xlab = "",ylab = "",bty = 'n')
  
  mtext(panel[i],side=2,las=1,adj=1,padj=-8.25,font=2)
  grid(nx=NULL,ny=NULL,col='gray')
  text(par("usr")[1]-70,-15,srt=90,adj = 0,labels = "Mean SOS Anomaly (days)",cex=1,xpd = TRUE)
  text(-100,par("usr")[3]-8,srt=0,adj = 0,labels = "Mean GDD Anomaly",cex=1,xpd = TRUE)
  text(-75, par("usr")[4] + 2,srt = 0,adj = c(0.25,NA),labels = L2_nam[i], cex = 1.1,xpd = TRUE, font = 1)
  
  lm1 <- lm(y~x)
  abline(lm1, col = 'gray60')
  text(100,15,labels = paste('Slope = \n',signif(lm1$coefficients[2],3)*100,'d/100 GDD',''))
}
dev.off()
  

#### Figure S1 ####

setwd('/projectnb/modislc/users/emelaas/scratch15/NASA_TE/R_files/')
load(file='meta_analysis_nodist_all7')

#Plot number of cloud-free observations at each site
cuts_nobs <- classIntervals(index$nobs,dataPrecision=3,style='fixed',
  fixedBreaks=c(100,200,300,400,500,600,700,800))
plotclr <- rev(brewer.pal(9,'PuOr'))
colcode_nobs <- findColours(cuts_nobs,plotclr)

cuts_pdec <- classIntervals(index$npan_spr/index$npan_tot,dataPrecision=3,style='fixed',
  fixedBreaks=c(0,0.2,0.4,0.6,0.8,1))
plotclr <- brewer.pal(9,'YlGn')
colcode_pdec <- findColours(cuts_pdec,plotclr)

cuts_nyrs <- classIntervals(index$nyrs_spr,dataPrecision=3,style='fixed',
  fixedBreaks=c(14,16,18,20,22,24))
plotclr <- rev(brewer.pal(9,'RdYlBu'))
colcode_nyrs <- findColours(cuts_nyrs,plotclr)

#x11(h=6,w=20)
pdf(h=6,w=20,'/projectnb/modislc/projects/te_phenology/figures/Storyboard_3/tmp2.pdf')
#pdf(h=6,w=20,'/projectnb/modislc/projects/te_phenology/figures/Storyboard_3/SM_figure1.pdf')
par(mfrow=c(1,3),mar=c(2.6,3,2,1))
plot(sidelaps,pch=21,cex=0.01,col=colcode_s,bg=colcode_s,
  main='Median No. Clear Observations')
lim <- par("usr")
rect(lim[1],lim[3],lim[2],lim[4],border="lightblue",col="lightblue")
plot(usa_shp_proj,axes=FALSE,add=TRUE,col='gray70',lwd=0.1,border='gray70')
plot(can_shp_proj,axes=FALSE,add=TRUE,col='gray70',lwd=0.1,border='gray70')
plot(mex_shp_proj,axes=FALSE,add=TRUE,col='gray70',lwd=0.1,border='gray70')
plot(ecoreg_shp,axes=FALSE,add=TRUE,col='white',lwd=0.2,border='gray40')
plot(sidelaps,pch=21,cex=1.6,col=colcode_nobs,bg=colcode_nobs,add=TRUE,lwd=0.5)
mtext('a',side=3,line=0.5,cex=1,adj=0,font=2)
legend("bottomright",title='',cex=1.4,legend=names(attr(colcode_nobs,'table')),
  fill=attr(colcode_nobs,'palette'),bty='n')
rect(grconvertX(0.025, from='npc'), grconvertY(0.025, from='npc'),
  grconvertX(0.325, from='npc'), grconvertY(0.55, from='npc'), 
  col="gray90", border=NA)

plot(sidelaps,pch=21,cex=0.01,col=colcode_s,bg=colcode_s,
  main='Fraction of 500 m Panels with Deciduous/Mixed Forest')
lim <- par("usr")
rect(lim[1],lim[3],lim[2],lim[4],border="lightblue",col="lightblue")
plot(usa_shp_proj,axes=FALSE,add=TRUE,col='gray70',lwd=0.1,border='gray70')
plot(can_shp_proj,axes=FALSE,add=TRUE,col='gray70',lwd=0.1,border='gray70')
plot(mex_shp_proj,axes=FALSE,add=TRUE,col='gray70',lwd=0.1,border='gray70')
plot(ecoreg_shp,axes=FALSE,add=TRUE,col='white',lwd=0.2,border='gray40')
plot(sidelaps,pch=21,cex=1.6,col=colcode_pdec,bg=colcode_pdec,add=TRUE,lwd=0.5)
mtext('b',side=3,line=0.5,cex=1,adj=0,font=2)
legend("bottomright",title='',cex=1.4,legend=names(attr(colcode_pdec,'table')),
  fill=attr(colcode_pdec,'palette'),bty='n')
rect(grconvertX(0.025, from='npc'), grconvertY(0.025, from='npc'),
  grconvertX(0.325, from='npc'), grconvertY(0.55, from='npc'), 
  col="gray90", border=NA)

plot(sidelaps,pch=21,cex=0.01,col=colcode_s,bg=colcode_s,
  main='Median No. Detected SOS Years')
lim <- par("usr")
rect(lim[1],lim[3],lim[2],lim[4],border="lightblue",col="lightblue")
plot(usa_shp_proj,axes=FALSE,add=TRUE,col='gray70',lwd=0.1,border='gray70')
plot(can_shp_proj,axes=FALSE,add=TRUE,col='gray70',lwd=0.1,border='gray70')
plot(mex_shp_proj,axes=FALSE,add=TRUE,col='gray70',lwd=0.1,border='gray70')
plot(ecoreg_shp,axes=FALSE,add=TRUE,col='white',lwd=0.2,border='gray40')
plot(sidelaps,pch=21,cex=1.6,col=colcode_nyrs,bg=colcode_nyrs,add=TRUE,lwd=0.5)
mtext('c',side=3,line=0.5,cex=1,adj=0,font=2)
legend("bottomright",title='',cex=1.4,legend=names(attr(colcode_nyrs,'table')),
  fill=attr(colcode_nyrs,'palette'),bty='n')
rect(grconvertX(0.025, from='npc'), grconvertY(0.025, from='npc'),
  grconvertX(0.325, from='npc'), grconvertY(0.55, from='npc'), 
  col="gray90", border=NA)

plotdim <- par("plt")
par(fig = c(plotdim[1]/2.5,plotdim[2]/8,plotdim[3]*1.2,plotdim[4]/1.7), new=TRUE, bty='n')
boxplot(index.table$nobs~L1_col,col='white', main='Level 1 Ecoregions')
abline(h=0,lty=2)

par(fig = c(plotdim[1]/2.5+(1/3),plotdim[2]/8+(1/3),plotdim[3]*1.2,plotdim[4]/1.7), new=TRUE, bty='n')
boxplot((index.table$npan_spr/index.table$npan_tot)~L1_col,col='white')

par(fig = c(plotdim[1]/2.5+(2/3),plotdim[2]/8+(2/3),plotdim[3]*1.2,plotdim[4]/1.7), new=TRUE, bty='n')
boxplot(index.table$nyrs_spr~L1_col,col='white')


dev.off() 


#### Figure S2 ####

setwd('/projectnb/modislc/users/emelaas/scratch15/NASA_TE/R_files/')
load(file="scene_panel_mats_all")

#Plot histograms of spring trends across each sidelap region
x11(h=8,w=11)
#pdf(h=8,w=11,'/projectnb/modislc/projects/te_phenology/figures/Storyboard_3/SM_figure8.pdf')
par(mfrow=c(9,9),mar=c(1,1,1,1))
for (i in 1:length(scenes)){
  ind <- as.numeric(rownames(index)[i])
  if (length(spr2_list[[ind]])>0){
    d <- density(spr2_list[[ind]])
    d$x <- c(0,d$x,0)
    d$y <- c(0,d$y,0)
    
    q1 <- quantile(spr2_list[[ind]],0.025)
    q2 <- quantile(spr2_list[[ind]],0.975)
    
    if (q1<(-1) | q2>1){
      plot(d,xlim=c(-2,2),xaxt='n',main=scenes[i],yaxt='n',cex.main=0.95)
      polygon(d$x[which(d$x<=0)],d$y[which(d$x<=0)],col='lightblue')
      polygon(d$x[which(d$x>=0)],d$y[which(d$x>=0)],col='indianred2')
      text(seq(-2,2,1),max(d$y)-0.05*max(d$y), seq(-2,2,1), cex = 0.7)
      #text(-1.7,seq(0,max(d$y),0.5), seq(0,max(d$y),0.5), cex = 0.8)
      grid(nx=NULL,ny=NA,col='black',lwd=0.4)
    } else {
      plot(d,xlim=c(-1,1),xaxt='n',main=scenes[i],yaxt='n',cex.main=0.95)
      polygon(d$x[which(d$x<=0)],d$y[which(d$x<=0)],col='lightblue')
      polygon(d$x[which(d$x>=0)],d$y[which(d$x>=0)],col='indianred2')
      text(seq(-1,1,0.5),max(d$y)-0.05*max(d$y), seq(-1,1,0.5), cex = 0.7)
      #text(-1.7,seq(0,max(d$y),0.5), seq(0,max(d$y),0.5), cex = 0.8)
      grid(nx=NULL,ny=NA,col='black',lwd=0.4)
    }
  }
}
dev.off()


#### Figure S3 ####

setwd('/projectnb/modislc/users/emelaas/scratch15/NASA_TE/R_files/')
load(file='meta_analysis_nodist_all7')

cuts_s1 <- classIntervals(index$perc_spr_sig_adv,dataPrecision=3,style='fixed',
  fixedBreaks=c(0,0.025,0.2,0.40,0.6,0.8,1))
cuts_s2 <- classIntervals(index$perc_spr_sig_del,dataPrecision=3,style='fixed',
  fixedBreaks=c(0,0.025,0.2,0.40,0.6,0.8,1))
plotclr <- brewer.pal(9,'Oranges')
colcode_s1 <- findColours(cuts_s1,plotclr)
colcode_s2 <- findColours(cuts_s2,plotclr)

#x11(h=6,w=40/3)
pdf(h=6,w=40/3,'/projectnb/modislc/projects/te_phenology/figures/Storyboard_3/tmp2.pdf')
par(mfrow=c(1,2),mar=c(2.6,3,2,1))
plot(sidelaps,pch=21,cex=0.01,col=colcode_s1,bg=colcode_s1,
  main='Fraction of 500 m Panels with Significant Advanced SOS')
lim <- par("usr")
rect(lim[1],lim[3],lim[2],lim[4],border="lightblue",col="lightblue")
plot(usa_shp_proj,axes=FALSE,add=TRUE,col='gray70',lwd=0.1,border='gray70')
plot(can_shp_proj,axes=FALSE,add=TRUE,col='gray70',lwd=0.1,border='gray70')
plot(mex_shp_proj,axes=FALSE,add=TRUE,col='gray70',lwd=0.1,border='gray70')
plot(ecoreg_shp,axes=FALSE,add=TRUE,col='white',lwd=0.2,border='gray40')
plot(sidelaps,pch=21,cex=1.6,col=colcode_s1,bg=colcode_s1,add=TRUE,lwd=0.5)
mtext('a',side=3,line=0.5,cex=1,adj=0,font=2)
legend("bottomright",cex=0.9,title='',
  legend=names(attr(colcode_s1,'table')),fill=attr(colcode_s1,'palette'),bty='n')
rect(grconvertX(0.025, from='npc'), grconvertY(0.025, from='npc'),
  grconvertX(0.325, from='npc'), grconvertY(0.55, from='npc'), 
  col="gray90", border=NA)

plot(sidelaps,pch=21,cex=0.01,col=colcode_s1,bg=colcode_s1,
  main='Fraction of 500 m Panels with Significant Delayed SOS')
lim <- par("usr")
rect(lim[1],lim[3],lim[2],lim[4],border="lightblue",col="lightblue")
plot(usa_shp_proj,axes=FALSE,add=TRUE,col='gray70',lwd=0.1,border='gray70')
plot(can_shp_proj,axes=FALSE,add=TRUE,col='gray70',lwd=0.1,border='gray70')
plot(mex_shp_proj,axes=FALSE,add=TRUE,col='gray70',lwd=0.1,border='gray70')
plot(ecoreg_shp,axes=FALSE,add=TRUE,col='white',lwd=0.2,border='gray40')
plot(sidelaps,pch=21,cex=1.6,col=colcode_s2,bg=colcode_s2,add=TRUE,lwd=0.5)
mtext('b',side=3,line=0.5,cex=1,adj=0,font=2)
legend("bottomright",cex=0.9,title='',
  legend=names(attr(colcode_s1,'table')),fill=attr(colcode_s1,'palette'),bty='n')
rect(grconvertX(0.025, from='npc'), grconvertY(0.025, from='npc'),
  grconvertX(0.325, from='npc'), grconvertY(0.55, from='npc'), 
  col="gray90", border=NA)

plotdim <- par("plt")
par(fig = c(plotdim[1]/2.3,plotdim[2]/4.7,plotdim[3]*1.2,plotdim[4]/1.65), new=TRUE, bty='n')
boxplot(index.table$perc_spr_sig_adv~L1_col,col='white', main='Level 1 Ecoregions',cex.main=0.8)

par(fig = c(plotdim[1]/2.3+(1/2),plotdim[2]/4.7+(1/2),plotdim[3]*1.2,plotdim[4]/1.65), new=TRUE, bty='n')
boxplot(index.table$perc_spr_sig_del~L1_col,col='white')

dev.off()       


#### Figure S4 ####

setwd('/projectnb/modislc/users/emelaas/scratch15/NASA_TE/R_files/')
load(file='meta_analysis_nodist_all7')
load('preseason_temp_cor_gdd')

#Plot correlation between mean preseason temperature and SOS
cuts_COR <- classIntervals(COR3^2,dataPrecision=3,style='fixed',
  fixedBreaks=c(0,0.4,0.5,0.6,0.7,0.8,1))
plotclr <- brewer.pal(6,'YlGnBu')
colcode_COR <- findColours(cuts_COR,plotclr)

x11(h=6,w=6)
#pdf(h=6,w=6,'/projectnb/modislc/projects/te_phenology/figures/Storyboard_3/FigureA4.pdf')
par(mfrow=c(1,1),mar=c(2.6,3,2,1))
plot(sidelaps,pch=21,cex=0.01,col=colcode_COR,bg=colcode_COR,
  main=expression(bold(R^2)),cex.main=1)
lim <- par("usr")
rect(lim[1],lim[3],lim[2],lim[4],border="lightblue",col="lightblue")
plot(usa_shp_proj,axes=FALSE,add=TRUE,col='gray70',lwd=0.1,border='gray70')
plot(can_shp_proj,axes=FALSE,add=TRUE,col='gray70',lwd=0.1,border='gray70')
plot(mex_shp_proj,axes=FALSE,add=TRUE,col='gray70',lwd=0.1,border='gray70')
plot(ecoreg_shp,axes=FALSE,add=TRUE,col='white',lwd=0.2,border='gray40')
plot(sidelaps,pch=21,cex=1.6,col=colcode_COR,bg=colcode_COR,add=TRUE,lwd=0.5)
legend("bottomright",title='',cex=0.9,legend=names(attr(colcode_COR,'table')),
  fill=attr(colcode_COR,'palette'),bty='n')
rect(grconvertX(0.025, from='npc'), grconvertY(0.025, from='npc'),
  grconvertX(0.325, from='npc'), grconvertY(0.55, from='npc'), 
  col="gray90", border=NA)

plotdim <- par("plt")
par(fig = c(plotdim[1]/1.2,plotdim[2]/2.5,plotdim[3]*1.2,plotdim[4]/1.65), new=TRUE, bty='n',cex.axis=0.7)
boxplot(COR3^2~L1_col,col='white',main='Level 1 Ecoregions',cex.main=0.8)

dev.off()


#### Figure S5 ####

setwd('/projectnb/modislc/users/emelaas/scratch15/NASA_TE/R_files/')
load(file='meta_analysis_nodist_all7')
load('preseason_temp_cor_gdd')

index$COR3 <- COR3^2

w <- apply(site_cor3,1,which.min)
cuts_Twindow <- classIntervals(w*5,dataPrecision=3,style='fixed',
  fixedBreaks=c(20,30,40,50,60,70,80,90))
plotclr <- rev(brewer.pal(7,'RdYlBu'))
colcode_Twindow <- findColours(cuts_Twindow,plotclr)

#x11(h=6,w=6)
pdf(h=6,w=6,'/projectnb/modislc/projects/te_phenology/figures/Storyboard_4/FigureS5.pdf')
par(mfrow=c(1,1),mar=c(2.6,3,2,1))
plot(sidelaps,pch=21,cex=0.01,col=colcode_Twindow,bg=colcode_Twindow,
  main='Optimal Preseason Length (days)',cex.main=0.9)
lim <- par("usr")
rect(lim[1],lim[3],lim[2],lim[4],border="lightblue",col="lightblue")
plot(usa_shp_proj,axes=FALSE,add=TRUE,col='gray70',lwd=0.1,border='gray70')
plot(can_shp_proj,axes=FALSE,add=TRUE,col='gray70',lwd=0.1,border='gray70')
plot(mex_shp_proj,axes=FALSE,add=TRUE,col='gray70',lwd=0.1,border='gray70')
plot(ecoreg_shp,axes=FALSE,add=TRUE,col='white',lwd=0.2,border='gray40')
plot(sidelaps,pch=21,cex=1.6,col=colcode_Twindow,bg=colcode_Twindow,add=TRUE,lwd=0.5)
legend("bottomright",title='',cex=0.9,legend=names(attr(colcode_Twindow,'table')),
  fill=attr(colcode_Twindow,'palette'),bty='n')
rect(grconvertX(0.025, from='npc'), grconvertY(0.025, from='npc'),
  grconvertX(0.325, from='npc'), grconvertY(0.55, from='npc'), 
  col="gray90", border=NA)

plotdim <- par("plt")
par(fig = c(plotdim[1]/1.2,plotdim[2]/2.5,plotdim[3]*1.2,plotdim[4]/1.65), new=TRUE, bty='n',cex.axis=0.7)
boxplot(w*5~L1_col,col='white', main='Level 1 Ecoregions',cex.main=0.8)
abline(h=0,lty=2)

dev.off() 



