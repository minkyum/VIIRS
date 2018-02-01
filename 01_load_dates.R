## This code generate is for VIIRS validation.
# 1. Load land cover type (MCD12Q1)
# 2. LSP dates
#     A. H12V04 tile (2012-2014)
#        - VIIRS
#        - MODIS
#     B. H08V05 and H11V05 tiles (2013)
#        - VIIRS
#        - MODIS
# 3. Overview via raster and histgram

rm(list = ls())
# Loading required libraries 
library(sp)
library(raster)
library(gdalUtils)
library(rgeos)
library(maptools)
library(rgdal)
library(RColorBrewer)

### 1. Load land cover type (MCD12Q1) 
# File path 
mcd12q1_path <- '/projectnb/modislc/data/mcd12_out/lc_out/c5_hdf/c5.1_deliv/bug_fix_081214/mcd12q1/'
# Info for specific time and tile 
for(i in 1:3){
  tile <- c('h12v04','h11v04','h08v05')
  type = c('LC1','LC2','LC3','LC4','LC5')
  years = seq(2001,2013)
  LC <- vector('list',length(years))
  lc <- vector('list',length(type))
  for(j in 1:length(years)){
    for(e in 1:length(type)){
      mcd12q1_full_path <- paste(mcd12q1_path,'hdf_',years[j],'/',type[e],'/',
                                 type[e],'.A',years[j],'001.',tile[i],'.hdf',sep='')
      lc[[e]] = raster(mcd12q1_full_path)
    }
    LC[[j]] <- lc
  }
  if(i==1){
    lct.1204 <- LC[[13]][[1]]  
  }else if(i==2){
    lct.1104 <- LC[[13]][[1]]  
  }else{
    lct.0805 <- LC[[13]][[1]]  
  }
}
ext1204 <- raster('/projectnb/modislc/users/mkmoon/VIIRS/NLCD_crosswalked_IGBP/nlcd_igbp_500m/h12v04.bsq')
ext1104 <- raster('/projectnb/modislc/users/mkmoon/VIIRS/NLCD_crosswalked_IGBP/nlcd_igbp_500m/h11v04.bsq')
ext0805 <- raster('/projectnb/modislc/users/mkmoon/VIIRS/NLCD_crosswalked_IGBP/nlcd_igbp_500m/h08v05.bsq')

lct1204 <- table(getValues(lct.1204))
lct1104 <- table(getValues(lct.1104))
lct0805 <- table(getValues(lct.0805))


### 2. LSP dates
## A. H12V04 tile
tile <- 'h12v04'
# VIIRS 
for(yy in 2012:2014){
  
  phe <- c('Increase','Maximum','Decrease','Minimum')
  phe_viirs <- vector('list',length(phe))
  str <- paste('H12V04_phenology_EVI2_',yy,'_04132017_Onset_Greenness_',sep='')
  rast <- vector('list',length(phe))
  for(i in 1:length(phe)){
    ncol <- 2400
    nrow <- 2400
    nbands <- 2
    cnt <- ncol*nrow*nbands
    data <- readBin(paste('/projectnb/modislc/users/mkmoon/VIIRS/tiles/h12v04/Dates_VIIRS/',str,phe[i],'.img',sep=''),what="integer",n=cnt,size=2,endian="little")
    data[data>30000] <- NA
    data <- data-(yy-2000)*366
    data2 <- array(data,c(nbands, ncol, nrow))
    data2 <- aperm(data2, c(3,2,1)) #for transposing
    rast[[i]] <- brick(data2,xmn=ext1204@extent@xmin,xmx=ext1204@extent@xmax,ymn=ext1204@extent@ymin,ymx=ext1204@extent@ymax,crs=ext1204@crs)
  }
  if(yy==2012){
    vv120412 <- rast
  }else if(yy==2013){
    vv120413 <- rast
  }else{
    vv120414 <- rast
  }
}
# C6
for(yy in 2012:2014){
  phe <- c('Greenup','Maturity','Senescence','Dormancy')  
  rast <- vector('list',length(phe))
  for(vari in 1:4){
    path <- paste('/projectnb/modislc/users/dsm/eval_modis_lc_061917/MCD12I6/',phe[vari],'/',phe[vari],'_',tile,'_',yy,sep='')
    ncol <- 2400
    nrow <- 2400
    nbands <- 2
    cnt <- ncol*nrow*nbands
    data <- readBin(path,what="integer",n=cnt,size=2,endian="little")
    data[data>30000] <- NA
    data1 <- data - as.numeric(as.Date(paste(yy-1,'-12-31',sep='')))
    data2 <- array(data1,c(nbands, ncol, nrow))
    data2 <- aperm(data2, c(3,2,1)) #for transposing
    aa <- brick(data2,xmn=ext1204@extent@xmin,xmx=ext1204@extent@xmax,ymn=ext1204@extent@ymin,ymx=ext1204@extent@ymax,crs=ext1204@crs)
    rast[[vari]] <- aa[[1]]  
  }
  if(yy==2012){
    mm120412 <- rast
  }else if(yy==2013){
    mm120413 <- rast
  }else{
    mm120414 <- rast
  }
}
 
## B. H08V05 and H11V05 tiles
# VIIRS 0805 and 1104
for(tt in 1:2){
  year <- 2013
  tile <- c('h08v05','h11v04')
  viirs_path <- paste('/projectnb/modislc/users/mkmoon/VIIRS/tiles/',tile[tt],'/Dates_VIIRS/',sep='')
  phe <- c('gri','gre','sei','see')
  phe_viirs <- vector('list',length(phe))
  
  if(tt==1){
    str <- paste('VIIRS_H08V05_',phe,'_EVI2_x2013',sep='')
  }else{
    str <- paste('VIIRS_H11V04_',phe,'_EVI2_x2013',sep='')  
  }
  rast <- vector('list',length(phe))
  for(i in 1:length(phe)){
    cnt <- ncol*nrow*nbands
    data <- readBin(paste(viirs_path,str[i],sep=''),what="integer",n=cnt,size=2,endian="little")
    data[data>30000] <- NA
    data2 <- array(data,c(nbands, ncol, nrow))
    data2 <- aperm(data2, c(3,2,1)) #for transposing
    if(tt==1){
      rast[[i]] <- brick(data2,xmn=ext0805@extent@xmin,xmx=ext0805@extent@xmax,ymn=ext0805@extent@ymin,ymx=ext0805@extent@ymax)  
    }else{
      rast[[i]] <- brick(data2,xmn=ext1104@extent@xmin,xmx=ext1104@extent@xmax,ymn=ext1104@extent@ymin,ymx=ext1104@extent@ymax)
    }
    
  }
  if(tt==1){
    vv0805 <- rast  
  }else{
    vv1104 <- rast  
  }
}
# MODIS 0805 and 1104
for(tt in 1:2){
  year <- 2013
  tile <- c('h08v05','h11v04')

  phe <- c('Greenup','Maturity','Senescence','Dormancy')  
  rast <- vector('list',length(phe))
  for(vari in 1:4){
    path <- paste('/projectnb/modislc/users/dsm/eval_modis_lc_061917/MCD12I6/',phe[vari],'/',phe[vari],'_',tile[tt],'_',year,sep='')
    ncol <- 2400
    nrow <- 2400
    nbands <- 2
    cnt <- ncol*nrow*nbands
    data <- readBin(path,what="integer",n=cnt,size=2,endian="little")
    data[data>30000] <- NA
    data1 <- data - as.numeric(as.Date(paste(year-1,'-12-31',sep='')))
    data2 <- array(data1,c(nbands, ncol, nrow))
    data2 <- aperm(data2, c(3,2,1)) #for transposing
    if(tt==1){
      aa <- brick(data2,xmn=ext0805@extent@xmin,xmx=ext0805@extent@xmax,ymn=ext0805@extent@ymin,ymx=ext0805@extent@ymax)  
    }else{
      aa <- brick(data2,xmn=ext1104@extent@xmin,xmx=ext1104@extent@xmax,ymn=ext1104@extent@ymin,ymx=ext1104@extent@ymax)
    }
    rast[[vari]] <- aa[[1]]  
  }
  if(tt==1){
    mm0805 <- rast  
  }else{
    mm1104 <- rast  
  }
}


### 3. Overview via raster and histgram
phe.plot <- function(rast1,rast2,tile,year,vari){
  
  par(mfrow=c(2,2),mgp=c(2,1,0),oma=c(1,1,1,1),mar=c(3,3,2.5,2.5)) 
  
  values(lct.1204)[values(lct.1204!=0)] <- NA
  values(lct.1104)[values(lct.1104!=0)] <- NA
  values(lct.0805)[values(lct.0805!=0)] <- NA
  
  if(as.numeric(substr(tile,2,3))==12){
    wm <- lct.1204 
    lp <- 2400*2400 - sum(values(wm)==0,na.rm=T)
    yylim=c(0,0.035)
  }else if(as.numeric(substr(tile,2,3))==11){
    wm <- lct.1104 
    lp <- 2400*2400 - sum(values(wm)==0,na.rm=T)
    yylim=c(0,0.035)
  }else{
    wm <- lct.0805 
    lp <- 2400*2400 - sum(values(wm)==0,na.rm=T)
    yylim=c(0,0.022)
  }
  
  bb <- rast1
  
  mycol <- brewer.pal(11,'Spectral')
  mycol <- rev(colorRampPalette(mycol)(1000))
  
  upp <- round(quantile(bb,0.98,na.rm=T))
  lwp <- round(quantile(bb,0.02,na.rm=T))
  leg.int <- round(seq(lwp,upp,length.out = 5))  
  values(bb)[values(bb)>=upp] <- upp
  values(bb)[values(bb)<=lwp] <- lwp
  
  plot(bb,zlim=c(lwp,upp),box=F,bty = "n",xaxt = "n", yaxt = "n",
       col=mycol,colNA='grey75',legend=F,main=paste("VIIRS_",tile,"_",year,"_",vari,sep=""))
  plot(wm,add=T,col='grey45',legend=F)
  
  bb <- rast2
  values(bb)[values(bb)>=upp] <- upp
  values(bb)[values(bb)<=lwp] <- lwp
  
  plot(bb,zlim=c(lwp,upp),box=F,bty = "n",xaxt = "n", yaxt = "n",
       col=mycol,colNA='grey75',legend=F,main=paste("MODIS_",tile,"_",year,"_",vari,sep=""))
  plot(wm,add=T,col='grey45',legend=F)
  plot(bb,legend.only=T,col=mycol,zlim=c(lwp,upp),
       legend.width=1.5,legend.shrink=0.6,
       smallplot=c(0.80,0.83,0.2,0.8),
       axis.args=list(at=leg.int,cex.axis=1.2,font=1,
                      labels=c(paste('<',lwp,sep=''),leg.int[2:4],paste('>',upp,sep='')))) 
  
  dd1 <- getValues(rast1)
  dd2 <- getValues(rast2)
  
  hist(dd1,xlim=c(-100,450),ylim=yylim,
       breaks = seq(-400,700,1),probability=T,lty="blank",col=rgb(1,0,0,0.5),
       main=NULL,cex.axis=1,xlab="Day of year",ylab="Density",cex.lab=1)
  hist(dd2,xlim=c(-100,450),breaks = seq(-400,700,1),probability=T,lty="blank",col=rgb(0,0,1,0.5),
       main=NULL,cex.axis=1,xlab="Day of year",ylab="Density",cex.lab=1,add=T)
  legend("topright",c("VIIRS","MODIS"),pch=22,pt.bg=c(rgb(1,0,0,0.5),rgb(0,0,1,0.5)),cex=1.2,bty="n")  
  
  pp1 <- 1 - sum(!is.na(dd1))/lp
  pp2 <- 1 - sum(!is.na(dd2))/lp
  
  barplot(c(pp1,pp2),ylim=c(0,1),name=c("VIIRS","MODIS"),ylab="NA fraction (land only)")
  text(0.7,pp1,round(pp1,3),pos=3)
  text(1.9,pp2,round(pp2,3),pos=3)
}

setwd('/projectnb/modislc/users/mkmoon/VIIRS/figures/')
pdf(file=paste('VIvsC6_diag.pdf',sep=''),width=9,height=7)

for(i in 1:4){
  phe.plot(vv120412[[i]][[1]],mm120412[[i]][[1]],"H12V04",2012,phe[i])
}
for(i in 1:4){
  phe.plot(vv120413[[i]][[1]],mm120413[[i]][[1]],"H12V04",2013,phe[i])
}
for(i in 1:4){
  phe.plot(vv120414[[i]][[1]],mm120414[[i]][[1]],"H12V04",2014,phe[i])
}

for(i in 1:4){
  phe.plot(vv1104[[i]][[1]],mm1104[[i]][[1]],"H11V04",2013,phe[i])
}
for(i in 1:4){
  phe.plot(vv0805[[i]][[1]],mm0805[[i]][[1]],"H08V05",2013,phe[i])
}

dev.off()





 
# 
# 
# #filter
# values(vv120412[[1]][[1]])[values(vv120412[[1]][[1]])<2] <- NA
# values(mm120412[[1]][[1]])[values(mm120412[[1]][[1]])<2] <- NA
# values(vv120413[[1]][[1]])[values(vv120413[[1]][[1]])<2] <- NA
# values(mm120413[[1]][[1]])[values(mm120413[[1]][[1]])<2] <- NA
# values(vv120414[[1]][[1]])[values(vv120414[[1]][[1]])<2] <- NA
# values(mm120414[[1]][[1]])[values(mm120414[[1]][[1]])<2] <- NA
# values(vv120412[[4]][[1]])[values(vv120412[[4]][[1]])<150] <- NA
# values(mm120412[[4]][[1]])[values(mm120412[[4]][[1]])<150] <- NA
# values(vv120413[[4]][[1]])[values(vv120413[[4]][[1]])<150] <- NA
# values(mm120413[[4]][[1]])[values(mm120413[[4]][[1]])<150] <- NA
# values(vv120414[[4]][[1]])[values(vv120414[[4]][[1]])<150] <- NA
# values(mm120414[[4]][[1]])[values(mm120414[[4]][[1]])<150] <- NA
# 
# par(mfrow=c(4,2),oma=c(1,1,1,1),mar=c(1,1,1,1),mgp=c(1,1,0),bty='n')
# plot(vv120413[[1]][[1]],colNA='grey',zlim=c(0,366),legend=F,axes=F)
# plot(mm120413[[1]][[1]],colNA='grey',zlim=c(0,366),legend=F,axes=F)
# plot(vv120413[[2]][[1]],colNA='grey',zlim=c(0,366),legend=F,axes=F)
# plot(mm120413[[2]][[1]],colNA='grey',zlim=c(0,366),legend=F,axes=F)
# plot(vv120413[[3]][[1]],colNA='grey',zlim=c(0,366),legend=F,axes=F)
# plot(mm120413[[3]][[1]],colNA='grey',zlim=c(0,366),legend=F,axes=F)
# plot(vv120413[[4]][[1]],colNA='grey',zlim=c(0,366),legend=F,axes=F)
# plot(mm120413[[4]][[1]],colNA='grey',zlim=c(0,366),legend=F,axes=F)
# hist(vv120413[[1]][[1]],breaks = seq(0,365,1))
# hist(mm120413[[1]][[1]],breaks = seq(0,365,1))
# hist(vv120413[[2]][[1]],breaks = seq(0,365,1))
# hist(mm120413[[2]][[1]],breaks = seq(0,365,1))
# hist(vv120413[[3]][[1]],breaks = seq(0,365,1))
# hist(mm120413[[3]][[1]],breaks = seq(0,365,1))
# hist(vv120413[[4]][[1]],breaks = seq(0,365,1))
# hist(mm120413[[4]][[1]],breaks = seq(0,365,1))
# 
# par(mfrow=c(4,2),oma=c(1,1,1,1),mar=c(1,1,1,1),mgp=c(1,1,0),bty='n')
# plot(vv0805[[1]][[1]],colNA='grey',zlim=c(0,366),legend=F,axes=F)
# plot(mm0805[[1]][[1]],colNA='grey',zlim=c(0,366),legend=F,axes=F)
# plot(vv0805[[2]][[1]],colNA='grey',zlim=c(0,366),legend=F,axes=F)
# plot(mm0805[[2]][[1]],colNA='grey',zlim=c(0,366),legend=F,axes=F)
# plot(vv0805[[3]][[1]],colNA='grey',zlim=c(0,366),legend=F,axes=F)
# plot(mm0805[[3]][[1]],colNA='grey',zlim=c(0,366),legend=F,axes=F)
# plot(vv0805[[4]][[1]],colNA='grey',zlim=c(0,366),legend=F,axes=F)
# plot(mm0805[[4]][[1]],colNA='grey',zlim=c(0,366),legend=F,axes=F)
# hist(vv0805[[1]][[1]],breaks = seq(0,365,1))
# hist(mm0805[[1]][[1]],breaks = seq(0,365,1))
# hist(vv0805[[2]][[1]],breaks = seq(0,365,1))
# hist(mm0805[[2]][[1]],breaks = seq(0,365,1))
# hist(vv0805[[3]][[1]],breaks = seq(0,365,1))
# hist(mm0805[[3]][[1]],breaks = seq(0,365,1))
# hist(vv0805[[4]][[1]],breaks = seq(0,365,1))
# hist(mm0805[[4]][[1]],breaks = seq(0,365,1))
# # filter
# values(vv0805[[1]][[1]])[values(vv0805[[1]][[1]])==1] <- NA
# values(mm0805[[1]][[1]])[values(mm0805[[1]][[1]])==1] <- NA
# values(vv0805[[2]][[1]])[values(vv0805[[2]][[1]])==1] <- NA
# values(mm0805[[2]][[1]])[values(mm0805[[2]][[1]])==1] <- NA
# 
# par(mfrow=c(4,2),oma=c(1,1,1,1),mar=c(1,1,1,1),mgp=c(1,1,0),bty='n')
# plot(vv1104[[1]][[1]],colNA='grey',zlim=c(0,366),legend=F,axes=F)
# plot(mm1104[[1]][[1]],colNA='grey',zlim=c(0,366),legend=F,axes=F)
# plot(vv1104[[2]][[1]],colNA='grey',zlim=c(0,366),legend=F,axes=F)
# plot(mm1104[[2]][[1]],colNA='grey',zlim=c(0,366),legend=F,axes=F)
# plot(vv1104[[3]][[1]],colNA='grey',zlim=c(0,366),legend=F,axes=F)
# plot(mm1104[[3]][[1]],colNA='grey',zlim=c(0,366),legend=F,axes=F)
# plot(vv1104[[4]][[1]],colNA='grey',zlim=c(0,366),legend=F,axes=F)
# plot(mm1104[[4]][[1]],colNA='grey',zlim=c(0,366),legend=F,axes=F)
# hist(vv1104[[1]][[1]],breaks = seq(0,365,1))
# hist(mm1104[[1]][[1]],breaks = seq(0,365,1))
# hist(vv1104[[2]][[1]],breaks = seq(0,365,1))
# hist(mm1104[[2]][[1]],breaks = seq(0,365,1))
# hist(vv1104[[3]][[1]],breaks = seq(0,365,1))
# hist(mm1104[[3]][[1]],breaks = seq(0,365,1))
# hist(vv1104[[4]][[1]],breaks = seq(0,365,1))
# hist(mm1104[[4]][[1]],breaks = seq(0,365,1))
# # filter
# values(vv1104[[1]][[1]])[values(vv1104[[1]][[1]])==1] <- NA
# values(mm1104[[1]][[1]])[values(mm1104[[1]][[1]])==1] <- NA
# values(vv1104[[4]][[1]])[values(vv1104[[4]][[1]])<150] <- NA
# values(mm1104[[4]][[1]])[values(mm1104[[4]][[1]])<150] <- NA
# 
# 
# ## Summary stats
# lctp1204 <- lct1204/(sum(lct1204[2:17]))*100
# lctp1104 <- lct1104/(sum(lct1104[2:17]))*100
# lctp0805 <- lct0805/(sum(lct0805[2:17]))*100
# 
# barplot(lctp1204) # 5 and 4 (Mixed and Deciduous forests)
# barplot(lctp1104) # 12 and 14 (Croplands and Cropland/natural VM)
# barplot(lctp0805) # 7 and 10 (Open shrublands and Grasslands)
# 
# lct1204[c(5+1,4+1)]
# lct1104[c(12+1,14+1)]
# lct0805[c(7+1,10+1)]
# 
# lctp1204[c(5+1,4+1)]
# lctp1104[c(12+1,14+1)]
# lctp0805[c(7+1,10+1)]
# 
# mean(vv120413[[1]][[1]][lct1204==5],narm=T)
# mean(vv120413[[2]][[1]][lct1204==5],narm=T)
# mean(vv120413[[3]][[1]][lct1204==5],narm=T)
# mean(vv120413[[4]][[1]][lct1204==5],narm=T)
# mean(mm120413[[1]][[1]][lct.1204==5],na.rm=T)
# mean(mm.1204.13[[2]][[1]][lct.1204==5],na.rm=T)
# mean(mm.1204.13[[3]][[1]][lct.1204==5],na.rm=T)
# mean(mm.1204.13[[4]][[1]][lct.1204==5],na.rm=T)
# 
# sd(vv.1204.13[[1]][[1]][lct.1204==5],na.rm=T)
# sd(vv.1204.13[[2]][[1]][lct.1204==5],na.rm=T)
# sd(vv.1204.13[[3]][[1]][lct.1204==5],na.rm=T)
# sd(vv.1204.13[[4]][[1]][lct.1204==5],na.rm=T)
# sd(mm.1204.13[[1]][[1]][lct.1204==5],na.rm=T)
# sd(mm.1204.13[[2]][[1]][lct.1204==5],na.rm=T)
# sd(mm.1204.13[[3]][[1]][lct.1204==5],na.rm=T)
# sd(mm.1204.13[[4]][[1]][lct.1204==5],na.rm=T)
# 
# 
# 
# 
# sum(getValues(vv.1204.13[[1]][[1]])==1,na.rm=T)
# sum(getValues(vv.1204.13[[2]][[1]])==1,na.rm=T)
# sum(getValues(vv.1204.13[[3]][[1]])==1,na.rm=T)
# sum(getValues(vv.1204.13[[4]][[1]])[lct1204==5]==1,na.rm=T)
# sum(getValues(mm.1204.13[[1]][[1]])[lct1204==5]==1,na.rm=T)
# sum(getValues(mm.1204.13[[2]][[1]])[lct1204==5]==1,na.rm=T)
# sum(getValues(mm.1204.13[[3]][[1]])[lct1204==5]==1,na.rm=T)
# sum(getValues(mm.1204.13[[4]][[1]])[lct1204==5]==1,na.rm=T)
# 
# sum(getValues(vv.1104[[1]][[1]])==1,na.rm=T)
# sum(getValues(vv.1104[[2]][[1]])==1,na.rm=T)
# sum(getValues(vv.1104[[3]][[1]])==1,na.rm=T)
# sum(getValues(vv.1104[[4]][[1]])==1,na.rm=T)
# sum(getValues(mm.1104[[1]][[1]])==1,na.rm=T)
# sum(getValues(mm.1104[[2]][[1]])==1,na.rm=T)
# sum(getValues(mm.1104[[3]][[1]])==1,na.rm=T)
# sum(getValues(mm.1104[[4]][[1]])==1,na.rm=T)
# 
# sum(getValues(vv.0805[[1]][[1]])==1,na.rm=T)
# sum(getValues(vv.0805[[2]][[1]])==1,na.rm=T)
# sum(getValues(vv.0805[[3]][[1]])==1,na.rm=T)
# sum(getValues(vv.0805[[4]][[1]])==1,na.rm=T)
# sum(getValues(mm.0805[[1]][[1]])==1,na.rm=T)
# sum(getValues(mm.0805[[2]][[1]])==1,na.rm=T)
# sum(getValues(mm.0805[[3]][[1]])==1,na.rm=T)
# sum(getValues(mm.0805[[4]][[1]])==1,na.rm=T)
