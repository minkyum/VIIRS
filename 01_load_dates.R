# This code generate is for VIIRS validation.
# 1. Load land cover type (MCD12Q1)

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

# File path
viirs_path <- '/projectnb/modislc/users/mkmoon/VIIRS/tiles/h12v04/Dates_VIIRS/'
modis_path6 <- '/projectnb/modislc/data/mcd12_out/phen_out/c6/'
modis_path5 <- '/projectnb/modislc/data/mcd12_out/phen_out/c5_hdf1/'
##  all tiles
path <- '/projectnb/modislc/users/dsm/eval_modis_lc_061917/MCD12I6/Dormancy/'
search_str <- paste('*2001',sep='')
files <- list.files(path=path,pattern=glob2rx(search_str),full.names=T,include.dirs=F)
tiles_c6 <- substr(files,78,83)
tiles <- tiles_c6

# One tile
phe_met <- c('Dormancy','EVI_Amplitude','EVI_Area','EVI_Minimum','Greenup','Maturity','MidGreendown','MidGreenup','NumCycles','Peak','QA_Detailed','QA_Overall','Senescence')
year <- seq(2001,2015)

## Land cover tye C5
# File path 
mcd12q1_path <- "/projectnb/modislc/data/mcd12_out/lc_out/c5_hdf/c5.1_deliv/bug_fix_081214/mcd12q1/"

lct <- vector("list",length(tiles))
for(i in 1:length(tiles)){
  mcd12q1_full_path <- paste(mcd12q1_path,"hdf_2013/LC1/LC1.A2013001.",tiles[i],".hdf",sep="")
  lct[[i]] <- raster(mcd12q1_full_path)
}

## C6 one tile multi-year
args <- commandArgs()
print(args)

tt <- as.numeric(args[3])
i <- tt
  
  setwd(paste('/projectnb/modislc/users/mkmoon/LCD_C6/v2/One tile/',sep=''))
  pdf(file=paste('C6_diag_',tiles[i],'.pdf',sep=''),width=12,height=7)
  
## C6
  # QA_overall
  vari <- 12
  mycol <- brewer.pal(9,'PuRd')
  par(mfcol=c(3,5),mgp=c(0,0,0),oma=c(0.2,0.2,0.2,0.2),mar=c(1.2,0,1.2,0))  
  for(j in 1:15){
    path <- paste('/projectnb/modislc/users/dsm/eval_modis_lc_061917/MCD12I6/',phe_met[vari],'/',phe_met[vari],'_',tiles[i],'_',year[j],sep='')
    ncol <- 2400
    nrow <- 2400
    nbands <- 2
    cnt <- ncol*nrow*nbands
    data <- readBin(path,what="integer",n=cnt,size=2,endian="little")
    data[data>30000] <- NA
    data <- data*0.0001
    data2 <- array(data,c(nbands, ncol, nrow))
    data2 <- aperm(data2, c(3,2,1)) #for transposing
    aa <- brick(data2)
    bb <- aa[[1]]
    cc <- getValues(bb)      
    
    if(sum(is.na(cc))>(2400*2400-10)){
      plot(raster(matrix(1,1,1)),box=F,bty = "n",xaxt = "n", yaxt = "n",
           col='grey45',legend=F)
      title(paste('C6_',tiles[i],'_',year[j],'_',phe_met[vari],'_NA',sep=''),line=0.25,cex.main=0.9)
    }else{
      if(j==1){
        plot(bb,box=F,bty = "n",xaxt = "n", yaxt = "n",
             col=mycol[c(2,4,6,8)],colNA='grey45',legend=F)
        title(paste('C6_',tiles[i],'_',year[j],'_',phe_met[vari],sep=''),line=0.25,cex.main=0.9)
        legend('bottomright',c('0','1','2','3'),pt.bg=mycol[c(2,4,6,8)],pch=22)
      }else{
        plot(bb,box=F,bty = "n",xaxt = "n", yaxt = "n",
             col=mycol[c(2,4,6,8)],colNA='grey45',legend=F)
        title(paste('C6_',tiles[i],'_',year[j],'_',phe_met[vari],sep=''),line=0.25,cex.main=0.9)
      }
    }
  } 

  # Num_Cycles
  vari <- 9
  mycol <- brewer.pal(9,'YlGnBu')
  par(mfcol=c(3,5),mgp=c(0,0,0),oma=c(0.2,0.2,0.2,0.2),mar=c(1.2,0,1.2,0))  
  for(j in 1:15){
    path <- paste('/projectnb/modislc/users/dsm/eval_modis_lc_061917/MCD12I6/',phe_met[vari],'/',phe_met[vari],'_',tiles[i],'_',year[j],sep='')
    ncol <- 2400
    nrow <- 2400
    nbands <- 1
    cnt <- ncol*nrow*nbands
    data <- readBin(path,what="integer",n=cnt,size=2,endian="little")
    data[data>30000] <- NA
    data <- data*0.0001
    data2 <- array(data,c(nbands, ncol, nrow))
    data2 <- aperm(data2, c(3,2,1)) #for transposing
    aa <- brick(data2)
    bb <- aa[[1]]
    cc <- getValues(bb)      
    
    if(sum(is.na(cc))>(2400*2400-10)){
      plot(raster(matrix(1,1,1)),box=F,bty = "n",xaxt = "n", yaxt = "n",
           col='grey45',legend=F)
      title(paste('C6_',tiles[i],'_',year[j],'_',phe_met[vari],'_NA',sep=''),line=0.25,cex.main=0.9)
    }else{
      if(j==1){
        plot(bb,box=F,bty = "n",xaxt = "n", yaxt = "n",
             col=mycol[c(2,4,6,8)],colNA='grey45',legend=F)
        title(paste('C6_',tiles[i],'_',year[j],'_',phe_met[vari],sep=''),line=0.25,cex.main=0.9)
        legend('bottomright',c('1','2','3','>3'),pt.bg=mycol[c(2,4,6,8)],pch=22)
      }else{
        plot(bb,box=F,bty = "n",xaxt = "n", yaxt = "n",
             col=mycol[c(2,4,6,8)],colNA='grey45',legend=F)
        title(paste('C6_',tiles[i],'_',year[j],'_',phe_met[vari],sep=''),line=0.25,cex.main=0.9)
      }
    }
  } 
  
  # EVI amplitude
  vari <- 2
  mycol <- brewer.pal(9,'BuGn')
  mycol <- colorRampPalette(mycol)(1000)
    par(mfcol=c(3,5),mgp=c(0,0,0),oma=c(0.2,0.2,0.2,0.2),mar=c(1.2,0,1.2,0))  
    for(j in 1:15){
      path <- paste('/projectnb/modislc/users/dsm/eval_modis_lc_061917/MCD12I6/',phe_met[vari],'/',phe_met[vari],'_',tiles[i],'_',year[j],sep='')
      ncol <- 2400
      nrow <- 2400
      nbands <- 2
      cnt <- ncol*nrow*nbands
      data <- readBin(path,what="integer",n=cnt,size=2,endian="little")
      data[data>30000] <- NA
      data <- data*0.0001
      data2 <- array(data,c(nbands, ncol, nrow))
      data2 <- aperm(data2, c(3,2,1)) #for transposing
      aa <- brick(data2)
      bb <- aa[[1]]
      cc <- getValues(bb)      
      
      if(sum(is.na(cc))>(2400*2400-10)){
        plot(raster(matrix(1,1,1)),box=F,bty = "n",xaxt = "n", yaxt = "n",
             col='grey45',legend=F)
        title(paste('C6_',tiles[i],'_',year[j],'_',phe_met[vari],'_NA',sep=''),line=0.25,cex.main=0.9)
      }else{
        if(j==1){
          upp <- round(quantile(bb,0.98,na.rm=T),3)
          lwp <- round(quantile(bb,0.02,na.rm=T),3) 
          leg.int <- round(seq(lwp,upp,length.out = 3),3)  
          values(bb)[values(bb)>=upp] <- upp
          values(bb)[values(bb)<=lwp] <- lwp
          plot(bb,zlim=c(lwp,upp),box=F,bty = "n",xaxt = "n", yaxt = "n",
               col=mycol,colNA='grey45',legend=F)
          title(paste('C6_',tiles[i],'_',year[j],'_',phe_met[vari],sep=''),line=0.25,cex.main=0.9)
          plot(bb,legend.only=T,col=mycol,zlim=c(lwp,upp),
               legend.width=1.5,legend.shrink=0.6,horiz=T,
               smallplot=c(0.15,0.8,0.11,0.14),
               axis.args=list(at=leg.int,cex.axis=0.8,font=2,
                              labels=c(paste('<',lwp,sep=''),leg.int[2],paste('>',upp,sep=''))))          
        }else{
          values(bb)[values(bb)>=upp] <- upp
          values(bb)[values(bb)<=lwp] <- lwp
          plot(bb,zlim=c(lwp,upp),box=F,bty = "n",xaxt = "n", yaxt = "n",
               col=mycol,colNA='grey45',legend=F)
          title(paste('C6_',tiles[i],'_',year[j],'_',phe_met[vari],sep=''),line=0.25,cex.main=0.9)
        }
      }  
    } 
  
  # EVI area
  vari <- 3
  mycol <- brewer.pal(9,'Blues')
  mycol <- colorRampPalette(mycol)(1000)
    par(mfcol=c(3,5),mgp=c(0,0,0),oma=c(0.2,0.2,0.2,0.2),mar=c(1.2,0,1.2,0))  
    for(j in 1:15){
    path <- paste('/projectnb/modislc/users/dsm/eval_modis_lc_061917/MCD12I6/',phe_met[vari],'/',phe_met[vari],'_',tiles[i],'_',year[j],sep='')
    ncol <- 2400
    nrow <- 2400
    nbands <- 2
    cnt <- ncol*nrow*nbands
    data <- readBin(path,what="integer",n=cnt,size=2,endian="little")
    data[data>30000] <- NA
    data <- data*0.1
    data2 <- array(data,c(nbands, ncol, nrow))
    data2 <- aperm(data2, c(3,2,1)) #for transposing
    aa <- brick(data2)
    bb <- aa[[1]]
    cc <- getValues(bb)      
    
    if(sum(is.na(cc))>(2400*2400-10)){
      plot(raster(matrix(1,1,1)),box=F,bty = "n",xaxt = "n", yaxt = "n",
           col='grey45',legend=F)
      title(paste('C6_',tiles[i],'_',year[j],'_',phe_met[vari],'_NA',sep=''),line=0.25,cex.main=0.9)
    }else{
      if(j==1){
        upp <- round(quantile(bb,0.98,na.rm=T),1)
        lwp <- round(quantile(bb,0.02,na.rm=T),1) 
        leg.int <- round(seq(lwp,upp,length.out = 5),1)  
        values(bb)[values(bb)>=upp] <- upp
        values(bb)[values(bb)<=lwp] <- lwp
        plot(bb,zlim=c(lwp,upp),box=F,bty = "n",xaxt = "n", yaxt = "n",
             col=mycol,colNA='grey45',legend=F)
        title(paste('C6_',tiles[i],'_',year[j],'_',phe_met[vari],sep=''),line=0.25,cex.main=0.9)
        plot(bb,legend.only=T,col=mycol,zlim=c(lwp,upp),
             legend.width=1.5,legend.shrink=0.6,horiz=T,
             smallplot=c(0.15,0.8,0.11,0.14),
             axis.args=list(at=leg.int,cex.axis=0.8,font=2,
                            labels=c(paste('<',lwp,sep=''),leg.int[2:4],paste('>',upp,sep=''))))          
      }else{
        values(bb)[values(bb)>=upp] <- upp
        values(bb)[values(bb)<=lwp] <- lwp
        plot(bb,zlim=c(lwp,upp),box=F,bty = "n",xaxt = "n", yaxt = "n",
             col=mycol,colNA='grey45',legend=F)
        title(paste('C6_',tiles[i],'_',year[j],'_',phe_met[vari],sep=''),line=0.25,cex.main=0.9)
      }
    }  
  } 
  
  # EVI minimum
  vari <- 4
  mycol <- brewer.pal(9,'Reds')
  mycol <- colorRampPalette(mycol)(1000)
    par(mfcol=c(3,5),mgp=c(0,0,0),oma=c(0.2,0.2,0.2,0.2),mar=c(1.2,0,1.2,0))  
    for(j in 1:15){
    path <- paste('/projectnb/modislc/users/dsm/eval_modis_lc_061917/MCD12I6/',phe_met[vari],'/',phe_met[vari],'_',tiles[i],'_',year[j],sep='')
    ncol <- 2400
    nrow <- 2400
    nbands <- 2
    cnt <- ncol*nrow*nbands
    data <- readBin(path,what="integer",n=cnt,size=2,endian="little")
    data[data>30000] <- NA
    data <- data*0.0001
    data2 <- array(data,c(nbands, ncol, nrow))
    data2 <- aperm(data2, c(3,2,1)) #for transposing
    aa <- brick(data2)
    bb <- aa[[1]]
    cc <- getValues(bb)      
    
    if(sum(is.na(cc))>(2400*2400-10)){
      plot(raster(matrix(1,1,1)),box=F,bty = "n",xaxt = "n", yaxt = "n",
           col='grey45',legend=F)
      title(paste('C6_',tiles[i],'_',year[j],'_',phe_met[vari],'_NA',sep=''),line=0.25,cex.main=0.9)
    }else{
      if(j==1){
        upp <- round(quantile(bb,0.98,na.rm=T),3)
        lwp <- round(quantile(bb,0.02,na.rm=T),3) 
        leg.int <- round(seq(lwp,upp,length.out = 3),3)  
        values(bb)[values(bb)>=upp] <- upp
        values(bb)[values(bb)<=lwp] <- lwp
        plot(bb,zlim=c(lwp,upp),box=F,bty = "n",xaxt = "n", yaxt = "n",
             col=mycol,colNA='grey45',legend=F)
        title(paste('C6_',tiles[i],'_',year[j],'_',phe_met[vari],sep=''),line=0.25,cex.main=0.9)
        plot(bb,legend.only=T,col=mycol,zlim=c(lwp,upp),
             legend.width=1.5,legend.shrink=0.6,horiz=T,
             smallplot=c(0.15,0.8,0.11,0.14),
             axis.args=list(at=leg.int,cex.axis=0.8,font=2,
                            labels=c(paste('<',lwp,sep=''),leg.int[2],paste('>',upp,sep=''))))          
      }else{
        values(bb)[values(bb)>=upp] <- upp
        values(bb)[values(bb)<=lwp] <- lwp
        plot(bb,zlim=c(lwp,upp),box=F,bty = "n",xaxt = "n", yaxt = "n",
             col=mycol,colNA='grey45',legend=F)
        title(paste('C6_',tiles[i],'_',year[j],'_',phe_met[vari],sep=''),line=0.25,cex.main=0.9)
      }
    }  
  } 

  # Peak
  vari <- 10
  mycol <- brewer.pal(11,'Spectral')
  mycol <- rev(colorRampPalette(mycol)(1000))
    par(mfcol=c(3,5),mgp=c(0,0,0),oma=c(0.2,0.2,0.2,0.2),mar=c(1.2,0,1.2,0))  
    for(j in 1:15){
      path <- paste('/projectnb/modislc/users/dsm/eval_modis_lc_061917/MCD12I6/',phe_met[vari],'/',phe_met[vari],'_',tiles[i],'_',year[j],sep='')
      ncol <- 2400
      nrow <- 2400
      nbands <- 2
      cnt <- ncol*nrow*nbands
      data <- readBin(path,what="integer",n=cnt,size=2,endian="little")    
      data[data>30000] <- NA
      data1 <- data - as.numeric(as.Date(paste(year[j]-1,'-12-31',sep='')))
      data2 <- array(data1,c(nbands, ncol, nrow))
      data2 <- aperm(data2, c(3,2,1)) #for transposing
      aa <- brick(data2)
      bb <- aa[[1]]
      cc <- getValues(bb)      
      
      if(sum(is.na(cc))>(2400*2400-10)){
        plot(raster(matrix(1,1,1)),box=F,bty = "n",xaxt = "n", yaxt = "n",
             col='grey45',legend=F)
        title(paste('C6_',tiles[i],'_',year[j],'_',phe_met[vari],'_NA',sep=''),line=0.25,cex.main=0.8)
      }else{
        if(j==1){
          upp <- round(quantile(bb,0.98,na.rm=T))
          lwp <- round(quantile(bb,0.02,na.rm=T)) 
          leg.int <- round(seq(lwp,upp,length.out = 5))  
          values(bb)[values(bb)>=upp] <- upp
          values(bb)[values(bb)<=lwp] <- lwp
          plot(bb,zlim=c(lwp,upp),box=F,bty = "n",xaxt = "n", yaxt = "n",
               col=mycol,colNA='grey45',legend=F)
          title(paste('C6_',tiles[i],'_',year[j],'_',phe_met[vari],sep=''),line=0.25,cex.main=0.9)
          plot(bb,legend.only=T,col=mycol,zlim=c(lwp,upp),
               legend.width=1.5,legend.shrink=0.6,horiz=T,
               smallplot=c(0.15,0.8,0.11,0.14),
               axis.args=list(at=leg.int,cex.axis=0.8,font=2,
                              labels=c(paste('<',lwp,sep=''),leg.int[2:4],paste('>',upp,sep=''))))          
        }else{
          values(bb)[values(bb)>=upp] <- upp
          values(bb)[values(bb)<=lwp] <- lwp
          plot(bb,zlim=c(lwp,upp),box=F,bty = "n",xaxt = "n", yaxt = "n",
               col=mycol,colNA='grey45',legend=F)
          title(paste('C6_',tiles[i],'_',year[j],'_',phe_met[vari],sep=''),line=0.25,cex.main=0.9)
        }
      }  
    } 
    
  # Pheno-metrics
  for(vari in c(5,8,6,13,7,1)){  
    par(mfcol=c(3,5),mgp=c(0,0,0),oma=c(0.2,0.2,0.2,0.2),mar=c(1,0,1,0))  
    for(j in 1:15){
        mycol <- brewer.pal(11,'Spectral')
        mycol <- rev(colorRampPalette(mycol)(1000))
        
        path <- paste('/projectnb/modislc/users/dsm/eval_modis_lc_061917/MCD12I6/',phe_met[vari],'/',phe_met[vari],'_',tiles[i],'_',year[j],sep='')
        ncol <- 2400
        nrow <- 2400
        nbands <- 2
        cnt <- ncol*nrow*nbands
        data <- readBin(path,what="integer",n=cnt,size=2,endian="little")
        data[data>30000] <- NA
        data1 <- data - as.numeric(as.Date(paste(year[j]-1,'-12-31',sep='')))
        data2 <- array(data1,c(nbands, ncol, nrow))
        data2 <- aperm(data2, c(3,2,1)) #for transposing
        aa <- brick(data2)
        bb <- aa[[1]]
        cc <- getValues(bb)      
        
        if(sum(is.na(cc))>(2400*2400-10)){
          plot(raster(matrix(1,1,1)),box=F,bty = "n",xaxt = "n", yaxt = "n",
               col='grey45',legend=F)
          title(paste('C6_',tiles[i],'_',year[j],'_',phe_met[vari],'_NA',sep=''),line=0.25,cex.main=0.8)
        }else{
          if(j==1){
            upp <- round(quantile(bb,0.98,na.rm=T))
            lwp <- round(quantile(bb,0.02,na.rm=T)) 
            leg.int <- round(seq(lwp,upp,length.out = 5))  
            values(bb)[values(bb)>=upp] <- upp
            values(bb)[values(bb)<=lwp] <- lwp
            plot(bb,zlim=c(lwp,upp),box=F,bty = "n",xaxt = "n", yaxt = "n",
                 col=mycol,colNA='grey45',legend=F)
            title(paste('C6_',tiles[i],'_',year[j],'_',phe_met[vari],sep=''),line=0.25,cex.main=0.9)
            plot(bb,legend.only=T,col=mycol,zlim=c(lwp,upp),
                 legend.width=1.5,legend.shrink=0.6,horiz=T,
                 smallplot=c(0.15,0.8,0.11,0.14),
                 axis.args=list(at=leg.int,cex.axis=0.8,font=2,
                                labels=c(paste('<',lwp,sep=''),leg.int[2:4],paste('>',upp,sep=''))))          
          }else{
            values(bb)[values(bb)>=upp] <- upp
            values(bb)[values(bb)<=lwp] <- lwp
            plot(bb,zlim=c(lwp,upp),box=F,bty = "n",xaxt = "n", yaxt = "n",
                 col=mycol,colNA='grey45',legend=F)
            title(paste('C6_',tiles[i],'_',year[j],'_',phe_met[vari],sep=''),line=0.25,cex.main=0.9)
          }
        }  
      } 
    
    # Boxplot/ NA fraction
    mycol <- brewer.pal(11,'Spectral')
    mycol <- colorRampPalette(mycol)(15)
    temp1 <- matrix(NA,2400*2400,15)
    nplan <- matrix(NA,1,15)
    ee <- getValues(lct[[i]])
    plan <- sum(table(ee)[c(-1,-16,-17)])
    
    par(mfcol=c(2,1),mgp=c(2,1,0),oma=c(1,1,1,1),mar=c(2,3,2,0))  
    for(j in 1:15){
      # C6
      path <- paste('/projectnb/modislc/users/dsm/eval_modis_lc_061917/MCD12I6/',phe_met[vari],'/',phe_met[vari],'_',tiles[i],'_',year[j],sep='')
      ncol <- 2400
      nrow <- 2400
      nbands <- 2
      cnt <- ncol*nrow*nbands
      data <- readBin(path,what="integer",n=cnt,size=2,endian="little")
      data[data>30000] <- NA
      data1 <- data - as.numeric(as.Date(paste(year[j]-1,'-12-31',sep='')))
      data2 <- array(data1,c(nbands, ncol, nrow))
      data2 <- aperm(data2, c(3,2,1)) #for transposing
      aa <- brick(data2)
      bb <- aa[[1]]
      values(bb)[values(lct[[i]]==0)] <- NA
      cc <- getValues(bb)  
      
      nplan[1,j] <- 1-sum(!is.na(cc))/plan
      temp1[,j] <- cc
    }  
    
    if(sum(!is.na(temp1))<10){
      plot.new()
      title(paste('C6_',tiles[i],'_',phe_met[vari],'_NA',sep=''),line=0.3,cex.main=0.9)
      
      barplot(nplan,ylim=c(0,1),name=year,cex.axis=0.8,cex.name=0.8)
      title(paste('C6_',tiles[i],'_',phe_met[vari],sep=''),line=0.3,cex.main=0.9)
      mtext('Fraction of NA (land-only)',2,line=2.5,cex=0.8)
      box()
    }else{
      boxplot(temp1,outline=F,col=rev(mycol),axes=F,frame.plot=T,las=2)
      axis(1,at=seq(1,15),year,cex.axis=0.8)
      axis(2,at=seq(-600,600,50),cex.axis=0.8)
      mtext('DOY',2,line=2.5,cex=0.8)
      title(paste('C6_',tiles[i],'_',phe_met[vari],sep=''),line=0.3,cex.main=0.9)
      
      barplot(nplan,ylim=c(0,1),name=year,cex.axis=0.8,cex.name=0.8)
      title(paste('C6_',tiles[i],'_',phe_met[vari],sep=''),line=0.3,cex.main=0.9)
      mtext('Fraction of NA (land-only)',2,line=2.5,cex=0.8)
      box()
    }
    #################################
  } 

  # Length of growing season
    par(mfcol=c(3,5),mgp=c(0,0,0),oma=c(0.2,0.2,0.2,0.2),mar=c(1,0,1,0))  
    for(j in 1:15){
      mycol <- brewer.pal(11,'Spectral')
      mycol <- rev(colorRampPalette(mycol)(1000))
      
      vari <- 1 # Dormancy
      path <- paste('/projectnb/modislc/users/dsm/eval_modis_lc_061917/MCD12I6/',phe_met[vari],'/',phe_met[vari],'_',tiles[i],'_',year[j],sep='')
      ncol <- 2400
      nrow <- 2400
      nbands <- 2
      cnt <- ncol*nrow*nbands
      data <- readBin(path,what="integer",n=cnt,size=2,endian="little")
      data[data>30000] <- NA
      data1 <- data - as.numeric(as.Date(paste(year[j]-1,'-12-31',sep='')))
      data2 <- array(data1,c(nbands, ncol, nrow))
      data2 <- aperm(data2, c(3,2,1)) #for transposing
      aa <- brick(data2)
      bb1 <- aa[[1]]    
      
      vari <- 5 # Greenup
      path <- paste('/projectnb/modislc/users/dsm/eval_modis_lc_061917/MCD12I6/',phe_met[vari],'/',phe_met[vari],'_',tiles[i],'_',year[j],sep='')
      ncol <- 2400
      nrow <- 2400
      nbands <- 2
      cnt <- ncol*nrow*nbands
      data <- readBin(path,what="integer",n=cnt,size=2,endian="little")
      data[data>30000] <- NA
      data1 <- data - as.numeric(as.Date(paste(year[j]-1,'-12-31',sep='')))
      data2 <- array(data1,c(nbands, ncol, nrow))
      data2 <- aperm(data2, c(3,2,1)) #for transposing
      aa <- brick(data2)
      bb2 <- aa[[1]] 
      
      bb <- bb1 - bb2
      cc <- getValues(bb)
      
      if(sum(is.na(cc))>(2400*2400-10)){
        plot(raster(matrix(1,1,1)),box=F,bty = "n",xaxt = "n", yaxt = "n",
             col='grey45',legend=F)
        title(paste('C6_',tiles[i],'_',year[j],'_LofGS_NA',sep=''),line=0.25,cex.main=0.8)
      }else{
        if(j==1){
          upp <- round(quantile(bb,0.98,na.rm=T))
          lwp <- round(quantile(bb,0.02,na.rm=T)) 
          leg.int <- round(seq(lwp,upp,length.out = 5))  
          values(bb)[values(bb)>=upp] <- upp
          values(bb)[values(bb)<=lwp] <- lwp
          plot(bb,zlim=c(lwp,upp),box=F,bty = "n",xaxt = "n", yaxt = "n",
               col=mycol,colNA='grey45',legend=F)
          title(paste('C6_',tiles[i],'_',year[j],'_LofGS',sep=''),line=0.25,cex.main=0.9)
          plot(bb,legend.only=T,col=mycol,zlim=c(lwp,upp),
               legend.width=1.5,legend.shrink=0.6,horiz=T,
               smallplot=c(0.15,0.8,0.11,0.14),
               axis.args=list(at=leg.int,cex.axis=0.8,font=2,
                              labels=c(paste('<',lwp,sep=''),leg.int[2:4],paste('>',upp,sep=''))))          
        }else{
          values(bb)[values(bb)>=upp] <- upp
          values(bb)[values(bb)<=lwp] <- lwp
          plot(bb,zlim=c(lwp,upp),box=F,bty = "n",xaxt = "n", yaxt = "n",
               col=mycol,colNA='grey45',legend=F)
          title(paste('C6_',tiles[i],'_',year[j],'_LofGS',sep=''),line=0.25,cex.main=0.9)
        }
      }  
    }  

    # Boxplot/ NA fraction
    mycol <- brewer.pal(11,'Spectral')
    mycol <- colorRampPalette(mycol)(15)
    temp1 <- matrix(NA,2400*2400,15)
    nplan <- matrix(NA,1,15)
    ee <- getValues(lct[[i]])
    plan <- sum(table(ee)[c(-1,-16,-17)])
    
    par(mfcol=c(2,1),mgp=c(2,1,0),oma=c(1,1,1,1),mar=c(2,3,2,0))  
    for(j in 1:15){
      vari <- 1 # Dormancy
      path <- paste('/projectnb/modislc/users/dsm/eval_modis_lc_061917/MCD12I6/',phe_met[vari],'/',phe_met[vari],'_',tiles[i],'_',year[j],sep='')
      ncol <- 2400
      nrow <- 2400
      nbands <- 2
      cnt <- ncol*nrow*nbands
      data <- readBin(path,what="integer",n=cnt,size=2,endian="little")
      data[data>30000] <- NA
      data1 <- data - as.numeric(as.Date(paste(year[j]-1,'-12-31',sep='')))
      data2 <- array(data1,c(nbands, ncol, nrow))
      data2 <- aperm(data2, c(3,2,1)) #for transposing
      aa <- brick(data2)
      bb1 <- aa[[1]]    
      
      vari <- 5 # Greenup
      path <- paste('/projectnb/modislc/users/dsm/eval_modis_lc_061917/MCD12I6/',phe_met[vari],'/',phe_met[vari],'_',tiles[i],'_',year[j],sep='')
      ncol <- 2400
      nrow <- 2400
      nbands <- 2
      cnt <- ncol*nrow*nbands
      data <- readBin(path,what="integer",n=cnt,size=2,endian="little")
      data[data>30000] <- NA
      data1 <- data - as.numeric(as.Date(paste(year[j]-1,'-12-31',sep='')))
      data2 <- array(data1,c(nbands, ncol, nrow))
      data2 <- aperm(data2, c(3,2,1)) #for transposing
      aa <- brick(data2)
      bb2 <- aa[[1]] 
      
      bb <- bb1 - bb2
      cc <- getValues(bb) 
      
      nplan[1,j] <- 1-sum(!is.na(cc))/plan
      temp1[,j] <- cc
    }  
    
    if(sum(!is.na(temp1))<10){
      plot.new()
      title(paste('C6_',tiles[i],'_LofGS_NA',sep=''),line=0.3,cex.main=0.9)
      
      barplot(nplan,ylim=c(0,1),name=year,cex.axis=0.8,cex.name=0.8)
      title(paste('C6_',tiles[i],'_',phe_met[vari],sep=''),line=0.3,cex.main=0.9)
      mtext('Fraction of NA (land-only)',2,line=2.5,cex=0.8)
      box()
    }else{
      boxplot(temp1,outline=F,col=rev(mycol),axes=F,frame.plot=T,las=2)
      axis(1,at=seq(1,15),year,cex.axis=0.8)
      axis(2,at=seq(-600,600,50),cex.axis=0.8)
      mtext('DOY',2,line=2.5,cex=0.8)
      title(paste('C6_',tiles[i],'_Length of growing season',sep=''),line=0.3,cex.main=0.9)
      
      barplot(nplan,ylim=c(0,1),name=year,cex.axis=0.8,cex.name=0.8)
      title(paste('C6_',tiles[i],'_Length of growing season',sep=''),line=0.3,cex.main=0.9)
      mtext('Fraction of NA (land-only)',2,line=2.5,cex=0.8)
      box()
    }

## C6 vs C5
  # Pheno-metrics - raster
  for(j in 1:14){
    
    mycol <- brewer.pal(11,'Spectral')
    mycol <- rev(colorRampPalette(mycol)(1000))
    par(mfcol=c(2,4),mgp=c(0,0,0),oma=c(0.2,0.2,0.2,0.2),mar=c(1,0,1,0))
    
    # C6
    for(vari in c(5,6,13,1)){  
      
      path <- paste('/projectnb/modislc/users/dsm/eval_modis_lc_061917/MCD12I6/',phe_met[vari],'/',phe_met[vari],'_',tiles[i],'_',year[j],sep='')
      ncol <- 2400
      nrow <- 2400
      nbands <- 2
      cnt <- ncol*nrow*nbands
      data <- readBin(path,what="integer",n=cnt,size=2,endian="little")
      data[data>30000] <- NA
      data1 <- data - as.numeric(as.Date(paste(year[j]-1,'-12-31',sep='')))
      data2 <- array(data1,c(nbands, ncol, nrow))
      data2 <- aperm(data2, c(3,2,1)) #for transposing
      aa <- brick(data2)
      bb <- aa[[1]]
      cc <- getValues(bb)      
      
      if(sum(is.na(cc))>(2400*2400-10)){
        plot(raster(matrix(1,1,1)),box=F,bty = "n",xaxt = "n", yaxt = "n",
             col='grey45',legend=F)
        title(paste('C6_',tiles[i],'_',year[j],'_',phe_met[vari],'_NA',sep=''),line=0.25,cex.main=0.8)
      }else{
        if(vari==5){
          upp <- 365
          lwp <- 1
          leg.int <- c(lwp,50,150,250,upp)  
          values(bb)[values(bb)>=upp] <- upp
          values(bb)[values(bb)<=lwp] <- lwp
          plot(bb,zlim=c(lwp,upp),box=F,bty = "n",xaxt = "n", yaxt = "n",
               col=mycol,colNA='grey45',legend=F)
          title(paste('C6_',tiles[i],'_',year[j],'_',phe_met[vari],sep=''),line=0.25,cex.main=0.9)
          plot(bb,legend.only=T,col=mycol,zlim=c(lwp,upp),
               legend.width=1.5,legend.shrink=0.6,horiz=T,
               smallplot=c(0.15,0.8,0.11,0.14),
               axis.args=list(at=leg.int,cex.axis=0.8,font=2,
                              labels=c(paste('<',lwp,sep=''),leg.int[2:4],paste('>',upp,sep=''))))          
        }else{
          values(bb)[values(bb)>=upp] <- upp
          values(bb)[values(bb)<=lwp] <- lwp
          plot(bb,zlim=c(lwp,upp),box=F,bty = "n",xaxt = "n", yaxt = "n",
               col=mycol,colNA='grey45',legend=F)
          title(paste('C6_',tiles[i],'_',year[j],'_',phe_met[vari],sep=''),line=0.25,cex.main=0.9)
        }
      }
    }
    
    # C5
    for(vari in c(5,6,13,1)){
      
      if(i>287 & j==2 & vari==13){
        plot(raster(matrix(1,1,1)),box=F,bty = "n",xaxt = "n", yaxt = "n",
             col='grey45',legend=F)
        title(paste('C5_',tiles[i],'_',year[j],'_',phe_met[vari],'_NA',sep=''),line=0.25,cex.main=0.8)
      }else{
        path <- paste('/projectnb/modislc/data/mcd12_out/phen_out/c5_hdf1/',year[j],'/001/',phe_met[vari],'/',sep='') 
        bb <- raster(paste(path,'Phe',substr(phe_met[vari],1,3),'.A',year[j],'001.',tiles[i],'.hdf',sep=''))
        cc <- getValues(bb)
        cc <- as.numeric(as.Date(cc,origin='2000-01-01')) - as.numeric(as.Date(paste(year[j]-1,'-12-31',sep=''),origin='2000-01-01'))
        bb <- raster(matrix(cc,2400,2400,byrow=T))
        values(bb)[values(lct[[i]]==0)] <- NA
        
        if(sum(is.na(cc))>(2400*2400-10)){
          plot(raster(matrix(1,1,1)),box=F,bty = "n",xaxt = "n", yaxt = "n",
               col='grey45',legend=F)
          title(paste('C5_',tiles[i],'_',year[j],'_',phe_met[vari],'_NA',sep=''),line=0.25,cex.main=0.8)
        }else{
          values(bb)[values(bb)>=upp] <- upp
          values(bb)[values(bb)<=lwp] <- lwp
          plot(bb,zlim=c(lwp,upp),box=F,bty = "n",xaxt = "n", yaxt = "n",
               col=mycol,colNA='grey45',legend=F)
          title(paste('C5_',tiles[i],'_',year[j],'_',phe_met[vari],sep=''),line=0.25,cex.main=0.9)
        } 
      } 
    }  
    
  }
  
  # Pheno-metrics - boxplot/ NA fraction
  for(vari in c(5,6,13,1)){
    
    # Boxplot/ NA fraction
    mycol <- brewer.pal(11,'Spectral')
    mycol <- colorRampPalette(mycol)(15)
    temp1 <- matrix(NA,2400*2400,15)
    temp2 <- matrix(NA,2400*2400,15)
    nplan <- matrix(NA,2,15)
    ee <- getValues(lct[[i]])
    plan <- sum(table(ee)[c(-1,-16,-17)])
    
    for(j in 1:14){      
      # C6
      path <- paste('/projectnb/modislc/users/dsm/eval_modis_lc_061917/MCD12I6/',phe_met[vari],'/',phe_met[vari],'_',tiles[i],'_',year[j],sep='')
      ncol <- 2400
      nrow <- 2400
      nbands <- 2
      cnt <- ncol*nrow*nbands
      data <- readBin(path,what="integer",n=cnt,size=2,endian="little")
      data[data>30000] <- NA
      data1 <- data - as.numeric(as.Date(paste(year[j]-1,'-12-31',sep='')))
      data2 <- array(data1,c(nbands, ncol, nrow))
      data2 <- aperm(data2, c(3,2,1)) #for transposing
      aa <- brick(data2)
      bb <- aa[[1]]
      values(bb)[values(lct[[i]]==0)] <- NA
      cc <- getValues(bb)  
      
      nplan[1,j] <- 1-sum(!is.na(cc))/plan
      temp1[,j] <- cc
      
      # C5
      if(i>287 & j==2 & vari==13){
        nplan[2,j] <- NA
        temp2[,j] <- NA
      }else{
        path <- paste('/projectnb/modislc/data/mcd12_out/phen_out/c5_hdf1/',year[j],'/001/',phe_met[vari],'/',sep='') 
        bb <- raster(paste(path,'Phe',substr(phe_met[vari],1,3),'.A',year[j],'001.',tiles[i],'.hdf',sep=''))
        cc <- getValues(bb)
        cc <- as.numeric(as.Date(cc,origin='2000-01-01')) - as.numeric(as.Date(paste(year[j]-1,'-12-31',sep=''),origin='2000-01-01'))
        bb <- raster(matrix(cc,2400,2400,byrow=T))
        values(bb)[values(lct[[i]]==0)] <- NA
        cc <- getValues(bb)
        
        nplan[2,j] <- 1-sum(!is.na(cc))/plan
        temp2[,j] <- cc
      }
    }
    
    par(mfcol=c(2,1),mgp=c(2,1,0),oma=c(1,1,1,1),mar=c(2,3,2,0)) 
    if(sum(!is.na(temp1))<10 & sum(!is.na(temp1))<10){
      plot.new()
      title(paste('C6vsC5_',tiles[i],'_',phe_met[vari],'_NA',sep=''),line=0.3,cex.main=0.9)
      
      barplot(nplan,ylim=c(0,1),name=year,cex.axis=0.8,cex.name=0.8)
      title(paste('C6vsC5_',tiles[i],'_',phe_met[vari],sep=''),line=0.3,cex.main=0.9)
      mtext('Fraction of NA (land-only)',2,line=2.5,cex=0.8)
      box()
    }else{
      boxplot(temp1[,1],temp2[,1],temp1[,2],temp2[,2],temp1[,3],temp2[,3],temp1[,4],temp2[,4],temp1[,5],temp2[,5],
              temp1[,6],temp2[,6],temp1[,7],temp2[,7],temp1[,8],temp2[,8],temp1[,9],temp2[,9],temp1[,10],temp2[,10],
              temp1[,11],temp2[,11],temp1[,12],temp2[,12],temp1[,13],temp2[,13],temp1[,14],temp2[,14],
              at=c(1,2,4,5,7,8,10,11,13,14,16,17,19,20,22,23,25,26,28,29,31,32,34,35,37,38,40,41),
              outline=F,col=rev(rep(mycol,each=2)),axes=F,frame.plot=T,las=2)
      axis(1,at=seq(1.3,40.5,3),year[1:14],cex.axis=0.8)
      axis(2,at=seq(-600,600,50),cex.axis=0.8)
      mtext('DOY',2,line=2.5,cex=0.8)
      title(paste('C6vsC5_',tiles[i],'_',phe_met[vari],sep=''),line=0.3,cex.main=0.9)
      
      barplot(nplan[,c(1:14)],beside=T,ylim=c(0,1),name=year[1:14],cex.axis=0.8,cex.name=0.8)
      title(paste('C6vsC5_',tiles[i],'_',phe_met[vari],sep=''),line=0.3,cex.main=0.9)
      mtext('Fraction of NA (land-only)',2,line=2.5,cex=0.8)
      box()
    }
  }  


  dev.off() 
}

 

