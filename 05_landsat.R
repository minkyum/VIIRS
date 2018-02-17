rm(list = ls())
# loading required libraries 
library(sp)
library(raster)
library(gdalUtils)

## Get dates
### 1. Load land cover type (MCD12Q1) 
ext1204 <- raster('/projectnb/modislc/users/mkmoon/VIIRS/NLCD_crosswalked_IGBP/nlcd_igbp_500m/h12v04.bsq')
ext1104 <- raster('/projectnb/modislc/users/mkmoon/VIIRS/NLCD_crosswalked_IGBP/nlcd_igbp_500m/h11v04.bsq')
ext0805 <- raster('/projectnb/modislc/users/mkmoon/VIIRS/NLCD_crosswalked_IGBP/nlcd_igbp_500m/h08v05.bsq')

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
      cor <- '+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'
      lc[[e]] = raster(mcd12q1_full_path,crs=ext1204@crs)
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
    rast[[i]] <- brick(data2,
                       xmn=ext1204@extent@xmin,
                       xmx=ext1204@extent@xmax,
                       ymn=ext1204@extent@ymin,
                       ymx=ext1204@extent@ymax,
                       crs=ext1204@crs)
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
    aa <- brick(data2,
                xmn=ext1204@extent@xmin,
                xmx=ext1204@extent@xmax,
                ymn=ext1204@extent@ymin,
                ymx=ext1204@extent@ymax,
                crs=ext1204@crs)
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
      rast[[i]] <- brick(data2,
                         xmn=ext0805@extent@xmin,
                         xmx=ext0805@extent@xmax,
                         ymn=ext0805@extent@ymin,
                         ymx=ext0805@extent@ymax,
                         crs=ext1204@crs)  
    }else{
      rast[[i]] <- brick(data2,
                         xmn=ext1104@extent@xmin,
                         xmx=ext1104@extent@xmax,
                         ymn=ext1104@extent@ymin,
                         ymx=ext1104@extent@ymax,
                         crs=ext1204@crs)
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
      aa <- brick(data2,
                  xmn=ext0805@extent@xmin,
                  xmx=ext0805@extent@xmax,
                  ymn=ext0805@extent@ymin,
                  ymx=ext0805@extent@ymax,
                  crs=ext1204@crs)  
    }else{
      aa <- brick(data2,
                  xmn=ext1104@extent@xmin,
                  xmx=ext1104@extent@xmax,
                  ymn=ext1104@extent@ymin,
                  ymx=ext1104@extent@ymax,
                  crs=ext1204@crs)
    }
    rast[[vari]] <- aa[[1]]  
  }
  if(tt==1){
    mm0805 <- rast  
  }else{
    mm1104 <- rast  
  }
}


### Landsat Phenology
args <- commandArgs()
print(args)

ss <- as.numeric(substr(args[3],1,1))
tt <- as.numeric(substr(args[3],2,3))

scene <- c('mwp_harvard','ah_bartlett','mwp_cary','ah_hubbard','mwp_proctor','mws_boundary_waters')
scene <- scene[ss]
load(paste('/projectnb/modislc/projects/te_phenology/landsat_stacks/',scene,'/phenology/pheno_nodist_mat',sep=''))  
lan_phe <- pheno_mat

# Extent for Landsat Scenes
setwd(paste('/projectnb/modislc/projects/te_phenology/landsat_stacks/',scene,'/NLCD',sep=''))  
nlcd_scene <- raster('nlcd_clip.bip')       
nlcd_scene_vals <- getValues(nlcd_scene)
tmp_vals <- matrix(NA,ncell(nlcd_scene),1)
tmp_vals[lan_phe[,1],1] <- lan_phe[,tt]
lan_phe_sp <- setValues(nlcd_scene,tmp_vals)
tmp_vals[lan_phe[,1],1] <- lan_phe[,73]
lan_phe_au <- setValues(nlcd_scene,tmp_vals)

# Get MODIS location
get_mod_xy <- function(in_coords) {
  ## so the actual ids we have been extracting are 1 indexed
  ## for this to work we need to subtract the base by 1 to make them 0 indexed
  temp_ids = in_coords - 1
  ## dimensions of our modis grid we are using to get the ids
  nrows = 12000
  ncols = 19200
  ## starting location of the upper left tile in the modis grid we are using
  tile_h_start = 8
  tile_v_start = 2
  ## convert the 0-indexed 1-d id array to x,y coordinates relative to the full 19200x12000 grid of pixels
  y_coords = floor(temp_ids/ncols)
  x_coords = temp_ids - (y_coords * ncols) 
  
  ## use the starting tileid to get the current tile and the pixels with a 0-indexed output
  tile_h = floor(x_coords/2400) + tile_h_start
  tile_v = floor(y_coords/2400) + tile_v_start
  pix_x = x_coords - ((tile_h - tile_h_start)*2400) 
  pix_y = y_coords - ((tile_v - tile_v_start)*2400) 
  
  out_dat = cbind(tile_h,tile_v,pix_x,pix_y)
  
  return(out_dat)
}

## Match MODIS and VIIRS to Landsat
# Get MODIS coordinates
if(ss==6){
  modcor <- data.frame(get_mod_xy(pheno_mat[,4]))
}else{
  modcor <- data.frame(get_mod_xy(pheno_mat[,4]))
  modcor$tile_h[modcor$tile_h!=12] <- NA
  modcor <- na.omit(modcor)  
}
# Get MODIS boundaries
modpix <- unique(modcor[,c('pix_x','pix_y')])
nsam <- 1000
# Extract LSPs for each MODIS pixel
lanbymod <- matrix(NA,nsam,250*9)
nlpix <- matrix(NA,nsam,1)
for(i in 1:nsam){
  
  temp <- rep(0)
  while(sum(!is.na(temp))<500){
    
    sam <- sample(1:nrow(modpix),1) 
    
    r1 <- which(modcor[,3]==(modpix[sam,1]-1) & modcor[,4]==(modpix[sam,2]-1))
    r2 <- which(modcor[,3]==(modpix[sam,1]-1) & modcor[,4]==(modpix[sam,2]-0))
    r3 <- which(modcor[,3]==(modpix[sam,1]-1) & modcor[,4]==(modpix[sam,2]+1))
    r4 <- which(modcor[,3]==(modpix[sam,1]-0) & modcor[,4]==(modpix[sam,2]-1))
    r5 <- which(modcor[,3]==(modpix[sam,1]-0) & modcor[,4]==(modpix[sam,2]-0))
    r6 <- which(modcor[,3]==(modpix[sam,1]-0) & modcor[,4]==(modpix[sam,2]+1))
    r7 <- which(modcor[,3]==(modpix[sam,1]+1) & modcor[,4]==(modpix[sam,2]-1))
    r8 <- which(modcor[,3]==(modpix[sam,1]+1) & modcor[,4]==(modpix[sam,2]-0))
    r9 <- which(modcor[,3]==(modpix[sam,1]+1) & modcor[,4]==(modpix[sam,2]+1))
    rr <- c(r1,r2,r3,r4,r5,r6,r7,r8,r9)
    
    temp <- pheno_mat[rr,tt] 
    temp[temp==0] <- NA
  }
  
  nlpix[i] <- sum(!is.na(temp))
  lanbymod[i,1:(length(rr))] <- temp
}

# Get Mid-Greenup and smoothing by 3 by 3 window
if(tt==41){
  if(ss==6){
    mmmid <- (mm1104[[2]]+mm1104[[1]])/2
    mmmid <- focal(mmmid,w=matrix(1,3,3),mean,na.rm=F)
    mmmid <- getValues(mmmid)
    vvmid <- (vv1104[[2]][[1]]+vv1104[[1]][[1]])/2
    vvmid <- focal(vvmid[[1]],w=matrix(1,3,3),mean,na.rm=F)
    vvmid <- getValues(vvmid[[1]])  
  }else{
    mmmid <- (mm120413[[2]]+mm120413[[1]])/2
    mmmid <- focal(mmmid,w=matrix(1,3,3),mean,na.rm=F)
    mmmid <- getValues(mmmid)
    vvmid <- (vv120413[[2]][[1]]+vv120413[[1]][[1]])/2
    vvmid <- focal(vvmid[[1]],w=matrix(1,3,3),mean,na.rm=F)
    vvmid <- getValues(vvmid[[1]])  
  }  
}else if(tt==40){
  if(ss==6){
    mmmid <- (mm1104[[2]]+mm1104[[1]])/2
    mmmid <- focal(mmmid,w=matrix(1,3,3),mean,na.rm=F)
    mmmid <- getValues(mmmid)
    vvmid <- (vv1104[[2]][[1]]+vv1104[[1]][[1]])/2
    vvmid <- focal(vvmid[[1]],w=matrix(1,3,3),mean,na.rm=F)
    vvmid <- getValues(vvmid[[1]])  
  }else{
    mmmid <- (mm120412[[2]]+mm120412[[1]])/2
    mmmid <- focal(mmmid,w=matrix(1,3,3),mean,na.rm=F)
    mmmid <- getValues(mmmid)
    vvmid <- (vv120412[[2]][[1]]+vv120412[[1]][[1]])/2
    vvmid <- focal(vvmid[[1]],w=matrix(1,3,3),mean,na.rm=F)
    vvmid <- getValues(vvmid[[1]])  
  } 
}else if(tt==73){
  if(ss==6){
    mmmid <- (mm1104[[4]]+mm1104[[3]])/2
    mmmid <- focal(mmmid,w=matrix(1,3,3),mean,na.rm=F)
    mmmid <- getValues(mmmid)
    vvmid <- (vv1104[[4]][[1]]+vv1104[[3]][[1]])/2
    vvmid <- focal(vvmid[[1]],w=matrix(1,3,3),mean,na.rm=F)
    vvmid <- getValues(vvmid[[1]])  
  }else{
    mmmid <- (mm120413[[4]]+mm120413[[3]])/2
    mmmid <- focal(mmmid,w=matrix(1,3,3),mean,na.rm=F)
    mmmid <- getValues(mmmid)
    vvmid <- (vv120413[[4]][[1]]+vv120413[[3]][[1]])/2
    vvmid <- focal(vvmid[[1]],w=matrix(1,3,3),mean,na.rm=F)
    vvmid <- getValues(vvmid[[1]])  
  }
}else{
  if(ss==6){
    mmmid <- (mm1104[[4]]+mm1104[[3]])/2
    mmmid <- focal(mmmid,w=matrix(1,3,3),mean,na.rm=F)
    mmmid <- getValues(mmmid)
    vvmid <- (vv1104[[4]][[1]]+vv1104[[3]][[1]])/2
    vvmid <- focal(vvmid[[1]],w=matrix(1,3,3),mean,na.rm=F)
    vvmid <- getValues(vvmid[[1]])  
  }else{
    mmmid <- (mm120412[[4]]+mm120412[[3]])/2
    mmmid <- focal(mmmid,w=matrix(1,3,3),mean,na.rm=F)
    mmmid <- getValues(mmmid)
    vvmid <- (vv120412[[4]][[1]]+vv120412[[3]][[1]])/2
    vvmid <- focal(vvmid[[1]],w=matrix(1,3,3),mean,na.rm=F)
    vvmid <- getValues(vvmid[[1]])  
  } 
}

pixnum <- modpix[sam,1]+2400*(modpix[sam,2]-1)

mmphe <- round(mmmid[pixnum])
vvphe <- round(vvmid[pixnum])

quant <- matrix(NA,nsam,2)
for(i in 1:nsam){
  if(sum(!is.na(lanbymod[i,]))<1000){
    quant[i,1] <- NA
    quant[i,2] <- NA
  }else{
    temp <- ecdf(lanbymod[i,])
    
    quant[i,1] <- temp(mmphe[i])
    quant[i,2] <- temp(vvphe[i])  
  }
}

setwd('/projectnb/modislc/users/mkmoon/VIIRS/R_data/')
if(tt==41){
  save(mmphe,vvphe,lanbymod,quant,nlpix,file=paste('landsat_c6_2013_sos_',scene,'.RData',sep=''))
}else if(tt==40){
  save(mmphe,vvphe,lanbymod,quant,nlpix,file=paste('landsat_c6_2012_sos_',scene,'.RData',sep=''))
}else if(tt==73){
  save(mmphe,vvphe,lanbymod,quant,nlpix,file=paste('landsat_c6_2013_eos_',scene,'.RData',sep=''))
}else{
  save(mmphe,vvphe,lanbymod,quant,nlpix,file=paste('landsat_c6_2012_eos_',scene,'.RData',sep=''))
}

