rm(list = ls())
# Loading required libraries 
library(sp)
library(raster)
library(gdalUtils)
library(h5)

args <- commandArgs()
print(args)

ss <- as.numeric(substr(args[3],1,1))
xx <- as.numeric(substr(args[3],2,2))
yy <- as.numeric(substr(args[3],3,3))


# Info for specific time and tile 
tile <- c('h12v04','h11v04','h08v05')
year <- c(2012,2013,2014)

tile <- tile[xx]
year <- year[yy]

if(ss==1){
  # VNP43I4 -VIIRS
  viirs_path <- paste('/projectnb/modislc/users/mkmoon/VIIRS/tiles/',tile,'/VIIRS/',year,sep='')
  
  if(xx==1 & yy==2){
    search_str <- paste('VNP43I4.A*.hdf',sep='')
  }else if(xx==1 & yy==3){
    search_str <- paste('VNP43I4.A*.hdf',sep='')
  }else{
    search_str <- paste('VNP43I4.A*.h5',sep='')  
  }
  
  files <- list.files(path=viirs_path,pattern=glob2rx(search_str),full.names=T,include.dirs=F)
  
  evi_vv <- matrix(NA,366,(2400*2400))
  for(i in c(1:length(files))){
      doy <- as.numeric(substr(files[i],76,78))
      if(year==2012 | xx==2 | xx==3){
        file <- h5file(name=files[i])
        dataset <- list.datasets(file)
        i1 <- readDataSet(file[dataset[[4]]])*0.0001
        i2 <- readDataSet(file[dataset[[5]]])*0.0001
        i1[i1>3] <- NA
        i2[i2>3] <- NA
        i1 <- raster(i1)
        i2 <- raster(i2)
      }else{
        sds <- get_subdatasets(files[i])
        i1 <- raster(sds[4])
        i2 <- raster(sds[5])
      }
      temp <- 2.5*((i2-i1)/(i2 + (2.4*i1) +1)) # calculate EVI2 
      temp <- getValues(temp)
      
      evi_vv[doy,] <- temp
    }
  setwd(paste('/projectnb/modislc/users/mkmoon/VIIRS/R_data/evi2_ts/VIIRS/',tile,'/',sep=''))
  for(dd in 1:240){
    temp <- evi_vv[,(dd*24000-23999):(dd*24000)]
    
    d_str <- sprintf('%03d',dd)
    write.csv(temp,file=paste('evi_vv_',tile,'_',year,'_',d_str,'.csv',sep=''))  
  }
  
}else{
  # MCD43A4 -MODIS
  modis_path <- '/projectnb/modislc/users/mkmoon/VIIRS/tiles/e4ftl01.cr.usgs.gov/MOTA/MCD43A4.006'
  search_str <- paste('*43A4.A',year,'*.',tile,'.006*',sep='')
  files <- list.files(path=modis_path,pattern=glob2rx(search_str),full.names=T,include.dirs=F,recursive=T)
  
  evi_mm <- matrix(NA,366,(2400*2400))
  for(i in 1:length(files)){
    doy <- as.numeric(substr(files[i],106,108))
    sds <- get_subdatasets(files[i])
      
    b1 <- raster(sds[[8]])
    b2 <- raster(sds[[9]])
    
    temp <- 2.5*((b2-b1)/(b2 + (2.4*b1) +1)) # calculate EVI2 
    temp <- getValues(temp)
        
    evi_mm[doy,] <- temp
  }
  setwd(paste('/projectnb/modislc/users/mkmoon/VIIRS/R_data/evi2_ts/MODIS/',tile,'/',sep=''))
  for(dd in 1:240){
    temp <- evi_mm[,(dd*24000-23999):(dd*24000)]
        
    d_str <- sprintf('%03d',dd)
    write.csv(evi_mm,file=paste('evi_mm_',tile,'_',year,'_',d_str,'.csv',sep=''))  
  }
}






