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
  # MCD43A4 -VIIRS
  viirs_path <- paste('/projectnb/modislc/users/mkmoon/VIIRS/tiles/',tile,'/VIIRS/',year,sep='')
  search_str <- paste('VNP43I4.A*.h5',sep='')
  files <- list.files(path=viirs_path,pattern=glob2rx(search_str),full.names=T,include.dirs=F)
  
  for(ii in 1:24){
    evi_vv <- matrix(NA,length(files),240000)

    for(i in c(1:length(files))){
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
      
      evi_vv[i,] <- temp[(ii*240000-239999):(ii*240000)]
    }
    
    setwd('/projectnb/modislc/users/mkmoon/VIIRS/R_data/evi2_ts/')
    d <- ii
    d_str <- sprintf('%03d',d)
    write.csv(evi_vv,file=paste('evi_vv_',tile,'_',year,'_',d_str,'.csv',sep=''))
  }
}else{
  if(xx==2 | xx==3){
    # MCD43A4 -MODIS
    modis_path <- '/projectnb/modislc/users/mkmoon/VIIRS/tiles/e4ftl01.cr.usgs.gov/MOTA/MCD43A4.006'
    files <- list.files(path=modis_path,full.names=T,include.dirs=F,recursive=T)
    
    
    for(ii in 1:24){
      evi_mm <- matrix(NA,365,240000)
      for(i in 1:365){
        d <- seq(1,365)
        d_str <- sprintf('%03d',d[i])
        
        dat <- as.Date(i,origin="2012-12-31")
        year <- substr(dat,1,4)
        mont <- substr(dat,6,7)
        days <- substr(dat,9,10)
        
        path <- paste(modis_path,year,'.',mont,'.',days,'/',sep='')
        s_str <- paste('*.A',year,d_str,'.',tile,'*.hdf',sep='')
        full_path <- list.files(path=path,pattern=glob2rx(s_str),full.names=T,include.dirs=F)
        sds <- get_subdatasets(full_path)
        
        b1 <- raster(sds[[8]])
        b2 <- raster(sds[[9]])
        
        temp <- 2.5*((b2-b1)/(b2 + (2.4*b1) +1)) # calculate EVI2 
        temp <- getValues(temp)
        
        evi_mm[i,] <- temp[(ii*240000-239999):(ii*240000)]
      }
      setwd('/projectnb/modislc/users/mkmoon/VIIRS/R_data/evi2_ts/')
      d <- ii
      d_str <- sprintf('%03d',d)
      write.csv(evi_mm,file=paste('evi_mm_',tile,'_',year,'_',d_str,'.csv',sep=''))
    }
  }else{
    # MCD43A4 - from Project cluster
    modis_path <- '/projectnb/modislc/data/mcd12_in/c6/mcd43a4'
    files <- list.files(path=modis_path,full.names=T,include.dirs=F,recursive=F)
    
    
    for(ii in 1:24){
      evi_mm <- matrix(NA,365,2400*100)
      for(i in 1:365){
        d <- seq(1,365)
        d_str <- sprintf('%03d',d[i])
        
        path <- paste(modis_path,year,'/',d_str,'/',sep='')
        s_str <- paste('*.A',year,d_str,'.',tile,'*.hdf',sep='')
        full_path <- list.files(path=path,pattern=glob2rx(s_str),full.names=T,include.dirs=F)
        sds <- get_subdatasets(full_path)
        
        b1 <- raster(sds[[8]])
        b2 <- raster(sds[[9]])
        
        temp <- 2.5*((b2-b1)/(b2 + (2.4*b1) +1)) # calculate EVI2 
        temp <- getValues(temp)
        
        evi_mm[i,] <- temp[(ii*240000-239999):(ii*240000)]
      }
      setwd('/projectnb/modislc/users/mkmoon/VIIRS/R_data/evi2_ts/')
      d <- ii
      d_str <- sprintf('%03d',d)
      write.csv(evi_mm,file=paste('evi_mm_',tile,'_',year,'_',d_str,'.csv',sep=''))  
    }
  }
}





