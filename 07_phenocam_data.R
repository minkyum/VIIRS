## First run 01 code

setwd('/projectnb/modislc/users/mkmoon/Phenocam/PhenoCam_V1_1511/data/')

# Extract site list
path <- '/projectnb/modislc/users/mkmoon/Phenocam/PhenoCam_V1_1511/data/'
search_str <- '*meta.txt'
files <- list.files(path=path,pattern=glob2rx(search_str),full.names=F,include.dirs=F)

# Get site coordiante
source('~/R/Codes/latlong2MODIStile.R', echo=TRUE)
site.name <- rep(0)
cor <- matrix(NA,length(files),6)
for(i in 1:length(files)){
  filename <- files[i]
  temp <- readChar(filename,file.info(filename)$size)
  
  lat <- sub('.*lat: ','',temp)
  lat <- sub('\n.*','',lat)
  
  lon <- sub('.*lon: ','',temp)
  lon <- sub('\n.*','',lon)
  
  cor[i,1] <- as.numeric(lat)
  cor[i,2] <- as.numeric(lon)
  
  site.name[i] <- sub('_.*','',files[i])
  
  # Get MODIS tile coordinates
  sam <- latlong2tile(lat,lon)
  cor[i,3:6] <- sam
}
sitecord <- data.frame(site.name,cor)  ## every sites from Dateset V1.0

# Extract sites located in 1204, 1104, and 0805 tiles only
sitecord1 <- subset(sitecord,X3==12&X4==4)
sitecord2 <- subset(sitecord,X3==11&X4==4)
sitecord3 <- subset(sitecord,X3==8&X4==5)
sitecord <- rbind(sitecord1,sitecord2,sitecord3)

# setwd('/projectnb/modislc/users/mkmoon/VIIRS/R_data/')
# write.csv(sitecord,file='phenocam_sites.csv')

# Extract 3 day transition data file list
search_str <- '*3day_transi*.csv'
files <- list.files(path=path,pattern=glob2rx(search_str),full.names=F,include.dirs=F)
sitefull <- sub('_3day.*','',files)
sites <- sub('_.*','',files)

# ROI vegetation type
vgt <- rep(0)
for(i in 1:length(sites)){
  vgt[i] <- sub(sites[i],'',files[i])
  vgt[i] <- substr(vgt[i],2,3)  
}

## Match site name and coordinate
cor <- as.matrix(sitecord[,2:7])
cord <- matrix(NA,length(sites),6)
diag <- rep(0)
for(i in 1:length(sites)){
  for(j in 1:nrow(sitecord)){
    if(sitecord$site.name[j]==sites[i]) cord[i,] <- cor[j,]  
  }
}
sitecord <- data.frame(sitefull,vgt,cord)
sitecord <- na.omit(sitecord)
colnames(sitecord) <- c('site','vegtype','lat','lon','htile','vtile','x','y')

## Extract 3 day transition data files
files <- rep(0)
for(i in 1:nrow(sitecord)){
  search_str <- paste(sitecord$site[i],'_3day_*nsition_dates.csv',sep='')
  files[i] <- list.files(path=path,pattern=glob2rx(search_str),full.names=F,include.dirs=F)
}
sitecord <- data.frame(files,sitecord)
colnames(sitecord) <- c('file','site','vegtype','lat','lon','htile','vtile','x','y')

## Extract 50% SOS and EOS dates
datesos <- matrix(NA,nrow(sitecord),16)
dateeos <- matrix(NA,nrow(sitecord),16)

for(ss in 1:nrow(sitecord)){
  dat <- read.csv(files[ss],skip=16,header=T)
  dat1 <- subset(dat,direction=='rising'&gcc_value=='gcc_90')
  dat2 <- subset(dat,direction=='falling'&gcc_value=='gcc_90')
  
  sos <- dat1$transition_50
  for(i in 1:length(sos)){
    yy <- as.numeric(substr(sos[(length(sos)-i+1)],1,4))
    cc <- yy-1999
    datesos[ss,cc] <- as.Date(sos[(length(sos)-i+1)])
  }
  
  eos <- dat2$transition_50
  for(i in 1:length(eos)){
    yy <- as.numeric(substr(eos[(length(eos)-i+1)],1,4))
    cc <- yy-1999
    dateeos[ss,cc] <- as.Date(eos[(length(eos)-i+1)])
  }
}

## Merge coordinate and dates
sos <- as.Date(datesos,origin='1970-1-1')
sos <- matrix(as.numeric(strftime(sos,format='%j')),nrow(sitecord),16)
colnames(sos) <- 2000:2015
phenosos <- cbind(sitecord,sos)
eos <- as.Date(dateeos,origin='1970-1-1')
eos <- matrix(as.numeric(strftime(eos,format='%j')),nrow(sitecord),16)
colnames(eos) <- 2000:2015
phenoeos <- cbind(sitecord,eos)


### MODIS data near Phenocam sites
get_lsp_sos <- function(sensor,year){
  phe <- matrix(NA,nrow(phenosos),9)
  
  for(i in 1:nrow(phenosos)){
    if(year==2013){
      if(sensor=='mm'){
        if(phenosos$htile[i]==12){
          rast <- (mm120413[[1]]+mm120413[[2]])/2
        }else if(phenosos$htile[i]==11){
          rast <- (mm1104[[1]]+mm1104[[2]])/2
        }else{
          rast <- (mm0805[[1]]+mm0805[[2]])/2
        }
      }else{
        if(phenosos$htile[i]==12){
          rast <- (vv120413[[1]][[1]]+vv120413[[2]][[1]])/2
        }else if(phenosos$htile[i]==11){
          rast <- (vv1104[[1]][[1]]+vv1104[[2]][[1]])/2
        }else{
          rast <- (vv0805[[1]][[1]]+vv0805[[2]][[1]])/2
        }
      }  
    }else if(year==2012){
      if(sensor=='mm'){
        if(phenosos$htile[i]==12){
          rast <- (mm120412[[1]]+mm120412[[2]])/2
        }else{
          rast <- NA
        }
      }else{
        if(phenosos$htile[i]==12){
          rast <- (vv120412[[1]][[1]]+vv120412[[2]][[1]])/2
        }else{
          rast <- NA
        }
      }
    }else{
      if(sensor=='mm'){
        if(phenosos$htile[i]==12){
          rast <- (mm120414[[1]]+mm120414[[2]])/2
        }else{
          rast <- NA
        }
      }else{
        if(phenosos$htile[i]==12){
          rast <- (vv120414[[1]][[1]]+vv120414[[2]][[1]])/2
        }else{
          rast <- NA
        }
      }
    }
    
    if(year==2013){
      phe[i,1] <- rast[round(phenosos$x[i])-1,round(phenosos$y[i])-1]
      phe[i,2] <- rast[round(phenosos$x[i])-1,round(phenosos$y[i])-0]
      phe[i,3] <- rast[round(phenosos$x[i])-1,round(phenosos$y[i])+1]
      phe[i,4] <- rast[round(phenosos$x[i])-0,round(phenosos$y[i])-1]
      phe[i,5] <- rast[round(phenosos$x[i])-0,round(phenosos$y[i])+0]
      phe[i,6] <- rast[round(phenosos$x[i])-0,round(phenosos$y[i])+1]
      phe[i,7] <- rast[round(phenosos$x[i])+1,round(phenosos$y[i])-1]
      phe[i,8] <- rast[round(phenosos$x[i])+1,round(phenosos$y[i])+0]
      phe[i,9] <- rast[round(phenosos$x[i])+1,round(phenosos$y[i])+1]
    }else{
      if(phenosos$htile[i]==11|phenosos$htile[i]==8){
        phe[i,] <- NA
      }else{
        phe[i,1] <- rast[round(phenosos$x[i])-1,round(phenosos$y[i])-1]
        phe[i,2] <- rast[round(phenosos$x[i])-1,round(phenosos$y[i])-0]
        phe[i,3] <- rast[round(phenosos$x[i])-1,round(phenosos$y[i])+1]
        phe[i,4] <- rast[round(phenosos$x[i])-0,round(phenosos$y[i])-1]
        phe[i,5] <- rast[round(phenosos$x[i])-0,round(phenosos$y[i])+0]
        phe[i,6] <- rast[round(phenosos$x[i])-0,round(phenosos$y[i])+1]
        phe[i,7] <- rast[round(phenosos$x[i])+1,round(phenosos$y[i])-1]
        phe[i,8] <- rast[round(phenosos$x[i])+1,round(phenosos$y[i])+0]
        phe[i,9] <- rast[round(phenosos$x[i])+1,round(phenosos$y[i])+1]
      }
    }
  }
  return(phe)
}
get_lsp_eos <- function(sensor,year){
  phe <- matrix(NA,nrow(phenosos),9)
  
  for(i in 1:nrow(phenosos)){
    if(year==2013){
      if(sensor=='mm'){
        if(phenosos$htile[i]==12){
          rast <- (mm120413[[3]]+mm120413[[4]])/2
        }else if(phenosos$htile[i]==11){
          rast <- (mm1104[[3]]+mm1104[[4]])/2
        }else{
          rast <- (mm0805[[3]]+mm0805[[4]])/2
        }
      }else{
        if(phenosos$htile[i]==12){
          rast <- (vv120413[[3]][[1]]+vv120413[[4]][[1]])/2
        }else if(phenosos$htile[i]==11){
          rast <- (vv1104[[3]][[1]]+vv1104[[4]][[1]])/2
        }else{
          rast <- (vv0805[[3]][[1]]+vv0805[[4]][[1]])/2
        }
      }  
    }else if(year==2012){
      if(sensor=='mm'){
        if(phenosos$htile[i]==12){
          rast <- (mm120412[[3]]+mm120412[[4]])/2
        }else{
          rast <- NA
        }
      }else{
        if(phenosos$htile[i]==12){
          rast <- (vv120412[[3]][[1]]+vv120412[[4]][[1]])/2
        }else{
          rast <- NA
        }
      }
    }else{
      if(sensor=='mm'){
        if(phenosos$htile[i]==12){
          rast <- (mm120414[[3]]+mm120414[[4]])/2
        }else{
          rast <- NA
        }
      }else{
        if(phenosos$htile[i]==12){
          rast <- (vv120414[[3]][[1]]+vv120414[[4]][[1]])/2
        }else{
          rast <- NA
        }
      }
    }
    
    if(year==2013){
      phe[i,1] <- rast[round(phenosos$x[i])-1,round(phenosos$y[i])-1]
      phe[i,2] <- rast[round(phenosos$x[i])-1,round(phenosos$y[i])-0]
      phe[i,3] <- rast[round(phenosos$x[i])-1,round(phenosos$y[i])+1]
      phe[i,4] <- rast[round(phenosos$x[i])-0,round(phenosos$y[i])-1]
      phe[i,5] <- rast[round(phenosos$x[i])-0,round(phenosos$y[i])+0]
      phe[i,6] <- rast[round(phenosos$x[i])-0,round(phenosos$y[i])+1]
      phe[i,7] <- rast[round(phenosos$x[i])+1,round(phenosos$y[i])-1]
      phe[i,8] <- rast[round(phenosos$x[i])+1,round(phenosos$y[i])+0]
      phe[i,9] <- rast[round(phenosos$x[i])+1,round(phenosos$y[i])+1]
    }else{
      if(phenosos$htile[i]==11|phenosos$htile[i]==8){
        phe[i,] <- NA
      }else{
        phe[i,1] <- rast[round(phenosos$x[i])-1,round(phenosos$y[i])-1]
        phe[i,2] <- rast[round(phenosos$x[i])-1,round(phenosos$y[i])-0]
        phe[i,3] <- rast[round(phenosos$x[i])-1,round(phenosos$y[i])+1]
        phe[i,4] <- rast[round(phenosos$x[i])-0,round(phenosos$y[i])-1]
        phe[i,5] <- rast[round(phenosos$x[i])-0,round(phenosos$y[i])+0]
        phe[i,6] <- rast[round(phenosos$x[i])-0,round(phenosos$y[i])+1]
        phe[i,7] <- rast[round(phenosos$x[i])+1,round(phenosos$y[i])-1]
        phe[i,8] <- rast[round(phenosos$x[i])+1,round(phenosos$y[i])+0]
        phe[i,9] <- rast[round(phenosos$x[i])+1,round(phenosos$y[i])+1]
      }
    }
  }
  return(phe)
}
get_lct <- function(phenosos){
  lct <- matrix(NA,nrow(phenosos),25)
  
  for(i in 1:nrow(phenosos)){
    if(phenosos$htile[i]==12){
      rast <- lct.1204
    }else if(phenosos$htile[i]==11){
      rast <- lct.1104
    }else{
      rast <- lct.0805
    }
    
    lct[i,1] <- rast[round(phenosos$x[i])-2,round(phenosos$y[i])-2]
    lct[i,2] <- rast[round(phenosos$x[i])-2,round(phenosos$y[i])-1]
    lct[i,3] <- rast[round(phenosos$x[i])-2,round(phenosos$y[i])]
    lct[i,4] <- rast[round(phenosos$x[i])-2,round(phenosos$y[i])+1]
    lct[i,5] <- rast[round(phenosos$x[i])-2,round(phenosos$y[i])+2]
    lct[i,6] <- rast[round(phenosos$x[i])-1,round(phenosos$y[i])-2]
    lct[i,7] <- rast[round(phenosos$x[i])-1,round(phenosos$y[i])-1]
    lct[i,8] <- rast[round(phenosos$x[i])-1,round(phenosos$y[i])]
    lct[i,9] <- rast[round(phenosos$x[i])-1,round(phenosos$y[i])+1]
    lct[i,10] <- rast[round(phenosos$x[i])-1,round(phenosos$y[i])+2]
    lct[i,11] <- rast[round(phenosos$x[i])+0,round(phenosos$y[i])-2]
    lct[i,12] <- rast[round(phenosos$x[i])+0,round(phenosos$y[i])-1]
    lct[i,13] <- rast[round(phenosos$x[i])+0,round(phenosos$y[i])]
    lct[i,14] <- rast[round(phenosos$x[i])+0,round(phenosos$y[i])+1]
    lct[i,15] <- rast[round(phenosos$x[i])+0,round(phenosos$y[i])+2]
    lct[i,16] <- rast[round(phenosos$x[i])+1,round(phenosos$y[i])-2]
    lct[i,17] <- rast[round(phenosos$x[i])+1,round(phenosos$y[i])-1]
    lct[i,18] <- rast[round(phenosos$x[i])+1,round(phenosos$y[i])]
    lct[i,19] <- rast[round(phenosos$x[i])+1,round(phenosos$y[i])+1]
    lct[i,20] <- rast[round(phenosos$x[i])+1,round(phenosos$y[i])+2]
    lct[i,21] <- rast[round(phenosos$x[i])+2,round(phenosos$y[i])-2]
    lct[i,22] <- rast[round(phenosos$x[i])+2,round(phenosos$y[i])-1]
    lct[i,23] <- rast[round(phenosos$x[i])+2,round(phenosos$y[i])]
    lct[i,24] <- rast[round(phenosos$x[i])+2,round(phenosos$y[i])+1]
    lct[i,25] <- rast[round(phenosos$x[i])+2,round(phenosos$y[i])+2]
  }
  return(lct)
}
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# ## Land cover type
# lct <- get_lct(phenosos)
# LC <- apply(lct,1,Mode) 
# LCT <- data.frame(phenosos$vegtype,lct[,c(7,8,9,12,13,14,17,18,19)])
# 
# lctGR <- LCT[LCT$phenosos.vegtype=='GR',c(2:10)]
# lctDB <- LCT[LCT$phenosos.vegtype=='DB',c(2:10)]
# lctEN <- LCT[LCT$phenosos.vegtype=='EN',c(2:10)]
# lctAG <- LCT[LCT$phenosos.vegtype=='AG',c(2:10)]

# ## Sites
# numsites <- phenosos[phenosos$vegtype=='GR'|phenosos$vegtype=='DB'|
#                        phenosos$vegtype=='EN'|phenosos$vegtype=='AG',]
# numsites <- unique(phenosos[,c('lat','lon')])
# setwd('/projectnb/modislc/users/mkmoon/VIIRS/R_data/')
# save(numsites,file='phenocam_site_coordi.rda')

## SOS
phevv12s <- get_lsp_sos('vv',2012)
phemm12s <- get_lsp_sos('mm',2012)
phevv13s <- get_lsp_sos('vv',2013)
phemm13s <- get_lsp_sos('mm',2013)
phevv14s <- get_lsp_sos('vv',2014)
phemm14s <- get_lsp_sos('mm',2014)
## EOS
phevv12e <- get_lsp_eos('vv',2012)
phemm12e <- get_lsp_eos('mm',2012)
phevv13e <- get_lsp_eos('vv',2013)
phemm13e <- get_lsp_eos('mm',2013)
phevv14e <- get_lsp_eos('vv',2014)
phemm14e <- get_lsp_eos('mm',2014)

# Get mean values for 3 by 3 window
## SOS
sosvv12 <- apply(phevv12s[,c(7,8,9,12,13,14,17,18,19)],1,mean,na.rm=T)
sosmm12 <- apply(phemm12s[,c(7,8,9,12,13,14,17,18,19)],1,mean,na.rm=T)
sosvv13 <- apply(phevv13s[,c(7,8,9,12,13,14,17,18,19)],1,mean,na.rm=T)
sosmm13 <- apply(phemm13s[,c(7,8,9,12,13,14,17,18,19)],1,mean,na.rm=T)
sosvv14 <- apply(phevv14s[,c(7,8,9,12,13,14,17,18,19)],1,mean,na.rm=T)
sosmm14 <- apply(phemm14s[,c(7,8,9,12,13,14,17,18,19)],1,mean,na.rm=T)

sosvv12s <- apply(phevv12s[,c(7,8,9,12,13,14,17,18,19)],1,sd,na.rm=T)
sosmm12s <- apply(phemm12s[,c(7,8,9,12,13,14,17,18,19)],1,sd,na.rm=T)
sosvv13s <- apply(phevv13s[,c(7,8,9,12,13,14,17,18,19)],1,sd,na.rm=T)
sosmm13s <- apply(phemm13s[,c(7,8,9,12,13,14,17,18,19)],1,sd,na.rm=T)
sosvv14s <- apply(phevv14s[,c(7,8,9,12,13,14,17,18,19)],1,sd,na.rm=T)
sosmm14s <- apply(phemm14s[,c(7,8,9,12,13,14,17,18,19)],1,sd,na.rm=T)
## EOS
eosvv12 <- apply(phevv12e[,c(7,8,9,12,13,14,17,18,19)],1,mean,na.rm=T)
eosmm12 <- apply(phemm12e[,c(7,8,9,12,13,14,17,18,19)],1,mean,na.rm=T)
eosvv13 <- apply(phevv13e[,c(7,8,9,12,13,14,17,18,19)],1,mean,na.rm=T)
eosmm13 <- apply(phemm13e[,c(7,8,9,12,13,14,17,18,19)],1,mean,na.rm=T)
eosvv14 <- apply(phevv14e[,c(7,8,9,12,13,14,17,18,19)],1,mean,na.rm=T)
eosmm14 <- apply(phemm14e[,c(7,8,9,12,13,14,17,18,19)],1,mean,na.rm=T)

eosvv12s <- apply(phevv12e[,c(7,8,9,12,13,14,17,18,19)],1,sd,na.rm=T)
eosmm12s <- apply(phemm12e[,c(7,8,9,12,13,14,17,18,19)],1,sd,na.rm=T)
eosvv13s <- apply(phevv13e[,c(7,8,9,12,13,14,17,18,19)],1,sd,na.rm=T)
eosmm13s <- apply(phemm13e[,c(7,8,9,12,13,14,17,18,19)],1,sd,na.rm=T)
eosvv14s <- apply(phevv14e[,c(7,8,9,12,13,14,17,18,19)],1,sd,na.rm=T)
eosmm14s <- apply(phemm14e[,c(7,8,9,12,13,14,17,18,19)],1,sd,na.rm=T)

# Get Phenocam transtision dates
## SOS
sospc12 <- phenosos$'2012'
sospc13 <- phenosos$'2013'
sospc14 <- phenosos$'2014'
## EOS
eospc12 <- phenoeos$'2012'
eospc13 <- phenoeos$'2013'
eospc14 <- phenoeos$'2014'

sospc <- cbind(sospc12,sospc13,sospc14)
eospc <- cbind(eospc12,eospc13,eospc14)

sospcDB <- sospc[phenosos$vegtype=='DB',]
sospcEN <- sospc[phenosos$vegtype=='EN',]
sospcAG <- sospc[phenosos$vegtype=='AG',]
sospcGR <- sospc[phenosos$vegtype=='GR'|phenosos$vegtype=='SH',]

eospcDB <- eospc[phenosos$vegtype=='DB',]
eospcEN <- eospc[phenosos$vegtype=='EN',]
eospcAG <- eospc[phenosos$vegtype=='AG',]
eospcGR <- eospc[phenosos$vegtype=='GR'|phenosos$vegtype=='SH',]

## Phenocam vs. VIIRS
par(mfrow=c(2,4))
# SOS
sosvv12DB <- sosvv12[phenosos$vegtype=='DB']
sosvv13DB <- sosvv13[phenosos$vegtype=='DB']
sosvv14DB <- sosvv14[phenosos$vegtype=='DB']

sosvv12EN <- sosvv12[phenosos$vegtype=='EN']
sosvv13EN <- sosvv13[phenosos$vegtype=='EN']
sosvv14EN <- sosvv14[phenosos$vegtype=='EN']

sosvv12AG <- sosvv12[phenosos$vegtype=='AG']
sosvv13AG <- sosvv13[phenosos$vegtype=='AG']
sosvv14AG <- sosvv14[phenosos$vegtype=='AG']

sosvv12GR <- sosvv12[phenosos$vegtype=='GR'|phenosos$vegtype=='SH']
sosvv13GR <- sosvv13[phenosos$vegtype=='GR'|phenosos$vegtype=='SH']
sosvv14GR <- sosvv14[phenosos$vegtype=='GR'|phenosos$vegtype=='SH']

# EOS
eosvv12DB <- eosvv12[phenoeos$vegtype=='DB']
eosvv13DB <- eosvv13[phenoeos$vegtype=='DB']
eosvv14DB <- eosvv14[phenoeos$vegtype=='DB']

eosvv12EN <- eosvv12[phenoeos$vegtype=='EN']
eosvv13EN <- eosvv13[phenoeos$vegtype=='EN']
eosvv14EN <- eosvv14[phenoeos$vegtype=='EN']

eosvv12AG <- eosvv12[phenoeos$vegtype=='AG']
eosvv13AG <- eosvv13[phenoeos$vegtype=='AG']
eosvv14AG <- eosvv14[phenoeos$vegtype=='AG']

eosvv12GR <- eosvv12[phenoeos$vegtype=='GR'|phenosos$vegtype=='SH']
eosvv13GR <- eosvv13[phenoeos$vegtype=='GR'|phenosos$vegtype=='SH']
eosvv14GR <- eosvv14[phenoeos$vegtype=='GR'|phenosos$vegtype=='SH']

## Phenocam vs. MODIS
# SOS
sosmm12DB <- sosmm12[phenosos$vegtype=='DB']
sosmm13DB <- sosmm13[phenosos$vegtype=='DB']
sosmm14DB <- sosmm14[phenosos$vegtype=='DB']

sosmm12EN <- sosmm12[phenosos$vegtype=='EN']
sosmm13EN <- sosmm13[phenosos$vegtype=='EN']
sosmm14EN <- sosmm14[phenosos$vegtype=='EN']

sosmm12AG <- sosmm12[phenosos$vegtype=='AG']
sosmm13AG <- sosmm13[phenosos$vegtype=='AG']
sosmm14AG <- sosmm14[phenosos$vegtype=='AG']

sosmm12GR <- sosmm12[phenosos$vegtype=='GR'|phenosos$vegtype=='SH']
sosmm13GR <- sosmm13[phenosos$vegtype=='GR'|phenosos$vegtype=='SH']
sosmm14GR <- sosmm14[phenosos$vegtype=='GR'|phenosos$vegtype=='SH']

# EOS
eosmm12DB <- eosmm12[phenoeos$vegtype=='DB']
eosmm13DB <- eosmm13[phenoeos$vegtype=='DB']
eosmm14DB <- eosmm14[phenoeos$vegtype=='DB']

eosmm12EN <- eosmm12[phenoeos$vegtype=='EN']
eosmm13EN <- eosmm13[phenoeos$vegtype=='EN']
eosmm14EN <- eosmm14[phenoeos$vegtype=='EN']

eosmm12AG <- eosmm12[phenoeos$vegtype=='AG']
eosmm13AG <- eosmm13[phenoeos$vegtype=='AG']
eosmm14AG <- eosmm14[phenoeos$vegtype=='AG']

eosmm12GR <- eosmm12[phenoeos$vegtype=='GR'|phenosos$vegtype=='SH']
eosmm13GR <- eosmm13[phenoeos$vegtype=='GR'|phenosos$vegtype=='SH']
eosmm14GR <- eosmm14[phenoeos$vegtype=='GR'|phenosos$vegtype=='SH']

### plots
# setwd('/projectnb/modislc/users/mkmoon/VIIRS/figures/')
# pdf(file='phenocam.pdf',width=10,height=8.5)
mycol <- brewer.pal(5,'Set1')

par(mfrow=c(2,2),oma=c(1,1,1,1),mar=c(4,4,1,1),mgp=c(2.5,1,0))
# SOS
plot(sospcDB[,1],sosvv12DB,xlim=c(60,190),ylim=c(60,190),bg=mycol[1],pch=21,cex=1.8,
     cex.lab=1.5,cex.axis=1.5,axe=F,
     xlab='GCC 50% amplitude during rising',ylab='VIIRS SOS')
axis(1,at=seq(0,300,30),cex.axis=1.5);axis(2,at=seq(0,300,30),cex.axis=1.5);box(lty=1)
points(sospcDB[,2],sosvv13DB,bg=mycol[1],pch=21,cex=1.8);points(sospcDB[,3],sosvv14DB,bg=mycol[1],pch=21,cex=1.8)
points(sospcEN[,1],sosvv12EN,bg=mycol[2],pch=21,cex=1.8);points(sospcEN[,2],sosvv13EN,bg=mycol[2],pch=21,cex=1.8)
points(sospcEN[,3],sosvv14EN,bg=mycol[2],pch=21,cex=1.8);points(sospcAG[,1],sosvv12AG,bg=mycol[3],pch=21,cex=1.8)
points(sospcAG[,2],sosvv13AG,bg=mycol[3],pch=21,cex=1.8);points(sospcAG[,3],sosvv14AG,bg=mycol[3],pch=21,cex=1.8)
points(sospcGR[,1],sosvv12GR,bg=mycol[4],pch=21,cex=1.8);points(sospcGR[,2],sosvv13GR,bg=mycol[4],pch=21,cex=1.8)
points(sospcGR[,3],sosvv14GR,bg=mycol[4],pch=21,cex=1.8);abline(0,1,lty=5)
legend('bottomright',c('DB','EN','AG','GR/SH'),pch=21,
       pt.bg=c(mycol[1],mycol[2],mycol[3],mycol[4]),bty='n',cex=1.5)
text(65,180,'(a)',cex=1.5)

plot(sospcDB[,1],sosmm12DB,xlim=c(60,190),ylim=c(60,190),bg=mycol[1],pch=21,cex=1.8,
     cex.lab=1.5,cex.axis=1.5,axe=F,
     xlab='GCC 50% amplitude during rising',ylab='MODIS SOS')
axis(1,at=seq(0,300,30),cex.axis=1.5);axis(2,at=seq(0,300,30),cex.axis=1.5);box(lty=1)
points(sospcDB[,2],sosmm13DB,bg=mycol[1],pch=21,cex=1.8);points(sospcDB[,3],sosmm14DB,bg=mycol[1],pch=21,cex=1.8)
points(sospcEN[,1],sosmm12EN,bg=mycol[2],pch=21,cex=1.8);points(sospcEN[,2],sosmm13EN,bg=mycol[2],pch=21,cex=1.8)
points(sospcEN[,3],sosmm14EN,bg=mycol[2],pch=21,cex=1.8);points(sospcAG[,1],sosmm12AG,bg=mycol[3],pch=21,cex=1.8)
points(sospcAG[,2],sosmm13AG,bg=mycol[3],pch=21,cex=1.8);points(sospcAG[,3],sosmm14AG,bg=mycol[3],pch=21,cex=1.8)
points(sospcGR[,1],sosmm12GR,bg=mycol[4],pch=21,cex=1.8);points(sospcGR[,2],sosmm13GR,bg=mycol[4],pch=21,cex=1.8)
points(sospcGR[,3],sosmm14GR,bg=mycol[4],pch=21,cex=1.8);abline(0,1,lty=5)
text(65,180,'(b)',cex=1.5)

# EOS
plot(eospcDB[,1],eosvv12DB,xlim=c(90,350),ylim=c(90,350),bg=mycol[1],pch=21,cex=1.8,
     cex.lab=1.5,cex.axis=1.5,axe=F,
     xlab='GCC 50% amplitude during falling',ylab='VIIRS EOS')
axis(1,at=seq(50,400,75),cex.axis=1.5);axis(2,at=seq(50,400,75),cex.axis=1.5);box(lty=1)
points(eospcDB[,2],eosvv13DB,bg=mycol[1],pch=21,cex=1.8);points(eospcDB[,3],eosvv14DB,bg=mycol[1],pch=21,cex=1.8)
points(eospcEN[,1],eosvv12EN,bg=mycol[2],pch=21,cex=1.8);points(eospcEN[,2],eosvv13EN,bg=mycol[2],pch=21,cex=1.8)
points(eospcEN[,3],eosvv14EN,bg=mycol[2],pch=21,cex=1.8);points(eospcAG[,1],eosvv12AG,bg=mycol[3],pch=21,cex=1.8)
points(eospcAG[,2],eosvv13AG,bg=mycol[3],pch=21,cex=1.8);points(eospcAG[,3],eosvv14AG,bg=mycol[3],pch=21,cex=1.8)
points(eospcGR[,1],eosvv12GR,bg=mycol[4],pch=21,cex=1.8);points(eospcGR[,2],eosvv13GR,bg=mycol[4],pch=21,cex=1.8)
points(eospcGR[,3],eosvv14GR,bg=mycol[4],pch=21,cex=1.8);abline(0,1,lty=5)
text(100,330,'(c)',cex=1.5)

plot(eospcDB[,1],eosmm12DB,xlim=c(90,350),ylim=c(90,350),bg=mycol[1],pch=21,cex=1.8,
     cex.lab=1.5,cex.axis=1.5,axe=F,
     xlab='GCC 50% amplitude during falling',ylab='MODIS EOS')
axis(1,at=seq(50,400,75),cex.axis=1.5);axis(2,at=seq(50,400,75),cex.axis=1.5);box(lty=1)
points(eospcDB[,2],eosmm13DB,bg=mycol[1],pch=21,cex=1.8);points(eospcDB[,3],eosmm14DB,bg=mycol[1],pch=21,cex=1.8)
points(eospcEN[,1],eosmm12EN,bg=mycol[2],pch=21,cex=1.8);points(eospcEN[,2],eosmm13EN,bg=mycol[2],pch=21,cex=1.8)
points(eospcEN[,3],eosmm14EN,bg=mycol[2],pch=21,cex=1.8);points(eospcAG[,1],eosmm12AG,bg=mycol[3],pch=21,cex=1.8)
points(eospcAG[,2],eosmm13AG,bg=mycol[3],pch=21,cex=1.8);points(eospcAG[,3],eosmm14AG,bg=mycol[3],pch=21,cex=1.8)
points(eospcGR[,1],eosmm12GR,bg=mycol[4],pch=21,cex=1.8);points(eospcGR[,2],eosmm13GR,bg=mycol[4],pch=21,cex=1.8)
points(eospcGR[,3],eosmm14GR,bg=mycol[4],pch=21,cex=1.8);abline(0,1,lty=5)
text(100,330,'(d)',cex=1.5)

# dev.off()

### Stat table
statmat <- matrix(NA,8,14)
for(i in c(1,3,5,7)){
  # SOS
  if(i==1){
    mat <- matrix(c(sospcDB[,1],sospcDB[,2],sospcDB[,3],
                    sosvv12DB,sosvv13DB,sosvv14DB),length(sosvv12DB)*3,2)  
  }else if(i==3){
    mat <- matrix(c(sospcEN[,1],sospcEN[,2],sospcEN[,3],
                    sosvv12EN,sosvv13EN,sosvv14EN),length(sosvv12EN)*3,2)  
  }else if(i==5){
    mat <- matrix(c(sospcAG[,1],sospcAG[,2],sospcAG[,3],
                    sosvv12AG,sosvv13AG,sosvv14AG),length(sosvv12AG)*3,2)  
  }else{
    mat <- matrix(c(sospcGR[,1],sospcGR[,2],sospcGR[,3],
                    sosvv12GR,sosvv13GR,sosvv14GR),length(sosvv12GR)*3,2)  
  }
  mat <- na.omit(mat)
  statmat[i,1] <- sum(!is.na(mat[,1]))
  statmat[i,2] <- mean(mat[,1])
  statmat[i,3] <- sd(mat[,1])
  statmat[i,4] <- mean(mat[,2])
  statmat[i,5] <- sd(mat[,2])
#   reg <- summary(lm(mat[,2]~mat[,1]))
#   statmat[i,6] <- reg$r.squared
  statmat[i,6] <- mean((mat[,2]-mat[,1]),na.rm=T)
  statmat[i,7] <- sqrt(mean((mat[,1]-mat[,2])^2,na.rm=T))
  
  # EOS
  if(i==1){
    mat <- matrix(c(eospcDB[,1],eospcDB[,2],eospcDB[,3],
                    eosvv12DB,eosvv13DB,eosvv14DB),length(eosvv12DB)*3,2)  
  }else if(i==3){
    mat <- matrix(c(eospcEN[,1],eospcEN[,2],eospcEN[,3],
                    eosvv12EN,eosvv13EN,eosvv14EN),length(eosvv12EN)*3,2)  
  }else if(i==5){
    mat <- matrix(c(eospcAG[,1],eospcAG[,2],eospcAG[,3],
                    eosvv12AG,eosvv13AG,eosvv14AG),length(eosvv12AG)*3,2)  
  }else{
    mat <- matrix(c(eospcGR[,1],eospcGR[,2],eospcGR[,3],
                    eosvv12GR,eosvv13GR,eosvv14GR),length(eosvv12GR)*3,2)  
  }
  mat <- na.omit(mat)
  statmat[i,8] <- sum(!is.na(mat[,1]))
  statmat[i,9] <- mean(mat[,1])
  statmat[i,10] <- sd(mat[,1])
  statmat[i,11] <- mean(mat[,2])
  statmat[i,12] <- sd(mat[,2])
#   reg <- summary(lm(mat[,2]~mat[,1]))
#   statmat[i,13] <- reg$r.squared
  statmat[i,13] <- mean((mat[,2]-mat[,1]),na.rm=T)
  statmat[i,14] <- sqrt(mean((mat[,1]-mat[,2])^2,na.rm=T))
  
}
for(i in c(1,3,5,7)){
  # SOS
  if(i==1){
    mat <- matrix(c(sospcDB[,1],sospcDB[,2],sospcDB[,3],
                    sosmm12DB,sosmm13DB,sosmm14DB),length(sosmm12DB)*3,2)  
  }else if(i==3){
    mat <- matrix(c(sospcEN[,1],sospcEN[,2],sospcEN[,3],
                    sosmm12EN,sosmm13EN,sosmm14EN),length(sosmm12EN)*3,2)  
  }else if(i==5){
    mat <- matrix(c(sospcAG[,1],sospcAG[,2],sospcAG[,3],
                    sosmm12AG,sosmm13AG,sosmm14AG),length(sosmm12AG)*3,2)  
  }else{
    mat <- matrix(c(sospcGR[,1],sospcGR[,2],sospcGR[,3],
                    sosmm12GR,sosmm13GR,sosmm14GR),length(sosmm12GR)*3,2)  
  }
  mat <- na.omit(mat)
  statmat[(i+1),1] <- sum(!is.na(mat[,1]))
  statmat[(i+1),2] <- mean(mat[,1])
  statmat[(i+1),3] <- sd(mat[,1])
  statmat[(i+1),4] <- mean(mat[,2])
  statmat[(i+1),5] <- sd(mat[,2])
#   reg <- summary(lm(mat[,2]~mat[,1]))
#   statmat[(i+1),6] <- reg$adj.r.squared
  statmat[(i+1),6] <- mean((mat[,2]-mat[,1]),na.rm=T)
  statmat[(i+1),7] <- sqrt(mean((mat[,1]-mat[,2])^2,na.rm=T))
  
  # EOS
  if(i==1){
    mat <- matrix(c(eospcDB[,1],eospcDB[,2],eospcDB[,3],
                    eosmm12DB,eosmm13DB,eosmm14DB),length(eosmm12DB)*3,2)  
  }else if(i==3){
    mat <- matrix(c(eospcEN[,1],eospcEN[,2],eospcEN[,3],
                    eosmm12EN,eosmm13EN,eosmm14EN),length(eosmm12EN)*3,2)  
  }else if(i==5){
    mat <- matrix(c(eospcAG[,1],eospcAG[,2],eospcAG[,3],
                    eosmm12AG,eosmm13AG,eosmm14AG),length(eosmm12AG)*3,2)  
  }else{
    mat <- matrix(c(eospcGR[,1],eospcGR[,2],eospcGR[,3],
                    eosmm12GR,eosmm13GR,eosmm14GR),length(eosmm12GR)*3,2)  
  }
  mat <- na.omit(mat)
  statmat[(i+1),8] <- sum(!is.na(mat[,1]))
  statmat[(i+1),9] <- mean(mat[,1])
  statmat[(i+1),10] <- sd(mat[,1])
  statmat[(i+1),11] <- mean(mat[,2])
  statmat[(i+1),12] <- sd(mat[,2])
#   reg <- summary(lm(mat[,2]~mat[,1]))
#   statmat[(i+1),13] <- reg$adj.r.squared
  statmat[(i+1),13] <- mean((mat[,2]-mat[,1]),na.rm=T)
  statmat[(i+1),14] <- sqrt(mean((mat[,1]-mat[,2])^2,na.rm=T))
  
}

setwd('/projectnb/modislc/users/mkmoon/VIIRS/R_data/')
write.csv(statmat,file='phenocam.csv')







# ### Plot
# plot_phenocam <- function(sensor){
#   
#   if(sensor=='vv'){
#     xx <- c(sosvv12,sosvv13,sosvv14,eosvv12,eosvv13,eosvv14)
#     yy <- c(sospc12,sospc13,sospc14,eospc12,eospc13,eospc14)
#     xsd <- c(sosvv12s,sosvv13s,sosvv14s,eosvv12s,eosvv13s,eosvv14s)
#     labx <- 'VIIRS pheno-metrics'  
#     
#     plot(sosvv13,sospc13,xlim=c(50,350),ylim=c(50,350),pch=1,
#          cex.lab=1.5,cex.axis=1.3,cex=1.5,
#          xlab=labx,ylab='Phenocam pheno-metrics')
#     points(sosvv12,sospc12,pch=1,cex=1.5)
#     points(sosvv14,sospc14,pch=1,cex=1.5)
#     points(eosvv12,eospc12,pch=0,cex=1.5)
#     points(eosvv13,eospc13,pch=0,cex=1.5)
#     points(eosvv14,eospc14,pch=0,cex=1.5)
#     abline(0,1,lty=5)
#     
#   }else{
#     xx <- c(sosmm12,sosmm13,sosmm14,eosmm12,eosmm13,eosmm14)
#     yy <- c(sospc12,sospc13,sospc14,eospc12,eospc13,eospc14)
#     xsd <- c(sosmm12s,sosmm13s,sosmm14s,eosmm12s,eosmm13s,eosmm14s)
#     labx <- 'MODIS pheno-metrics'
#     
#     plot(sosmm13,sospc13,xlim=c(50,350),ylim=c(50,350),pch=1,
#          cex.lab=1.5,cex.axis=1.3,cex=1.5,
#          xlab=labx,ylab='Phenocam pheno-metrics')
#     points(sosmm12,sospc12,pch=1,cex=1.5)
#     points(sosmm14,sospc14,pch=1,cex=1.5)
#     points(eosmm12,eospc12,pch=0,cex=1.5)
#     points(eosmm13,eospc13,pch=0,cex=1.5)
#     points(eosmm14,eospc14,pch=0,cex=1.5)
#     abline(0,1,lty=5)
#   }
#   
# #   arrows(xx-xsd,yy,xx+xsd,yy,length=0.05,angle=0)
#   
#   abline(lm(yy~xx))
#   reg <- summary(lm(yy~xx))
#   temp <- na.omit(yy-xx)
#   rmse <- round(sqrt(mean((temp)^2)),1)
#   rseq <- round(reg$r.squared,2)
#   b1 <- round(reg$coefficients[2],1)
#   b1c1 <- round(b1+reg$coefficients[4],1)
#   b1c2 <- round(b1-reg$coefficients[4],1)
#   b0 <- round(reg$coefficients[1],1)
#   b0c1 <- round(b0+reg$coefficients[3],1)
#   b0c2 <- round(b0-reg$coefficients[3],1)
#   text(200,120,expression(paste(R^2,'=')),pos=4,cex=1)
#   text(230,120,rseq,pos=4,cex=1)
#   if(reg$coefficient[8]<0.01){
#     text(260,120,substitute(paste('(', italic(p),' < 0.01)',sep='')),pos=4,cex=1)
#   }else{
#     text(260,120,substitute(paste('(', italic(p),' = ',sep='')),pos=4,cex=1)
#     text(280,120,paste(round(reg$coefficient[8],2),' )',sep=''),pos=4,cex=1)
#   }
#   text(200,100,paste('RMSE =',rmse),pos=4,cex=1)
#   text(200,80,expression(paste(beta[1],'=')),pos=4,cex=1)
#   text(230,80,paste(b1,' (',b1c1,',',b1c2,')',sep=''),pos=4,cex=1)
#   text(200,60,expression(paste(beta[0],'=')),pos=4,cex=1)
#   text(230,60,paste(b0,' (',b0c1,',',b0c2,')',sep=''),pos=4,cex=1)  
# }
# 
# # setwd('/projectnb/modislc/users/mkmoon/VIIRS/figures/')
# # pdf(file='phenocam1.pdf',width=12,height=6)  
# 
# par(mfrow=c(1,2),oma=c(1,1,1,1),mar=c(4,4,4,4),mgp=c(2.5,1,0))
# plot_phenocam('vv')
# plot_phenocam('mm')
# 
# # dev.off()
# 
# par(mfrow=c(2,2),oma=c(1,1,1,1),mar=c(4,4,4,4),mgp=c(2.5,1,0))
# 
# vgt <- rep(phenosos$vegtype,3)
# year <- c(rep(2012,nrow(phenosos)),rep(2013,nrow(phenosos)),rep(2014,nrow(phenosos)))
# xx1 <- c(sosvv12,sosvv13,sosvv14)
# xx2 <- c(sosmm12,sosmm13,sosmm14)
# yy <- c(sospc12,sospc13,sospc14)
# dat <- data.frame(xx1,xx2,yy,vgt,year,lct)
# dat <- na.omit(dat)
# dat <- subset(dat,vgt=='DB')
# temp <- na.omit(dat$yy-dat$xx1)
# rmse <- round(sqrt(mean((temp)^2)),1)
# temp <- na.omit(dat$yy-dat$xx2)
# rmse <- round(sqrt(mean((temp)^2)),1)
# 
# plot(dat[,1],dat[,3],xlim=c(100,180),ylim=c(100,180),pch=19,
#      xlab='VIIRS SOS',ylab='Phenocam SOS')
# legend('topleft',c('DB','EN'),pch=c(19,1))
# abline(lm(dat[,3]~dat[,1]))
# summary(lm(dat[,3]~dat[,1]))
# dat <- data.frame(xx1,xx2,yy,vgt,year,lct)
# dat <- na.omit(dat)
# dat <- subset(dat,vgt=='EN')
# points(dat[,1],dat[,3])
# abline(0,1,lty=5)
# abline(lm(dat[,3]~dat[,1]))
# summary(lm(dat[,3]~dat[,1]))
# 
# dat <- data.frame(xx1,xx2,yy,vgt,year,lct)
# dat <- na.omit(dat)
# dat <- subset(dat,vgt=='DB')
# plot(dat[,2],dat[,3],xlim=c(100,180),ylim=c(100,180),pch=19,
#      xlab='MODIS SOS',ylab='Phenocam SOS')
# legend('topleft',c('DB','EN'),pch=c(19,1))
# abline(lm(dat[,3]~dat[,2]))
# summary(lm(dat[,3]~dat[,2]))
# dat <- data.frame(xx1,xx2,yy,vgt,year,lct)
# dat <- na.omit(dat)
# dat <- subset(dat,vgt=='EN')
# points(dat[,2],dat[,3])
# abline(0,1,lty=5)
# abline(lm(dat[,3]~dat[,2]))
# summary(lm(dat[,3]~dat[,2]))
# 
# xx1 <- c(eosvv12,eosvv13,eosvv14)
# xx2 <- c(eosmm12,eosmm13,eosmm14)
# yy <- c(eospc12,eospc13,eospc14)
# dat <- data.frame(xx1,xx2,yy,vgt,year,lct)
# dat <- na.omit(dat)
# dat <- subset(dat,vgt=='DB')
# temp <- na.omit(dat$yy-dat$xx1)
# rmse <- round(sqrt(mean((temp)^2)),1)
# temp <- na.omit(dat$yy-dat$xx2)
# rmse <- round(sqrt(mean((temp)^2)),1)
# 
# plot(dat[,1],dat[,3],xlim=c(210,350),ylim=c(210,350),pch=19,
#      xlab='VIIRS EOS',ylab='Phenocam EOS')
# legend('topleft',c('DB','EN'),pch=c(19,1))
# abline(lm(dat[,3]~dat[,1]))
# summary(lm(dat[,3]~dat[,1]))
# dat <- data.frame(xx1,xx2,yy,vgt,year,lct)
# dat <- na.omit(dat)
# dat <- subset(dat,vgt=='EN')
# points(dat[,1],dat[,3])
# abline(0,1,lty=5)
# abline(lm(dat[,3]~dat[,1]))
# summary(lm(dat[,3]~dat[,1]))
# 
# dat <- data.frame(xx1,xx2,yy,vgt,year,lct)
# dat <- na.omit(dat)
# dat <- subset(dat,vgt=='DB')
# plot(dat[,2],dat[,3],xlim=c(210,350),ylim=c(210,350),pch=19,
#      xlab='MODIS EOS',ylab='Phenocam EOS')
# legend('topleft',c('DB','EN'),pch=c(19,1))
# abline(lm(dat[,3]~dat[,2]))
# summary(lm(dat[,3]~dat[,2]))
# dat <- data.frame(xx1,xx2,yy,vgt,year,lct)
# dat <- na.omit(dat)
# dat <- subset(dat,vgt=='EN')
# points(dat[,2],dat[,3])
# abline(0,1,lty=5)
# abline(lm(dat[,3]~dat[,2]))
# summary(lm(dat[,3]~dat[,2]))