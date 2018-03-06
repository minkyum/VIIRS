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
sitecord <- data.frame(site.name,cor)

# Extract 1204, 1104, and 0805 tiles
sitecord1 <- subset(sitecord,X3==12&X4==4)
sitecord2 <- subset(sitecord,X3==11&X4==4)
sitecord3 <- subset(sitecord,X3==8&X4==5)
sitecord <- rbind(sitecord1,sitecord2,sitecord3)

# Extract 3 day transition data
search_str <- '*3day_transi*.csv'
files <- list.files(path=path,pattern=glob2rx(search_str),full.names=F,include.dirs=F)
sites <- sub('_.*','',files)

## Match site name and coordinate
cor <- as.matrix(sitecord[,2:7])
cord <- matrix(NA,length(sites),6)
for(i in 1:length(sites)){
  for(j in 1:nrow(sitecord)){
    if(sitecord$site.name[j]==sites[i]) cord[i,] <- cor[j,]  
  }
}
sitecord <- data.frame(sites,cord)
sitecord <- na.omit(sitecord)
colnames(sitecord) <- c('site','lat','lon','htile','vtile','x','y')

## Extract 50% SOS and EOS dates
datesos <- matrix(NA,nrow(sitecord),16)
dateeos <- matrix(NA,nrow(sitecord),16)
colnames(datesos) <- 2001:2020
colnames(dateeos) <- 2001:2020

for(ss in 1:nrow(sitecord)){
  setwd('/projectnb/modislc/users/mkmoon/NEphenology/phenocam/data/')
  dat <- read.csv(files[ss],skip=16,header=T)
  dat1 <- subset(dat,direction=='rising'&gcc_value=='gcc_90',)
  dat2 <- subset(dat,direction=='falling'&gcc_value=='gcc_90',)
  
  sos <- dat1$transition_50
  for(i in 1:length(sos)){
    yy <- as.numeric(substr(sos[i],1,4))
    cc <- yy-2000
    datesos[ss,cc] <- as.Date(sos[i])
  }
  
  eos <- dat2$transition_50
  for(i in 1:length(eos)){
    yy <- as.numeric(substr(eos[i],1,4))
    cc <- yy-2000
    dateeos[ss,cc] <- as.Date(eos[i])
  }
}

## Merge coordinate and dates
phenosos <- cbind(sitecord,datesos)
phenoeos <- cbind(sitecord,dateeos)

phenosos <- subset(phenosos,vg.type=='AG'|vg.type=='DB'|vg.type=='DN'|vg.type=='EN'|vg.type=='GR'|vg.type=='SH')
phenoeos <- subset(phenoeos,vg.type=='AG'|vg.type=='DB'|vg.type=='DN'|vg.type=='EN'|vg.type=='GR'|vg.type=='SH')

### MODIS data near Phenocam sites
get_lsp_sos <- function(sensor,year){
  phe <- matrix(NA,nrow(phenosos),25)
  
  for(i in 1:nrow(phenosos)){
    if(year==2013){
      if(sensor=='mm'){
        if(phenosos[i,3]==12){
          rast <- (mm120413[[1]]+mm120413[[2]])/2
        }else if(phenosos[i,3]==11){
          rast <- (mm1104[[1]]+mm1104[[2]])/2
        }else{
          rast <- (mm0805[[1]]+mm0805[[2]])/2
        }
      }else{
        if(phenosos[i,3]==12){
          rast <- (vv120413[[1]][[1]]+vv120413[[2]][[1]])/2
        }else if(phenosos[i,3]==11){
          rast <- (vv1104[[1]][[1]]+vv1104[[2]][[1]])/2
        }else{
          rast <- (vv0805[[1]][[1]]+vv0805[[2]][[1]])/2
        }
      }  
    }else if(year==2012){
      if(sensor=='mm'){
        if(phenosos[i,3]==12){
          rast <- (mm120412[[1]]+mm120412[[2]])/2
        }else{
          rast <- NA
        }
      }else{
        if(phenosos[i,3]==12){
          rast <- (vv120412[[1]][[1]]+vv120412[[2]][[1]])/2
        }else{
          rast <- NA
        }
      }
    }else{
      if(sensor=='mm'){
        if(phenosos[i,3]==12){
          rast <- (mm120414[[1]]+mm120414[[2]])/2
        }else{
          rast <- NA
        }
      }else{
        if(phenosos[i,3]==12){
          rast <- (vv120414[[1]][[1]]+vv120414[[2]][[1]])/2
        }else{
          rast <- NA
        }
      }
    }
    
    if(year==2013){
      phe[i,1] <- rast[round(phenosos[i,5])-2,round(phenosos[i,6])-2]
      phe[i,2] <- rast[round(phenosos[i,5])-2,round(phenosos[i,6])-1]
      phe[i,3] <- rast[round(phenosos[i,5])-2,round(phenosos[i,6])]
      phe[i,4] <- rast[round(phenosos[i,5])-2,round(phenosos[i,6])+1]
      phe[i,5] <- rast[round(phenosos[i,5])-2,round(phenosos[i,6])+2]
      phe[i,6] <- rast[round(phenosos[i,5])-1,round(phenosos[i,6])-2]
      phe[i,7] <- rast[round(phenosos[i,5])-1,round(phenosos[i,6])-1]
      phe[i,8] <- rast[round(phenosos[i,5])-1,round(phenosos[i,6])]
      phe[i,9] <- rast[round(phenosos[i,5])-1,round(phenosos[i,6])+1]
      phe[i,10] <- rast[round(phenosos[i,5])-1,round(phenosos[i,6])+2]
      phe[i,11] <- rast[round(phenosos[i,5])+0,round(phenosos[i,6])-2]
      phe[i,12] <- rast[round(phenosos[i,5])+0,round(phenosos[i,6])-1]
      phe[i,13] <- rast[round(phenosos[i,5])+0,round(phenosos[i,6])]
      phe[i,14] <- rast[round(phenosos[i,5])+0,round(phenosos[i,6])+1]
      phe[i,15] <- rast[round(phenosos[i,5])+0,round(phenosos[i,6])+2]
      phe[i,16] <- rast[round(phenosos[i,5])+1,round(phenosos[i,6])-2]
      phe[i,17] <- rast[round(phenosos[i,5])+1,round(phenosos[i,6])-1]
      phe[i,18] <- rast[round(phenosos[i,5])+1,round(phenosos[i,6])]
      phe[i,19] <- rast[round(phenosos[i,5])+1,round(phenosos[i,6])+1]
      phe[i,20] <- rast[round(phenosos[i,5])+1,round(phenosos[i,6])+2]
      phe[i,21] <- rast[round(phenosos[i,5])+2,round(phenosos[i,6])-2]
      phe[i,22] <- rast[round(phenosos[i,5])+2,round(phenosos[i,6])-1]
      phe[i,23] <- rast[round(phenosos[i,5])+2,round(phenosos[i,6])]
      phe[i,24] <- rast[round(phenosos[i,5])+2,round(phenosos[i,6])+1]
      phe[i,25] <- rast[round(phenosos[i,5])+2,round(phenosos[i,6])+2]
    }else{
      if(phenosos[i,3]==11|phenosos[i,3]==8){
        phe[i,] <- NA
      }else{
        phe[i,1] <- rast[round(phenosos[i,5])-2,round(phenosos[i,6])-2]
        phe[i,2] <- rast[round(phenosos[i,5])-2,round(phenosos[i,6])-1]
        phe[i,3] <- rast[round(phenosos[i,5])-2,round(phenosos[i,6])]
        phe[i,4] <- rast[round(phenosos[i,5])-2,round(phenosos[i,6])+1]
        phe[i,5] <- rast[round(phenosos[i,5])-2,round(phenosos[i,6])+2]
        phe[i,6] <- rast[round(phenosos[i,5])-1,round(phenosos[i,6])-2]
        phe[i,7] <- rast[round(phenosos[i,5])-1,round(phenosos[i,6])-1]
        phe[i,8] <- rast[round(phenosos[i,5])-1,round(phenosos[i,6])]
        phe[i,9] <- rast[round(phenosos[i,5])-1,round(phenosos[i,6])+1]
        phe[i,10] <- rast[round(phenosos[i,5])-1,round(phenosos[i,6])+2]
        phe[i,11] <- rast[round(phenosos[i,5])+0,round(phenosos[i,6])-2]
        phe[i,12] <- rast[round(phenosos[i,5])+0,round(phenosos[i,6])-1]
        phe[i,13] <- rast[round(phenosos[i,5])+0,round(phenosos[i,6])]
        phe[i,14] <- rast[round(phenosos[i,5])+0,round(phenosos[i,6])+1]
        phe[i,15] <- rast[round(phenosos[i,5])+0,round(phenosos[i,6])+2]
        phe[i,16] <- rast[round(phenosos[i,5])+1,round(phenosos[i,6])-2]
        phe[i,17] <- rast[round(phenosos[i,5])+1,round(phenosos[i,6])-1]
        phe[i,18] <- rast[round(phenosos[i,5])+1,round(phenosos[i,6])]
        phe[i,19] <- rast[round(phenosos[i,5])+1,round(phenosos[i,6])+1]
        phe[i,20] <- rast[round(phenosos[i,5])+1,round(phenosos[i,6])+2]
        phe[i,21] <- rast[round(phenosos[i,5])+2,round(phenosos[i,6])-2]
        phe[i,22] <- rast[round(phenosos[i,5])+2,round(phenosos[i,6])-1]
        phe[i,23] <- rast[round(phenosos[i,5])+2,round(phenosos[i,6])]
        phe[i,24] <- rast[round(phenosos[i,5])+2,round(phenosos[i,6])+1]
        phe[i,25] <- rast[round(phenosos[i,5])+2,round(phenosos[i,6])+2]
      }
    }
  }
  return(phe)
}
get_lsp_eos <- function(sensor,year){
  phe <- matrix(NA,nrow(phenosos),25)
  
  for(i in 1:nrow(phenosos)){
    if(year==2013){
      if(sensor=='mm'){
        if(phenosos[i,3]==12){
          rast <- (mm120413[[3]]+mm120413[[4]])/2
        }else if(phenosos[i,3]==11){
          rast <- (mm1104[[3]]+mm1104[[4]])/2
        }else{
          rast <- (mm0805[[3]]+mm0805[[4]])/2
        }
      }else{
        if(phenosos[i,3]==12){
          rast <- (vv120413[[3]][[1]]+vv120413[[4]][[1]])/2
        }else if(phenosos[i,3]==11){
          rast <- (vv1104[[3]][[1]]+vv1104[[4]][[1]])/2
        }else{
          rast <- (vv0805[[3]][[1]]+vv0805[[4]][[1]])/2
        }
      }  
    }else if(year==2012){
      if(sensor=='mm'){
        if(phenosos[i,3]==12){
          rast <- (mm120412[[3]]+mm120412[[4]])/2
        }else{
          rast <- NA
        }
      }else{
        if(phenosos[i,3]==12){
          rast <- (vv120412[[3]][[1]]+vv120412[[4]][[1]])/2
        }else{
          rast <- NA
        }
      }
    }else{
      if(sensor=='mm'){
        if(phenosos[i,3]==12){
          rast <- (mm120414[[3]]+mm120414[[4]])/2
        }else{
          rast <- NA
        }
      }else{
        if(phenosos[i,3]==12){
          rast <- (vv120414[[3]][[1]]+vv120414[[4]][[1]])/2
        }else{
          rast <- NA
        }
      }
    }
    
    if(year==2013){
      phe[i,1] <- rast[round(phenosos[i,5])-2,round(phenosos[i,6])-2]
      phe[i,2] <- rast[round(phenosos[i,5])-2,round(phenosos[i,6])-1]
      phe[i,3] <- rast[round(phenosos[i,5])-2,round(phenosos[i,6])]
      phe[i,4] <- rast[round(phenosos[i,5])-2,round(phenosos[i,6])+1]
      phe[i,5] <- rast[round(phenosos[i,5])-2,round(phenosos[i,6])+2]
      phe[i,6] <- rast[round(phenosos[i,5])-1,round(phenosos[i,6])-2]
      phe[i,7] <- rast[round(phenosos[i,5])-1,round(phenosos[i,6])-1]
      phe[i,8] <- rast[round(phenosos[i,5])-1,round(phenosos[i,6])]
      phe[i,9] <- rast[round(phenosos[i,5])-1,round(phenosos[i,6])+1]
      phe[i,10] <- rast[round(phenosos[i,5])-1,round(phenosos[i,6])+2]
      phe[i,11] <- rast[round(phenosos[i,5])+0,round(phenosos[i,6])-2]
      phe[i,12] <- rast[round(phenosos[i,5])+0,round(phenosos[i,6])-1]
      phe[i,13] <- rast[round(phenosos[i,5])+0,round(phenosos[i,6])]
      phe[i,14] <- rast[round(phenosos[i,5])+0,round(phenosos[i,6])+1]
      phe[i,15] <- rast[round(phenosos[i,5])+0,round(phenosos[i,6])+2]
      phe[i,16] <- rast[round(phenosos[i,5])+1,round(phenosos[i,6])-2]
      phe[i,17] <- rast[round(phenosos[i,5])+1,round(phenosos[i,6])-1]
      phe[i,18] <- rast[round(phenosos[i,5])+1,round(phenosos[i,6])]
      phe[i,19] <- rast[round(phenosos[i,5])+1,round(phenosos[i,6])+1]
      phe[i,20] <- rast[round(phenosos[i,5])+1,round(phenosos[i,6])+2]
      phe[i,21] <- rast[round(phenosos[i,5])+2,round(phenosos[i,6])-2]
      phe[i,22] <- rast[round(phenosos[i,5])+2,round(phenosos[i,6])-1]
      phe[i,23] <- rast[round(phenosos[i,5])+2,round(phenosos[i,6])]
      phe[i,24] <- rast[round(phenosos[i,5])+2,round(phenosos[i,6])+1]
      phe[i,25] <- rast[round(phenosos[i,5])+2,round(phenosos[i,6])+2]
    }else{
      if(phenosos[i,3]==11|phenosos[i,3]==8){
        phe[i,] <- NA
      }else{
        phe[i,1] <- rast[round(phenosos[i,5])-2,round(phenosos[i,6])-2]
        phe[i,2] <- rast[round(phenosos[i,5])-2,round(phenosos[i,6])-1]
        phe[i,3] <- rast[round(phenosos[i,5])-2,round(phenosos[i,6])]
        phe[i,4] <- rast[round(phenosos[i,5])-2,round(phenosos[i,6])+1]
        phe[i,5] <- rast[round(phenosos[i,5])-2,round(phenosos[i,6])+2]
        phe[i,6] <- rast[round(phenosos[i,5])-1,round(phenosos[i,6])-2]
        phe[i,7] <- rast[round(phenosos[i,5])-1,round(phenosos[i,6])-1]
        phe[i,8] <- rast[round(phenosos[i,5])-1,round(phenosos[i,6])]
        phe[i,9] <- rast[round(phenosos[i,5])-1,round(phenosos[i,6])+1]
        phe[i,10] <- rast[round(phenosos[i,5])-1,round(phenosos[i,6])+2]
        phe[i,11] <- rast[round(phenosos[i,5])+0,round(phenosos[i,6])-2]
        phe[i,12] <- rast[round(phenosos[i,5])+0,round(phenosos[i,6])-1]
        phe[i,13] <- rast[round(phenosos[i,5])+0,round(phenosos[i,6])]
        phe[i,14] <- rast[round(phenosos[i,5])+0,round(phenosos[i,6])+1]
        phe[i,15] <- rast[round(phenosos[i,5])+0,round(phenosos[i,6])+2]
        phe[i,16] <- rast[round(phenosos[i,5])+1,round(phenosos[i,6])-2]
        phe[i,17] <- rast[round(phenosos[i,5])+1,round(phenosos[i,6])-1]
        phe[i,18] <- rast[round(phenosos[i,5])+1,round(phenosos[i,6])]
        phe[i,19] <- rast[round(phenosos[i,5])+1,round(phenosos[i,6])+1]
        phe[i,20] <- rast[round(phenosos[i,5])+1,round(phenosos[i,6])+2]
        phe[i,21] <- rast[round(phenosos[i,5])+2,round(phenosos[i,6])-2]
        phe[i,22] <- rast[round(phenosos[i,5])+2,round(phenosos[i,6])-1]
        phe[i,23] <- rast[round(phenosos[i,5])+2,round(phenosos[i,6])]
        phe[i,24] <- rast[round(phenosos[i,5])+2,round(phenosos[i,6])+1]
        phe[i,25] <- rast[round(phenosos[i,5])+2,round(phenosos[i,6])+2]
      }
    }
  }
  return(phe)
}

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

# Get mean values for 5 by 5 window
## SOS
sosvv12 <- apply(phevv12s,1,mean)
sosmm12 <- apply(phemm12s,1,mean)
sosvv13 <- apply(phevv13s,1,mean)
sosmm13 <- apply(phemm13s,1,mean)
sosvv14 <- apply(phevv14s,1,mean)
sosmm14 <- apply(phemm14s,1,mean)

sosvv12s <- apply(phevv12s,1,sd)
sosmm12s <- apply(phemm12s,1,sd)
sosvv13s <- apply(phevv13s,1,sd)
sosmm13s <- apply(phemm13s,1,sd)
sosvv14s <- apply(phevv14s,1,sd)
sosmm14s <- apply(phemm14s,1,sd)
## EOS
eosvv12 <- apply(phevv12e,1,mean)
eosmm12 <- apply(phemm12e,1,mean)
eosvv13 <- apply(phevv13e,1,mean)
eosmm13 <- apply(phemm13e,1,mean)
eosvv14 <- apply(phevv14e,1,mean)
eosmm14 <- apply(phemm14e,1,mean)

eosvv12s <- apply(phevv12e,1,sd)
eosmm12s <- apply(phemm12e,1,sd)
eosvv13s <- apply(phevv13e,1,sd)
eosmm13s <- apply(phemm13e,1,sd)
eosvv14s <- apply(phevv14e,1,sd)
eosmm14s <- apply(phemm14e,1,sd)


# Get Phenocam transtision dates
## SOS
sospc12 <- as.Date(phenosos[,18],origin='1970-1-1')
sospc12 <- as.numeric(strftime(sospc12,format='%j'))
sospc13 <- as.Date(phenosos[,19],origin='1970-1-1')
sospc13 <- as.numeric(strftime(sospc13,format='%j'))
sospc14 <- as.Date(phenosos[,20],origin='1970-1-1')
sospc14 <- as.numeric(strftime(sospc14,format='%j'))
## EOS
eospc12 <- as.Date(phenoeos[,18],origin='1970-1-1')
eospc12 <- as.numeric(strftime(eospc12,format='%j'))
eospc13 <- as.Date(phenoeos[,19],origin='1970-1-1')
eospc13 <- as.numeric(strftime(eospc13,format='%j'))
eospc14 <- as.Date(phenoeos[,20],origin='1970-1-1')
eospc14 <- as.numeric(strftime(eospc14,format='%j'))


### Plot
plot_phenocam <- function(sensor){
  
  if(sensor=='vv'){
    xx <- c(sosvv12,sosvv13,sosvv14,eosvv12,eosvv13,eosvv14)
    yy <- c(sospc12,sospc13,sospc14,eospc12,eospc13,eospc14)
    xsd <- c(sosvv12s,sosvv13s,sosvv14s,eosvv12s,eosvv13s,eosvv14s)
    labx <- 'VIIRS pheno-metrics'  
    
    plot(sosvv13,sospc13,xlim=c(50,350),ylim=c(50,350),pch=1,
         cex.lab=1.5,cex.axis=1.3,cex=1.5,
         xlab=labx,ylab='Phenocam pheno-metrics')
    points(sosvv12,sospc12,pch=1,cex=1.5)
    points(sosvv14,sospc14,pch=1,cex=1.5)
    points(eosvv12,eospc12,pch=0,cex=1.5)
    points(eosvv13,eospc13,pch=0,cex=1.5)
    points(eosvv14,eospc14,pch=0,cex=1.5)
    abline(0,1,lty=5)
    
  }else{
    xx <- c(sosmm12,sosmm13,sosmm14,eosmm12,eosmm13,eosmm14)
    yy <- c(sospc12,sospc13,sospc14,eospc12,eospc13,eospc14)
    xsd <- c(sosmm12s,sosmm13s,sosmm14s,eosmm12s,eosmm13s,eosmm14s)
    labx <- 'MODIS pheno-metrics'
    
    plot(sosmm13,sospc13,xlim=c(50,350),ylim=c(50,350),pch=1,
         cex.lab=1.5,cex.axis=1.3,cex=1.5,
         xlab=labx,ylab='Phenocam pheno-metrics')
    points(sosmm12,sospc12,pch=1,cex=1.5)
    points(sosmm14,sospc14,pch=1,cex=1.5)
    points(eosmm12,eospc12,pch=0,cex=1.5)
    points(eosmm13,eospc13,pch=0,cex=1.5)
    points(eosmm14,eospc14,pch=0,cex=1.5)
    abline(0,1,lty=5)
  }
  
#   arrows(xx-xsd,yy,xx+xsd,yy,length=0.05,angle=0)
  
  abline(lm(yy~xx))
  reg <- summary(lm(yy~xx))
  temp <- na.omit(yy-xx)
  rmse <- round(sqrt(mean((temp)^2)),1)
  rseq <- round(reg$r.squared,2)
  b1 <- round(reg$coefficients[2],1)
  b1c1 <- round(b1+reg$coefficients[4],1)
  b1c2 <- round(b1-reg$coefficients[4],1)
  b0 <- round(reg$coefficients[1],1)
  b0c1 <- round(b0+reg$coefficients[3],1)
  b0c2 <- round(b0-reg$coefficients[3],1)
  text(200,120,expression(paste(R^2,'=')),pos=4,cex=1)
  text(230,120,rseq,pos=4,cex=1)
  if(reg$coefficient[8]<0.01){
    text(260,120,substitute(paste('(', italic(p),' < 0.01)',sep='')),pos=4,cex=1)
  }else{
    text(260,120,substitute(paste('(', italic(p),' = ',sep='')),pos=4,cex=1)
    text(280,120,paste(round(reg$coefficient[8],2),' )',sep=''),pos=4,cex=1)
  }
  text(200,100,paste('RMSE =',rmse),pos=4,cex=1)
  text(200,80,expression(paste(beta[1],'=')),pos=4,cex=1)
  text(230,80,paste(b1,' (',b1c1,',',b1c2,')',sep=''),pos=4,cex=1)
  text(200,60,expression(paste(beta[0],'=')),pos=4,cex=1)
  text(230,60,paste(b0,' (',b0c1,',',b0c2,')',sep=''),pos=4,cex=1)  
}

# setwd('/projectnb/modislc/users/mkmoon/VIIRS/figures/')
# pdf(file='phenocam1.pdf',width=12,height=6)  

par(mfrow=c(1,2),oma=c(1,1,1,1),mar=c(4,4,4,4),mgp=c(2.5,1,0))
plot_phenocam('vv')
plot_phenocam('mm')

# dev.off()


xx <- c(sosvv12,sosvv13,sosvv14,eosvv12,eosvv13,eosvv14)
yy <- c(sospc12,sospc13,sospc14,eospc12,eospc13,eospc14)

xx <- c(sosmm12,sosmm13,sosmm14,eosmm12,eosmm13,eosmm14)
yy <- c(sospc12,sospc13,sospc14,eospc12,eospc13,eospc14)

xx <- c(sosvv12,sosvv13,sosvv14,eosvv12,eosvv13,eosvv14)
yy <- c(sospc12,sospc13,sospc14,eospc12,eospc13,eospc14)

xx <- c(sosmm12,sosmm13,sosmm14,eosmm12,eosmm13,eosmm14)
yy <- c(sospc12,sospc13,sospc14,eospc12,eospc13,eospc14)
