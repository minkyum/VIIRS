## First run 01 code

# Site location
setwd('/projectnb/modislc/users/mkmoon/NEphenology/phenocam/')
sites <- read.csv('site_modistile.csv')

site.1204 <- subset(sites,h_tile==12&v_tile==4)
site.1104 <- subset(sites,h_tile==11&v_tile==4)
site.0805 <- subset(sites,h_tile==8&v_tile==5)

path <- '/projectnb/modislc/users/mkmoon/NEphenology/phenocam/data/'
search_str <- paste('*transi*.csv',sep='')
files <- list.files(path=path,pattern=glob2rx(search_str),full.names=F,include.dirs=F)

## Get only for DB ROI files
site.name <- sub('_.*','',files)

vg.type <- rep(0)
for(i in 1:length(site.name)){
  vg.type[i] <- sub(site.name[i],'',files[i])  
  vg.type[i] <- substr(vg.type[i],2,3)
}


## Match site name and coordinate
sitecord <- matrix(NA,length(site.name),4)
# rownames(sitecord) <- site.name
cor <- as.matrix(sites[,3:6])
for(i in 1:length(site.name)){
  for(j in 1:nrow(cor)){
    if(sites$sites[j]==site.name[i]) sitecord[i,] <- cor[j,]  
  }
}

## Extract 50% SOS and EOS dates
datesos <- matrix(NA,length(site.name),20)
dateeos <- matrix(NA,length(site.name),20)
colnames(datesos) <- 2001:2020
colnames(dateeos) <- 2001:2020

for(ss in 1:length(site.name)){
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


### MODIS data near Phenocam sites
get_lsp_sos <- function(sensor,year){
  phe <- matrix(NA,length(site.name),25)
  
  for(i in 1:length(site.name)){
    if(year==2013){
      if(sensor=='mm'){
        if(phenosos[i,1]==12){
          rast <- (mm120413[[1]]+mm120413[[2]])/2
        }else if(phenosos[i,1]==11){
          rast <- (mm1104[[1]]+mm1104[[2]])/2
        }else{
          rast <- (mm0805[[1]]+mm0805[[2]])/2
        }
      }else{
        if(phenosos[i,1]==12){
          rast <- (vv120413[[1]][[1]]+vv120413[[2]][[1]])/2
        }else if(phenosos[i,1]==11){
          rast <- (vv1104[[1]][[1]]+vv1104[[2]][[1]])/2
        }else{
          rast <- (vv0805[[1]][[1]]+vv0805[[2]][[1]])/2
        }
      }  
    }else if(year==2012){
      if(sensor=='mm'){
        if(phenosos[i,1]==12){
          rast <- (mm120412[[1]]+mm120412[[2]])/2
        }else{
          rast <- NA
        }
      }else{
        if(phenosos[i,1]==12){
          rast <- (vv120412[[1]][[1]]+vv120412[[2]][[1]])/2
        }else{
          rast <- NA
        }
      }
    }else{
      if(sensor=='mm'){
        if(phenosos[i,1]==12){
          rast <- (mm120414[[1]]+mm120414[[2]])/2
        }else{
          rast <- NA
        }
      }else{
        if(phenosos[i,1]==12){
          rast <- (vv120414[[1]][[1]]+vv120414[[2]][[1]])/2
        }else{
          rast <- NA
        }
      }
    }
    
    if(year==2013){
      phe[i,1] <- rast[round(phenosos[i,3])-2,round(phenosos[i,4])-2]
      phe[i,2] <- rast[round(phenosos[i,3])-2,round(phenosos[i,4])-1]
      phe[i,3] <- rast[round(phenosos[i,3])-2,round(phenosos[i,4])]
      phe[i,4] <- rast[round(phenosos[i,3])-2,round(phenosos[i,4])+1]
      phe[i,5] <- rast[round(phenosos[i,3])-2,round(phenosos[i,4])+2]
      phe[i,6] <- rast[round(phenosos[i,3])-1,round(phenosos[i,4])-2]
      phe[i,7] <- rast[round(phenosos[i,3])-1,round(phenosos[i,4])-1]
      phe[i,8] <- rast[round(phenosos[i,3])-1,round(phenosos[i,4])]
      phe[i,9] <- rast[round(phenosos[i,3])-1,round(phenosos[i,4])+1]
      phe[i,10] <- rast[round(phenosos[i,3])-1,round(phenosos[i,4])+2]
      phe[i,11] <- rast[round(phenosos[i,3])+0,round(phenosos[i,4])-2]
      phe[i,12] <- rast[round(phenosos[i,3])+0,round(phenosos[i,4])-1]
      phe[i,13] <- rast[round(phenosos[i,3])+0,round(phenosos[i,4])]
      phe[i,14] <- rast[round(phenosos[i,3])+0,round(phenosos[i,4])+1]
      phe[i,15] <- rast[round(phenosos[i,3])+0,round(phenosos[i,4])+2]
      phe[i,16] <- rast[round(phenosos[i,3])+1,round(phenosos[i,4])-2]
      phe[i,17] <- rast[round(phenosos[i,3])+1,round(phenosos[i,4])-1]
      phe[i,18] <- rast[round(phenosos[i,3])+1,round(phenosos[i,4])]
      phe[i,19] <- rast[round(phenosos[i,3])+1,round(phenosos[i,4])+1]
      phe[i,20] <- rast[round(phenosos[i,3])+1,round(phenosos[i,4])+2]
      phe[i,21] <- rast[round(phenosos[i,3])+2,round(phenosos[i,4])-2]
      phe[i,22] <- rast[round(phenosos[i,3])+2,round(phenosos[i,4])-1]
      phe[i,23] <- rast[round(phenosos[i,3])+2,round(phenosos[i,4])]
      phe[i,24] <- rast[round(phenosos[i,3])+2,round(phenosos[i,4])+1]
      phe[i,25] <- rast[round(phenosos[i,3])+2,round(phenosos[i,4])+2]
    }else{
      if(phenosos[i,1]==11|phenosos[i,1]==8){
        phe[i,] <- NA
      }else{
        phe[i,1] <- rast[round(phenosos[i,3])-2,round(phenosos[i,4])-2]
        phe[i,2] <- rast[round(phenosos[i,3])-2,round(phenosos[i,4])-1]
        phe[i,3] <- rast[round(phenosos[i,3])-2,round(phenosos[i,4])]
        phe[i,4] <- rast[round(phenosos[i,3])-2,round(phenosos[i,4])+1]
        phe[i,5] <- rast[round(phenosos[i,3])-2,round(phenosos[i,4])+2]
        phe[i,6] <- rast[round(phenosos[i,3])-1,round(phenosos[i,4])-2]
        phe[i,7] <- rast[round(phenosos[i,3])-1,round(phenosos[i,4])-1]
        phe[i,8] <- rast[round(phenosos[i,3])-1,round(phenosos[i,4])]
        phe[i,9] <- rast[round(phenosos[i,3])-1,round(phenosos[i,4])+1]
        phe[i,10] <- rast[round(phenosos[i,3])-1,round(phenosos[i,4])+2]
        phe[i,11] <- rast[round(phenosos[i,3])+0,round(phenosos[i,4])-2]
        phe[i,12] <- rast[round(phenosos[i,3])+0,round(phenosos[i,4])-1]
        phe[i,13] <- rast[round(phenosos[i,3])+0,round(phenosos[i,4])]
        phe[i,14] <- rast[round(phenosos[i,3])+0,round(phenosos[i,4])+1]
        phe[i,15] <- rast[round(phenosos[i,3])+0,round(phenosos[i,4])+2]
        phe[i,16] <- rast[round(phenosos[i,3])+1,round(phenosos[i,4])-2]
        phe[i,17] <- rast[round(phenosos[i,3])+1,round(phenosos[i,4])-1]
        phe[i,18] <- rast[round(phenosos[i,3])+1,round(phenosos[i,4])]
        phe[i,19] <- rast[round(phenosos[i,3])+1,round(phenosos[i,4])+1]
        phe[i,20] <- rast[round(phenosos[i,3])+1,round(phenosos[i,4])+2]
        phe[i,21] <- rast[round(phenosos[i,3])+2,round(phenosos[i,4])-2]
        phe[i,22] <- rast[round(phenosos[i,3])+2,round(phenosos[i,4])-1]
        phe[i,23] <- rast[round(phenosos[i,3])+2,round(phenosos[i,4])]
        phe[i,24] <- rast[round(phenosos[i,3])+2,round(phenosos[i,4])+1]
        phe[i,25] <- rast[round(phenosos[i,3])+2,round(phenosos[i,4])+2]
      }
    }
      
    
  }
  
  return(phe)
}

phevv12 <- get_lsp_sos('vv',2012)
phemm12 <- get_lsp_sos('mm',2012)
phevv13 <- get_lsp_sos('vv',2013)
phemm13 <- get_lsp_sos('mm',2013)
phevv14 <- get_lsp_sos('vv',2014)
phemm14 <- get_lsp_sos('mm',2014)


sosvv12 <- apply(phevv12,1,mean)
sosmm12 <- apply(phemm12,1,mean)
sosvv13 <- apply(phevv13,1,mean)
sosmm13 <- apply(phemm13,1,mean)
sosvv14 <- apply(phevv14,1,mean)
sosmm14 <- apply(phemm14,1,mean)

sospc12 <- as.Date(phenosos[,16],origin='1970-1-1')
sospc12 <- as.numeric(strftime(sospc12,format='%j'))
sospc13 <- as.Date(phenosos[,17],origin='1970-1-1')
sospc13 <- as.numeric(strftime(sospc13,format='%j'))
sospc14 <- as.Date(phenosos[,18],origin='1970-1-1')
sospc14 <- as.numeric(strftime(sospc14,format='%j'))

par(mfrow=c(1,2))
plot(sosvv13,sosmm13,xlim=c(100,180),ylim=c(100,180),pch=19)
points(sosvv12,sosmm12,pch=18)
points(sosvv14,sosmm14,pch=17)
abline(0,1)

plot(sosvv13,sospc13,xlim=c(100,180),ylim=c(100,180),pch=19,col='red')
points(sosmm13,sospc13,pch=19,col='blue')
points(sosvv12,sospc12,pch=18,col='red')
points(sosmm12,sospc12,pch=18,col='blue')
points(sosvv14,sospc14,pch=17,col='red')
points(sosmm14,sospc14,pch=17,col='blue')
abline(0,1)


par(mfrow=c(1,2))
plot(sosvv13,sosmm13,xlim=c(100,180),ylim=c(100,180),pch=19)
abline(0,1)
plot(sosvv13,sospc13,xlim=c(100,180),ylim=c(100,180),pch=vg.type,col='red')
points(sosmm13,sospc13,pch=vg.type,col='blue')
abline(0,1)
pcsos <- phenosos[,c(1,2,3,4,17)]
pceos <- phenoeos[,c(1,2,3,4,17)]
pcsos <- cbind(vg.type,pcsos,sosvv13,sosmm13)
pceos <- cbind(vg.type,pceos,eosvv13,eosmm13)
pcsos <- na.omit(pcsos)
pceos <- na.omit(pceos)
table(pcsos[,1])
table(pceos[,1])
