## First run 01 code

# Site location
setwd('/projectnb/modislc/users/mkmoon/NEphenology/phenocam/')
sites <- read.csv('site_modistile.csv')

site.1204 <- subset(sites,h_tile==12&v_tile==4)
site.1104 <- subset(sites,h_tile==11&v_tile==4)
site.0805 <- subset(sites,h_tile==8&v_tile==5)

path <- '/projectnb/modislc/users/mkmoon/NEphenology/phenocam/data/'
search_str <- paste('*DB*transi*.csv',sep='')
files <- list.files(path=path,pattern=glob2rx(search_str),full.names=F,include.dirs=F)

## Get only for DB ROI files
site.name <- sub('_DB_.*','',files)

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
get_lsp_sos <- function(sensor){
  phe <- matrix(NA,length(site.name),25)
  
  for(i in 1:length(site.name)){
    
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
  
  return(phe)
}

phevv <- get_lsp_sos('vv')
phemm <- get_lsp_sos('mm')


sosvv <- apply(phevv,1,mean)
sosmm <- apply(phemm,1,mean)
sospc <- as.Date(phenosos[,17],origin='1970-1-1')
sospc <- as.numeric(strftime(sospc,format='%j'))

plot(sosvv,sospc,xlim=c(90,170),ylim=c(90,170),pch=19,col='red')
points(sosmm,sospc,pch=19,col='blue')
abline(0,1)


