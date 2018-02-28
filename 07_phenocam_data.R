# Site location
setwd('/projectnb/modislc/users/mkmoon/NEphenology/phenocam/')
sites <- read.csv('site_modistile.csv')

site.1204 <- subset(sites,h_tile==12&v_tile==4)
site.1104 <- subset(sites,h_tile==11&v_tile==4)
site.0805 <- subset(sites,h_tile==8&v_tile==5)

path <- '/projectnb/modislc/users/mkmoon/NEphenology/phenocam/data/'
search_str <- paste('*DB*transi*.csv',sep='')
files <- list.files(path=path,pattern=glob2rx(search_str),full.names=F,include.dirs=F)

site.name <- sub('_DB_.*','',files)

sitecord <- matrix(NA,length(site.name),4)
rownames(sitecord) <- site.name
cor <- as.matrix(sites[,3:6])
for(i in 1:length(site.name)){
  for(j in 1:nrow(cor)){
    if(sites$sites[j]==site.name[i]) sitecord[i,] <- cor[j,]  
  }
}

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
    datesos[ss,cc] <- (sos[i])
  }
  
  eos <- dat2$transition_50
  for(i in 1:length(eos)){
    yy <- as.numeric(substr(eos[i],1,4))
    cc <- yy-2000
    dateeos[ss,cc] <- (eos[i])
  }
}

