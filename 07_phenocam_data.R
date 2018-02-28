# if(!require(devtools)){install.package(devtools)}
# devtools::install_github("khufkens/daymetr")
# devtools::install_github("khufkens/phenocamr")
# devtools::install_github("khufkens/phenor")

# Site location
setwd('/projectnb/modislc/users/mkmoon/NEphenology/phenocam/')
sites <- read.csv('site_modistile.csv')

sam <- as.matrix(sites[,c('h_tile','v_tile','line','samp')])

sam.1204 <- sam[which(sam[,1]==12 & sam[,2]==4),]
sam.1104 <- sam[which(sam[,1]==11 & sam[,2]==4),]
sam.0805 <- sam[which(sam[,1]==8 & sam[,2]==5),]

site.1204 <- as.character(sites[which(sam[,1]==12 & sam[,2]==4),2])
site.1104 <- as.character(sites[which(sam[,1]==11 & sam[,2]==4),2])
site.0805 <- as.character(sites[which(sam[,1]==08 & sam[,2]==5),2])


## Download data
setwd('/projectnb/modislc/users/mkmoon/NEphenology/phenocam/data')
for(i in 1:length(site.1204)){
  download_phenocam(site = site.1204[i],
                    frequency = 3,
                    phenophase = TRUE)  
}



