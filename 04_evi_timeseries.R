



evitime <- function(sensor,year,tile,pixels){
  
  ifelse(sensor==1,sensor <- 'MODIS',sensor <- 'VIIRS')
  setwd(paste('/projectnb/modislc/users/mkmoon/VIIRS/R_data/evi2_ts/',sensor,'/',tile,sep=''))
  
  pixels <- pixels - 1
  rr <- pixels%/%240000 + 1
  rr <- sprintf('%03d',rr)
  
  ifelse(sensor=='MODIS',sensor <- 'mm',sensor <- 'vv')
  aa <- read.csv(paste('evi_',sensor,'_',tile,'_',year,'_',rr,'.csv',sep=''),header=T)
  
  
  return(evi)
}

evitime('m',2013,'h11v04',4:6)

setwd('/projectnb/modislc/users/mkmoon/VIIRS/figures/')
pdf(file=paste('evi2time_test.pdf',sep=''),width=15,height=7)

plot(evitime)

def.off()