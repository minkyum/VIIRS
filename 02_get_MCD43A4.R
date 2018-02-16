# Get NBAR (MCD43A4 C6) data from USGS
args <- commandArgs()
print(args)

tile <- substr(args[3],1,6)
year <- as.numeric(substr(args[3],7,10))

# MCD43A4
setwd('/projectnb/modislc/users/mkmoon/VIIRS/tiles/')
mod_dates <- seq(1,365)
for(i in 1:365){
  d = mod_dates[i]
  d_str = sprintf("%03d",d)
  dates <- as.Date(d,format='%Y.%m.%d',origin=paste((year-1),'.12.31',sep=''))
  mm <- substr(dates,6,7)
  dd <- substr(dates,9,10)
  
  url <- paste('https://e4ftl01.cr.usgs.gov/MOTA/MCD43A4.006/',year,'.',mm,'.',dd,'/',sep='')
  
  system(paste('wget --user=mkmoon --password=M159k258! -l1 -r --no-parent -A "MCD43A4.A',year,d_str,'.',tile,'.006.*.hdf" ',url,sep=''))    
}