# Get NBAR (MCD43A4 C6) data from USGS

setwd('/projectnb/modislc/users/mkmoon/VIIRS/tiles/')
tile <- c('h12v04','h11v04','h08v05')

# MCD43A4
for(i in 1:1096){
  d = mod_dates[i]
  d_str = sprintf("%03d",d)
  dates <- as.Date(d,format='%Y.%m.%d',origin='2011.12.31')
  year <- substr(dates,1,4)
  month <- substr(dates,6,7)
  day <- substr(dates,9,10)
  
  url <- paste('https://e4ftl01.cr.usgs.gov/MOTA/MCD43A4.006/',year,'.',month,'.',day,'/',sep='')
  
#   for(j in 1:length(tile)){
#     system(paste('wget --user=mkmoon --password=M159k258! -l1 -r --no-parent -A "MCD43A4.A',year,d_str,'.',tile[j],'.006.*.hdf" ',url,sep=''))    
#   }
  system(paste('wget --user=mkmoon --password=M159k258! -l1 -r --no-parent -A "MCD43A4.A',year,d_str,'.',tile[j],'.006.*.hdf" ',url,sep=''))    
}