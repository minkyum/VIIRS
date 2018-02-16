library(RColorBrewer)

### Load data
setwd('/projectnb/modislc/users/mkmoon/VIIRS/R_data/')
scene <- c('mwp_harvard','ah_bartlett','mwp_cary','ah_hubbard','mwp_proctor','mws_boundary_waters')
lphe <- vector('list',12)
mphe <- vector('list',12)
vphe <- vector('list',12)
qphe <- vector('list',12)

### SOS
for(i in 1:12){
  if(i<=6){
    load(paste('landsat_',scene[i],'_41.RData',sep=''))  
  }else{
    load(paste('landsat_',scene[i-6],'_40.RData',sep=''))  
  }
  lphe[[i]] <- lanbymod
  mphe[[i]] <- mmphe
  vphe[[i]] <- vvphe
  qphe[[i]] <- quant  
}

# Take a mean for each panel (i.e. 3 by 3 MODIS pixel)
pixmean <- vector('list',12)
for(i in 1:12){
  pixmean[[i]] <- apply(lphe[[i]],1,mean,na.rm=T)  
}

mesd <- matrix(NA,11,6)
for(i in 1:11){
  mesd[i,1] <- mean(mphe[[i]],na.rm=T)
  mesd[i,2] <- mean(vphe[[i]],na.rm=T)
  mesd[i,3] <- mean(lphe[[i]],na.rm=T)
  mesd[i,4] <- sd(mphe[[i]],na.rm=T)
  mesd[i,5] <- sd(vphe[[i]],na.rm=T)
  mesd[i,6] <- sd(lphe[[i]],na.rm=T)
}

## 1:1 plot
# all points
for(i in 1:11){
  if(i==1){
    plot(mphe[[i]],pixmean[[i]],pch=i,
         xlim=c(100,180),ylim=c(100,180),
         main=scene[i])  
    points(vphe[[i]],pixmean[[i]],col='red',pch=i)    
    abline(0,1,lty=5)
#     abline(lm(pixmean[[i]]~mphe[[i]]))
  }else{
    plot(mphe[[i]],pixmean[[i]],pch=i,
         xlim=c(100,180),ylim=c(100,180),
         main=scene[i])    
    points(vphe[[i]],pixmean[[i]],col='red',pch=i)    
    abline(0,1)
#     abline(lm(pixmean[[i]]~mphe[[i]]))
  } 
}

# Mean with error bars
plot_1to1 <- function(xmean,ymean,xsd,ysd,xlabel,ylabel,start){
  mycol <- brewer.pal(6,'Set1')
  plot(xmean,ymean,
       xlim=c(start,(start+60)),ylim=c(start,(start+60)),
       xlab=xlabel,
       ylab=ylabel,
       cex.lab=1.1,cex.axis=1.0)
  arrows(xmean,ymean-ysd,xmean,ymean+ysd,
         length=0.05,angle=0)
  arrows(xmean-xsd,ymean,xmean+xsd,ymean,
         length=0.05,angle=0)
  abline(0,1)
  abline(lm(ymean~xmean),lwd=1.3,lty=5)
  points(xmean,ymean,bg=mycol,pch=c(rep(21,6),rep(23,5)),cex=1.7)
  legend('topleft',c('Harvard','Bartlett','Cary','Hubbard','Proctor','Boundary waters'),
         bty='n',pch=21,pt.bg=mycol,pt.cex=1.5)
  reg <- summary(lm(ymean~xmean))
  rseq <- round(reg$r.squared,2)
  b1 <- round(reg$coefficients[2],1)
  b1c1 <- b1+round(reg$coefficients[4],1)
  b1c2 <- b1-round(reg$coefficients[4],1)
  b0 <- round(reg$coefficients[1],1)
  b0c1 <- b0+round(reg$coefficients[3],1)
  b0c2 <- b0-round(reg$coefficients[3],1)
  rmse <- round(sqrt(mean((ymean-xmean)^2)),1)
  text((start+33),(start+13),expression(paste(R^2,'=')),pos=4,cex=0.8)
  text((start+38),(start+13),rseq,pos=4,cex=0.8)
  text((start+45),(start+13),substitute(paste('(', italic(p),' < 0.01)',sep='')),pos=4,cex=0.8)
  text((start+33),(start+9),paste('RMSE =',rmse),pos=4,cex=0.8)
  text((start+33),(start+5),expression(paste(beta[1],'=')),pos=4,cex=0.8)
  text((start+38),(start+5),paste(b1,' (',b1c1,',',b1c2,')',sep=''),pos=4,cex=0.8)
  text((start+33),(start+1),expression(paste(beta[0],'=')),pos=4,cex=0.8)
  text((start+38),(start+1),paste(b0,' (',b0c1,',',b0c2,')',sep=''),pos=4,cex=0.8)
}

par(mfrow=c(1,3),oma=c(1,1,1,1),mar=c(4,4,4,4),mgp=c(2.5,1,0))
plot_1to1(mesd[,1],mesd[,2],mesd[,4],mesd[,5],'MODIS SOS','VIIRS SOS',110)
plot_1to1(mesd[,1],mesd[,3],mesd[,4],mesd[,6],'MODIS SOS','Landsat SOS',110)
plot_1to1(mesd[,2],mesd[,3],mesd[,5],mesd[,6],'VIIRS SOS','Landsat SOS',110)


### EOS
for(i in 1:12){
  if(i<=6){
    load(paste('landsat_',scene[i],'_73.RData',sep=''))  
  }else{
    load(paste('landsat_',scene[i-6],'_72.RData',sep=''))  
  }
  lphe[[i]] <- lanbymod
  mphe[[i]] <- mmphe
  vphe[[i]] <- vvphe
  qphe[[i]] <- quant  
}

# Take a mean for each panel (i.e. 3 by 3 MODIS pixel)
pixmean <- vector('list',12)
for(i in 1:12){
  pixmean[[i]] <- apply(lphe[[i]],1,mean,na.rm=T)  
}

mesd <- matrix(NA,11,6)
for(i in 1:11){
  mesd[i,1] <- mean(mphe[[i]],na.rm=T)
  mesd[i,2] <- mean(vphe[[i]],na.rm=T)
  mesd[i,3] <- mean(lphe[[i]],na.rm=T)
  mesd[i,4] <- sd(mphe[[i]],na.rm=T)
  mesd[i,5] <- sd(vphe[[i]],na.rm=T)
  mesd[i,6] <- sd(lphe[[i]],na.rm=T)
}

## 1:1 plot
# all points
for(i in 1:11){
  if(i==1){
    plot(mphe[[i]],pixmean[[i]],pch=i,
         xlim=c(220,320),ylim=c(220,320),
         main=scene[i])  
    points(vphe[[i]],pixmean[[i]],col='red',pch=i)    
    abline(0,1,lty=5)
    #     abline(lm(pixmean[[i]]~mphe[[i]]))
  }else{
    plot(mphe[[i]],pixmean[[i]],pch=i,
         xlim=c(220,320),ylim=c(220,320),
         main=scene[i])    
    points(vphe[[i]],pixmean[[i]],col='red',pch=i)    
    abline(0,1)
    #     abline(lm(pixmean[[i]]~mphe[[i]]))
  } 
}

# Mean with error bars
plot_1to1 <- function(xmean,ymean,xsd,ysd,xlabel,ylabel,start){
  mycol <- brewer.pal(6,'Set1')
  plot(xmean,ymean,
       xlim=c(start,(start+60)),ylim=c(start,(start+60)),
       xlab=xlabel,
       ylab=ylabel,
       cex.lab=1.1,cex.axis=1.0)
  arrows(xmean,ymean-ysd,xmean,ymean+ysd,
         length=0.05,angle=0)
  arrows(xmean-xsd,ymean,xmean+xsd,ymean,
         length=0.05,angle=0)
  abline(0,1)
  abline(lm(ymean~xmean),lwd=1.3,lty=5)
  points(xmean,ymean,bg=mycol,pch=c(rep(21,6),rep(23,5)),cex=1.7)
  legend('topleft',c('Harvard','Bartlett','Cary','Hubbard','Proctor','Boundary waters'),
         bty='n',pch=21,pt.bg=mycol,pt.cex=1.5)
  reg <- summary(lm(ymean~xmean))
  rseq <- round(reg$r.squared,2)
  b1 <- round(reg$coefficients[2],1)
  b1c1 <- b1+round(reg$coefficients[4],1)
  b1c2 <- b1-round(reg$coefficients[4],1)
  b0 <- round(reg$coefficients[1],1)
  b0c1 <- b0+round(reg$coefficients[3],1)
  b0c2 <- b0-round(reg$coefficients[3],1)
  rmse <- round(sqrt(mean((ymean-xmean)^2)),1)
  text((start+33),(start+13),expression(paste(R^2,'=')),pos=4,cex=0.8)
  text((start+38),(start+13),rseq,pos=4,cex=0.8)
  text((start+45),(start+13),substitute(paste('(', italic(p),' < 0.01)',sep='')),pos=4,cex=0.8)
  text((start+33),(start+9),paste('RMSE =',rmse),pos=4,cex=0.8)
  text((start+33),(start+5),expression(paste(beta[1],'=')),pos=4,cex=0.8)
  text((start+38),(start+5),paste(b1,' (',b1c1,',',b1c2,')',sep=''),pos=4,cex=0.8)
  text((start+33),(start+1),expression(paste(beta[0],'=')),pos=4,cex=0.8)
  text((start+38),(start+1),paste(b0,' (',b0c1,',',b0c2,')',sep=''),pos=4,cex=0.8)
}

par(mfrow=c(1,3),oma=c(1,1,1,1),mar=c(4,4,4,4),mgp=c(2.5,1,0))
plot_1to1(mesd[,1],mesd[,2],mesd[,4],mesd[,5],'MODIS EOS','VIIRS EOS',220)
plot_1to1(mesd[,1],mesd[,3],mesd[,4],mesd[,6],'MODIS EOS','Landsat EOS',220)
plot_1to1(mesd[,2],mesd[,3],mesd[,5],mesd[,6],'VIIRS EOS','Landsat EOS',220)








# ### Percentile for LSP
# for(i in c(1,4)){
#   difmm <- lphe[[i]] - mphe[[i]]
#   difvv <- lphe[[i]] - vphe[[i]]
#   
#   if(su)
#   ecmm <- ecdf(difmm[1,])
#   ecvv <- ecdf(difvv)
#   
#   ranbi <- -20:20
#   plot(ecmm(ranbi)*100,ranbi,type='l',col='blue',
#        main=scene[i],
#        xlab='Percentile for LSP (%)',
#        ylab='Bias (Days)')
#   legend('topleft',c('VIIRS','MODIS'),col=c('red','blue'),lty=1)
#   lines(ecvv(ranbi)*100,ranbi,col='red')
#   abline(v=ecmm(0)*100,col='blue')
#   abline(v=ecvv(0)*100,col='red')
#   abline(h=0,lty=5)  
# }
