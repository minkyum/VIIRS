library(RColorBrewer)

## Select collection (6 for C6 and 5 for C5)
clt <- 6

## Save fitures as PDF
setwd('/projectnb/modislc/users/mkmoon/VIIRS/figures/')
pdf(file=paste('1to1_percentile_c',clt,'.pdf',sep=''),width=12,height=8)  

### Load data
setwd('/projectnb/modislc/users/mkmoon/VIIRS/R_data/landsat/')
scene <- c('mwp_harvard','ah_bartlett','mwp_cary','ah_hubbard','mwp_proctor','mws_boundary_waters')
lphe <- vector('list',12)
mphe <- vector('list',12)
vphe <- vector('list',12)
npix <- vector('list',12)

### SOS
if(clt==6){
  for(i in 1:12){
    if(i<=6){
      load(paste('landsat_c6_2013_sos_',scene[i],'.RData',sep=''))  
    }else{
      load(paste('landsat_c6_2012_sos_',scene[i-6],'.RData',sep=''))  
    }
    lphe[[i]] <- lanbymod
    mphe[[i]] <- mmphe
    vphe[[i]] <- vvphe
    npix[[i]] <- nlpix  
  }  
}else{
  for(i in 1:12){
    if(i<=6){
      load(paste('landsat_c5_2013_sos_',scene[i],'.RData',sep=''))  
    }else{
      load(paste('landsat_c5_2012_sos_',scene[i-6],'.RData',sep=''))  
    }
    lphe[[i]] <- lanbymod
    mphe[[i]] <- mmphe
    vphe[[i]] <- vvphe
    npix[[i]] <- nlpix  
  }
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
# for(i in 1:11){
#   if(i==1){
#     plot(mphe[[i]],pixmean[[i]],pch=i,
#          xlim=c(100,180),ylim=c(100,180),
#          main=scene[i])  
#     points(vphe[[i]],pixmean[[i]],col='red',pch=i)    
#     abline(0,1,lty=5)
# #     abline(lm(pixmean[[i]]~mphe[[i]]))
#   }else{
#     plot(mphe[[i]],pixmean[[i]],pch=i,
#          xlim=c(100,180),ylim=c(100,180),
#          main=scene[i])    
#     points(vphe[[i]],pixmean[[i]],col='red',pch=i)    
#     abline(0,1)
# #     abline(lm(pixmean[[i]]~mphe[[i]]))
#   } 
# }

# Mean with error bars
plot_1to1 <- function(xmean,ymean,xsd,ysd,xlabel,ylabel,start){
  mycol <- brewer.pal(6,'Set1')
  plot(xmean,ymean,
       xlim=c(start,(start+60)),ylim=c(start,(start+60)),
       xlab=xlabel,
       ylab=ylabel,
       cex.lab=1.5,cex.axis=1.3)
  arrows(xmean,ymean-ysd,xmean,ymean+ysd,
         length=0.05,angle=0)
  arrows(xmean-xsd,ymean,xmean+xsd,ymean,
         length=0.05,angle=0)
  abline(0,1,lty=2)
  abline(lm(ymean~xmean),lwd=1)
  points(xmean,ymean,bg=mycol,pch=c(rep(21,6),rep(23,5)),cex=1.7)
  legend('topleft',c('Harvard','Bartlett','Cary','Hubbard','Proctor','Boundary waters'),
         bty='n',pch=21,pt.bg=mycol,pt.cex=1.5)
  reg <- summary(lm(ymean~xmean))
  rseq <- round(reg$r.squared,2)
  b1 <- round(reg$coefficients[2],1)
  b1c1 <- round(b1+reg$coefficients[4],1)
  b1c2 <- round(b1-reg$coefficients[4],1)
  b0 <- round(reg$coefficients[1],1)
  b0c1 <- round(b0+reg$coefficients[3],1)
  b0c2 <- round(b0-reg$coefficients[3],1)
  rmse <- round(sqrt(mean((ymean-xmean)^2)),1)
  text((start+33-1),(start+13),expression(paste(R^2,'=')),pos=4,cex=1)
  text((start+38-1),(start+13),rseq,pos=4,cex=1)
  if(reg$coefficient[8]<0.01){
    text((start+45-1),(start+13),substitute(paste('(', italic(p),' < 0.01)',sep='')),pos=4,cex=1)
  }else{
    text((start+45-1),(start+13),substitute(paste('(', italic(p),' = ',sep='')),pos=4,cex=1)
    text((start+50-1),(start+13),paste(round(reg$coefficient[8],2),' )',sep=''),pos=4,cex=1)
  }
  text((start+33-1),(start+9),paste('RMSE =',rmse),pos=4,cex=1)
  text((start+33-1),(start+5),expression(paste(beta[1],'=')),pos=4,cex=1)
  text((start+38-1),(start+5),paste(b1,' (',b1c1,',',b1c2,')',sep=''),pos=4,cex=1)
  text((start+33-1),(start+1),expression(paste(beta[0],'=')),pos=4,cex=1)
  text((start+38-1),(start+1),paste(b0,' (',b0c1,',',b0c2,')',sep=''),pos=4,cex=1)

  return(reg$coefficient[8])
}

par(mfrow=c(2,3),oma=c(1,1,1,1),mar=c(4,4,4,4),mgp=c(2.5,1,0))
plot_1to1(mesd[,1],mesd[,2],mesd[,4],mesd[,5],'MODIS SOS','VIIRS SOS',110)
plot_1to1(mesd[,1],mesd[,3],mesd[,4],mesd[,6],'MODIS SOS','Landsat SOS',110)
plot_1to1(mesd[,2],mesd[,3],mesd[,5],mesd[,6],'VIIRS SOS','Landsat SOS',110)


### EOS
if(clt==6){
  for(i in 1:12){
    if(i<=6){
      load(paste('landsat_c6_2013_eos_',scene[i],'.RData',sep=''))  
    }else{
      load(paste('landsat_c6_2012_eos_',scene[i-6],'.RData',sep=''))  
    }
    lphe[[i]] <- lanbymod
    mphe[[i]] <- mmphe
    vphe[[i]] <- vvphe
    npix[[i]] <- nlpix  
  }
}else{
  for(i in 1:12){
    if(i<=6){
      load(paste('landsat_c5_2013_eos_',scene[i],'.RData',sep=''))  
    }else{
      load(paste('landsat_c5_2012_eos_',scene[i-6],'.RData',sep=''))  
    }
    lphe[[i]] <- lanbymod
    mphe[[i]] <- mmphe
    vphe[[i]] <- vvphe
    npix[[i]] <- nlpix  
  }  
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
# for(i in 1:11){
#   if(i==1){
#     plot(mphe[[i]],pixmean[[i]],pch=i,
#          xlim=c(220,320),ylim=c(220,320),
#          main=scene[i])  
#     points(vphe[[i]],pixmean[[i]],col='red',pch=i)    
#     abline(0,1,lty=5)
#     #     abline(lm(pixmean[[i]]~mphe[[i]]))
#   }else{
#     plot(mphe[[i]],pixmean[[i]],pch=i,
#          xlim=c(220,320),ylim=c(220,320),
#          main=scene[i])    
#     points(vphe[[i]],pixmean[[i]],col='red',pch=i)    
#     abline(0,1)
#     #     abline(lm(pixmean[[i]]~mphe[[i]]))
#   } 
# }

# Mean with error bars
plot_1to1(mesd[,1],mesd[,2],mesd[,4],mesd[,5],'MODIS EOS','VIIRS EOS',240)
plot_1to1(mesd[,1],mesd[,3],mesd[,4],mesd[,6],'MODIS EOS','Landsat EOS',240)
plot_1to1(mesd[,2],mesd[,3],mesd[,5],mesd[,6],'VIIRS EOS','Landsat EOS',240)



### Percentile for LSP
library(vioplot)
mycol <- brewer.pal(5,'Set1')
mvio <- function (x, ..., range = 1.5, h = NULL, ylim = NULL, names = NULL, 
                  horizontal = FALSE, col = "magenta", border = "black", lty = 1, 
                  lwd = 1, rectCol = "black", colMed = "white", pchMed = 19, 
                  at, add = FALSE, wex = 1, drawRect = TRUE){
  datas <- list(x, ...)
  n <- length(datas)
  if (missing(at)) 
    at <- 1:n
  upper <- vector(mode = "numeric", length = n)
  lower <- vector(mode = "numeric", length = n)
  q1 <- vector(mode = "numeric", length = n)
  q3 <- vector(mode = "numeric", length = n)
  med <- vector(mode = "numeric", length = n)
  base <- vector(mode = "list", length = n)
  height <- vector(mode = "list", length = n)
  baserange <- c(Inf, -Inf)
  args <- list(display = "none")
  if (!(is.null(h))) 
    args <- c(args, h = h)
  for (i in 1:n) {
    data <- datas[[i]]
    data.min <- min(data)
    data.max <- max(data)
    q1[i] <- quantile(data, 0.25)
    q3[i] <- quantile(data, 0.75)
    med[i] <- median(data)
    iqd <- q3[i] - q1[i]
    upper[i] <- min(q3[i] + range * iqd, data.max)
    lower[i] <- max(q1[i] - range * iqd, data.min)
    est.xlim <- c(min(lower[i], data.min), max(upper[i], 
                                               data.max))
    smout <- do.call("sm.density", c(list(data, xlim = est.xlim), 
                                     args))
    hscale <- 0.4/max(smout$estimate) * wex
    base[[i]] <- smout$eval.points
    height[[i]] <- smout$estimate * hscale
    t <- range(base[[i]])
    baserange[1] <- min(baserange[1], t[1])
    baserange[2] <- max(baserange[2], t[2])
  }
  if (!add) {
    xlim <- if (n == 1) 
      at + c(-0.5, 0.5)
    else range(at) + min(diff(at))/2 * c(-1, 1)
    if (is.null(ylim)) {
      ylim <- baserange
    }
  }
  if (is.null(names)) {
    label <- 1:n
  }
  else {
    label <- names
  }
  boxwidth <- 0.05 * wex
  if (!add) 
    plot.new()
  if (!horizontal) {
    if (!add) {
      plot.window(xlim = xlim, ylim = ylim)
      axis(2)
      axis(1, at = at, label = label)
    }
    box()
    for (i in 1:n) {
      polygon(c(at[i] - height[[i]], rev(at[i] + height[[i]])), 
              c(base[[i]], rev(base[[i]])), col = col[i], border = border, 
              lty = lty, lwd = lwd)
      if (drawRect) {
        lines(at[c(i, i)], c(lower[i], upper[i]), lwd = lwd, 
              lty = lty)
        rect(at[i] - boxwidth/2, q1[i], at[i] + boxwidth/2, 
             q3[i], col = rectCol)
        points(at[i], med[i], pch = pchMed, col = colMed)
      }
    }
  }
  else {
    if (!add) {
      plot.window(xlim = ylim, ylim = xlim)
      axis(1)
      axis(2, at = at, label = label)
    }
    box()
    for (i in 1:n) {
      polygon(c(base[[i]], rev(base[[i]])), c(at[i] - height[[i]], 
                                              rev(at[i] + height[[i]])), col = col[i], border = border, 
              lty = lty, lwd = lwd)
      if (drawRect) {
        lines(c(lower[i], upper[i]), at[c(i, i)], lwd = lwd, 
              lty = lty)
        rect(q1[i], at[i] - boxwidth/2, q3[i], at[i] + 
               boxwidth/2, col = rectCol)
        points(med[i], at[i], pch = pchMed, col = colMed)
      }
    }
  }
  invisible(list(upper = upper, lower = lower, median = med, 
                 q1 = q1, q3 = q3))
}
plot_vio <- function(phe){
  # Get quantile
  mq <- vector('list',11)
  vq <- vector('list',11)
  for(i in 1:11){
    for(j in 1:1000){
      temp <- ecdf(lphe[[i]][j,])
      mq[[i]][j] <- temp(mphe[[i]][j])
      vq[[i]][j] <- temp(vphe[[i]][j])  
    }
    mq[[i]] <- na.omit(mq[[i]])
    vq[[i]] <- na.omit(vq[[i]])
  }
  
  par(mfrow=c(2,1),oma=c(1,1,1,1),mar=c(3,3,3,3),mgp=c(2.5,1,0))
  plot(0:1,0:1,xlim=c(0,14.5),ylim=c(0,1),type='n',axes=F,ann=F)
  mvio(mq[[1]],vq[[1]],mq[[2]],vq[[2]],mq[[3]],vq[[3]],
       mq[[4]],vq[[4]],mq[[5]],vq[[5]],
       at=c(1,2,4,5,7,8,10,11,13,14),col=rep(mycol,each=2),add=T)
  axis(1,at=c(1,2,4,5,7,8,10,11,13,14),
       c('M.Ha','V.Ha','M.Ba','V.Ba','M.Ca','V.Ca','M.Hu','V.Hu','M.Pr','V.Pr'),
       cex.axis=1.1)
  axis(2,at=seq(0,1,0.2),seq(0,100,20),cex.axis=0.9)
  mtext('Percentile for LSP (%)',2,line=2.5,cex=1.1)
  text(0,0.9,phe,cex=0.8)
  text(0,0.8,2013,cex=0.8)
  
  plot(0:1,0:1,xlim=c(0,14.5),ylim=c(0,1),type='n',axes=F,ann=F)
  mvio(mq[[6]],vq[[6]],mq[[7]],vq[[7]],mq[[8]],vq[[8]],
       mq[[9]],vq[[9]],mq[[10]],vq[[10]],
       at=c(1,2,4,5,7,8,10,11,13,14),col=rep(mycol,each=2),add=T)
  axis(1,at=c(1,2,4,5,7,8,10,11,13,14),
       c('M.Ha','V.Ha','M.Ba','V.Ba','M.Ca','V.Ca','M.Hu','V.Hu','M.Pr','V.Pr'),
       cex.axis=1.1)
  axis(2,at=seq(0,1,0.2),seq(0,100,20),cex.axis=0.9)
  mtext('Percentile for LSP (%)',2,line=2.5,cex=1.1)
  text(0,0.9,phe,cex=0.8)
  text(0,0.8,2012,cex=0.8)
}

## SOS
if(clt==6){
  for(i in 1:12){
    if(i<=6){
      load(paste('landsat_c6_2013_sos_',scene[i],'.RData',sep=''))  
    }else{
      load(paste('landsat_c6_2012_sos_',scene[i-6],'.RData',sep=''))  
    }
    lphe[[i]] <- lanbymod
    mphe[[i]] <- mmphe
    vphe[[i]] <- vvphe
    npix[[i]] <- nlpix  
  }  
}else{
  for(i in 1:12){
    if(i<=6){
      load(paste('landsat_c5_2013_sos_',scene[i],'.RData',sep=''))  
    }else{
      load(paste('landsat_c5_2012_sos_',scene[i-6],'.RData',sep=''))  
    }
    lphe[[i]] <- lanbymod
    mphe[[i]] <- mmphe
    vphe[[i]] <- vvphe
    npix[[i]] <- nlpix  
  }
}
plot_vio('SOS')

## EOS
if(clt==6){
  for(i in 1:12){
    if(i<=6){
      load(paste('landsat_c6_2013_eos_',scene[i],'.RData',sep=''))  
    }else{
      load(paste('landsat_c6_2012_eos_',scene[i-6],'.RData',sep=''))  
    }
    lphe[[i]] <- lanbymod
    mphe[[i]] <- mmphe
    vphe[[i]] <- vvphe
    npix[[i]] <- nlpix  
  }
}else{
  for(i in 1:12){
    if(i<=6){
      load(paste('landsat_c5_2013_eos_',scene[i],'.RData',sep=''))  
    }else{
      load(paste('landsat_c5_2012_eos_',scene[i-6],'.RData',sep=''))  
    }
    lphe[[i]] <- lanbymod
    mphe[[i]] <- mmphe
    vphe[[i]] <- vvphe
    npix[[i]] <- nlpix  
  }  
}
plot_vio('EOS')

dev.off()
