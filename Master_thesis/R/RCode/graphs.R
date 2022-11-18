setwd("~/Desktop/Thesis/R/test29_allall")
options(scipen=999)

library(readxl)
library(matrixcalc)
library(optimx)

load(file="../RData/Switz/results_KF/estimated_var.RData")
load(file="../RData/Switz/results_KF/estimated_DN_var.RData")

data <- read_xlsx("../../Data/data_inf_var_switz.xlsx",sheet=1)
data <- data.frame(data)
data<-merge(data,data.frame(year=1870:2020),by="year",all=TRUE)

obs=3
nb=(dim(data)[2]-1)/obs

countries_DN <- c("CAN","FRA","DEU","ITA","JPN","CHE","GBR","USA")
colors=c("orangered1","gold","royalblue3","seagreen3","plum4","hotpink","turquoise","peachpuff4")

CIs <- function(j,xi.est,P.est){
  up_95<-matrix(0)
  down_95<-matrix(0)
  up_68<-matrix(0)
  down_68<-matrix(0)
  T <- dim(xi.est)[1]
  r <- sqrt(dim(P.est)[2])
  for (t in 1:T) {
    if (!is.na(xi.est[t,j])){
      up_95[t] <- xi.est[t,j]+1.96*sqrt(matrix(P.est[t,],r,r)[j,j])
      down_95[t] <- xi.est[t,j]-1.96*sqrt(matrix(P.est[t,],r,r)[j,j])
      up_68[t] <- xi.est[t,j]+sqrt(matrix(P.est[t,],r,r)[j,j])
      down_68[t] <- xi.est[t,j]-sqrt(matrix(P.est[t,],r,r)[j,j])
    }
  }
  return(cbind(cbind(down_95,up_95),cbind(down_68,up_68)))
}

# ts double 
{
  j = 1+nb+1
  ts_data <- data.frame(data$year)
  ts_data$tsCanada <- data$ltrateCanada-data$stirCanada
  ts_data$tsFrance <- data$ltrateFrance-data$stirFrance
  ts_data$tsGermany <- data$ltrateGermany-data$stirGermany
  ts_data$tsItaly <- data$ltrateItaly-data$stirItaly
  ts_data$tsJapan <- data$ltrateJapan-data$stirJapan
  ts_data$tsSwitzerland <- data$ltrateSwitzerland-data$stirSwitzerland
  ts_data$tsUK <- data$ltrateUK-data$stirUK
  ts_data$tsUSA <- data$ltrateUSA-data$stirUSA
  
  FILE="/../Plots/Switz/All_sigma/RLR_ts.pdf"
  pdf(file=paste(getwd(),FILE,sep=""),pointsize=10, width=10, height=4)
  par(mfrow=c(1,2), oma = c(2, 0, 0, 0), mar = c(4.5,2,1,1))
  plot(data$year,ts_world$tsw,
       xaxt="non",type="l",xlab="",ylab="",lty=2,lwd=2,ylim=c(-2,4),xaxs="i")
  axis(1,seq(1880,2020,20))
  abline(h=0,col=adjustcolor(("black"),alpha.f=.3))
  for (i in 1:nb) {
    lines(data$year,ts_data[,i+1],type="l",xlab="",ylab="",lty="dotted",col=colors[i])
  }
  CI <- CIs(j,xi.est.smoothed,P.est.smoothed)
  polygon(c(data$year[1:length(which(!is.na(xi.est.smoothed[,j])))],
            rev(data$year[1:length(which(!is.na(xi.est.smoothed[,j])))])),
          c(CI[,1],
            rev(CI[,2])),
          col=adjustcolor("black", alpha.f = 0.2), border=NA)
  polygon(c(data$year[1:length(which(!is.na(xi.est.smoothed[,j])))],
            rev(data$year[1:length(which(!is.na(xi.est.smoothed[,j])))])),
          c(CI[,3],
            rev(CI[,4])),
          col=adjustcolor("black", alpha.f = 0.3), border=NA)
  plot(data$year,ts_world$tsw,
       xaxt="non",type="l",xlab="",ylab="",lty=2,lwd=2,ylim=c(-2,4),xaxs="i")
  axis(1,seq(1880,2020,20))
  abline(h=0,col=adjustcolor(("black"),alpha.f=.3))
  for (i in 1:nb) {
    lines(data$year,ts_i[,1+i],type="l",xlab="",ylab="",lty="dotted",lwd=2,col=colors[i])
  }
  par(fig = c(0, 1, 0, 1), oma = c(1, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
  legend("bottom",legend = c(countries_DN), bty = 'n',
         col=c(colors), lwd=1.5, cex=1.3, horiz = TRUE,xpd=TRUE)
  dev.off()
}

# ts + pi
{
  FILE="/../Plots/Switz/All_sigma/pi+ts.pdf"
  pdf(file=paste(getwd(),FILE,sep=""),pointsize=10, width=6, height=4)
  par(mfrow=c(2,2), oma=c(5, 0, 0, 0), mar=c(2,2,1,1))
  #pi_w and pi
  {
    j = 1+nb+1+nb+1
    plot(pi_world$year,pi_world$piw,xaxt="non",type="l",xlab="",ylab="",lty=2,lwd=2,ylim=c(-2,15),xaxs="i")
    axis(1,seq(1880,2020,20))
    abline(h=0,col=adjustcolor(("black"),alpha.f=.3))
    for (i in 1:nb) {
      lines(data$year,data[,(1+nb+nb+i)],type="l",lty="dotted",col=colors[i])
    }
    CI <- CIs(j,xi.est.smoothed,P.est.smoothed)
    polygon(c(data$year[1:length(which(!is.na(xi.est.smoothed[,j])))],
              rev(data$year[1:length(which(!is.na(xi.est.smoothed[,j])))])),
            c(CI[,1],
              rev(CI[,2])),
            col=adjustcolor("black", alpha.f = 0.2), border=NA)
    polygon(c(data$year[1:length(which(!is.na(xi.est.smoothed[,j])))],
              rev(data$year[1:length(which(!is.na(xi.est.smoothed[,j])))])),
            c(CI[,3],
              rev(CI[,4])),
            col=adjustcolor("black", alpha.f = 0.3), border=NA)
  }
  
  # pi_w and pi_i
  {
    plot(data$year,pi_world$piw,xaxt="non",type="l",xlab="",ylab="",lty=2,lwd=2,ylim=c(-2,15),xaxs="i")
    axis(1,seq(1880,2020,20))
    abline(h=0,col=adjustcolor(("black"),alpha.f=.3))
    for (i in 1:nb) {
      lines(data$year,pi_i[,i+1],type="l",lty="dotted",lwd=2,col=colors[i])
    }
  }
  
  #ts_w and RL-R
  {
    j = 1+nb+1
    ts_data <- data.frame(data$year)
    ts_data$tsCanada <- data$ltrateCanada-data$stirCanada
    ts_data$tsFrance <- data$ltrateFrance-data$stirFrance
    ts_data$tsGermany <- data$ltrateGermany-data$stirGermany
    ts_data$tsItaly <- data$ltrateItaly-data$stirItaly
    ts_data$tsJapan <- data$ltrateJapan-data$stirJapan
    ts_data$tsSwitzerland <- data$ltrateSwitzerland-data$stirSwitzerland
    ts_data$tsUK <- data$ltrateUK-data$stirUK
    ts_data$tsUSA <- data$ltrateUSA-data$stirUSA
    
    plot(data$year,ts_world$tsw,
         xaxt="non",type="l",xlab="",ylab="",lty=2,lwd=2,ylim=c(-2,4),xaxs="i")
    axis(1,seq(1880,2020,20))
    abline(h=0,col=adjustcolor(("black"),alpha.f=.3))
    for (i in 1:nb) {
      lines(data$year,ts_data[,i+1],type="l",xlab="",ylab="",lty="dotted",col=colors[i])
    }
    CI <- CIs(j,xi.est.smoothed,P.est.smoothed)
    polygon(c(data$year[1:length(which(!is.na(xi.est.smoothed[,j])))],
              rev(data$year[1:length(which(!is.na(xi.est.smoothed[,j])))])),
            c(CI[,1],
              rev(CI[,2])),
            col=adjustcolor("black", alpha.f = 0.2), border=NA)
    polygon(c(data$year[1:length(which(!is.na(xi.est.smoothed[,j])))],
              rev(data$year[1:length(which(!is.na(xi.est.smoothed[,j])))])),
            c(CI[,3],
              rev(CI[,4])),
            col=adjustcolor("black", alpha.f = 0.3), border=NA)
  }
  
  #ts_w and ts_i
  {
    plot(data$year,ts_world$tsw,
         xaxt="non",type="l",xlab="",ylab="",lty=2,lwd=2,ylim=c(-2,4),xaxs="i")
    axis(1,seq(1880,2020,20))
    abline(h=0,col=adjustcolor(("black"),alpha.f=.3))
    for (i in 1:nb) {
      lines(data$year,ts_i[,1+i],type="l",xlab="",ylab="",lty="dotted",lwd=2,col=colors[i])
    }
  }
  par(fig = c(0, 1, 0, 1), oma = c(2, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
  legend("bottom",legend = c(countries_DN), bty = 'n',
         col=c(colors), lwd=1.5, cex=.95, horiz = TRUE,xpd=TRUE)
  dev.off()
}

# r_w and R-pi
{
  j = 1
  real_data <- data.frame(year=data$year)
  real_data$realCanada <- data$stirCanada-data$infCanada
  real_data$realFrance <- data$stirFrance - data$infFrance
  real_data$realGermany <- data$stirGermany - data$infGermany
  real_data$realItaly <- data$stirItaly - data$infItaly
  real_data$realJapan <- data$stirJapan - data$infJapan
  real_data$realSwitzerland <- data$stirSwitzerland- data$infSwitzerland
  real_data$realUK <- data$stirUK - data$infUK
  real_data$realUSA <- data$stirUSA - data$infUSA
  
  par(mfrow=c(1,1))
  FILE="/../Plots/Switz/All_sigma/Rpi.pdf"
  pdf(file=paste(getwd(),FILE,sep=""),pointsize=10, width=6, height=4)
  par(oma = c(2, 0, 0, 0), mar = c(4.5,2,1,1))
  plot(data$year,r_world[,2],
       xaxt="non",type="l",xlab="",ylab="",lty=2,lwd=2,ylim=c(-5,11),xaxs="i")
  axis(1,seq(1880,2020,20))
  abline(h=0,col=adjustcolor(("black"),alpha.f=.3))
  for (i in 1:nb) {
    lines(real_data$year,real_data[,i+1],type="l",lty=3,col=colors[i])
  }
  
  CI <- CIs(j,xi.est.smoothed,P.est.smoothed)
  polygon(c(data$year[1:length(which(!is.na(xi.est.smoothed[,j])))],
            rev(data$year[1:length(which(!is.na(xi.est.smoothed[,j])))])),
          c(CI[,1],
            rev(CI[,2])),
          col=adjustcolor("black", alpha.f = 0.2), border=NA)
  polygon(c(data$year[1:length(which(!is.na(xi.est.smoothed[,j])))],
            rev(data$year[1:length(which(!is.na(xi.est.smoothed[,j])))])),
          c(CI[,3],
            rev(CI[,4])),
          col=adjustcolor("black", alpha.f = 0.3), border=NA)

   legend(1870,-8.5,legend=c(NULL,countries_DN),
         col=c(NULL,colors),
         lty=c(NULL,c(matrix(3,nb))),lwd=c(NULL,matrix(2,nb)),
         bty="n",xpd=TRUE,horiz=TRUE,cex=.85)
  dev.off()
}

#r_w and r_i
{
  
  FILE="/../Plots/Switz/All_sigma/r_i.pdf"
  pdf(file=paste(getwd(),FILE,sep=""),pointsize=10, width=6, height=4)
  par(oma = c(2, 0, 0, 0), mar = c(4.5,2,1,1))
  plot(r_world$year,r_world$rw,xaxt="non",type="l",xlab="",ylab="",lty=2,lwd=2,ylim=c(-3,6),xaxs="i")
  axis(1,seq(1880,2000,20))
  abline(h=0,col=adjustcolor(("black"),alpha.f=.3))
  for (i in 1:nb) {
    lines(r_world$year,r_i[,i+1],type="l",lty="dotted",lwd=2,col=colors[i])
  }
  legend(1870,-5,legend=c(NULL,countries_DN),
         col=c(NULL,colors),
         lty=c(NULL,matrix(3,nb)),lwd=c(NULL,matrix(2,nb)),
         bty="n",xpd=TRUE,horiz=TRUE,cex=.85
  )
  dev.off()
}

#pi_w and pi
{
  FILE="/../Plots/Switz/All_sigma/Pi.pdf"
  pdf(file=paste(getwd(),FILE,sep=""),pointsize=10, width=6, height=4)
  par(oma = c(2, 0, 0, 0), mar = c(4.5,2,1,1))
  # par(oma = c(2, 0, 0, 0), mar = c(1,2,1,1))
  j = 1+nb+1+nb+1
  plot(pi_world$year,pi_world$piw,xaxt="non",type="l",xlab="",ylab="",lty=2,lwd=2,ylim=c(-2,15),xaxs="i")
  axis(1,seq(1880,2020,20))
  abline(h=0,col=adjustcolor(("black"),alpha.f=.3))
  for (i in 1:nb) {
    lines(data$year,data[,(1+nb+nb+i)],type="l",lty="dotted",col=colors[i])
  }
  CI <- CIs(j,xi.est.smoothed,P.est.smoothed)
  polygon(c(data$year[1:length(which(!is.na(xi.est.smoothed[,j])))],
            rev(data$year[1:length(which(!is.na(xi.est.smoothed[,j])))])),
          c(CI[,1],
            rev(CI[,2])),
          col=adjustcolor("black", alpha.f = 0.2), border=NA)
  polygon(c(data$year[1:length(which(!is.na(xi.est.smoothed[,j])))],
            rev(data$year[1:length(which(!is.na(xi.est.smoothed[,j])))])),
          c(CI[,3],
            rev(CI[,4])),
          col=adjustcolor("black", alpha.f = 0.3), border=NA)
  legend(1870,-6,legend=c(NULL,countries_DN),
         col=c(NULL,colors),
         lty=c(NULL,matrix(3,nb)),lwd=c(NULL,matrix(2,nb)),
         bty="n",xpd=TRUE,horiz=TRUE,cex=.85
  )
  dev.off()
}

# pi_w and pi_i
{
  FILE="/../Plots/Switz/All_sigma/Pi_i.pdf"
  pdf(file=paste(getwd(),FILE,sep=""),pointsize=10, width=6, height=4)
  par(oma = c(2, 0, 0, 0), mar = c(4.5,2,1,1))
  # par(oma = c(2, 0, 0, 0), mar = c(1,2,1,1))
  plot(data$year,pi_world$piw,xaxt="non",type="l",xlab="",ylab="",lty=2,lwd=2,ylim=c(-2,15),xaxs="i")
  axis(1,seq(1880,2020,20))
  abline(h=0,col=adjustcolor(("black"),alpha.f=.3))
  for (i in 1:nb) {
    lines(data$year,pi_i[,i+1],type="l",lty="dotted",lwd=2,col=colors[i])
  }
  legend(1870,-6,legend=c(NULL,countries_DN),
         col=c(NULL,colors),
         lty=c(NULL,matrix(3,nb)),lwd=c(NULL,matrix(2,nb)),
         bty="n",xpd=TRUE,horiz=TRUE,cex=.85
  )
  dev.off()
}

#ts_w and RL-R
{
  j = 1+nb+1
  ts_data <- data.frame(data$year)
  ts_data$tsCanada <- data$ltrateCanada-data$stirCanada
  ts_data$tsFrance <- data$ltrateFrance-data$stirFrance
  ts_data$tsGermany <- data$ltrateGermany-data$stirGermany
  ts_data$tsItaly <- data$ltrateItaly-data$stirItaly
  ts_data$tsJapan <- data$ltrateJapan-data$stirJapan
  ts_data$tsSwitzerland <- data$ltrateSwitzerland-data$stirSwitzerland
  ts_data$tsUK <- data$ltrateUK-data$stirUK
  ts_data$tsUSA <- data$ltrateUSA-data$stirUSA
  
  FILE="/../Plots/Switz/All_sigma/RLR.pdf"
  pdf(file=paste(getwd(),FILE,sep=""),pointsize=10, width=6, height=4)
  par(oma = c(2, 0, 0, 0), mar = c(4.5,2,1,1))
  plot(data$year,ts_world$tsw,
       xaxt="non",type="l",xlab="",ylab="",lty=2,lwd=2,ylim=c(-2,4),xaxs="i")
  axis(1,seq(1880,2020,20))
  abline(h=0,col=adjustcolor(("black"),alpha.f=.3))
  for (i in 1:nb) {
    lines(data$year,ts_data[,i+1],type="l",xlab="",ylab="",lty="dotted",col=colors[i])
  }
  CI <- CIs(j,xi.est.smoothed,P.est.smoothed)
  polygon(c(data$year[1:length(which(!is.na(xi.est.smoothed[,j])))],
            rev(data$year[1:length(which(!is.na(xi.est.smoothed[,j])))])),
          c(CI[,1],
            rev(CI[,2])),
          col=adjustcolor("black", alpha.f = 0.2), border=NA)
  polygon(c(data$year[1:length(which(!is.na(xi.est.smoothed[,j])))],
            rev(data$year[1:length(which(!is.na(xi.est.smoothed[,j])))])),
          c(CI[,3],
            rev(CI[,4])),
          col=adjustcolor("black", alpha.f = 0.3), border=NA)
  legend(1870,-3.5,legend=c(NULL,countries_DN),
         col=c(NULL,colors),
         lty=c(NULL,matrix(3,nb)),lwd=c(NULL,matrix(2,nb)),
         bty="n",xpd=TRUE,horiz=TRUE,cex=.85
  )
  dev.off()
}

#ts_w and ts_i
{
  FILE="/../Plots/Switz/All_sigma/ts_i.pdf"
  pdf(file=paste(getwd(),FILE,sep=""),pointsize=10, width=6, height=4)
  par(oma = c(2, 0, 0, 0), mar = c(4.5,2,1,1))
  plot(data$year,ts_world$tsw,
       xaxt="non",type="l",xlab="",ylab="",lty=2,lwd=2,ylim=c(-2,4),xaxs="i")
  axis(1,seq(1880,2020,20))
  abline(h=0,col=adjustcolor(("black"),alpha.f=.3))
  for (i in 1:nb) {
    lines(data$year,ts_i[,1+i],type="l",xlab="",ylab="",lty="dotted",lwd=2,col=colors[i])
  }
  legend(1870,-3.5,legend=c(NULL,countries_DN),
         col=c(NULL,colors),
         lty=c(NULL,matrix(3,nb)),lwd=c(NULL,matrix(2,nb)),
         bty="n",xpd=TRUE,horiz=TRUE,cex=.85
  )
  dev.off()
}

#1x3 r_w, pi_w + ts_w
{
  FILE="/../Plots/Switz/All_sigma/world_trends.pdf"
  pdf(file=paste(getwd(),FILE,sep=""),pointsize=10, width=6, height=3)
  par(mfrow=c(1,3),oma = c(2, 0, 0, 0),mar = c(5.5,2,1,1))
  # r_w
  j=1
  plot(data$year,r_world$rw,
       col="firebrick4",xaxt="non",type="l",xlab="",ylab="",lty=2,lwd=1.5,ylim=c(-1,8),xaxs="i")
  axis(1,seq(1880,2020,20))
  abline(h=0,col=adjustcolor(("black"),alpha.f=.3))
  CI <- CIs(j,xi.est.smoothed,P.est.smoothed)
  polygon(c(data$year[1:length(which(!is.na(xi.est.smoothed[,j])))],
            rev(data$year[1:length(which(!is.na(xi.est.smoothed[,j])))])),
          c(CIs(j,xi.est.smoothed,P.est.smoothed)[,1], rev(CIs(j,xi.est.smoothed,P.est.smoothed)[,2])),
          col=adjustcolor("firebrick4", alpha.f = 0.2), border=NA)
  polygon(c(data$year[1:length(which(!is.na(xi.est.smoothed[,j])))],
            rev(data$year[1:length(which(!is.na(xi.est.smoothed[,j])))])),
          c(CIs(j,xi.est.smoothed,P.est.smoothed)[,3], rev(CIs(j,xi.est.smoothed,P.est.smoothed)[,4])),
          col=adjustcolor("firebrick4", alpha.f = 0.3), border=NA)
  legend(1920,-3,legend=expression({bar(r)[t]}^w),col="firebrick4",lty="dashed",lwd=1.5,
         bty="n",xpd=TRUE,horiz=TRUE,cex=1.2)
  # pi_w
  j=1+nb+1+nb+1
  plot(data$year,pi_world$piw,xaxt="non",type="l",xlab="",ylab="",lty=2,lwd=1.5,ylim=c(-1,8),xaxs="i")
  axis(1,seq(1880,2000,20))
  abline(h=0,col=adjustcolor(("steelblue4"),alpha.f=.3))
  CI <- CIs(j,xi.est.smoothed,P.est.smoothed)
  polygon(c(data$year[1:length(which(!is.na(xi.est.smoothed[,j])))],
            rev(data$year[1:length(which(!is.na(xi.est.smoothed[,j])))])),
          c(CIs(j,xi.est.smoothed,P.est.smoothed)[,1], rev(CIs(j,xi.est.smoothed,P.est.smoothed)[,2])),
          col=adjustcolor("steelblue4", alpha.f = 0.2), border=NA)
  polygon(c(data$year[1:length(which(!is.na(xi.est.smoothed[,j])))],
            rev(data$year[1:length(which(!is.na(xi.est.smoothed[,j])))])),
          c(CIs(j,xi.est.smoothed,P.est.smoothed)[,3], rev(CIs(j,xi.est.smoothed,P.est.smoothed)[,4])),
          col=adjustcolor("steelblue4", alpha.f = 0.3), border=NA)
  legend(1920,-3,legend=expression({bar(pi)[t]}^w),col="steelblue4",lty="dashed",lwd=1.5,
         bty="n",xpd=TRUE,horiz=TRUE,cex=1.2)
  # ts_w
  j=1+nb+1
  plot(data$year,ts_world$tsw,
       xaxt="non",type="l",xlab="",ylab="",lty=2,lwd=1.5,ylim=c(-1,8),xaxs="i")
  axis(1,seq(1880,2000,20))
  abline(h=0,col=adjustcolor(("darkolivegreen"),alpha.f=.3))
  CI <- CIs(j,xi.est.smoothed,P.est.smoothed)
  polygon(c(data$year[1:length(which(!is.na(xi.est.smoothed[,j])))],
            rev(data$year[1:length(which(!is.na(xi.est.smoothed[,j])))])),
          c(CIs(j,xi.est.smoothed,P.est.smoothed)[,1], rev(CIs(j,xi.est.smoothed,P.est.smoothed)[,2])),
          col=adjustcolor("darkolivegreen", alpha.f = 0.2), border=NA)
  polygon(c(data$year[1:length(which(!is.na(xi.est.smoothed[,j])))],
            rev(data$year[1:length(which(!is.na(xi.est.smoothed[,j])))])),
          c(CIs(j,xi.est.smoothed,P.est.smoothed)[,3], rev(CIs(j,xi.est.smoothed,P.est.smoothed)[,4])),
          col=adjustcolor("darkolivegreen", alpha.f = 0.3), border=NA)
  legend(1920,-3,legend=expression({bar(ts)[t]}^w),col="darkolivegreen",lty="dashed",lwd=1.5,
         bty="n",xpd=TRUE,horiz=TRUE,cex=1.2)
  
  dev.off()
}

# 2x4 r_w, r_i, R-pi
{
  FILE="/../Plots/Switz/All_sigma/real_countries_DN.pdf"
  pdf(file=paste(getwd(),FILE,sep=""),pointsize=10, width=6, height=4)
  par(mfrow=c(2,4),oma = c(4, 0, 0, 0),mar = c(2,2,1,1))
  for (j in 1:nb){
    plot(data$year,r_world$rw,
         xaxt="non",type="l",xlab="",ylab="",lty=2,lwd=1.5,ylim=c(-5,10),xaxs="i")
    lines(data$year,r_i[,j+1],type="l",lty=1,col=colors[j])
    lines(data$year,real_data[,j+1],type="l",lty=3,col=colors[j])
    axis(1,seq(1880,2020,20))
    abline(h=0,col=adjustcolor(("black"),alpha.f=.3))
  }
  par(fig = c(0, 1, 0, 1), oma = c(2, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
  legend("bottom",legend = c(countries_DN), bty = 'n',
         col=c(colors), lwd=1.5, cex=1, horiz = TRUE,xpd=TRUE)
  dev.off()
}

# 2x4 pi_w, pi_i,pi
{
  FILE="/../Plots/Switz/All_sigma/pi_countries_DN.pdf"
  pdf(file=paste(getwd(),FILE,sep=""),pointsize=10, width=6, height=4)
  par(mfrow=c(2,4),oma = c(4, 0, 0, 0),mar = c(2,2,1,1))
  # r_w
  for (j in 1:nb){
    plot(data$year,pi_world$piw,
         xaxt="non",type="l",xlab="",ylab="",lty=2,lwd=1.5,ylim=c(-3,15),xaxs="i")
    lines(data$year,pi_i[,j+1],type="l",lty=1,col=colors[j])
    lines(data$year,data[,nb+nb+j+1],type="l",lty=3,col=colors[j])
    axis(1,seq(1880,2020,20))
    abline(h=0,col=adjustcolor(("black"),alpha.f=.3))
  }
  par(fig = c(0, 1, 0, 1), oma = c(2, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
  legend("bottom",legend = c(countries_DN), bty = 'n',
         col=c(colors), lwd=1.5, cex=1, horiz = TRUE,xpd=TRUE)
  dev.off()
}

# 2x4  ts_w, ts_i,RL-R
{
  FILE="/../Plots/Switz/All_sigma/ts_countries_DN.pdf"
  pdf(file=paste(getwd(),FILE,sep=""),pointsize=10, width=6, height=4)
  par(mfrow=c(2,4),oma = c(6, 0, 0, 0),mar = c(2,2,1,1))
  # ts_w
  for (j in 1:nb){
    plot(data$year,ts_world$tsw,
         xaxt="non",type="l",xlab="",ylab="",lty=2,lwd=1.5,ylim=c(-2,4),xaxs="i")
    lines(data$year,ts_i[,j+1],type="l",lty=1,col=colors[j])
    lines(data$year,ts_data[,j+1],type="l",lty=3,col=colors[j])
    axis(1,seq(1880,2020,20))
    abline(h=0,col=adjustcolor(("black"),alpha.f=.3))
  }
  par(fig = c(0, 1, 0, 1), oma = c(2, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
  legend("bottom",legend = c(countries_DN), bty = 'n',
         col=c(colors), lwd=1.5, cex=1, horiz = TRUE,xpd=TRUE)
  dev.off()
}

# r_i
{
  FILE="/../Plots/Switz/All_sigma/r_i_r_i.pdf"
  pdf(file=paste(getwd(),FILE,sep=""),pointsize=10, width=6, height=4)
  par(mfrow=c(2,4),oma = c(4, 0, 0, 0),mar = c(2,2,1,1))
  for (i in 1:nb) {
    plot(r_world$year,r_i[,i+1],xaxt="non",type="l",xlab="",ylab="",
         lty=2,lwd=1,ylim=c(-2,6),xaxs="i",col=colors[i])
    lines(r_world$year,xi.est.smoothed[,i+1],col=colors[i])
    axis(1,seq(1880,2000,20))
    abline(h=0,col=adjustcolor(("black"),alpha.f=.3))
    abline(v=1980,col=adjustcolor(("black"),alpha.f=.3))
  }
  par(fig = c(0, 1, 0, 1), oma = c(2, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
  legend(-1,-0.97,legend = c(countries_DN), bty = 'n',
         col=c(colors), lwd=1.6, cex=1.3, horiz = TRUE,xpd=TRUE)
  dev.off()
}


dev.off()
