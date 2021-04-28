library(plotrix)
library(mvtnorm)
library(dplyr)

#Data
if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')
setwd(handl_OneDrive('Analyses\\Growth\\Greynurse recapture'))
  #tag and recapture info
TL.rel=188
TL.rec=246
Time.liberty=22.16  #years

  #male growth pars from Goldman et al 2006
Linf=249.5
K=0.16
to=-3.4

  #eye-balled data from Goldman et al 2006
rough.dat=data.frame(age=c(0,0,1,1,2,2,3,4,4,4,4,5,5,6,7,7,8,8,9,10,11,11,11,
                           12,12,12,12,13,13,14,15,16),
                     TL=c(100,110,125,130,140,145,150,165,167,170,193,165,202,190,
                          200,215,205,215,210,220,225,230,232,205,220,223,245,214,
                          230,227,245,245))


#1. create var covar matrix and take sample
set.seed(666)
Start.pars=list(Linf.est=jitter(Linf,1),K.est=jitter(K,1),to.est=jitter(to,1))  
fit=nls(TL~Linf.est*(1-exp(-K.est*(age-to.est))),data=rough.dat,start=Start.pars)
SIGMA=vcov(fit)
N.sim=1000
age=seq(0,40,.1)
Sim.tl=matrix(nrow=length(age),ncol=N.sim)

VBGF=function(t,linf,k,t0) linf*(1-exp(-k*(t-t0)))

x <- 1
repeat {
  par.samp=rmvnorm(1,mean=c(Linf,K,to),sigma=SIGMA)
  dummy=VBGF(age,par.samp[1],par.samp[2],par.samp[3])
  Age.dummy=age[which.min(abs(dummy - TL.rec))]
  if( TL.rec<max(dummy))
  {
    Sim.tl[,x]=dummy
    x = x+1
  }
  if (x == N.sim+1){
    break
  }
}  #make sure rec length/age is meaningful

#plot simulated lengths
plot(rough.dat$age,rough.dat$TL,col='transparent',xlim=c(0,40),ylim=c(0,300))
LT=VBGF(age,Linf,K,to)
for(i in 1:N.sim) lines(age,Sim.tl[,i],col=rgb(.1,.5,.1,alpha=.1))
lines(age,LT,lwd=2,col=2)
points(rough.dat$age,rough.dat$TL,pch=19)
abline(h=TL.rel,col=4)
abline(h=TL.rec,col=4)


#extract age from closet predicted TL to TL.rel
Age.rel=numeric(N.sim)
for(i in 1:N.sim) Age.rel[i]=age[which.min(abs(Sim.tl[,i] - TL.rel))] 
Smry.age.rel=round(as.data.frame(t(as.matrix(summary(Age.rel)))),1)

  
#derive age at recapture from time at liberty
Age.rec=Age.rel+Time.liberty
Smry.age.rec=round(as.data.frame(t(as.matrix(summary(Age.rec)))),1)

#derive age at recapture from predicted TL to TL.rec
Age.rec.TL=numeric(N.sim)
for(i in 1:N.sim) Age.rec.TL[i]=age[which.min(abs(Sim.tl[,i] - TL.rec))] 
Age.rec.TL=Age.rec.TL[Age.rec.TL>=(0.5*min(Age.rel)+Time.liberty)]

Smry.age.rec.TL=round(as.data.frame(t(as.matrix(summary(Age.rec.TL)))),1)

#figure 2
add.Tabl=F
Ymax=300
Xmax=max(age)

jpeg(file='figure2.jpg',width=2400,height=2000,units="px",res=300)
par(mar=c(1,1,1,.1),oma=c(2,3,.1,.5),las=1,mgp=c(1,.5,0),cex.axis=1.25)
plot(1,1,col='transparent',xlim=c(0,Xmax),ylim=c(0,Ymax), yaxs="i",xaxs="i",
     ylab="",xlab="")
for(i in 1:N.sim) lines(age,Sim.tl[,i],col=rgb(.1,.1,.3,alpha=.02))
mtext("Total length (cm)",2,las=3,line=2.5,cex=1.5)
mtext("Age",1,line=1.5,cex=1.5)

#add distribution
CL.dis=rgb(.1,.1,.1,alpha=.45)
den=density(Age.rel,adjust=2)
X=den$x
Y=rescale(den$y,c(TL.rel,TL.rel*1.2))
polygon(X,Y,col=CL.dis,border=CL.dis)
if(add.Tabl)addtable2plot(Smry.age.rel$Median ,TL.rel-20,Smry.age.rel,bty="o",xjust=0.2,
              display.rownames=F,hlines=F,vlines=F,cex=1)
den=density(Age.rec.TL,adjust=2)
X=den$x
Y=rescale(den$y,c(TL.rec,TL.rec*1.2))
polygon(X,Y,col=CL.dis,border=CL.dis)
if(add.Tabl)addtable2plot(Smry.age.rec.TL$Median ,TL.rec-20,Smry.age.rec.TL,bty="o",xjust=0.4,
              display.rownames=F,hlines=F,vlines=F,cex=1)
lines(age,LT,lwd=2,col="black")
abline(h=TL.rel,lwd=1.5,lty=2)
text(-.5,TL.rel+5,'Release',pos=4)
abline(h=TL.rec,lwd=1.5,lty=2)
text(-.5,TL.rec+5,'Recapture',pos=4)
dev.off()


# Quants=apply(Sim.tl, 1, quantile, probs = c(0,0.025,0.975,1),  na.rm = TRUE)
# jpeg(file='figure2.jpg',width=2400,height=2000,units="px",res=300)
# par(mar=c(1,1,1,.1),oma=c(2,3,.1,.5),las=1,mgp=c(1,.5,0),cex.axis=1.25)
# plot(1,1,col='transparent',xlim=c(0,Xmax),ylim=c(0,Ymax), yaxs="i",xaxs="i",
#      ylab="",xlab="")
# CL.dis=rgb(.1,.1,.3,alpha=.5)
# X=c(age,rev(age))
# Y=c(Quants[1,],rev(Quants[4,]))
# polygon(X,Y,col=CL.dis,border=CL.dis)
# 
# CL.dis=rgb(.1,.1,.3,alpha=.2)
# Y=c(Quants[2,],rev(Quants[3,]))
# polygon(X,Y,col=CL.dis,border=CL.dis)
# 
# 
# mtext("Total length (cm)",2,las=3,line=2.5,cex=1.5)
# mtext("Age",1,line=1.5,cex=1.5)
# 
# #add distribution
# CL.dis=rgb(.1,.1,.3,alpha=.5)
# den=density(Age.rel,adjust=2)
# X=den$x
# Y=rescale(den$y,c(TL.rel,TL.rel*1.2))
# polygon(X,Y,col=CL.dis,border=CL.dis)
# den=density(Age.rec.TL,adjust=2)
# X=den$x
# Y=rescale(den$y,c(TL.rec,TL.rec*1.2))
# polygon(X,Y,col=CL.dis,border=CL.dis)
# lines(age,LT,lwd=2,col="black")
# abline(h=TL.rel,lwd=1.5,lty=2)
# text(-.5,TL.rel+5,'Release',pos=4)
# abline(h=TL.rec,lwd=1.5,lty=2)
# text(-.5,TL.rec+5,'Recapture',pos=4)
# dev.off()



#-------------Figure 1------------------
library(PBSmapping)
library(lubridate)
library(geosphere)
library(raster) #for scalebar
SMN=read.csv("Acous.tags_SMN.csv",stringsAsFactors = F, sep = ";")
PerthIs=read.table("WAislandsPointsNew.txt", header=T)
Rottnest.Is=subset(PerthIs,ID%in%c("ROTT1"))
Garden.Is=subset(PerthIs,ID%in%c("ROTT3"))
STATIONS=read.csv('Receivers.csv',stringsAsFactors = F)
#For SCALER
# fuür SMN
d=read.csv("Acous.tags_SMN.csv", sep = ";")

SMN=SMN%>%mutate(DateTime.local=as.POSIXct(DateTime.local,format="%d/%m/%y %H:%M"),
                 Date.local=as.Date(SMN$Date.local,format="%d/%m/%y"),
                 Time.local=hms(Time.local),
                 Year=year(Date.local),
                 Month=month(Date.local),
                 Hour=hour(Time.local),
                 Mins=minute(Time.local))%>%
                arrange(DateTime.local)

d=d%>%mutate(DateTime.local=as.POSIXct(DateTime.local,format="%d/%m/%y %H:%M"),
             Date.local=as.Date(SMN$Date.local,format="%d/%m/%y"),
             Time.local=hms(Time.local),
             Year=year(Date.local),
             Month=month(Date.local),
             Hour=hour(Time.local),
             Mins=minute(Time.local))%>%
            arrange(DateTime.local)


d$ReleaseLatitude=-abs(d$ReleaseLatitude)
d$Latitude=-abs(d$Latitude)
d=subset(SMN,Sex=="F")
k=subset(SMN,Sex=="M")

REls=SMN[!duplicated(SMN$TagCode),match(c("TagCode","ReleaseLatitude","ReleaseLongitude"),names(SMN))]
REls=rbind(REls,data.frame(TagCode=29445,ReleaseLatitude=34.351,ReleaseLongitude=115.22))
REls$ReleaseLatitude=-abs(REls$ReleaseLatitude)
m=SMN
#female
d$Station=with(d,paste(Longitude,Latitude))
Tab=table(d$Station)
Tab=data.frame(Tab)
names(Tab)=c('Station',"detections")
Tab1=data.frame(do.call('rbind', strsplit(as.character(Tab$Station),' ',fixed=TRUE)))
names(Tab1)=c("Long","Lat")
Tab=cbind(Tab,Tab1)
Tab$Long=as.numeric(as.character(Tab$Long))
Tab$Lat=as.numeric(as.character(Tab$Lat))
#male
k$Station=with(k,paste(Longitude,Latitude))
Tab3=table(k$Station)
Tab3=data.frame(Tab3)
names(Tab3)=c('Station',"detections")
Tab4=data.frame(do.call('rbind', strsplit(as.character(Tab3$Station),' ',fixed=TRUE)))
names(Tab4)=c("Long","Lat")
Tab3=cbind(Tab3,Tab4)
Tab3$Long=as.numeric(as.character(Tab3$Long))
Tab3$Lat=as.numeric(as.character(Tab3$Lat))
#for scale 
m$Station=with(m,paste(Longitude,Latitude))
Tab5=table(m$Station)
Tab5=data.frame(Tab5)
names(Tab5)=c('Station',"detections")
Tab6=data.frame(do.call('rbind', strsplit(as.character(Tab5$Station),' ',fixed=TRUE)))
names(Tab6)=c("Long","Lat")
Tab5=cbind(Tab5,Tab6)
Tab5$Long=as.numeric(as.character(Tab5$Long))
Tab5$Lat=as.numeric(as.character(Tab5$Lat))

#fn.scale=function(x,scaler) ((x/max(x))^0.5)*scaler
fn.scale=function(x,scaler) scaler*log(x)/max(log(x))
x=c(10,100,250,500,1000)
scaler=5
#CEX=mapply(fn.scale,x,max(x),scaler)
#CEX=3
#plot(1:length(x),cex=CEX)

#MAPPING


fn.plt.mp=function(PlotlonG,PlotlatT,Add.depth,add.closure)
{
  plotMap(worldLLhigh, xlim=PlotlonG,ylim=PlotlatT,plt = c(.001, 1, 0.075, 1),
          col=COLOR,tck = 0.025, tckMinor = 0.0125, xlab="",ylab="",axes=F)
  #if(add.closure=="YES") plot(WA_Northern_Shark_2,add=T,col="grey85",border=1)
  if(Add.depth=="YES")
  {
    contour(xbat, ybat, reshaped[,2:ncol(reshaped)],ylim=YLIM,xlim=XLIM, zlim=c(-1,-200),
            nlevels = 3,labcex=.8,lty = 1,col=c("darkgrey","darkgrey","darkgrey","darkgrey","transparent"),add=T)
  }
  
}
fn.poli=function(XLIM,YLIM,CL,LW)
{
  polygon(x=c(XLIM,rev(XLIM)),y=c(YLIM[1],YLIM[1],YLIM[2],YLIM[2]),lwd=LW,border=CL)
}
fnX=function(x=0.9)points(STATIONS$longitude,STATIONS$latitude,pch=21,bg="white",col="grey30",cex=x)

data(worldLLhigh)
COLOR="grey75"

tiff(file="Figure 1.tiff",width = 3200, height = 3000,   
     units = "px", res = 300,compression = "lzw")
par(mar = c(3.2, 5, 0.95, 0),oma=c(1,.1,.1,1),mgp=c(1,.75,0))
layout(matrix(cbind(rep(c(3,1,2),1),c(4,4,4)),3,2))
layout.show(4)

# Perth array 
EE1=c(115.17,116)
FF1=c(-32.3,-31.8)
LONGSEQ=EE1[1]:EE1[2]
LATSEQ=seq(ceiling(FF1[1]),ceiling(FF1[2]),by=1)
LATSEQ2=seq(ceiling(FF1[1]),ceiling(FF1[2]),by=2)
fn.plt.mp(EE1,FF1,"NO",'YES')
box(lwd=2,col="forestgreen")
polygon(x=Rottnest.Is$Longitude,y=Rottnest.Is$Latitude,col="grey")  #add missing islands
polygon(x=Garden.Is$Longitude,y=Garden.Is$Latitude,col="grey")
fnX()
points(Tab5$Long,Tab5$Lat, pch=21,cex=fn.scale(Tab5$detections,scaler=10),col="pink",
       bg=rgb(.5,.1,.1,alpha=.2))
points(Tab3$Long,Tab3$Lat, pch=21,cex=fn.scale(Tab5$detections,scaler=10),col="blue",
       bg=rgb(.1,.1,.5,alpha=.2))

text(115.75,-32.09,"Cockburn",cex=1.5,col='grey10',font=2,pos=4)
text(115.8,-32.15,"sound",cex=1.5,col='grey10',font=2,pos=4)
legend('topright',"2",bty="n",cex=3)

#Southern Array
EE=c(114.5305,116.6)
FF=c(-35.39,-34)
LONGSEQ=EE[1]:EE[2]
LATSEQ=seq(ceiling(FF[1]),ceiling(FF[2]),by=1)
LATSEQ2=seq(ceiling(FF[1]),ceiling(FF[2]),by=2)
fn.plt.mp(EE,FF,"NO",'YES')
box(lwd=2,col="forestgreen")
fnX()
points(Tab$Long,Tab$Lat, pch=21,cex=fn.scale(Tab5$detections,scaler=10),col="pink",
       bg=rgb(.5,.1,.1,alpha=.2))
points(Tab3$Long,Tab3$Lat, pch=21,cex=fn.scale(Tab5$detections,scaler=10),col="blue",
       bg=rgb(.1,.1,.5,alpha=.2))
legend('topright',"3",bty="n",cex=3)

#Ningaloo array
CC=c(113.5, 114.4)
DD=c(-23.29, -21.7)
LONGSEQ=CC[1]:CC[2]
LATSEQ=seq(ceiling(DD[1]),ceiling(DD[2]),by=1)
LATSEQ2=seq(ceiling(DD[1]),ceiling(DD[2]),by=2)
fn.plt.mp(CC,DD,"NO",'YES')
box(lwd=2,col="forestgreen")
fnX()
points(Tab$Long,Tab$Lat, pch=21,cex=fn.scale(Tab5$detections,scaler=10),col="pink",
       bg=rgb(.5,.1,.1,alpha=.2))
points(Tab3$Long,Tab3$Lat, pch=21,cex=fn.scale(Tab5$detections,scaler=10),col="blue",
       bg=rgb(.1,.1,.5,alpha=.2))
legend('topright',"1",bty="n",cex=3)



# Whole WA
XX=c(112,119)
YY=c(-35.5,-21.5)
LONGSEQ=XX[1]:XX[2]
LATSEQ=seq(ceiling(YY[1]),ceiling(YY[2]),by=1)
LATSEQ2=seq(ceiling(YY[1]),ceiling(YY[2]),by=2)
fn.plt.mp(XX,YY,"NO",'YES')
axis(side = 1, at =LONGSEQ, labels = paste(LONGSEQ,"ºE",sep=""), tcl = .5,las=1,cex.axis=1.5,padj=-.5)
axis(side = 4, at = LATSEQ, labels =F,tcl = .5,las=2,cex.axis=1.5,hadj=.15)
axis(side = 4, at = LATSEQ2, labels =paste(-LATSEQ2,"ºS",sep=""),tcl = .5,las=2,cex.axis=1.5,hadj=.15)
box(lwd=1.5)
polygon(x=Rottnest.Is$Longitude,y=Rottnest.Is$Latitude,col="grey")  #add missing islands
polygon(x=Garden.Is$Longitude,y=Garden.Is$Latitude,col="grey")
text(117,-31.5,"Perth array",cex=2,col='grey10',font=2)
text(117.5,-34,"Southern",cex=2,col='grey10',font=2)
text(117.5,-34.5,"array",cex=2,col='grey10',font=2)
text(115.5,-22.5,"Ningaloo",cex=2,col='grey10',font=2)
text(115.5,-23,"array",cex=2,col='grey10',font=2)
fn.poli(EE,FF,'darkgreen',2)
text(113.1268,-21.91353,"1",cex=3,font=2)
fn.poli(EE1,FF1,'darkgreen',2)
text(114.8443,-32.04497,"2",cex=3,font=2)
fn.poli(CC, DD,'darkgreen',2)
text(114.2326,-34.42276,"3",cex=3,font=2)
fnX(x=.7)
Shark.bay=cbind(113,-25.497)
Cape.Leeuwin=cbind(115.18,-34.17)

points(Cape.Leeuwin[1],Cape.Leeuwin[2],pch=21,col=1,bg=rgb(.1,.9,.3,alpha=.5),cex=3.5)
points(Shark.bay[1],Shark.bay[2],pch=21,col=1,bg=rgb(.1,.9,.3,alpha=.5),cex=3.5)
#legend(115.5,-26,c("Releases","Receiver"),pch= c(4,1),col=c("black","grey40"),cex=2.2,bty='n')
LEG=round(quantile(Tab5$detections,probs=c(.15,.5,.95)))
scalebar(100,cex=1.5)
legend("right",paste(LEG),pch=21,col="black",cex = 2, bg=rgb(.1,.1,.1,alpha=.2),
       pt.cex=fn.scale(LEG,scaler=10),bty='n',
       y.intersp = c(1,1,1.5), x.intersp = 1.3)
text(117.5,-27.5,"No. of detections",cex=2)
points(REls$ReleaseLongitude,REls$ReleaseLatitude,pch=24,cex=2,bg="orange")

dev.off()



#Stand alone Perth array
tiff(file="Figure 1_appendix_Perth.Array.tiff",width = 3000, height = 2000,   
     units = "px", res = 300,compression = "lzw")
par(mar = c(3.2, 5, 0.95, 1),oma=c(1,.1,.1,1),mgp=c(1,.75,0))
EE1=c(115.17,116)
FF1=c(-32.3,-31.8)
LONGSEQ=EE1[1]:EE1[2]
LATSEQ=seq(ceiling(FF1[1]),ceiling(FF1[2]),by=1)
LATSEQ2=seq(ceiling(FF1[1]),ceiling(FF1[2]),by=2)
fn.plt.mp(EE1,FF1,"NO",'YES')
box(lwd=2,col="black")
polygon(x=Rottnest.Is$Longitude,y=Rottnest.Is$Latitude,col="grey")  #add missing islands
polygon(x=Garden.Is$Longitude,y=Garden.Is$Latitude,col="grey")
fnX()
points(Tab5$Long,Tab5$Lat, pch=21,cex=fn.scale(Tab5$detections,scaler=10),col="pink",
       bg=rgb(.5,.1,.1,alpha=.2))
points(Tab3$Long,Tab3$Lat, pch=21,cex=fn.scale(Tab5$detections,scaler=10),col="blue",
       bg=rgb(.1,.1,.5,alpha=.2))

text(115.82,-32.09,"Cockburn",cex=1.5,col='grey10',font=2)
text(115.82,-32.12,"sound",cex=1.5,col='grey10',font=2)
legend('topright',"2",bty="n",cex=3)
LEG=round(quantile(Tab5$detections,probs=c(.15,.5,.95)))
scalebar(25,cex=1.5)
legend("topleft",paste(LEG),pch=21,col="black",cex = 2, bg=rgb(.1,.1,.1,alpha=.2),
       pt.cex=fn.scale(LEG,scaler=10),bty='n',title='No. of detections',
       y.intersp = c(1,1,1.25), x.intersp = 1.3)
axis(side = 1, at =c(115.35,115.7), labels = paste(c(115.35,115.7),"ºE",sep=""), tcl = .5,las=1,cex.axis=1.5,padj=-.5)
axis(side = 4, at = c(-32,-32.2), labels =paste(-c(-32,-32.2),"ºS",sep=""),tcl = .5,las=2,cex.axis=1.5,hadj=.15)
dev.off()




#Figure S3
#Stand alone Perth array
EE1=c(115.17,116)
FF1=c(-32.3,-31.8)
LONGSEQ=EE1[1]:EE1[2]
LATSEQ=seq(ceiling(FF1[1]),ceiling(FF1[2]),by=1)
LATSEQ2=seq(ceiling(FF1[1]),ceiling(FF1[2]),by=2)

This.mn=c(10:12,1:5)
This.yr=c(rep(2014,3),rep(2015,5))
fn.scale=function(x,scaler) scaler*log(x)/max(log(Tab5$detections))
#Tab1=table(m$Station)

tiff(file="Figure 3_appendix.tiff",width = 1400, height = 2000,   
     units = "px", res = 300,compression = "lzw")
par(mfcol=c(4,2),mar = c(3.2, 5, 0.95, 2),oma=c(1,4,.1,1),mgp=c(1,.75,0))
for(i in 1:length(This.mn))
{
  xx=subset(m,Year==This.yr[i] & Month==This.mn[i])
  Tab5=table(xx$Station)
  Tab5=data.frame(Tab5)
  names(Tab5)=c('Station',"detections")
  Tab6=data.frame(do.call('rbind', strsplit(as.character(Tab5$Station),' ',fixed=TRUE)))
  names(Tab6)=c("Long","Lat")
  Tab5=cbind(Tab5,Tab6)
  Tab5$Long=as.numeric(as.character(Tab5$Long))
  Tab5$Lat=as.numeric(as.character(Tab5$Lat))
  
  fn.plt.mp(EE1,FF1,"NO",'YES')
  box(lwd=2,col="black")
  polygon(x=Rottnest.Is$Longitude,y=Rottnest.Is$Latitude,col="grey")  #add missing islands
  polygon(x=Garden.Is$Longitude,y=Garden.Is$Latitude,col="grey")
  fnX()
  points(Tab5$Long,Tab5$Lat, pch=21,cex=fn.scale(Tab5$detections,scaler=10),col="pink",
         bg=rgb(.5,.1,.1,alpha=.2))
  legend('bottomleft',paste(as.character(month(xx$Month[1],label = TRUE, abbr = FALSE)),
                             This.yr[i]),bty="n",cex=1.1)
  LEG=round(quantile(Tab5$detections,probs=c(.5)))
  if(i==8)scalebar(25,xy=c(115.75,-32.26),cex=1.5)
  legend("topleft",paste(LEG),pch=21,col="black",cex = .9, bg=rgb(.1,.1,.1,alpha=.2),
         pt.cex=fn.scale(LEG,scaler=10),bty='n',title='No. of detections',
         y.intersp = c(1.8), x.intersp = c(-.5))
  axis(side = 1, at =c(115.35,115.7), labels = F, tcl = .5,las=1,cex.axis=1.5,padj=-.5)
  axis(side = 2, at = c(-32,-32.2), labels =F,tcl = .5,las=2,cex.axis=1.5,hadj=.8)
if(i%in%c(4,8))  axis(side = 1, at =c(115.35,115.7), labels = paste(c(115.35,115.7),"ºE",sep=""), tcl = .5,las=1,cex.axis=1.5,padj=-.5)
 if(i%in%1:4)  axis(side = 2, at = c(-32,-32.2), labels =paste(-c(-32,-32.2),"ºS",sep=""),tcl = .5,las=2,cex.axis=1.5,hadj=.8)
  
  
}
dev.off()


#Residency 
This.mn=c(10:12,1)
This.yr=c(rep(2014,3),rep(2015,1))

tiff(file="Figure 5_appendix.tiff",width = 2000, height = 2000,   
     units = "px", res = 300,compression = "lzw")
par(mfcol=c(2,2),mar = c(2, .5, 0.95, .5),oma=c(1,.1,.1,.1),mgp=c(1,.75,0))
for(i in 1:length(This.mn))
{
  xx=subset(SMN,Sex=="F" & Year==This.yr[i] & Month==This.mn[i])%>%
        mutate(day=day(Date.local),
               day.hour=day+(Hour/24))
    
  xx=xx[!duplicated(xx$day.hour),]
  plot(xx$day.hour,rep(1,nrow(xx)),xlim=c(0,24),xlab="",ylab="",yaxt='n',
       main=as.character(month(xx$Month[1],label = TRUE, abbr = FALSE)),
       pch=".",cex=2,col=2)
}
mtext("Day-hour",1,outer=T,cex=1.5)
dev.off()


#Max distance and time
xx=subset(SMN,Sex=="F")
xx=xx[-1,]%>%
        arrange(DateTime.local)%>%
        mutate(SerialNumber.prev=lag(SerialNumber,1),
               Latitude.prev=lag(Latitude,1),
               Longitude.prev=lag(Longitude,1),
               DateTime.local.prev=lag(DateTime.local,1),
               Distance=NA,
               Time=ifelse(!SerialNumber.prev==SerialNumber,DateTime.local-DateTime.local.prev,NA))
for(i in 1:nrow(xx))
{
  if(!xx$SerialNumber.prev[i]==xx$SerialNumber[i] & !is.na(xx$Latitude.prev[i]))
  {
    xx$Distance[i]=with(xx,distCosine(c(Longitude[i],Latitude[i]),
                                      c(Longitude.prev[i],Latitude.prev[i]))/1000)
  }

}
xx[which.max(xx$Distance),]


#test
##Daily movement 
This.mn=c(10:12,1)
This.yr=c(rep(2014,3),rep(2015,1))

tiff(file="Figure 4_appendix.tiff",width = 2000, height = 2000,   
     units = "px", res = 300,compression = "lzw")
par(mfcol=c(2,2),mar = c(2, 2, 0.95, .5),oma=c(1,2,.1,.1),mgp=c(1,.75,0),las=1)
for(i in 1:length(This.mn))
{
  xx=subset(SMN,Sex=="F" & Year==This.yr[i] & Month==This.mn[i])%>%
    arrange(DateTime.local)%>%
    mutate(SerialNumber.prev=lag(SerialNumber,1),
           same=ifelse(SerialNumber.prev==SerialNumber,"YES","NO"))%>%
    filter(same=="NO")

  hist(xx$Hour,breaks=24,col="grey",ylab="",xlab="",
       main=paste(as.character(month(xx$Month[1],label = TRUE, abbr = FALSE)),
                  This.yr[i]))
  box()
}
mtext("Frequency",2,outer=T,cex=1.5,las=3)
mtext("Hour",1,outer=T,cex=1.5)
dev.off()
