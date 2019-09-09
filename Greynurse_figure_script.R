library(plotrix)
library(mvtnorm)

#Data
setwd('C:\\Matias\\Analyses\\Growth\\Greynurse recapture')
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