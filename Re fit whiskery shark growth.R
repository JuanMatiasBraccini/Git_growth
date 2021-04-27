# Re fit whiskery shark growth model
handl_OneDrive=function(x)paste('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias',x,sep='/')

DATA=read.csv(handl_OneDrive("Data/Age and growth/Whiskery.csv"))

#Add half a year to Stevens' data s done by simpfendorfer et al 2000
DATA$Counts=with(DATA,ifelse(Source=="Stevens 1990",Counts+0.5,Counts))

#Convert TL to FL
a.w=1.0044
b.w=13.171

TL_to_FL=function(TL,a,b) FL=(TL-b)/a
DATA$FL=with(DATA,ifelse(Source=="Stevens 1990",TL_to_FL(TL,a.w,b.w),FL))

Stevens.data=subset(DATA,Source=="Stevens 1990")
Simpfen.data=subset(DATA,Source=="Simpfendorfer et al 2000")

write.csv(Simpfen.data,handl_OneDrive("Data/Age and growth/Simpfen.data.csv"),row.names=F)

  #Size at birth
#Lo=33.5 #LT
Lo=22   #FL (Simpfendorfer et al 2000)
age=0:16
Linf=120
k=0.4

#calculate parameters

#optim approach
pars=c(k=log(k),Linf=log(Linf),SD=log(5))

fit.fn=function(pars)
{
  Par=exp(pars)
  FL=Lo+(Par[2]-Lo)*(1-exp(-Par[1]*dat$Counts))
  epsilon = dat$FL-FL  #residuals
  nloglike=-1.0*sum(dnorm(epsilon,0,Par[3]))
  return(list(nloglike=nloglike,FL=FL,epsilon=epsilon))
  
}
fn=function(pars) fit.fn(pars)$nloglike

dat=Simpfen.data
fit = optim(pars,fn,method="BFGS",hessian=T)
V=solve(fit$hessian)  #getting the variance covariance matrix
std=sqrt(diag(V))		#Standard deviation of parameters
R=V/(std%o%std)			#Parameter correlation, the V divided by the outer product of std
R=round(R,3)			#round R to 3 decimals

REs=fit.fn(fit$par)$epsilon
VAR=sum(REs^2)/(nrow(dat)-1)
SDev=(VAR)^.5

#nls approach
Von.B <- FL~Linf*(1-exp(-k*(DATA$Counts-to)))
Von.B.2par=FL~Lo+(Linf-Lo)*(1-exp(-k*DATA$Counts))  #modified version


  #standard VonB
Simple.Theta=c(Linf=Linf,k=k,to=-1)  #initial par values
fit.Von.B <- nls(Von.B,data=DATA,start=as.list(Simple.Theta))
summary(fit.Von.B)  

Fits=coef(fit.Von.B)
Fits.vonB=Fits
plot(Stevens.data$Counts,Stevens.data$FL,xlim=c(0,17),ylim=c(0,140),pch=19,col=2,xaxt='n')
points(Simpfen.data$Counts,Simpfen.data$FL,pch=19,col=3)
lines(age,Fits[1]*(1-exp(-Fits[2]*(age-Fits[3]))),lwd=2,col=4)
axis(1,age,age,cex.axis=0.9)


  #2 par von B
Simple.Theta=c(Linf=Linf,k=k)  #initial par values
fit.Von.B.2par <- nls(Von.B.2par,data=dat,start=as.list(Simple.Theta))
summary(fit.Von.B.2par)  

Fits=coef(fit.Von.B.2par)
Fits.vonB.2par=Fits
lines(age,Lo+(Fits[1]-Lo)*(1-exp(-Fits[2]*age)),lwd=2,col=1)

Fits.vonB
Fits.vonB.2par
