#Age and growth analysis of 

library(fishmethods)

#DATA SECTION
setwd("D:/Growth")
Dat=read.csv("Data.csv",stringsAsFactors=F)

Specs=unique(Dat$Species)


#PARAMETER SECTION


#PROCEDURE SECTION
  #Calculate growth parameters
fn.growth=function(sp)
{
  #select species of interest
  a=subset(Dat,Species==sp,select=c(Length,Age))
  
  #fit growth model
  model.fit<-nls(Length ~ Linf * (1 - exp(-K*(Age-t0))),
                 data=a, start = list (Linf = max(a$Length), K = 0.2, t0 = 0))

 return(model.fit)
}

LISTA=vector('list',length(Specs))
names(LISTA)=Specs
for(i in 1:length(Specs))
{
  LISTA[[i]]=fn.growth(sp=Specs[i])
}
  
  #Do likelihood ratio tests
fun.LRT=function(spa,spb)
{
  Kimura=subset(Dat,Species%in%c(spa,spb),select=c(Length,Age,Species))
  LRT=vblrt(len=Kimura$Length,age=Kimura$Age,group=Kimura$Species,error=2,select=1,plottype=1)
}

Spec.combo=combn(Specs, 2)

LRT.lis=vector('list',ncol(Spec.combo))
names(LRT.lis)=paste(Spec.combo[1,],Spec.combo[2,],sep="_vs_")
for(s in 1:ncol(Spec.combo))
{
  LRT.lis[[s]]=fun.LRT(spa=Spec.combo[1,s],spb=Spec.combo[2,s])
}




#REPORT SECTION
  #plots
fn.plt.curve=function(model,SPE,MAIN)
{
  a=subset(Dat,Species==SPE,select=c(Length,Age))
  NewData=data.frame(Age=seq(1,max(a$Age)))
  Preds=predict(model,newdata=NewData,type='response')
  SE=summary(model)$coefficients
  MLE=SE[,1]
  SE=SE[,2]
  
  Preds.upSE=(MLE[1]+1.96*SE[1]) * (1 - exp(-(MLE[2]+1.96*SE[2])*(NewData$Age-(MLE[3]+1.96*SE[3]))))
  Preds.lowSE=(MLE[1]-1.96*SE[1]) * (1 - exp(-(MLE[2]-1.96*SE[2])*(NewData$Age-(MLE[3]-1.96*SE[3]))))
  
  plot(NewData$Age,Preds,ylim=c(0,max(Preds.upSE)),main=MAIN)
  lines(NewData$Age,Preds.upSE)
  lines(NewData$Age,Preds.lowSE)

}
SP.nams=c("Eibli","Flavissima","Hybrids","Vrolikii")
for(i in 1:length(Specs))
{
  fn.plt.curve(model=LISTA[[i]],SPE=Specs[i],MAIN=SP.nams[i])
}


#Export LRT to excel
for(s in 1:ncol(Spec.combo))
{
  write.csv(LRT.lis[[s]]$results,paste("LRT_",names(LRT.lis[s]),".csv",sep=""),row.names=F)
}
