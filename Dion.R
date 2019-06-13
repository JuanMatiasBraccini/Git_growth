#bring in packages
library(fishmethods)
library(FSA)

#Set working directory
setwd('C:\\Matias\\Analyses\\Dion')

#Read in data
dat=read.csv('Data.csv',stringsAsFactors = F)

#Define vector of regions
Reg=c("Kimberley","Pilbara","GC")
n.reg=length(Reg)


#initial par values
Linf=379
k=0.34
to=-0.45


#First. Select initial parameter values
fn.see.pin=function(DAT,Linf,k,to)
{
  svTall <- vbStarts(TL~decAge,data=DAT,type="typical") 
  plot(DAT$decAge,DAT$TL,ylab="TL",xlab="Age")
  AgE=sort(unique(DAT$decAge))
  lines(AgE,svTall$Linf*(1-exp(-svTall$K*(AgE-svTall$t0))),col=2,lwd=3)
  lines(AgE,Linf*(1-exp(-k*(AgE- to))),col=3,lwd=3)
  legend("topleft",c("model.based","eyeball"),lty=1,col=2:3,lwd=3,bty='n')
  legend("bottomright",paste(unique(DAT$Region)),bty='n')
  return(svTall)
}
Store.init.par=vector('list',n.reg)
names(Store.init.par)=Reg
par(mfcol=c(3,1),mar=c(3,4,.1,.1))
for(r in 1:n.reg)Store.init.par[[r]]=fn.see.pin(DAT=subset(dat,Region==Reg[r]),
                                                Linf=Linf,k=k,to=to)
PARS=Store.init.par
PARS$Pilbara=list(Linf=Linf,K=k,t0=to) #model based intial pars not good, used eyeballed


#Second. Fit model to data and predict new data
  #Generic function for estimating nls parameters
fn.fit.mod=function(Data,Pars,Formula,newD)
{
  mod <- nls(Formula, data=Data,start=Pars,
             nls.control(maxiter = 100, tol = 1e-06, minFactor = 1/1024,printEval = FALSE, warnOnly = T))
  PRED=predict(mod,newdata=newD)
  TABL=as.data.frame(summary(mod)$coefficients[,c(1:length(Pars),4)])
  return(list(mod=mod,ANOVA=TABL, PREDS=cbind(newD,PRED)))
}

  #function for prediction Confidence intervals
nlsint <-function(object, newdata = eval(object$call$data), interval = c("confidence", "prediction"), level = 0.95) 
{
  library(numDeriv)
  if (class(object) != "nls") 
    stop("nls objects only")
  type <- match.arg(interval)
  p.names <- names(coef(object))
  x.names <- names(object$dataClasses)
  for (i in 1:length(x.names)) {
    assign(x.names[i], newdata[[x.names[i]]])
  }
  f <- function(theta) {
    for (i in 1:length(p.names)) {
      assign(p.names[i], theta[i])
    }
    eval(parse(text = as.character(formula(object))[3]))
  }
  dmat <- jacobian(f, coef(object))
  vmat <- vcov(object)
  df <- summary(object)$df[2]
  yh <- f(coef(object))
  va <- diag(dmat %*% vmat %*% t(dmat))
  se <- switch(type, confidence = sqrt(va), prediction = sqrt(va + summary(object)$sigma^2))
  lw <- yh - qt(level + (1 - level)/2, df) * se
  up <- yh + qt(level + (1 - level)/2, df) * se
  out <- data.frame(fit = yh, se = se, lwr = lw, upr = up)
  rownames(out) <- NULL
  return(out)
}

  #Von B Formula
FORMULA=as.formula("TL~Linf*(1-exp(-K*(decAge-t0)))")

  #wrapper function
fun.fit=function(Data,PARs)
{
  Pred.data=data.frame(decAge=floor(seq(min(Data$decAge),max(Data$decAge),1)))
  
  Growth=fn.fit.mod(Data=Data,Pars=PARs,Formula=FORMULA,newD=Pred.data)
  Growth.CIs=nlsint(object=Growth$mod, newdata = Pred.data)
  return(list(fit=Growth,CI=Growth.CIs,Data=Data))
}
Store.out=Store.init.par
for(r in 1:n.reg)Store.out[[r]]=fun.fit(Data=subset(dat,Region==Reg[r]),PARs=PARS[[r]])


#Third. Output stuff
  #Anovas
for(r in 1:n.reg)  write.csv(Store.out[[r]]$fit$ANOVA,paste("ANOVA_",Reg[r],".csv",sep=""))

  #Plots
tiff(file=paste("Compare.models.tiff",sep=""),width = 1400, height = 2400,units = "px", res = 300, compression = "lzw")
par(mfcol=c(3,1),mar=c(2.5,3.5,.1,.1),oma=c(1.25,1.25,.1,.1),las=1,mgp=c(1,.75,0),cex.axis=1.25)
for(r in 1:n.reg) 
{
  with(Store.out[[r]]$Data,plot(decAge,TL,ylab="",xlab="",ylim=c(0,600),xlim=c(0,20),pch=21,cex=1.25,bg="grey70"))
  with(Store.out[[r]]$fit$PREDS,lines(decAge,PRED,lwd=1.5))
  lines(Store.out[[r]]$fit$PREDS$decAge,Store.out[[r]]$CI$lwr,col="red",lwd=1.5)
  lines(Store.out[[r]]$fit$PREDS$decAge,Store.out[[r]]$CI$upr,col="red",lwd=1.5)
  legend("bottomright",paste(Reg[r]),cex=1.5,bty='n')
}
mtext("Age",1,outer=T,line=-.5,cex=1.5)
mtext("Total length (mm)",2,outer=T,las=3,line=-.5,cex=1.5)
dev.off()  


#Fourth. Compare curves thru Likelihood ratio tests
  #source: Kimura (1980)
#error the error variance assumption: 
#           1= constant variance for all lijs;
#           2= constant variance for all mean lengths at age;
#           3=var of lij varies with age. See methods a-c in Kimura (1980: pp. 766).
LRT=with(subset(dat,!Region=="GC"),vblrt(len=TL,age=decAge,group=Region,
          error=1,select=1,plottype=1))
  
LRT$results
LRT$`model Ho`  #the general model
LRT$`model H1`  #model H1 (i.e. same Linf)
LRT$`model H2`  #model H2 (i.e. same k)
LRT$`model H3`  #model H3 (i.e. same to)
LRT$`model H4`  #model H4 (i.e. same Linf, k and to)
LRT$rss
head(LRT$residuals)


#Add GC to comparison: note, it doesn't converge, cannot find gradient...
LRT1=with(dat,vblrt(len=TL,age=decAge,group=Region,
                    error=1,select=2,plottype=1,
                    Linf = c(PARS$Kimberley$Linf,PARS$Pilbara$Linf,PARS$GC$Linf),
                    K = c(PARS$Kimberley$K,PARS$Pilbara$K,PARS$GC$K),
                    t0 = c(PARS$Kimberley$t0,PARS$Pilbara$t0,PARS$GC$t0)))





