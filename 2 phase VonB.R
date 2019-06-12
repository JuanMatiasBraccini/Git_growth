t.sim=0:20
Linf=300
k=0.2
to=-1
h=.5
th=5

#Matias
Pars=c(Linf=Linf,k=k,h=h,th=th,to=to)
fun.2PVonB=function(PARS) PARS[1]*(1-exp(-PARS[2]*(1-PARS[3]/((t.sim-PARS[4])^2+1))*(t.sim-PARS[5])))
Lt=fun.2PVonB(Pars)

#Marina's
Loo=Linf
E=t.sim
Lt_marina= Loo*(1-exp(-k* (1-(h/(((E-th)^2)+1)))  *(E-to)))


#compare
plot(t.sim,Lt,ylim=c(0,Linf*1.1))
lines(E,Lt_marina,col=2)


#Simulate data
JIT=200
t.sim=rep(seq(0:20),10)
Lt.sim=jitter(fun.2PVonB(Pars),JIT)
plot(t.sim,Lt.sim,ylim=c(0,Linf*1.1))


  #simulation evaluation
n=1000
Par.list=vector('list',n)
for(i in 1:n)
{
  Lt.sim=jitter(fun.2PVonB(Pars),JIT)
  FIT=nls(Lt.sim~Linf*(1-exp(-k*(1-h/((t.sim-th)^2+1))*(t.sim-to))),start=Pars)
  Par.list[[i]]=coef(FIT)
}

Par.list=do.call(rbind, Par.list)
Init.pars=matrix(Pars,nrow=nrow(Par.list),ncol=length(Pars),byrow=T)
boxplot(Par.list/Init.pars)
abline(h=1,col=2)
#El modelo que implementaste esta bien, deben ser tus datos
#fijate que en la simulacion se recupera el valor de los parameters usados
#para simular los datos. Si aumentas el jitter, ahi si tenes problemas de estimacion.
#Fijate que si cambias el valor de JIT de 200 a 1000, ahi ya empieza a ver problemas 
#con los parameter estimates