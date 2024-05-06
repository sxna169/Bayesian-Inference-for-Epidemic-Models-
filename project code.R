#Github code 

install.packages("MASS")
install.packages("ggplot2")
install.packages("coda")
install.packages("tidyverse")
install.packages("dplyr")
install.packages("gridExtra")
install.packages("mvtnorm")

library(coda)
library(MASS)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(gridExtra)
library(mvtnorm)

par(mfrow=c(1,1))

#Solving deterministic ODEs using eulers method

values<-function(theta,dt,t,X0){
  X<-matrix(0,nrow=(t/dt)+1,ncol=2)
  time<-seq(from=0,to=t,by=dt)
  N<-X0[1]+X0[2]
  X[1,]<-X0
  for (i in 1:(t/dt)){
    X[i+1,]<-c(X[i,1]+dt*(-theta[1]*(X[i,1])*(X[i,2])),
               X[i,2]+dt*((theta[1])*(X[i,1])*(X[i,2])-theta[2]*(X[i,2])))}
  X<-data.frame(X) 
  X$time<-time
  X$Rt<-763-X$X1-X$X2
  names(X)[1]<-"St"
  names(X)[2]<-"It"
  X<-relocate(X,St,It,Rt,time)
  return(X)
}

#Figure 2.2
ODEdata<-values(c(exp(-6),0.5),0.01,14,c(762,1))
ggplot()+geom_line(data=ODEdata,aes(x=time,y=St),color="red", linewidth=0.8)+geom_line(data=ODEdata,aes(x=time,y=It),color="green",linewidth=0.8)+geom_line(data=ODEdata,aes(x=time,y=Rt),color="blue",linewidth=0.8) + xlab("Time") + ylab("Number of individuals") + ggtitle("St, It, Rt simulated using the ODE model")

#Using Gillespies algorithm to generate simulations
#if the simlation shows that the epidemic dies early it gives us an error message and we need to resimulate the values  
GAvalues<-function(x,N,T,c,S){
  t<-0
  data<-data.frame(St=x[1],It=x[2],Rt=(N-x[1]-x[2]),t=0)
  while (t< T){
    h1<-c[1]*x[1]*x[2]
    h2<-c[2]*x[2]
    h0<-h1+h2
    dt<-rexp(1,rate=h0)
    t<-t+dt
    j<-sample(c(1,2), size=1, prob=c(h1/h0, h2/h0))
    x[1]<-x[1]+S[1,j]
    x[2]<-x[2]+S[2,j]
    data<-rbind(data,c(x[1],x[2],(N-x[1]-x[2]),t))
  }
  return(data) 
}

#Figure 2.3
GAdata<-GAvalues(c(762,1),763,14,c(exp(-6),0.5),matrix(c(-1,0,1,-1),nrow=2,ncol=2,byrow=TRUE))
ggplot()+geom_line(data=GAdata,aes(x=t,y=St),color="red", linewidth=0.8)+geom_line(data=GAdata,aes(x=t,y=It),color="green", linewidth=0.8)+geom_line(data=GAdata,aes(x=t,y=Rt),color="blue", linewidth=0.8)+ xlab("Time") +ylab("Number of individuals") +ggtitle("St, It, Rt generated using the Gillespie Algorithm")

#poisson leap process
PLval<-function(x,N,T,c,S,tau){
  t<-0
  data<-data.frame(St=x[1],It=x[2],Rt=(N-x[1]-x[2]),t=0)
  while (t< T-tau){
    h1<-c[1]*x[1]*x[2]
    h2<-c[2]*x[2]
    r1<-rpois(1,lambda = h1*tau)
    r2<-rpois(1,lambda = h2*tau)
    t<-t+tau
    x[1]<-x[1]-r1
    x[2]<-x[2]+r1-r2
    while(x[1]<0){
      x[1]<-0
    }
    while(x[2]<0){
      x[2]<-0
    }
    data<-rbind(data,c(x[1],x[2],(N-x[1]-x[2]),t))
  }
  return(data)
}

# plvalues returns poisson leap values without the possibility of having a simulation that goes wrong 
PLvalues<-function(x,N,T,c,S,tau){
  max<-1
  while(max<=10){
    PLval<-PLval(x,N,T,c,S,tau) 
    max<-max(PLval$It)
  }
  return(PLval)
}

#Figure 2.4
data1<-PLvalues(c(762,1),763,14,c(exp(-6),0.5),matrix(c(-1,0,1,-1),nrow=2,ncol=2,byrow=TRUE),0.01)
data2<-PLvalues(c(762,1),763,14,c(exp(-6),0.5),matrix(c(-1,0,1,-1),nrow=2,ncol=2,byrow=TRUE),0.1)
PLdata<-data1
ggplot()+geom_line(data=data1,aes(x=t,y=St),color="red", linewidth=0.8)+geom_line(data=data1,aes(x=t,y=It),color="green", linewidth=0.8)+geom_line(data=data1,aes(x=t,y=Rt),color="blue", linewidth=0.8)+geom_line(data=data2,aes(x=t,y=St),color="purple", linewidth=0.8)+geom_line(data=data2,aes(x=t,y=It),color="pink", linewidth=0.8)+geom_line(data=data2,aes(x=t,y=Rt),color="orange", linewidth=0.8)+geom_line(data=GAdata,aes(x=t,y=St),color="black", linewidth=0.8)+geom_line(data=GAdata,aes(x=t,y=It),color="black", linewidth=0.8)+geom_line(data=GAdata,aes(x=t,y=Rt),color="black", linewidth=0.8)+ xlab("Time") + ylab("Number of individuals") + ggtitle("St, It, Rt simulations using poisson leap with tau=0.1 and 0.01 and gillespies algorithm")

#function that generates X number of PL simulations of It and their mean

PLIt<- function(x,N,T,c,S,tau){
  data<-PLvalues(x,N,T,c,S,tau)
  return(data$It)
}

PLXdata<-function(x,N,T,c,S,tau,X){
  p<-PLIt(x,N,T,c,S,tau)
  n<-length(p)
  data<-matrix(rep(0,n*X),ncol=X,nrow=n)
  data<-data.frame(data)
  i=1
  for(i in 1:X){
    data[,i]<-PLIt(x,N,T,c,S,tau)
    i=i+1
  }
  p<-PLval(x,N,T,c,S,tau)
  time<-p$t
  data$mean<-rowMeans(data)
  data$time<-time
  return(data)
}
data<-PLXdata(c(762,1),763,14,c(exp(-6),0.5),matrix(c(-1,0,1,-1),nrow=2,ncol=2,byrow=TRUE),0.05,100)

#Figure 2.5 plots 100 Pl simulations of It 
PLXplot<-function(x,N,T,c,S,tau,X){
  data<-PLXdata(x,N,T,c,S,tau,X)
  plot(ts(data[,1], start=0, end=14,deltat=0.05), xlab="Time",ylab=" Number of infected",col=rgb(runif(1),runif(1),runif(1)),main="100 simulations of It using the poisson leap method overlayed with the mean of all simulations",ylim=c(0,350))
  for(i in 2:X){
    lines(y=data[,i],x=data[,X+2],type="l", col=rgb(runif(1),runif(1),runif(1)))
    i=i+1
  }
  lines(y=data[,X+1],x=data[,X+2],type="l", col="black",lwd=2)
}
PLXplot(c(762,1),763,14,c(exp(-6),0.5),matrix(c(-1,0,1,-1),nrow=2,ncol=2,byrow=TRUE),0.05,100)

#calculating mean peak number of infected and mean of the day that they peak at 
#given PLXdata
meanPeak<-function(Data){
  X<-ncol(Data) - 2
  data<-Data[,1:X]
  col_max<-apply(data,2,max)
  mean<-mean(col_max)
  max_rows<-apply(data,2,which.max)
  meanrow<-mean(max_rows)
  meanrow<-round(meanrow)
  meanday<-Data[121,X+2]
  m<-c(mean,meanday)
  return(m)
}

#Figure 2.6 summary graph 
x<-c(762,1)
N<-763
c<-c(exp(-6),0.5)
S<-matrix(c(-1,0,1,-1),nrow=2,ncol=2,byrow=TRUE)
X<-100

Data<-PLXdata(c(762,1),763,14,c(exp(-6),0.5),matrix(c(-1,0,1,-1),nrow=2,ncol=2,byrow=TRUE),0.05,100) 
data<-Data[,1:(X-2)]
meandata<-Data[,(X+1)]
m<-meanPeak(Data)
SD<-apply(data,1,sd)
lower_quantile<-apply(data,1,function(x) quantile(x,probs=0.25))
upper_quantile<-apply(data,1,function(x) quantile(x,probs=0.75))
upper_sd_bound<-meandata+SD
lower_sd_bound<-meandata-SD
time<-Data[,(X+2)]
plot(ts(meandata, start=0, end=14,deltat=0.05), xlab="Time",ylab=" Number of infected",col="red",lwd=2,main="Summary statistics of 100 poisson leap simulations",ylim=c(0,350))
lines(x=time, y=upper_quantile,col="black")
lines(x=time, y=lower_quantile,col="black")
lines(x=time, y=upper_sd_bound,col="grey")
lines(x=time, y=lower_sd_bound,col="grey")
lines(x=m[2],y=m[1],type="p",col="blue")

#Generating simulations using 
CLE<-function(x,N,T,c,dt){
  t<-0
  data<-data.frame(St=x[1],It=x[2],Rt=(N-x[1]-x[2]),t=0)
  while (t< T-dt){
    h1<-c[1]*x[1]*x[2]
    h2<-c[2]*x[2]
    dW1<- rnorm(1,0,sqrt(dt))
    dW2<-rnorm(1,0,sqrt(dt))
    t<-t+dt
    x[1]<-x[1] - c[1]*x[1]*x[2]*dt - sqrt(c[1]*x[1]*x[2])*dW1
    x[2]<-x[2]+ c[1]*x[1]*x[2]*dt - c[2]*x[2]*dt + sqrt(c[1]*x[1]*x[2])*dW1 - sqrt(c[2]*x[2])*dW2
    data<-rbind(data,c(x[1],x[2],(N-x[1]-x[2]),t))
  }
  return(data)  
  
}

#Figure 2.7 
CLEdata<-CLE(c(762,1),763,14,c(exp(-6),0.5),0.01)
ggplot()+geom_line(data=CLEdata,aes(x=t,y=St),color="red", linewidth=0.8)+geom_line(data=CLEdata,aes(x=t,y=It),color="green",linewidth=0.8)+geom_line(data=CLEdata,aes(x=t,y=Rt),color="blue",linewidth=0.8)+xlab("Time after the first case in Days")+ylab("Number of individuals ")+ggtitle("Simulations of St, It, Rt using CLE solution")

#lna first algorithm 

LNA<-function(a,b,X0,V0,dt,T,N){
  X<-X0
  Z<-X0
  M<-X0-Z
  t<-0
  S<-Z[1,1]
  I<-Z[2,1]
  R<-N-S-I
  V<-V0
  dataLNA<-data.frame(St=S, It=I, Rt=R, t=0)
  while (t <T){
    Z<-Z+matrix(data=c(-a*Z[1,1]*Z[2,1],a*Z[1,1]*Z[2,1] - (b*Z[2,1])), nrow=2,ncol=1)*dt
    M<-M+matrix(data=c(-a*Z[2,1]*M[1,1]-(a*Z[1,1]*M[2,1]),a*Z[2,1]*M[1,1]+a*Z[1,1]*M[2,1]-b*M[2,1]),nrow=2,ncol=1)*dt
    V<-V+matrix(data=c(a*Z[1,1]*Z[2,1],-a*Z[1,1]*Z[2,1],-a*Z[1,1]*Z[2,1],a*Z[1,1]*Z[2,1]+b*Z[2,1]),nrow=2,ncol=2,byrow = TRUE)*dt
    t<-t+dt
    X[2,1]<-0
    while(X[2,1]<=0){
      X<-mvrnorm(1,mu=Z+M,Sigma=V)
      X<-as.matrix(X)
      # print(X)
    }
    M<-X-Z
    V<-V0
    dataLNA<-rbind(dataLNA,c(X[1,1],X[2,1],(N-X[1,1]-X[2,1]),t))
  }
  return(dataLNA)
}
LNAdata=LNA(exp(-6),0.5,X0<- matrix(data=c(762,1),nrow=2,ncol=1),V0<-matrix(data=rep(0,4),nrow=2,ncol=2,byrow=TRUE),0.01,14,763)

#lna second algorithm

LNA2<-function(a,b,X0,V0,dt,T,N){
  X<-X0
  Z<-X0
  M<-X0-Z
  t<-0
  S<-Z[1,1]
  I<-Z[2,1]
  R<-N-S-I
  V<-V0
  LNAdata<-data.frame(St=S, It=I, Rt=R, t=0)
  while (t <T){
    Z<-Z+matrix(data=c(-a*Z[1,1]*Z[2,1],a*Z[1,1]*Z[2,1] - (b*Z[2,1])), nrow=2,ncol=1)*dt
    V<-matrix(data=c(a*Z[1,1]*Z[2,1],-a*Z[1,1]*Z[2,1],-a*Z[1,1]*Z[2,1],a*Z[1,1]*Z[2,1]+b*Z[2,1]),nrow=2,ncol=2,byrow = TRUE)*dt
    t<-t+dt
    X[2,1]<-0
    while(X[2,1]<=0){
      X<-mvrnorm(1,mu=Z,Sigma=V)
      X<-as.matrix(X)
      # print(X)
    }
    V<-V0
    Z<-X
    LNAdata<-rbind(LNAdata,c(X[1,1],X[2,1],(N-X[1,1]-X[2,1]),t))
  }
  return(LNAdata)
}
LNA2data=LNA2(exp(-6),0.5,X0<- matrix(data=c(762,1),nrow=2,ncol=1),V0<-matrix(data=rep(0,4),nrow=2,ncol=2,byrow=TRUE),0.01,14,763)

# Figure 2.8
p1<-ggplot()+geom_line(data=LNAdata,aes(x=t,y=St),color="red")+geom_line(data=LNAdata,aes(x=t,y=It),color="green")+geom_line(data=LNAdata,aes(x=t,y=Rt),color="blue")+xlab("Time after the first case in days")+ylab("Number of individuals")+ggtitle("Simulations of St, It, Rt using LNA algorithm 1")
p2<-ggplot()+geom_line(data=LNA2data,aes(x=t,y=St),color="red")+geom_line(data=LNA2data,aes(x=t,y=It),color="green")+geom_line(data=LNA2data,aes(x=t,y=Rt),color="blue")+xlab("Time after the first case in days")+ylab("Number of individuals")+ggtitle("Simulations of St, It, Rt using LNA algorithm 2")
grid.arrange(p1, p2, ncol=1)

#Figure 3.1
y<-c(1,3,6,25,73,221,294,257,236,189,125,67,26,10,3)
t<-c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14)
#graph of observed
par(mfrow=c(1,1))
plot(y=y,x=t,type="p",xlab="Time",ylab=" Number of infected",col="black",lwd=2,main="Observed data")

#Figure 3.2
#graph of all simulations
par(mfrow=c(3,2))

#ODEdata<-values(c(exp(-6),0.5),0.01,14,c(762,1))
plot(ts(ODEdata[,2], start=0, end=14,deltat=0.01), xlab="Time",ylab=" Number of infected",col="orange",main="ODE simulations compared with observed data")
lines(x=t,y=y,type="p")

#GAdata<-GAvalues(c(762,1),763,14,c(exp(-6),0.5),matrix(c(-1,0,1,-1),nrow=2,ncol=2,byrow=TRUE))
plot(y=GAdata[,2],x=GAdata[,4],type="l",xlab="Time",ylab=" Number of infected",col="green",main="GA simulations compared with observed data")
lines(x=t,y=y,type="p")

#PLdata<-PLvalues(c(762,1),763,14,c(exp(-6),0.5),matrix(c(-1,0,1,-1),nrow=2,ncol=2,byrow=TRUE),0.01)
plot(ts(PLdata[,2], start=0, end=14,deltat=0.01), xlab="Time",ylab=" Number of infected",col="purple",main="Poisson Leap simulations compared with observed data")
lines(x=t,y=y,type="p")

#CLEdata<-CLE(c(762,1),763,14,c(exp(-6),0.5),0.01)
plot(ts(CLEdata[,2], start=0, end=14,deltat=0.01), xlab="Time",ylab=" Number of infected",col="red",main="CLE simulations compared with observed data")
lines(x=t,y=y,type="p")

#LNAdata<-LNA(exp(-6),0.5,matrix(data=c(762,1),nrow=2,ncol=1),matrix(data=rep(0,4),nrow=2,ncol=2,byrow=TRUE),0.01,14,763)
plot(ts(LNAdata[,2], start=0, end=14,deltat=0.01), xlab="Time",ylab=" Number of infected",col="blue",main="LNA compared with observed data")
lines(x=t,y=y,type="p")

plot(ts(ODEdata[,2], start=0, end=14,deltat=0.01), xlab="Time",ylab=" Number of infected",col="orange",main="Comparison of all simulation methods with observed data")
lines(y=GAdata[,2],x=GAdata[,4],type="l", col="green")
lines(y=PLdata[,2],x=PLdata[,4],type="l", col="purple")
lines(y=CLEdata[,2],x=CLEdata[,4],type="l", col="red")
lines(y=LNAdata[,2],x=LNAdata[,4],type="l", col="blue")
lines(x=t,y=y,type="p")

#Function that calculates MLE for generated noisey data through the ODE model 

loglikeNeg<-function(param, Ydata, dt, intert, X0)
{
  n<-length(Ydata)
  T<-(n-1)*intert 
  incr<-intert/dt
  out1<-values(t=T, dt=dt, X0=X0, theta=exp(param))
  out2<-out1[1+(0:(n-1))*(incr),2] 
  -1*sum(dnorm(Ydata, out2, exp(param[3]), log=TRUE)) 
}

MLE<-function(param,dt,T,X0,intert){
  t=T
  ODEdata<-values(param,dt,t=T,X0)
  set.seed(1)
  ODEIdata<-ODEdata[1+(0:(T/intert))*(intert/dt),2]   
  Ydata<-ODEIdata+ rnorm(length(ODEIdata),0,param[3])
  init<-log(param-(param/5))
  out<-optim(init,loglikeNeg,Ydata=Ydata,dt=dt,intert=intert,X0=X0)
  MLE<-exp(out$par)
  return(MLE)
}

#Code for values in table 3.1
MLE(param=c(exp(-6),0.5,0.5),dt=0.01,T=14,X0=c(762,1),intert=1)

#Metropolis hastings algorithm for sir models 
#plot(ts(y,start=0,deltat=1),main="",xlab="Time (days)",ylab="Yt")
alpha<-exp(-6)
beta<-0.5
param<-c(alpha,beta)
dt<-0.01
X0<-c(762,1)
t<-14

ODEdata<-values(param,dt,t,X0)
sigma<-0.5
theta<-c(log(alpha),log(beta),log(sigma)) ## c(theta1,theta2,theta3)

## we take the prior of all the log parameters to be standard normal 
## therefore say the prior for all parameters are g(theta1), g(theta2), 
## g(theta3)
## the posterior density is therefore proportional to all three priors and the 
## likelihood 
## g(theta|y) is proportional to g(theta1)*g(theta2)*g(theta3)*f(y|theta)

## what is log f(y|theta)?
## this function returns the loglikelihood of It which is the mproposed model
## y which are our observed data set is out target model 
theta<-c(log(alpha),log(beta),log(sigma)) ## c(theta1,theta2,theta3)

## data=y which is our observed data 
loglike<-function(theta,data,dt,X0){
  n<-length(data)
  incr<-1/dt
  out<-values(exp(theta),dt,n-1,X0) ##returns It values at dt inytervals 
  out<-out[1+(0:(n-1))*incr,2] ##returns It values at uit time intervals 
  return(sum(dnorm(data,out,exp(theta[3]),log=TRUE)))
}

##we then find the log prior log(g(theta))

logprior<-function(theta){
  return(sum(dnorm(theta,log=TRUE)))
}

## we us the posterior = g(theta|y) as our target density and a random walk 
## proposal = theta*|theta ~ N(theta, SigmaMatrix)

## metropolis hastings algorithm 

## N = number of itterations
## init= initial values for theta 
##SigmaMatrix is our tuning matrix 

MH=function(N,data,dt,SigmaMatrix,init,X0)
{
  theta.mat = matrix(0,nrow=N,ncol=3)  ##create a matrix of 0s with N rows 
  theta.mat[1,] = init ## first row of our theta matrix has our initial values 
  curr = init 
  count = 0 
  for (j in 2:N) 
  {
    can = mvrnorm(1,curr,SigmaMatrix) #proposal 
    #Evaluate log acceptance probability
    laprob = logprior(can)-logprior(curr)+
      loglike(can,data,dt,X0)-loglike(curr,data,dt,X0)
    if (log(runif(1)) < laprob) #Accept with the correct prob.
    {
      curr = can; #chain moves
      count = count + 1 #add one to acceptance counter
    }
    theta.mat[j,] = curr
  }
  print(count/(N-1)) # empirical acceptance prob
  return(theta.mat)
}

## N=1000
## dt=0.01
## tuning matrix SigmaMatrix<- diag(0.001,0.001,0.001)
## initial values of theta = (-6, log(0.5),2.5)
## X0= (762,1)

#Figure 3.3
MHdata<-MH(5000,y,0.01,diag(rep(0.001,3)),c(-6,log(0.5),2.5),c(762,1))
plot(ts(exp(MHdata)),xlab="Iteration count", main="Trace plots for log(alpha), log(beta) and log(sigma)")


## we get a decent trace plot for alpha and beta but not sigma 
## we choose our tuning matrix to be (((2.38)^2)/3)* var(theta|y)
#figure 3.4
newSigma<- ((2.38^2)/3)*var(MHdata)
system.time(MHdata2<-MH(5000,y,0.01,newSigma,c(-6,log(0.5),2.5),c(762,1)))
effectiveSize(MHdata2)
plot(ts(exp(MHdata2)),xlab="Iteration count", main="Trace plots for log(alpha), log(beta) and log(sigma) using optimal tuning matrix")
##plot with 10000 iterations 
#figure 3.5
system.time(MHdata3<-MH(10000,y,0.01,newSigma,c(-6,log(0.5),2.5),c(762,1)))
effectiveSize(MHdata3)
plot(ts(exp(MHdata3)),xlab="Iteration count", main="Trace plots for log(alpha), log(beta) and log(sigma) using optimal tuning matrix and 10000 iterations")
## as theta is log of all parameters we transform it back to original form 

MHdata4<-exp(MHdata3)

#figure 3.6 plotting 20 runs of 5000 iters on the same trace plot
par(mfrow=c(2,1))
ALPHA<-MHdata2[,1]
plot(ts(ALPHA),ylim=c(-6.2,-6),ylab="log(alpha)", xlab="Iteration count", main="20 trace plots for log(alpha) created using different initial values")
for(i in 1:20)
{
  a<-runif(1,-6.25,-5.75)
  b<-runif(1,-1,-0.5)
  sigma<-runif(1,2,3.5)
  param<-c(a,b,sigma)
  #print(param)
  data<-MH(5000,y,0.01, newSigma,init=c(a,b,sigma),c(762,1))
  log_alpha<-data[,1]
  lines(ts(log_alpha),col=(i+1))
}

BETA<-MHdata2[,2]
plot(ts(BETA),ylim=c(-1,-0.5),ylab="log(beta)", xlab="Iteration count", main="20 trace plots for log(beta) created using different initial values")
for(i in 1:20)
{
  a<-runif(1,-6.25,-5.75)
  b<-runif(1,-1,-0.5)
  sigma<-runif(1,2,3.5)
  param<-c(a,b,sigma)
  #print(param)
  data<-MH(5000,y,0.01, newSigma,init=c(a,b,sigma),c(762,1))
  log_beta<-data[,2]
  lines(ts(log_beta),col=(i+1))
  
}

#after removing 500 iterations 
par(mfrow=c(2,1))
ALPHA<-MHdata2[,1]
plot(ts(ALPHA[501:5000]),ylim=c(-6.2,-6),ylab="log(alpha)", xlab="Iteration count", main="20 trace plots for log(alpha) created using different initial values")
for(i in 1:20)
{
  a<-runif(1,-6.25,-5.75)
  b<-runif(1,-0.9,-0.6)
  sigma<-runif(1,2,3.5)
  param<-c(a,b,sigma)
  print(param)
  data<-MH(5000,y,0.01, newSigma,init=c(a,b,sigma),c(762,1))
  log_alpha<-data[,1]
  lines(ts(log_alpha[501:5000]),col=(i+1))
}

BETA<-MHdata2[,2]

plot(ts(BETA[501:5000]),ylim=c(-0.9,-0.6),ylab="log(beta)", xlab="Iteration count", main="20 trace plots for log(beta) created using different initial values")
for(i in 1:20)
{
  a<-runif(1,-6.25,-5.75)
  b<-runif(1,-0.9,-0.6)
  sigma<-runif(1,2,3.5)
  param<-c(a,b,sigma)
  print(param)
  data<-MH(5000,y,0.01, newSigma,init=c(a,b,sigma),c(762,1))
  log_beta<-data[,2]
  lines(ts(log_beta[501:5000]),col=(i+1))
  
}

##figure 3.7 we plot the kernal density estimates 
par(mfrow=c(1,2))
alpha<-MHdata4[,1]
beta<-MHdata4[,2]
sigma<-MHdata4[,3]

plot(density(alpha),xlab="alpha", xlim=c(0.002,0.0024),main="Kernel density plot for alpha")
plot(density(beta),xlab="beta",main="Kernel density plot for beta")
#plot(density(sigma),xlab="sigma")

#table 3.3
## summaries of alpha(infection rate)

mean(MHdata4[,1])
sd(MHdata4[,1])
quantile(MHdata4[,1],c(0.025,0.975))

## summaries of beta

mean(MHdata4[,2])
sd(MHdata4[,2])
quantile(MHdata4[,2],c(0.025,0.975))



#figure 3.8
##histograms for theta
hist(alpha,xlab="alpha", breaks=15,freq=FALSE)
hist(beta,xlab="beta", breaks=15,freq=FALSE)
#hist(sigma,xlab="sigma", breaks=15, freq=FALSE)

hist(alpha,xlab="alpha", breaks=25,freq=FALSE,xlim=c(0.002,0.0024), ylim=c(0,13000) ,main="Density histogram for alpha")
lines(density(alpha),xlab="alpha",col="red")
hist(beta,xlab="beta", breaks=25,freq=FALSE,ylim=c(0,30),main="Density histogram for beta")
lines(density(beta),xlab="beta", col="red")
#hist(sigma,xlab="sigma", breaks=15, freq=FALSE)
#lines(density(sigma),xlab="sigma", col="red")

#predictive posterior distributions
#need values function 
#already generated MH data 
par(mfrow=c(1,1))
N<-10000
dt=0.1
endT<-14
S.mat = matrix(0,nrow=endT/dt+1,ncol=N)
I.mat = matrix(0,nrow=endT/dt+1,ncol=N)
R.mat = matrix(0,nrow=endT/dt+1,ncol=N)
for(i in 1:N)
{
  data<-values(theta=exp(MHdata3[i,1:2]),dt=dt,t=endT,X0=c(762,1))
  data$Rt<-763-data$X1-data$X2
  S.mat[,i]<-data[,1]
  I.mat[,i]<-data[,2]
  R.mat[,i]<-data[,4]
}

S.mean <- apply(S.mat,1,mean)
S.lq <- apply(S.mat,1,quantile,0.025)
S.uq <- apply(S.mat,1,quantile,0.975)
I.mean <- apply(I.mat,1,mean)
I.lq <- apply(I.mat,1,quantile,0.025)
I.uq <- apply(I.mat,1,quantile,0.975)
R.mean <- apply(R.mat,1,mean)
R.lq <- apply(R.mat,1,quantile,0.025)
R.uq <- apply(R.mat,1,quantile,0.975)

#figure 3.9 PREDICTIVE POSTERIOR DISTRIBUTION GRAPH 
plot(ts(S.mean,start=0,deltat=0.1),ylim=c(0,800),lwd=2,col="red",xlab="Time", ylab="Number of individuals", main="Predictive posterior distribution")
lines(ts(S.lq,start=0,deltat=0.1),col=2)
lines(ts(S.uq,start=0,deltat=0.1),col=2)
lines(ts(I.mean,start=0,deltat=0.1),lwd=2,col="green")
lines(ts(I.lq,start=0,deltat=0.1),col=3)
lines(ts(I.uq,start=0,deltat=0.1),col=3)
lines(ts(R.mean,start=0,deltat=0.1),lwd=2,col="blue")
lines(ts(R.lq,start=0,deltat=0.1),col=4)
lines(ts(R.uq,start=0,deltat=0.1),col=4)

#figure 3.10 COMPARISON GRAPH 
plot(ts(I.mean,start=0,deltat=0.1),ylim=c(0,350),lwd=2,col="green",xlab="Time", ylab="Number of infected individuals")
ODEdata<-values(c(exp(-6),0.5),0.1,14,c(762,1))
lines(ts(ODEdata[,2],start=0,deltat=0.1),col="black")
y<-c(1,3,6,25,73,221,294,257,236,189,125,67,26,10,3)
t<-c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14)
lines(x=t,y=y,type="p")

## basic reproduction number R0

##figure 3.11
R0=MHdata4[,1]/MHdata4[,2]*763
hist(R0,xlab="R0", main="Basic reproduction number for each parameter output of MH MCMC")
quantile(R0,c(0.025,0.975))

#rejection ABC using euclidean distance between points 
#observed data 

y<-c(1,3,6,25,73,221,294,257,236,189,125,67,26,10,3)
alpha<-exp(-6)
beta<-0.5
sigma<-0.5
param<-c(alpha,beta,sigma)
tau<-0.1
x<-c(762,1)
T<-14
N<-763
iters=1000
S=matrix(c(-1,0,1,-1),nrow=2,ncol=2,byrow=TRUE)
dt=0.1
# euclidean distance 

euclid<-function(data1,data2){
  vec<-data1-data2
  sqrt(sum(vec^2))
}


ABCrejection<-function(iters,epsilon,y,x,N,S,dt){
  T<-length(y)-1
  theta.mat<-matrix(rep(0,2*iters),ncol=2,nrow=iters)
  j<-1
  overall<-0
  while(j<=iters){
    overall<-overall+1
    logalpha<-rnorm(1,-6,sqrt(0.1))
    logbeta<-rnorm(1,log(0.5),sqrt(0.01))
    alpha<-exp(logalpha)
    beta<-exp(logbeta)
    sim<-PLvalues(x,N,T,c(alpha,beta),S,dt)
    sim<-sim$It
    sim<-sim[1+(0:T)*(1/dt)]
    distance<-euclid(y, sim)
    print(distance)
    if(distance<=epsilon){
      theta.mat[j,1]<-alpha
      theta.mat[j,2]<-beta
      j<-j+1
      
    } 
    
  }
  
  print(overall)
  return(theta.mat)
}

#Figure 4.1
system.time(out<-ABCrejection(1000,150, y,c(762,1),763,matrix(c(-1,0,1,-1),nrow=2,ncol=2,byrow=TRUE),0.1))
effectiveSize(out)

par(mfrow=c(2,2))      
hist(out[,1],freq=FALSE,main="Posterior distribution of alpha using epsilon=150", xlab="value",breaks=15)
lines(density(out[,1]),col="red")
hist(out[,2],freq=FALSE,main="Posterior distribution of beta using epsilon=150", xlab="value",breaks=15)
lines(density(out[,2]),col="red")

system.time(out2<-ABCrejection(1000,100, y,c(762,1),763,matrix(c(-1,0,1,-1),nrow=2,ncol=2,byrow=TRUE),0.1))
effectiveSize(out2)

hist(out2[,1],freq=FALSE,main="Posterior distribution of alpha using epsilon=100", xlab="value",breaks=15)
lines(density(out2[,1]),col="red")
hist(out2[,2],freq=FALSE,main="Posterior distribution of beta using epsilon=100", xlab="value",breaks=15)
lines(density(out2[,2]),col="red")
#Table 4.1
effectiveSize(out)
mean(out[,1])
sd(out[,1])
quantile(out[,1],c(0.025,0.975))

## summaries of beta

mean(out[,2])
sd(out[,2])
quantile(out[,2],c(0.025,0.975))

#Table 4.2
effectiveSize(out2)
mean(out2[,1])
sd(out2[,1])
quantile(out2[,1],c(0.025,0.975))

## summaries of beta

mean(out2[,2])
sd(out2[,2])
quantile(out2[,2],c(0.025,0.975))

#abc rejection using summary statistics 
stats<-function(data){
  data<-as.data.frame(data)
  infected<-data[,2]
  s1<-max(infected)
  row<-data[data$It==s1,]
  row<-as.data.frame(row)
  s2<-row[1,4]
  daily<-infected[1+(0:14)*(1/0.1)]
  s3<-sum(daily)
  s4<-daily[15]
  s<-c(s1,s2,s3,s4)
  return(s)
}

ABCrejection2<-function(iters,epsilon,y,x,N,S,dt){
  T<-length(y)-1
  theta.mat<-matrix(rep(0,2*iters),ncol=2,nrow=iters)
  summary<-matrix(rep(0,4*iters),ncol=4,nrow=iters)
  j<-1
  overall<-0
  while(j<=iters){
    overall<-overall+1
    logalpha<-rnorm(1,-6,sqrt(0.1))
    logbeta<-rnorm(1,log(0.5),sqrt(0.01))
    alpha<-exp(logalpha)
    beta<-exp(logbeta)
    data<-PLvalues(x,N,T,c(alpha,beta),S,dt)
    summarySim<-stats(data)
    summaryY<-c(294,6,1536,3)
    distance<-euclid(summaryY, summarySim)
    #print(distance)
    if(distance<=epsilon){
      theta.mat[j,1]<-alpha
      theta.mat[j,2]<-beta
      #print(summarySim)
      #print(distance)
      j<-j+1
    } 
  }
  print(overall)
  return(theta.mat)
}

#Figure 4.2
system.time(out3<-ABCrejection2(1000,100, y,c(762,1),763,matrix(c(-1,0,1,-1),nrow=2,ncol=2,byrow=TRUE),0.1))
effectiveSize(out3)
system.time(out4<-ABCrejection2(1000,40, y,c(762,1),763,matrix(c(-1,0,1,-1),nrow=2,ncol=2,byrow=TRUE),0.1))
effectiveSize(out4)

par(mfrow=c(2,2))      
hist(out3[,1],freq=FALSE,main="Posterior distribution of alpha using epsilon=100", xlab="value",breaks=15)
lines(density(out3[,1]),col="red")
hist(out3[,2],freq=FALSE,main="Posterior distribution of beta using epsilon=100", xlab="value",breaks=15)
lines(density(out3[,2]),col="red")

hist(out4[,1],freq=FALSE,main="Posterior distribution of alpha using epsilon=50", xlab="value",breaks=15)
lines(density(out4[,1]),col="red")
hist(out4[,2],freq=FALSE,main="Posterior distribution of beta using epsilon=50", xlab="value",breaks=15)
lines(density(out4[,2]),col="red")
#Table 4.3
effectiveSize(out3)

mean(out3[,1])
sd(out3[,1])
quantile(out3[,1],c(0.025,0.975))

## summaries of beta

mean(out3[,2])
sd(out3[,2])
quantile(out3[,2],c(0.025,0.975))

#Table 4.4
effectiveSize(out4)
mean(out4[,1])
sd(out4[,1])
quantile(out4[,1],c(0.025,0.975))

## summaries of beta

mean(out4[,2])
sd(out4[,2])
quantile(out4[,2],c(0.025,0.975))


#figure 4.3
summary<-function(iters,T,x,N,S,dt){
  summary<-matrix(rep(0,4*iters),ncol=4,nrow=iters)
  j<-1
  overall<-0
  while(j<=iters){
    overall<-overall+1
    logalpha<-rnorm(1,-6,sqrt(0.1))
    logbeta<-rnorm(1,log(0.5),sqrt(0.01))
    alpha<-exp(logalpha)
    beta<-exp(logbeta)
    data<-PLvalues(x,N,T,c(alpha,beta),S,dt)
    summarySim<-stats(data)
    summary[j,1]<-summarySim[1]
    summary[j,2]<-summarySim[2]
    summary[j,3]<-summarySim[3]
    summary[j,4]<-summarySim[4]
    j<-j+1
  } 
  return(summary)
}

summarystat<-summary(1000,14,c(762,1),763,matrix(c(-1,0,1,-1),nrow=2,ncol=2,byrow=TRUE),0.1)
maxinfected<-summarystat[,1]
maxday<-summarystat[,2]
totalinfected<-summarystat[,3]
endinfected<-summarystat[,4]


sd1<-sd(maxinfected)
sd2<-sd(maxday)
sd3<-sd(totalinfected)
sd4<-sd(endinfected)
sd<-c(sd1,sd2,sd3,sd4)
sd

NormalisedDist<-function(summarySim, summaryY,sd){
  x<-summarySim-summaryY
  x<-c(x[1]/sd[1],x[2]/sd[2],x[3]/sd[3],x[4]/sd[4])
  x<-x*x
  x<-sum(x)
  x<-sqrt(x) 
  return(x)
}


ABCnorm<-function(iters,epsilon){
  T<-14
  theta.mat<-matrix(rep(0,2*iters),ncol=2,nrow=iters)
  j<-1
  overall<-0
  summaryY<-c(294,6,1536,3)
  while(j<=iters){
    overall<-overall+1
    logalpha<-rnorm(1,-6,sqrt(0.1))
    logbeta<-rnorm(1,log(0.5),sqrt(0.01))
    alpha<-exp(logalpha)
    beta<-exp(logbeta)
    data<-PLvalues(c(762,1),763,14,c(alpha,beta),matrix(c(-1,0,1,-1),nrow=2,ncol=2,byrow=TRUE),0.1)
    s<-stats(data)
    distance<-NormalisedDist(s, summaryY,sd)
    #print(distance)
    if(distance<=epsilon){
      theta.mat[j,1]<-alpha
      theta.mat[j,2]<-beta
      #print(theta.mat[j,])
      j<-j+1
    } 
  }
  print(overall)
  return(theta.mat)
}

system.time(norm<-ABCnorm(1000,0.75))
effectiveSize(norm)


par(mfrow=c(2,2))      
hist(norm[,1],freq=FALSE,main="Posterior distribution of alpha using epsilon=0.75", xlab="value",breaks=15)
lines(density(norm[,1]),col="red")
hist(norm[,2],freq=FALSE,main="Posterior distribution of beta using epsilon=0.75", xlab="value",breaks=15)
lines(density(norm[,2]),col="red")


effectiveSize(norm)
mean(norm[,1])
sd(norm[,1])
quantile(norm[,1],c(0.025,0.975))

## summaries of beta

mean(norm[,2])
sd(norm[,2])
quantile(norm[,2],c(0.025,0.975))

system.time(norm2<-ABCnorm(1000,0.5))
effectiveSize(norm2)

hist(norm2[,1],freq=FALSE,main="Posterior distribution of alpha using epsilon=0.5", xlab="value",breaks=15)
lines(density(norm2[,1]),col="red")
hist(norm2[,2],freq=FALSE,main="Posterior distribution of beta using epsilon=0.5", xlab="value",breaks=15)
lines(density(norm2[,2]),col="red")


effectiveSize(norm2)

mean(norm2[,1])
sd(norm2[,1])
quantile(norm2[,1],c(0.025,0.975))


## summaries of beta

mean(norm2[,2])
sd(norm2[,2])
quantile(norm2[,2],c(0.025,0.975))


#ABC MCMC

#theta<-c(log(alpha),log(beta))

##we then find the log prior log(g(theta))

#Make consistent with rejection sampler
logprior2<-function(theta){
  return(sum(dnorm(theta,c(-6,log(0.5)),c(0.1,0.1),log=TRUE)))
}


#install.packages("MASS")
#library(MASS)

ABCMCMC=function(iters,y,N,T,dt,SigmaMatrix,init,x,S,epsilon)
{
  theta.mat = matrix(0,nrow=iters,ncol=2) 
  theta.mat[1,] = init 
  curr = init 
  count = 0
  summaryY<-c(294,6,1536,3)
  for (j in 2:iters) 
  {
    can = mvrnorm(1,curr,SigmaMatrix)
    data<-PLvalues(x,N,T,exp(can),S,dt)
    s<-stats(data)
    distance<-NormalisedDist(summaryY,s,sd)
    if (distance<=epsilon){
      laprob = logprior2(can)-logprior2(curr)
      if (log(runif(1)) < laprob) 
      {
        curr = can; 
        count = count + 1 
      }}
    theta.mat[j,] = curr
  }
  print(count/(iters-1)) 
  return(theta.mat)
}
#run 1 using inital variance matrix and epsilon is 0.5 for normalised euclidean distance 

SIG1<-diag(c(0.001,0.001))
SIG1
system.time(DATA1<-ABCMCMC(25000,y,763,14,0.1,SIG1,c(-6,log(0.5)),c(762,1),matrix(c(-1,0,1,-1),nrow=2,ncol=2,byrow=TRUE),0.5))
plot(ts(exp(DATA1)),xlab="Iteration count", main="Trace plots for log(alpha) and log(beta) using initial variance matrix ")
effectiveSize(DATA1)


#run 2 using optimal variance matrix and epsilon is 0.5 for normalised euclidean distance 
SIG2<- ((2.38^2)/2)*var(DATA1)
SIG2
system.time(DATA2<-ABCMCMC(50000,y,763,14,0.1,SIG2,c(-6,log(0.5)),c(762,1),matrix(c(-1,0,1,-1),nrow=2,ncol=2,byrow=TRUE),0.5))
plot(ts(exp(DATA2)),xlab="Iteration count", main="Trace plots for log(alpha) and log(beta) using optimal tuning matrix and 50000 runs")
effectiveSize(DATA2)


DATA1<-exp(DATA1)

mean(DATA1[,1])
sd(DATA1[,1])
quantile(DATA1[,1],c(0.025,0.975))


## summaries of beta

mean(DATA1[,2])
sd(DATA1[,2])
quantile(DATA1[,2],c(0.025,0.975))

DATA2<-exp(DATA2)
mean(DATA2[,1])
sd(DATA2[,1])
quantile(DATA2[,1],c(0.025,0.975))


## summaries of beta

mean(DATA2[,2])
sd(DATA2[,2])
quantile(DATA2[,2],c(0.025,0.975))

#abc mcmc that iteratively calculates sd to normalise

ABCMCMC2=function(iters,y,N,T,dt,SigmaMatrix,init,x,S,epsilon){
  theta.mat<-matrix(0, nrow=iters, ncol=2)
  theta.mat[1,]<-init
  curr<-init
  count<-0
  summaryY<-c(294, 6, 1536, 3)
  data<-PLvalues(x,N,T,exp(curr),S,dt)
  s<-stats(data)
  all_summaries<-matrix(0, nrow=iters, ncol=4)
  all_summaries[1,]<-s
  sd<-rep(1, length(s))
  for(j in 2:iters){
    can<-mvrnorm(1,curr,SigmaMatrix)
    data<-PLvalues(x,N,T,exp(can),S,dt)
    s<-stats(data)
    all_summaries[j,]<-s 
    sd<-apply(all_summaries[1:j,,drop=FALSE],2,sd,na.rm=TRUE)
    distance<-NormalisedDist(summaryY,s,sd)
    #print(distance)
    if(distance<=epsilon){
      laprob<-logprior2(can)-logprior2(curr)
      if(log(runif(1))<laprob){
        curr<-can
        count<-count+1
      }
    }
    theta.mat[j,]<-curr
  }
  print(count/iters)
  return(theta.mat)
}


SIG3<-diag(c(0.001,0.001))
SIG3
system.time(DATA3<-ABCMCMC2(25000,y,763,14,0.1,SIG3,c(-6,log(0.5)),c(762,1),matrix(c(-1,0,1,-1),nrow=2,ncol=2,byrow=TRUE),1))
plot(ts(exp(DATA3)),xlab="Iteration count", main="Trace plots for log(alpha) and log(beta) using initial variance matrix ")
effectiveSize(DATA3)

#run 2 using optimal variance matrix and epsilon is 1 for normalised euclidean distance 
SIG4<-((2.38^2)/2)*var(DATA3)
SIG4
system.time(DATA4<-ABCMCMC2(50000,y,763,14,0.1,SIG4,c(-6,log(0.5)),c(762,1),matrix(c(-1,0,1,-1),nrow=2,ncol=2,byrow=TRUE),1))
plot(ts(exp(DATA4)),xlab="Iteration count", main="Trace plots for log(alpha) and log(beta) using optimal tuning matrix and 50000 runs")
effectiveSize(DATA4)


DATA3<-exp(DATA3)

mean(DATA3[,1])
sd(DATA3[,1])
quantile(DATA3[,1],c(0.025,0.975))


## summaries of beta

mean(DATA3[,2])
sd(DATA3[,2])
quantile(DATA3[,2],c(0.025,0.975))

DATA4<-exp(DATA4)
mean(DATA4[,1])
sd(DATA4[,1])
quantile(DATA4[,1],c(0.025,0.975))


## summaries of beta

mean(DATA4[,2])
sd(DATA4[,2])
quantile(DATA4[,2],c(0.025,0.975))

#comparing the three algorithms for the same epsilon 

system.time(norm2<-ABCnorm(1000,1))
effectiveSize(norm2)
## summaries for alpha
mean(norm2[,1])
sd(norm2[,1])
quantile(norm2[,1],c(0.025,0.975))

## summaries of beta

mean(norm2[,2])
sd(norm2[,2])
quantile(norm2[,2],c(0.025,0.975))

#ABC MCMC predetermined sd 
SIG1<-diag(c(0.001,0.001))
SIG1
system.time(DATA1<-ABCMCMC(25000,y,763,14,0.1,SIG1,c(-6,log(0.5)),c(762,1),matrix(c(-1,0,1,-1),nrow=2,ncol=2,byrow=TRUE),1))
plot(ts(exp(DATA1)),xlab="Iteration count", main="Trace plots for log(alpha) and log(beta) using initial variance matrix ")
effectiveSize(DATA1)


#run 2 using optimal variance matrix and epsilon is 0.5 for normalised euclidean distance 
SIG2<- ((2.38^2)/2)*var(DATA1)
SIG2
system.time(DATA2<-ABCMCMC(50000,y,763,14,0.1,SIG2,c(-6,log(0.5)),c(762,1),matrix(c(-1,0,1,-1),nrow=2,ncol=2,byrow=TRUE),1))
plot(ts(exp(DATA2)),xlab="Iteration count", main="Trace plots for log(alpha) and log(beta) using optimal tuning matrix and 50000 runs")
effectiveSize(DATA2)

DATA2<-exp(DATA2)
mean(DATA2[,1])
sd(DATA2[,1])
quantile(DATA2[,1],c(0.025,0.975))


## summaries of beta

mean(DATA2[,2])
sd(DATA2[,2])
quantile(DATA2[,2],c(0.025,0.975))

#posterior plot comparing three abc methods 
par(mfrow=c(1,1))
y<-c(1,3,6,25,73,221,294,257,236,189,125,67,26,10,3)
t<-c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14)


posteriorIt<-function(simulation,colour,iterations){
  I.mat = matrix(0,nrow=14/0.1+1,ncol=iterations)
  for(i in 1:iterations){
    data<-PLvalues(c(762,1),763,14,simulation[i,1:2],matrix(c(-1,0,1,-1),nrow=2,ncol=2,byrow=TRUE),0.1)
    I.mat[,i]<-data[,2]
  }
  I.mean <- apply(I.mat,1,mean)
  I.lq <- apply(I.mat,1,quantile,0.025)
  I.uq <- apply(I.mat,1,quantile,0.975)
  lines(ts(I.mean,start=0,deltat=0.1),lwd=2,col=colour)
  lines(ts(I.lq,start=0,deltat=0.1),col=colour)
  lines(ts(I.uq,start=0,deltat=0.1),col=colour)
}
plot(x=t,y=y,type="p",ylim = c(0,400))
posteriorIt(norm2,1,1000)
posteriorIt(DATA2,2,50000)
posteriorIt(DATA4,3,50000)
