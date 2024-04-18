#### PROJECT: Mimulus cardinalis demography 2010-2014
#### PURPOSE: Create vital rate functions and build integral projection model to calculate lambda, population growth rate
#### AUTHOR: Seema Sheth
#### DATE LAST MODIFIED: 20171110

#*******************************************************************************
### 1. Build vital rate functions
#*******************************************************************************

### Survival probability function ###

s.x=function(x,params1) {
  u=exp(params1$surv.int+params1$surv.slope*x)
  return(u/(1+u))
}

### Growth function ###

g.yx=function(xp,x,params1) {
  dnorm(xp,mean=params1$growth.int+params1$growth.slope*x,sd=params1$growth.sd)
}

### Fecundity function ###

#flowering probability function
p.flower.x=function(x,params1) {
  u=exp(params1$flowering.int+params1$flowering.slope*x)
  return(u/(1+u))
}

f.yx=function(xp,x,params1) {
  p.flower.x(x,params1)*
    exp(params1$fruit.int+params1$fruit.slope*x)*
    params1$seed.ct*
    params1$establishment.prob*
    dnorm(xp,mean=params1$recruit.size.mean,sd=params1$recruit.size.sd) }

#*******************************************************************************
### 2. Combine vital rate functions to build discretized IPM kernel 
#*******************************************************************************

min.logSize=0.9*min(c(data1$logSize,data1$logSizeNext),na.rm=T)
max.logSize=1.1*max(c(data1$logSize,data1$logSizeNext),na.rm=T) # NOTE: This includes full size range across entire dataset; alternatively, I could use site-specific size range
n=100 # number of cells in matrix
tol = 1.e-8 # tolerance for iterations 
b=min.logSize+c(0:n)*(max.logSize-min.logSize)/n # boundary points 
y=0.5*(b[1:n]+b[2:(n+1)]) # mesh points 
h=y[2]-y[1] # step logSize

### Create IPM matrix ###

G<-array(0,dim=c(n,n))
G=h*outer(y,y,g.yx,params=params1) # growth matrix 

#create IPM matrix
S=s.x(y,params=params1) # survival
F=h*outer(y,y,f.yx,params=params1) # reproduction matrix 
P<-array(0,dim=c(n,n))
P=G # placeholder; redefine P on the next line 
for(i in 1:n) P[,i]=G[,i]*S[i] # growth/survival matrix 

#fix eviction of offpsring
for(i in 1:(n/2)) { G[1,i]<-G[1,i]+1-sum(G[,i]) 
P[,i]<-G[,i]*S[i] } 

#fix eviction of large adults 
for(i in (n/2+1):n) { G[n,i]<-G[n,i]+1-sum(G[,i]) 
P[,i]<-G[,i]*S[i] }

K=P+F # full matrix