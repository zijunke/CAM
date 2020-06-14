library(MASS)
library(coda)
library(Matrix)
library(metafor)
library(hbmem)
library(truncnorm)

wd = '/home/ZijunKe/Research/MetaR/Simulation4/'
#wd <- "E:/MetaR/Simulation3/"
openbugs.d = NULL
source(paste(wd,'RCode/RFuncs_Sim.R',sep=''))

nsim = 500
n.iter = 60000
n.burnin = 30000
Nstudy.all = c(28,60,100)
mu.N.all = c(150,350)

L.All <- cbind(rep(.6,4),rep(.8,4))
rho0.All <- c(0,0.3,-0.5)
Phi.All <- c(0.1,0.2)
sdL.All <- c(0,0.06)
sdPhi.All <- c(0,0.1)
ind.All <- list(L = matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,2,1,1,2,1,2,1,1,1,1,2,1,1,1,2,1,
	5,5,5,5,5,4,4,4,6,8,7,3,3,3,3,3,4,5,4,5,4,3,4,3,6,6,4,3),28,2),
	Phi = c(3,3,3,3,3,2,2,2,0,0,0,1,0,1,1,0,2,0,2,3,2,1,0,1,0,0,0,1))
resample.s = vector('list',length(Nstudy.all))
resample.s[[1]] = 1:28
set.seed(1037302)
for(i in 2:length(Nstudy.all)){
	resample.s[[i]] = sample(1:nrow(ind.All$L),Nstudy.all[i],replace = T)
}

sdPhii = 1
sdLi = sdPhii
Ri = 2
SRi = 1
RRi = 2
NSi = 1
Nbari = 2

Nstudy = Nstudy.all[NSi]
mu.N = mu.N.all[Nbari]
indL = ind.All$L[resample.s[[NSi]],]
indP = ind.All$Phi[resample.s[[NSi]]]

T.Values = list(
   T.Values1 = list(rho0 = rho0.All[Ri],V.rho=0.04),
   T.Values2 = list(rho0 = rho0.All[Ri],V.rho=0.04),
   T.Values3 = list(rho0 = rho0.All[Ri],V.rho=0.04,
   sd.rho = 0.2,Phi = rep(Phi.All[SRi],3)))

prm <- list(
   prm1 = c('rho0','V.rho'),
   prm2 = c('rho0','V.rho'),
   prm3 = c('rho0','V.rho','sd.rho','Phi'))

senario.n <- paste('VPhi',sdPhii-1,'VL',sdLi-1,sep='')
cond.n <- paste('R',Ri,'SR',SRi,'RR',RRi,sep='')
print(paste('Senario.Name=',senario.n,sep=''))
print(paste('Cond.Name=',cond.n,sep=''))
newf = paste('mkdir ',wd,'Results/',senario.n,'/',sep='')
try(system(newf))
newf = paste('mkdir ',wd,'Results/',senario.n,'/',cond.n,'/',sep='')
try(system(newf))

work.d <- paste(wd,'Results/',senario.n,'/',cond.n,'/',Nstudy,mu.N,'/',sep='')
newf = paste('mkdir ',work.d,sep='')
try(system(newf))

data.fn = c(
   paste(wd,'Data/',senario.n,'/',cond.n,'/',Nstudy,mu.N,'.dat',sep=''),
   paste(wd,'Data/',senario.n,'/',cond.n,'/',Nstudy,mu.N,'.N.dat',sep=''),
   paste(wd,'Data/',senario.n,'/',cond.n,'/',Nstudy,mu.N,'.rr.dat',sep=''),
   paste(wd,'Data/',senario.n,'/',cond.n,'/',Nstudy,mu.N,'.valL.dat',sep='') )
out.fn = c(
   paste(work.d,Nstudy,mu.N,'.M1',sep=''),
   paste(work.d,Nstudy,mu.N,'.M2',sep=''),
   paste(work.d,Nstudy,mu.N,'.M3',sep='') )

mR <- read.table(data.fn[1])[,-1]
mN <- read.table(data.fn[2])[,-1]
mrr <- read.table(data.fn[3])[,-1]
mVrr <- read.table(data.fn[4])[,-1]


for(simi in 1:nsim){
print(paste('simi =',simi,sep=''))
if(is.na(mrr[simi,1])==0){
   r = as.numeric(mR[simi,]) # ind study corrs
   N = as.numeric(mN[simi,]) # ind study sample size
   vi = ((1-r^2)^2)/(N-1)   # Sampling variance of correlaiton
   mr <- sum(r/vi)/sum(1/vi)	   # weighted average mean corrs
	
   rr.obs = as.numeric(mrr[simi,]) # reliabilities
   Linfo = v2m.Lambda(as.numeric(mVrr[simi,]),8)
   vL.obs = Linfo$v
   V.L = Linfo$V

   # Data analysis
   # bare-bone + HS
   res.meta4 = vector('list',2)
   res.meta4[[1]] =  rma(r,vi,weights = N,method  = 'HS' )
   res.meta4[[2]] =  rma(r/rr.obs,vi/(rr.obs^2),weights = N*(rr.obs^2),method  = 'HS' )
   for(mi in 1:2){	org.res.meta4(simi,res.meta4[[mi]],out.fn[mi],T.Values[[mi]])	}
	
   # Bayesian methods
   inits = list(rho0 = 0,V.rho = 0.09,Phi = rep(0,3),rhoi = r)
   prior = list(Phi = list(mu = rep(0,3),sigma = rep(100,3)))
   res = try(wrbugs(r,N,vL.obs,V.L,indL,indP,inits,nburnin=1,niter=10000,nthin=1,prior))
   if(inherits(res,'try-error')==0){
	org.res(simi,res,out.fn[3],T.Values[[3]],prm[[3]])
   }else{
	org.resNA(simi,out.fn[3],T.Values[[3]],prm[[3]])
   }	
}else{
   for(mi in 1:3){
	org.resNA(simi,out.fn[mi],T.Values[[mi]],prm[[mi]])
   }
}
}

print(paste('Nstudy = ',Nstudy,'Nbar = ', mu.N,sep=''))
for(mi in 1:3){
	summary.s(mi,out.fn[mi],T.Values[[mi]],prm[[mi]])
}



