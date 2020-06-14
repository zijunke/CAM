library(ZIM) # for generate sample sizes within primary studies
library(metafor) # standard meta-analyses
library(xlsx) # read in .xlsx files
library(lavaan) # for CFA analysis
library(MASS) # for bayesian sampling
library(hbmem) # for bayesian sampling
library(truncnorm) # for bayesian sampling
library(coda) # for bayesian sampling

#wd <- "E:/MetaR/Simulation3/"
wd = '/home/ZijunKe/Research/MetaR/Simulation4/'
source(paste(wd,'RCode/RFuncs_Sim.R',sep=''))

# Simulation settings
nsim = 500
SZ = list(Nstudy = c(28,60,100),muN= c(150,350))
rho = list(v = c(0,0.3,-0.5),sdv = 0.2)
Phi = list(v = c(0.1,0.2), sdv = c(0,0.1),
	ind = c(3,3,3,3,3,2,2,2,4,4,4,1,4,1,1,4,2,4,2,3,2,1,4,1,4,4,4,1))
L = list(v = cbind(rep(0.6,8),rep(0.8,8)),sdv = c(0,0.06),
	ind = matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,2,1,1,2,1,2,1,1,1,1,2,1,1,1,2,1,
		5,5,5,5,5,4,4,4,6,8,7,3,3,3,3,3,4,5,4,5,4,3,4,3,6,6,4,3),28,2))
repStudy = vector('list',length(SZ$Nstudy))
repStudy[[1]] = 1:28
set.seed(1037302)
for(i in 2:length(SZ$Nstudy)){
	repStudy[[i]] = sample(1:nrow(L$ind),SZ$Nstudy[i],replace = T)
}

# Values for reliability
# Setting for CFA study
Nrr = 500
Frho.list = list(v=rep(0.7,3),sdv = rep(0.1,3),
	M.ind = matrix(c(4,1,2,1,4,3,2,3,4),3,3,byrow = F))
L.list = list(v=rep(NA,6),sdv = rep(NA,6),
	M.ind = matrix(c(1,7,7,7,2,7,7,7,3,4,7,7,7,5,7,7,7,6),6,3,byrow = T))
Phi.list = list(v=c(0.25,0.14,0,0.23,0.17,0.21),sdv = rep(0.05,6),
	M.ind = matrix(c(8,1,2,7,7,7,1,9,3,7,7,7,2,3,10,7,7,7,
	  7,7,7,11,4,5,7,7,7,4,12,6,7,7,7,5,6,13),6,6,byrow = T))
rr.model.SWB <- list(Frho = Frho.list,L = L.list,Phi = Phi.list)

Frho.list = list(v=0.06,sdv = 0.1,
	M.ind = matrix(c(2,1,1,2),2,2,byrow = F))
L.list = list(v=rep(NA,8),sdv = rep(NA,8),
	M.ind = matrix(c(1,9,9,2,3,9,9,4,5,9,9,6,7,9,9,8),8,2,byrow = T))
tmp = matrix(5,8,8); diag(tmp) = 6:13
tmp[1,2] = tmp[2,1] = 1;tmp[3,4] = tmp[4,3] = 2
tmp[5,6] = tmp[6,5] = 3;tmp[7,8] = tmp[8,7] = 4
Phi.list = list(v=c(0.12,0.13,0.2,0.14),sdv = rep(0.05,4),
	M.ind = tmp)
rr.model.Personality <- list(Frho = Frho.list,L = L.list,Phi = Phi.list)
rr.model = list(Personality = rr.model.Personality,SWB = rr.model.SWB)
rr.model$Nrr = Nrr

CFAM.Personality <- '
	E =~ L1*V1 + L2*V3 + L3*V5 + L4*V7
	Nn =~ L5*V2 + L6*V4 + L7*V6 + L8*V8
	V1 ~~ V2;V3 ~~ V4;V5 ~~ V6;V7 ~~ V8;
'
CFAM.SWB <- '
	LS =~ L1*V1 + L4*V4
	PA =~ L2*V2 + L5*V5
	NAn =~ L3*V3 + L6*V6
	V1 ~~ V2; V1 ~~ V3;V2 ~~ V3
	V4 ~~ V5; V4 ~~ V6;V5 ~~ V6
	V1 ~~ VE1*V1; V4 ~~ VE1*V4
	V2 ~~ VE2*V2; V5 ~~ VE2*V5
	V3 ~~ VE3*V3; V6 ~~ VE3*V6
'
CFAModel <- list(Personality = CFAM.Personality,SWB = CFAM.SWB)

#sdPhii = sdLi = sdLi = Ri = SRi = rri = NSi = Nbari = 1
#---------------------------------------------
# Simulation: Data Generation
# 2(Variance of Relibility)*3(Effect Size)*2(Reliability)*2(Nstudy)*2(mu.N)
#---------------------------------------------------------------------
for(sdPhii in 1:2){#1:2
sdLi = sdPhii
wdsubi = paste(wd,'Data/VPhi',sdPhii-1,'VL',sdLi-1,'/',sep='')
try(system(paste('mkdir ',wdsubi,sep='')))

for(Ri in 1:3){#1:3
for(SRi in 1:1){
for(rri in 1:2){#1:2
   dir = paste(wdsubi,'R',Ri,'SR',SRi,'RR',rri,'/',sep='')
   newf = paste('mkdir ',dir,sep='')
   try(system(newf))
   setwd(dir)
   
   rr.model$Personality$L$v = rep(L$v[1,rri],8)
   rr.model$Personality$L$sdv = rep(L$sdv[sdLi],8)
   rr.model$SWB$L$v = L$v[3:8,rri]
   rr.model$SWB$L$sdv = rep(L$sdv[sdLi],6)
   Prr = list(Personality = getPrr(rr.model$Personality),SWB = getPrr(rr.model$SWB))
   rr.model$Prr = Prr
		
   for(NSi in 1:3){ #1:3
	indL = L$ind[repStudy[[NSi]],]
	indPhi = Phi$ind[repStudy[[NSi]]]

	for(Nbari in 1:2){ #1:2
	   filename = c(paste(SZ$Nstudy[NSi],SZ$muN[Nbari],'.dat',sep=''),
		paste(SZ$Nstudy[NSi],SZ$muN[Nbari],'.N.dat',sep=''),
		paste(SZ$Nstudy[NSi],SZ$muN[Nbari],'.rr.dat',sep=''),
		paste(SZ$Nstudy[NSi],SZ$muN[Nbari],'.valL.dat',sep='')	)

		meta.model = list(SZ = list(Nstudy = SZ$Nstudy[NSi],muN = SZ$muN[Nbari]),
			Phi = list(v = rep(Phi$v[SRi],3),sdv = rep(Phi$sdv[sdPhii],3),ind = indPhi),
			rho = list(v = rho$v[Ri],sdv = rho$sdv),indL = indL)
		
	   for(simi in 1:nsim){
		print(paste('simi = ',simi,sep=''))				
		GenData(simi,meta.model,rr.model,CFAModel,filename)
	   }
	}
   }
}
}	
}
}
