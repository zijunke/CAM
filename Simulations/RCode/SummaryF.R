wd = '/home/ZijunKe/Research/MetaR/Simulation4/'
source(paste(wd,'RCode/RFuncs_Sim.R',sep=''))
nsim = 500

Nstudy.all = c(28,60,100)
mu.N.all = c(150,350)
L.All = cbind(rep(.6,8),rep(.8,8))
rho0.All = c(0,0.3,-0.5)
Phi.All = c(.1,.2)
sdL.All = c(0,0.06)
sdPhi.All = c(0,0.1)

allres.fname = paste(wd,'Results/AllRes.out',sep='')

SRi = 1
FID = 0
for(sdPhii in 1:1){
   sdLi = sdPhii
   for(RRi in 1:2){
   for(Nbari in 1:2){
   for(NSi in 1:3){
   for(Ri in 1:3){
   FID = FID +1 

   Nstudy = Nstudy.all[NSi]
   mu.N = mu.N.all[Nbari]
   senario.n = paste('VPhi',sdPhii-1,'VL',sdLi-1,sep='')
   cond.n = paste('R',Ri,'SR',SRi,'RR',RRi,sep='')
   work.d = paste(wd,'Results/',senario.n,'/',cond.n,'/',Nstudy,mu.N,'/',sep='')
   print(paste(cond.n,'Nstudy',Nstudy,'mu.N',mu.N,sep=''))
 

   T.Values = list(
   T.Values1 = list(rho0 = rho0.All[Ri],V.rho=0.04),
   T.Values2 = list(rho0 = rho0.All[Ri],V.rho=0.04),
   T.Values3 = list(rho0 = rho0.All[Ri],V.rho=0.04,sd.rho = 0.2, Phi = rep(Phi.All[SRi],3))	)

   prm <- list(prm0 = c('rho0','b','V.rho'),
   prm1 = c('rho0','V.rho'),
   prm2 = c('rho0','V.rho'),
   prm3 = c('rho0','V.rho','sd.rho','Phi')	)

   out.fn = c(
   paste(work.d,Nstudy,mu.N,'.M1',sep=''),
   paste(work.d,Nstudy,mu.N,'.M2',sep=''),
   paste(work.d,Nstudy,mu.N,'.M3',sep='') )

   AllRes = vector('list',3)
   for(mi in 1:3){ AllRes[[mi]] = summary.s(mi,out.fn[[mi]],T.Values[[mi]],prm[[mi]])	}
   AllRes2 = cbind(t(AllRes[[1]]),t(AllRes[[2]]),t(AllRes[[3]][1:3,]))
   print(round(AllRes2,3))
   write.table(round(AllRes2,3),file = allres.fname,row.names=FALSE,col.names = FALSE,append = TRUE)

   }
   }
   }
   }	

}
