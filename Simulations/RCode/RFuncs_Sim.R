# -------------------------------------
# Implement the proposed method CAM
# the wrapper function
wrbugs <- function(ri,ni,mu.L,V.L,ind.L,gphi,initials,nburnin=1,niter=30001,nthin=5,prior){
	# ri: vector of observed correlations from primary studies
	# mu.L: estimated vector of factor loadings for various measurement methods
	# V.L: sampling covariance matrix of the estimated loading matrix
	# ind.L: matrix of indices indicating which two measurement methods are used in the primary study
	# gphi: vector of indices indicating nonzero residual correlations
	# ni: within-study sample sizes
	# initials: list of initial values
	# nburnin: number of iterations for the burnin period
	# niter: number of iterations
	# nthin: thinning interval

	#ri = r;mu.L=vL.obs; ind.L=indL;gphi=indP;ni=N;initials=inits;nburnin=1;niter = 30001;nthin=5
	tryi = 1
	fit = cam.bugs(ri,ni,mu.L,V.L,ind.L,gphi,initials,niter,prior)
	fit = cam.bugs(ri,ni,mu.L,V.L,ind.L,gphi,fit$initials,niter,prior)
	for(tryi in 3:20){
		fit.coda <- mcmc(data = fit$mcmc.chain,start = nburnin+1,end = niter,thin = nthin)
		conv = geweke.diag(fit.coda)[[1]]
		if(sum(abs(conv)>1.96)==0){
			break;
		}else{
			#print(fit$initials)
			fit = cam.bugs(ri,ni,mu.L,V.L,ind.L,gphi,fit$initials,niter,prior)
		}
	}
	return(list(tryi=tryi,fit.coda = fit.coda,fit = fit,conv = conv))
}

# The core function that implement the bayesian sampling algorithm for the proposed method
cam.bugs <- function(ri,ni,mu.L,V.L,ind.L,gphi,initials,niter,prior){
   # uninformative prior is used
   # ri: vector of observed correlations from primary studies
   # mu.L: estimated vector of factor loadings for measurement methods
   # V.L: sampling covariance matrix of the estimated loading matrix
   # ind.L: matrix of indices indicating which measurement methods are used for the two factors
   # gphi: vector of indices indicating nonzero residual correlations
   # ni: within-study sample sizes
   # initials: list of initial values
   # niter: number of iterations

   #ri = r;mu.L=vL.obs; ind.L=indL;gphi=indP;ni = N;initials = inits;niter = 30001
   Nstudy = length(ri)
   gphi.value = sort(unique(gphi))[-1]
   nphi = length(gphi.value)
   Iphii = matrix(0,Nstudy,nphi)
   for(i in 1:nphi){
	Iphii[which(gphi==gphi.value[i]),i] = 1
   }
	
   rho0 = c(initials$rho0,rep(NA,niter))
   V.rho = c(initials$V.rho,rep(NA,niter))
   sd.rho = c(sqrt(initials$V.rho),rep(NA,niter))
   Phi = cbind(initials$Phi,matrix(NA,nphi,niter))
   rhoi = cbind(initials$rhoi,matrix(NA,Nstudy,niter))
   V.ri = (1-(ri^2))^2/(ni-1)

   mup.Phi = prior$Phi$mu
   Vp.Phi = prior$Phi$sigma
	
   for(bi in 1:niter){
	vL = my.mvrnorm(mu.L,V.L)
	vL = apply(matrix(vL,ncol=1),1,function(x) max(min(x,1),-1))
	rri = vL[ind.L[,1]]*vL[ind.L[,2]]
		
	# update V.rho
	# truncated gamma 
	a = Nstudy/2 - 1
	b = sum((rhoi[,bi]-rho0[bi])^2)/2
 	V.rho[bi+1] = max(0.00001,1/rtgamma(1,shape = a,scale = 1/b,a = 1,b=10000))
	sd.rho[bi+1] = sqrt(V.rho[bi+1]) 
		
	# update rho0 
	# truncated normal 
	rho0[bi+1] = rtruncnorm(1,a=-1,b=1,mean=mean(rhoi[,bi]),sd=sqrt(V.rho[bi+1]/Nstudy))
		
	# update Phi
	# truncated normal
	sig = 1/(t(Iphii)%*%(1/V.ri) + 1/Vp.Phi)
	mu = (t(Iphii)%*%((ri-rri*rhoi[,bi])/V.ri) + mup.Phi/Vp.Phi)*sig
	Phi[,bi+1] = apply(cbind(mu,sqrt(sig)),1,function(x) rtruncnorm(1,a=-1,b=1,mean=x[1],sd=x[2]))
		
	# update rhoi
	# unconstrained
	numerator = rho0[bi+1]/V.rho[bi+1] + rri*(ri-Iphii%*%Phi[,bi+1])/V.ri
	denominator =  1/V.rho[bi+1] + (rri^2)/V.ri
	mu = numerator/denominator
	sig = 1/denominator
	rhoi[,bi+1] = apply(cbind(mu,sqrt(sig)),1,function(x) rnorm(1,mean=x[1],sd=x[2]))
   }
	
   simiseq = data.frame(rho0=rho0[-1],V.rho=V.rho[-1],sd.rho=sd.rho[-1],Phi=t(Phi[,-1]))
   last.values = list(rho0=rho0[niter+1],V.rho = V.rho[niter+1],sd.rho = sd.rho[niter+1],
			Phi = Phi[,niter+1],rhoi = rhoi[,niter+1]) 
   return(list(mcmc.chain=simiseq,initials = last.values))
}

# Function that handles the square root of reliability esitmates
# if we are 100% confident about our reliability estimates
# we can set their variances to be zero
my.mvrnorm <- function(m,V){
	s = m
	sel0 = which(diag(V)==0)
	n0var = length(sel0)
	if(n0var==0){ # all reliability variances > 0
		s = mvrnorm(1,m,V)	
	}else if( n0var < length(s) ){ # some reliability variances > 0
		seln0 = which(diag(V)>0)
		s[seln0] = mvrnorm(1,m[-sel0],V[-sel0,-sel0])	
	}
	return(s)
}

#--------------------------------------------------
# Data generation for the simulation

# Transform the mean and sd of the normal distribution for true effect sizes (correlation)
# to those of a truncated normal
# The resulting mean and sd can be quite different from the original values when 
# the true average correlaiton is large and the heterogeneity is substantial
# Formulus are from wiki https://en.wikipedia.org/wiki/Truncated_normal_distribution
msd.Norm2Trunc <- function(pars,ll,ul){
   # ll: lower limit of the truncated normal
   # ul: upper limit of the truncated normal
   mu <- pars[1] # mean of the normal distribution
   sd <- pars[2] # sd of the normal distribution

   Zll <- (ll-mu)/sd
   Zul <- (ul-mu)/sd
   Z <- max(pnorm(Zul)-pnorm(Zll),.00001)
	
   mu.trunc <- mu + sd*(dnorm(Zll)-dnorm(Zul))/Z
   sd.trunc <- sd*sqrt(1+(Zll*dnorm(Zll)-Zul*dnorm(Zul))/Z-((dnorm(Zll)-dnorm(Zul))/Z)^2)
	
   return(c(mu.trunc,sd.trunc))
}

# From vector to matrix 
v2m <- function(v,Mind){
	M = Mind
	for(i in 1:ncol(Mind)){	M[,i] = v[Mind[,i]] }
	return(M)
}

# Simulate Population parameter matrix
# v: vector of parameter values
# sdv: vector of between-study heterogeneity of the parameter
# M: Parameter matrix
# M.ind: indicator matrix to transform v to M
# ex.values: extra values in the parameter matrix; usually 0 or 1
simPar <- function(model.list,ex.values=0){
   v = model.list$v
   sdv = model.list$sdv
   M.ind = model.list$M.ind
   seln0 = which(sdv>0)
   n.n0 = length(seln0)
   if(n.n0>0){
	if(n.n0==1){
	   vs <- msd.Norm2Trunc(c(v,sdv),-1,1)
	   v = rtruncnorm(1,a=-1,b=1,mean=vs[1],sd=vs[2])
	}else if(n.n0 == length(v)){
	   ms = t(apply(cbind(v,sdv),1,msd.Norm2Trunc,ll=-1,ul=1))
	   v = apply(ms,1,function(x) rtruncnorm(1,a=-1,b=1,mean=x[1],sd=x[2]))
	}else{
	   ms.sub = t(apply(cbind(v[seln0],sdv[seln0]),1,msd.Norm2Trunc,ll=-1,ul=1))
	   vsub = apply(ms.sub,1,function(x) rtruncnorm(1,a=-1,b=1,mean=x[1],sd=x[2]))
	   v[seln0] = vsub
	}
   }
   if(is.null(M.ind) == 0){
	M = v2m(c(v,ex.values),M.ind)
   }else{M = NULL}
   return(list(M=M,v=v))
}

# To obtain population correlation matrix for reliability estimation
# vFrho: vector of factor correlations
# sdFrho: between-study heterogeneity of factor correlations
# Frho: Factor correlation matrix
# Frho.ind: indicator matrix to transform vFrho to mFrho
# vL: vector of factor loadings
# sdL: between-study heterogeneity of factor loadings
# L: Factor loading matrix
# L.ind: indicator matrix to transform vL to L
# vPhi: vector of factor loadings
# sdL: between-study heterogeneity of factor loadings
# L: Factor loading matrix
# L.ind: indicator matrix to transform vL to L
getPrr <- function(rr.model){

   Frho = simPar(rr.model$Frho,1)$M # generate factor correlation matrix
   L = simPar(rr.model$L,0)$M # generate factor loading matrix
   P <- L%*%Frho%*%t(L)
   varP = diag(P)

   for(tryi in 1:50){ # generate residual correlation matrix
	Phi = simPar(rr.model$Phi,c(0,1-varP))$M
	if(sum(eigen(Phi)$values>0)==nrow(Phi)){break}
   }

   P.rr <- P + Phi # population correlation matrix
   P.rr = as.matrix(nearPD(P.rr)$mat)
   return(P.rr)
}

# Generate observed squre root of reliabilities (factor loadings)
estL <- function(P.rr,myModel,N=500,extract.names){

   p = nrow(P.rr)
   rr.data = as.data.frame(mvrnorm(N,rep(0,p),P.rr))
   colnames(rr.data) = paste('V',1:p,sep='')
   nL = length(extract.names)

   # model with residual correlations; 	
   cfa.res = seq.CFA(rr.data,myModel,extract.names)
   est 	   = coef(cfa.res$fit)
   Vest    = try(vcov(cfa.res$fit))
   if(inherits(Vest,'try-error')==0){
	Lhat   = est[extract.names]
	VL     = as.matrix(nearPD(Vest[extract.names,extract.names])$mat)
	converged = 1
   }else{
	Lhat = rep(NA,nL)
	VL = matrix(N,nL,nL)
	converged = 0
   }
   return(list(Lhat = Lhat, VL = VL,converged = converged))
}

seq.CFA <- function(d,myM,extract.names){
   n.M <- length(myM)
   for(mi in 1:n.M){
	fit <- cfa(model = myM[[mi]],data = d,std.lv = TRUE)
	crit1 = fit@optim$converged
	Lval <- c(coef(fit)[extract.names])
	crit2 = sum(abs(Lval)>1)
	Vest <- try(vcov(fit))
	crit3 = abs(inherits(Vest,'try-error')-1)
	crit4 = (sum(eigen(inspect(fit,"theta"))$values<0)==0)
	crit.all = (crit1 == 1)*(crit2 == 0)*(crit3==1)*(crit4==1)
	if(crit.all == 1){break}
   }
   crit12 = (crit1 == 1)*(crit2 == 0) 
   if(crit12==0){mi = mi+1}
   return(list(mi = mi,fit = fit))
}

# meta.model: simulation setting for meta-analysis
# rr.model: model specification for reliability generation
GenData <- function(simi,meta.model,rr.model,CFAModel,filename){

   Nstudy = meta.model$SZ$Nstudy
   mu.N = meta.model$SZ$muN
   indL = meta.model$indL
   indPhi = meta.model$Phi$ind
   Nrr = rr.model$N
   Prr = rr.model$Prr

   # Generate sample sizes per study
   N <- rzinb(n = Nstudy, k = 0.4, lambda = mu.N*0.7, omega = 0)
   N <- N + mu.N*0.3

   # Generate true individual study reliability
   vL.Per = simPar(rr.model$Personality$L,0)$v[c(1,3)]
   vL.SWB = simPar(rr.model$SWB$L,0)$v
   vL = c(vL.Per,vL.SWB)
   rr = vL[indL[,1]]*vL[indL[,2]]
   
   # Generate effect sizes
   vtmp = msd.Norm2Trunc(c(meta.model$rho$v,meta.model$rho$sdv),-1,1)
   rhoi = rtruncnorm(Nstudy,-1,1,mean=vtmp[1],sd=vtmp[2]) # individual study true correlations
   rhoi.rr = rhoi*rr	# attenuated population correlations

   Phi.values = c(simPar(meta.model$Phi)$v,0)
   vPhi = Phi.values[indPhi] # individual study residual correlations
   
   rho.manifest <- rhoi.rr + vPhi # population correlation between manifest variables
   # observed correlations for individual primary studies
   r.data = rep(NA,Nstudy)
   for(si in 1:Nstudy){
   	Sigma = matrix(rho.manifest[si],2,2)
	diag(Sigma) = 1
	Sigma = as.matrix(nearPD(Sigma)$mat)
	data = mvrnorm(N[si],rep(0,2),Sigma) # generated data for primary studies
	r.data[si] = cor(data)[2,1] # sampling error
   }
   r.data = round(r.data,3)
  
   # Generate reliability estimates
   L.personality = estL(Prr$Personality,CFAModel$Personality,Nrr,c('L1','L3'))
   L.SWB = estL(Prr$SWB,CFAModel$SWB,Nrr,paste('L',1:6,sep=''))
   vL.est = c(L.personality$Lhat,L.SWB$Lhat)
   VL.est = rbind(cbind(L.personality$VL,matrix(0,2,6)),cbind(matrix(0,6,2),L.SWB$VL))
   rr.est = round(vL.est[indL[,1]]*vL.est[indL[,2]],3)
   LvalV = c(vL.est,VL.est[lower.tri(VL.est,diag = T)])

   write.table(t(c(simi,r.data)),file = filename[1],append=T,quote=F,row.names=F,col.names=F)
   write.table(t(c(simi,N)),file = filename[2],append=T,quote=F,row.names=F,col.names=F)
   write.table(t(c(simi,rr.est)),file=filename[3],append=T,quote=F,row.names=F,col.names=F)
   write.table(t(c(simi,LvalV)),file = filename[4],append=T,quote=F,row.names=F,col.names=F)
}

v2m.Lambda <- function(vVal,p){
	vL <- vVal[1:p]
	vV.L <- vVal[-c(1:p)]
	mV.L <- matrix(0,p,p)
	mV.L[lower.tri(mV.L,diag=T)] <- vV.L
	mV.L <- mV.L + t(mV.L)
	diag(mV.L) <- diag(mV.L)/2
	mV.L <- as.matrix(nearPD(mV.L)$mat)
	return(list(v = vL,V = mV.L))
}

org.res <- function(simi,res,out.fn,T.Values,prm){
	fit  = summary(res$fit.coda)
	tryi = res$tryi
	conv = res$conv

	est <- c(simi,unlist(fit$quantiles[,'50%']))
	sd <- c(simi,unlist(fit$statistics[,'SD']))
	np = length(unlist(T.Values[prm]))
	CI = cbind(HPDinterval(res$fit.coda,prob = .95),unlist(T.Values[prm]))
	cover = c(simi,apply(CI,1,function(x) (x[1]< x[3]) * (x[3] < x[2])))
	CI[,3] = rep(0,np)
	sig = c(simi,1-apply(CI,1,function(x) (x[1]< x[3]) * (x[3] < x[2])))
	conv.ex <- c(simi,tryi,conv)
	
	write.table(t(est),paste(out.fn,'.Est.out',sep=''),append = T,
			quote = F,row.names = F,col.names = F)
	write.table(t(sd),paste(out.fn,'.SD.out',sep=''),append = T,
			quote = F,row.names = F,col.names = F)
	write.table(t(cover),paste(out.fn,'.CI.out',sep=''),append = T,
			quote = F,row.names = F,col.names = F)
	write.table(t(CI[,1]),paste(out.fn,'.CIL.out',sep=''),append = T,
			quote = F,row.names = F,col.names = F)
	write.table(t(CI[,2]),paste(out.fn,'.CIU.out',sep=''),append = T,
			quote = F,row.names = F,col.names = F)
	write.table(t(sig),paste(out.fn,'.Sig.out',sep=''),append = T,
			quote = F,row.names = F,col.names = F)
	write.table(t(conv.ex),paste(out.fn,'.Conv.out',sep=''),append = T,
			quote = F,row.names = F,col.names = F)
	
}

org.res.meta4 <- function(simi,res,out.fn,T.Values){
	np = sum(unlist(lapply(T.Values,length)))
	est = c(simi,coef(res),res$tau2)
	se =  c(simi,res$se,res$se.tau2)
	CI = cbind(res$ci.lb,res$ci.ub)
	CI = rbind(CI,c(res$tau2-1.96*res$se.tau2,res$tau2+1.96*res$se.tau2))
	CI = cbind(CI,unlist(T.Values))
	cover = c(simi,apply(CI,1,function(x) (x[1]< x[3]) * (x[3] < x[2])))
	CI[,3] = rep(0,np)
	sig = c(simi,1-apply(CI,1,function(x) (x[1]< x[3]) * (x[3] < x[2])))
	
	write.table(t(est),paste(out.fn,'.Est.out',sep=''),append = T,
			quote = F,row.names = F,col.names = F)
	write.table(t(se),paste(out.fn,'.SD.out',sep=''),append = T,
			quote = F,row.names = F,col.names = F)
	write.table(t(cover),paste(out.fn,'.CI.out',sep=''),append = T,
			quote = F,row.names = F,col.names = F)
	write.table(t(CI[,1]),paste(out.fn,'.CIL.out',sep=''),append = T,
			quote = F,row.names = F,col.names = F)
	write.table(t(CI[,2]),paste(out.fn,'.CIU.out',sep=''),append = T,
			quote = F,row.names = F,col.names = F)
	write.table(t(sig),paste(out.fn,'.Sig.out',sep=''),append = T,
			quote = F,row.names = F,col.names = F)
}

org.resNA <- function(simi,out.fn,T.Values,prm){
	np = sum(unlist(lapply(T.Values[prm],length)))
	write.table(t(rep(NA,np+1)),paste(out.fn,'.Est.out',sep=''),append = T,
			quote = F,row.names = F,col.names = F)
	write.table(t(rep(NA,np+1)),paste(out.fn,'.SD.out',sep=''),append = T,
			quote = F,row.names = F,col.names = F)
	write.table(t(rep(NA,np+1)),paste(out.fn,'.CI.out',sep=''),append = T,
			quote = F,row.names = F,col.names = F)
	write.table(t(rep(NA,np+1)),paste(out.fn,'.CIL.out',sep=''),append = T,
			quote = F,row.names = F,col.names = F)
	write.table(t(rep(NA,np+1)),paste(out.fn,'.CIU.out',sep=''),append = T,
			quote = F,row.names = F,col.names = F)
	write.table(t(rep(NA,np+1)),paste(out.fn,'.Sig.out',sep=''),append = T,
			quote = F,row.names = F,col.names = F)
	write.table(t(c(simi,20,rep(NA,np))),paste(out.fn,'.Conv.out',sep=''),append = T,
			quote = F,row.names = F,col.names = F)
}

summary.s <- function(mi,out.fn,T.Values,prm){
	vTV = unlist(T.Values)

	res <- vector('list',4)
	names(res) <- c('est','sd','CI','conv')
	res$est <- as.matrix(read.table(paste(out.fn,'.Est.out',sep='')))
	res$sd <- as.matrix(read.table(paste(out.fn,'.SD.out',sep='')))
	res$CI <- as.matrix(read.table(paste(out.fn,'.CI.out',sep='')))
	res$sig <- as.matrix(read.table(paste(out.fn,'.Sig.out',sep='')))

	id.nconv = rep(0,nrow(res$est))
	if(mi > 2){
		res$conv <- as.matrix(read.table(paste(out.fn,'.Conv.out',sep='')))
		id.nconv <- which(res$conv[,2]>9)
	}
	if(sum(id.nconv)==0){
		#print(id.nconv)
		summary.r <- cbind(vTV,
			apply(res$est[,-1],2,mean,na.rm = T),
			apply(res$est[,-1],2,sd,na.rm = T),
			apply(res$sd[,-1],2,mean,na.rm = T),
			apply(res$CI[,-1],2,mean,na.rm = T),
			apply(res$sig[,-1],2,mean,na.rm = T))
	}else{
		summary.r <- cbind(unlist(T.Values),
			apply(res$est[-id.nconv,-1],2,mean,na.rm = T),
			apply(res$est[-id.nconv,-1],2,sd,na.rm = T),
			apply(res$sd[-id.nconv,-1],2,mean,na.rm = T),
			apply(res$CI[-id.nconv,-1],2,mean,na.rm = T),
			apply(res$sig[-id.nconv,-1],2,mean,na.rm = T))
	}
	if(mi==3){
		summary.r[2,6] = mean(abs(res$est[,3]/res$sd[,3])>1.96,na.rm=T)
		summary.r[3,6] = mean(abs(res$est[,4]/res$sd[,4])>1.96,na.rm=T)
	}
	colnames(summary.r) <- c('TrueValues','Est','ESE','ASE','CR','Sig')
	rownames(summary.r) <- names(vTV)
	print(summary.r)
	return(summary.r)
}

