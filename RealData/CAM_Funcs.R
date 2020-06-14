# A wrapper of the function cam
# Stop every niter iterations and test whether convergence is achieved
# The maximum number of iterations is niter*20
# The minimum number of iterations is niter*2
# Convergence is examined using the Geweke method
cam <- function(ri,ni,mu.L,V.L,ind.L,gphi,
	initials,nburnin=1,niter=30001,nthin=5,prior){
	# ri: observed correlations
	# ni: within-study sample sizes
	# mu.L: vector of factor loadings for various measures
	# V.L: sampling covariance matrix of the estimated loadings
	# ind.L: (Nstudy*2) matrix of indices indicating which two 
	#        measures are used in primary studies
	# gphi: vector of indices indicating nonzero residual correlations
	# initials: list of initial values
	# nburnin: number of iterations for the burnin period
	# niter: number of iterations
	# nthin: thinning interval
	# prior: prior for residual correlation parameters

	tryi = 1
	fit = cam.bugs(ri,ni,mu.L,V.L,ind.L,gphi,initials,niter,prior)
	fit = cam.bugs(ri,ni,mu.L,V.L,ind.L,gphi,fit$initials,niter,prior)
	for(tryi in 2:20){
		fit.coda <- mcmc(data = fit$mcmc.chain,start = nburnin+1,
			end = niter,thin = nthin)
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

# Main function that implements CAM
cam.bugs <- function(ri,ni,mu.L,V.L,ind.L,gphi,initials,niter,prior){

   Nstudy = length(ri)
   gphi.value = sort(unique(gphi))
   if(gphi.value[1]==0){gphi.value = gphi.value[-1]}
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
   
   # Gibbs sampling
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
	
   if(nphi == 1){Phis = Phi[-1];Phi.initial = Phi[niter+1] }
   if(nphi > 1){Phis = t(Phi[,-1]);Phi.initial = Phi[,niter+1] } 
   simiseq = data.frame(rho0=rho0[-1],V.rho=V.rho[-1],sd.rho=sd.rho[-1],Phi=Phis)
   last.values = list(rho0=rho0[niter+1],V.rho = V.rho[niter+1],sd.rho = sd.rho[niter+1],
		Phi = Phi.initial,rhoi = rhoi[,niter+1])
   return(list(mcmc.chain=simiseq,initials = last.values))
}


# Generate data from a multivariate normal with mean of m and covariance matrix of V
# V can be a zero matrix
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

# Orgnize results return by cam
out.res <- function(res,CI.alpha = 0.95){
	fit  = summary(res$fit.coda)
	tryi = res$tryi
	conv = res$conv

	est <- unlist(fit$quantiles[,'50%'])
	sd <- unlist(fit$statistics[,'SD'])
	CI = HPDinterval(res$fit.coda,prob = CI.alpha)
	CI = cbind(CI,rep(0,nrow(CI)))
	sig = 1-apply(CI,1,function(x) (x[1]< x[3]) * (x[3] < x[2]))
	conv.ex <- c(tryi,conv)
	return(list(Estimate= round(est,3),PosteriorSD= round(sd,3),
		CI = round(CI[,-3],3),Sig=sig,Convergence = round(conv.ex,3)))
}

# Compute the standard error and confidence interval for tau for results returned by metafor
seCI.metafor <- function(fit){
	se = fit.tau = sqrt(1/fit$tau2/4)*fit$se.tau2
	res = list(tau2_CI = c(fit$tau2-1.96*fit$se.tau2,fit$tau2+1.96*fit$se.tau2),
		tau_se = se,tau_CI = c(sqrt(fit$tau2)-1.96*se,sqrt(fit$tau2)+1.96*se))
	return(res)
}

# Impute missing cronbach's alpha with group mean
# alpha: vector of alpha with missing values
# type: vector indicating the type of measures used in the primary study
alpha.impute <- function(alpha,type){
   rr.m = aggregate(alpha, list(type), mean,na.rm = T)
   sel = match(type[is.na(alpha)],rr.m[,1])
   alpha[is.na(alpha)] = rr.m[sel,2]
   return(alpha)
}
