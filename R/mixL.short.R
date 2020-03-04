
#' denLp : xxx.
#' @description xxx
#' @param rr The parameter.
#' @param sig The parameter.
#' @return Return value.
denLp=function(rr,sig){
        exp(-abs(sqrt(2)*(rr))/sig)/(sqrt(2)*sig)
}


#' mixLp_one function.
#' @description One step for the mixLp mixture regression. For more details, please see the 'mixLp' function.
#' @param formula A symbolic description of the model to be fit.
#' @param data A data frame containing the variables in the model.
#' @param nc The number of component. Default value is 2.
#' @description The robust EM algorithm to fit the mixture of linear regression based on bisquare function. Bai, X., Yao, W._, and Boyer, J. E. (2012). Robust fitting of mixture regression models.
#' @return xxx.
################
mixLp_one<-function(formula,data,nc=2){
	nx=ncol(data)-1; n = nrow(data)
	p=nx+1; n1=2*p;
	bet= matrix(rep(0,nc*p),nrow=nc);
	sig=0;
	xx=as.matrix(data[,1:nx,drop=FALSE])
	X=cbind(rep(1,n),xx);
	y=data[,nx+1]
	for(j in seq(nc)){
		ind1= sample(1:n,n1);
		tmp_mod=lm(formula,data=data[ind1,])
		bet[j,]=coefficients(tmp_mod)
		sig=sig+sum((resid(tmp_mod))^2)
	}
	pr=rep(1/nc,nc);sig=rep(sig/n1/nc,nc);#sig
sig0=sig;pr0=pr;bet0=bet;
sig=sig0;pr=pr0;bet=bet0
	r=matrix(rep(0,nc*n),nrow=n);
	pk=dk=r;
        for(j in seq(nc)){
		r[,j]=y-X%*%bet[j,]
        }

	run=0;acc=10^(-4)*max(abs(c(bet,sig,pr)));
#E-steps
	repeat{
		#print(run);	print(bet);
		prest=c(sig,bet,pr);run=run+1;
	        for(j in seq(nc)){
	                pk[,j]=pr[j]*denLp(r[,j],sig[j]);
	                dk[,j]=sig[j]/(sqrt(2)*abs(r[,j]));
	        }
		pk[pk<(10^(-6))]=10^(-6)
		pk=pk/matrix(rep(apply(pk,1,sum),nc),nrow=n);
	        dk[dk==Inf]=10^6;
#M-step
	        pr=apply(pk,2,mean);
	        for(j in seq(nc)){
			w = pk[,j]*dk[,j]
#			tmp_mod=lm(formula,data=data,weights=w)
#			bet[j,] = tmp_mod$coefficients;
#			r[,j]=resid(tmp_mod)
#	                if(any(is.na(bet[j,]))){
        	                W=diag(w);
                	        bet[j,]=ginv(t(X)%*%W%*%X)%*%t(X)%*%W%*%y
				r[,j]=y-X%*%bet[j,]
#	                }
	                sig[j]=sqrt(  2*sum(w*(r[,j])^2)/sum(pk[,j])    );
	        }

		dif=max(abs(c(sig,bet,pr)-prest))
		if(dif<acc|run>500){break}
	}
	theta=matrix(c(bet,sig,pr),nrow=nc);
	est=list(theta=theta,difpar=dif,run=run)
	return(est)
}


#' Method mixLp.
#' @name mixLp
#' @description mixLp estimates the mixture regression parameters robustly using bisquare function based on multiple initial values. The solution is found by the modal solution.
#' @rdname mixLp-methods
#' @exportMethod mixLp.
#' @param formula A symbolic description of the model to be fit.
#' @param data A data frame containing the variables in the model.
#' @param nc An optional number of clusters.
#' @param numini The minimal value of nu.
#' @usage mixLp(formula, data, nc=2, numini=200)
#' @return xxxxxx
##########################################################################################
##### setGeneric function mixLp
##########################################################################################
setGeneric("mixLp",
	function(formula,data, nc=2,numini=200)
	standardGeneric("mixLp"))

#' @rdname mixLp-methods
#' @aliases CTLE,formula,ANY,numeric,numeric-method
### #######################################################################################
### #######################################################################################
### #######################################################################################
### ## setMethod mixLp
### #######################################################################################
### if called with existing mixLp object (as result=), then restart estimate...
### #######################################################################################
setMethod("mixLp",
	signature(formula="formula",data="ANY",nc="numeric",numini="numeric"),
	function(formula,data, nc=2,numini=20)
	{
		mycall = match.call();

		nx=ncol(data)-1; n = nrow(data)
		p=nx+1; n1=2*p;
		perm=permutations(nc,nc);
		est=mixLp_one(formula,data,nc);
		lenpar=length(c(est$theta));
		theta= matrix(rep(0,lenpar*(numini)),ncol=lenpar);
		theta[1,]=c(est$theta);
		trimbet=matrix(theta[1,1:(p*nc)],nrow=nc);
		trimbet=matrix(rep(matrix(t(trimbet),ncol=p*nc,byrow=T),gamma(nc+1)),ncol=p*nc,byrow=T);
		ind=matrix(rep(0,numini),nrow=1);ind[1]=1;numsol=1;solindex=1; sol=matrix(theta[1,],nrow=1);
		for(i in 2:numini){
			est=mixLp_one(formula,data,nc)
			theta[i,]=est$theta;
			temp= matrix(theta[i,1:(p*nc)],nrow=nc);
			temp=matrix(t(temp[t(perm),]),ncol=p*nc,byrow=T);
			dif=apply((trimbet-temp )^2,1,sum);
			temp1=which(dif==min(dif));
			theta[i,]=c(c(matrix(temp[temp1[1],],nrow=nc,byrow=T)),theta[i,p*nc+perm[temp1[1],]],theta[i,p*nc+nc+perm[temp1[1],]]);
			dif= apply((matrix(rep(theta[i,1:(p*nc)],numsol),nrow=numsol,byrow=T)-sol[,1:(p*nc)])^2,1,sum);
			if(min(dif)>0.1){
				sol=rbind(sol,theta[i,]);
				numsol=numsol+1;
				solindex=c(solindex,i);
				ind=rbind(ind,rep(0,numini));
				ind[numsol,i]=1
			}else{
				ind1=which(dif==min(dif));
				ind[ind1,i]=1;
			}
		}
		num=apply(ind,1,sum);
		ind1=order(-num);
		bestindex=ind1;
		for(j in seq(numsol)){
			if(min(sol[ind1[j], (p*nc+nc+1):(p*nc+2*nc)])>0.05){
				index=1;
				est=matrix(sol[ind1[j],],nrow=nc);
				for(l in seq(nc-1)){
					temp=matrix(rep(est[l,1:p],nc-l),nrow=nc-l,byrow=T)-est[(l+1):nc,1:p];
					temp=matrix(temp,nrow=nc-l);
					dif=apply(temp^2,1,sum);
					if(min(dif)<0.1){
						index=0;break
					}
				}
				if(index==1){
					bestindex=ind1[j];break
				}
			}
		}
		est= sol[bestindex[1],];
		coffs_final = t(matrix(est,nrow=nc))
		cates = NA

                result = new("RobMixReg", inds_in=NA,indout=NA,ctleclusters=cates,compcoef=coffs_final,comppvals=NA)
                result@call <- mycall
                result
	}
)
