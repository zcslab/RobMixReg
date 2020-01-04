
library(robust)
library(flexmix)
library(robustbase)
flexmix_2<-function(formula,data1,k,mprior){
        liks=NULL
        mix_list=list()
        niter=10
	ii=1
	flag=0
	while(ii<niter & flag==0){
		tmp_mix = try(flexmix(formula,data1,k=k,control = list(minprior = mprior)),silent=TRUE)
                if(!inherits(tmp_mix, "try-error")){
			if(tmp_mix@k==k){
				flag=1
			}
                }
		ii=ii+1
        }
        return(tmp_mix)
}


setClass("RobMixReg",
	representation(inds_in="numeric",
	indout="ANY",
	ctleclusters="ANY",
	compcoef="ANY",
	comppvals="ANY",
	call="call"))

##########################################################################################
##### setGeneric function TLE
##########################################################################################
setGeneric("TLE",
	function(formula,data, nc=2,tRatio,MaxIt=200)
	standardGeneric("TLE"))


### #######################################################################################
### #######################################################################################
### #######################################################################################
### ## setMethod TLE
### #######################################################################################
### if called with existing TLE object (as result=), then restart estimate...
### #######################################################################################
setMethod("TLE",
	signature(formula="formula",data="ANY",nc="numeric",tRatio="numeric",MaxIt="numeric"),
	function(formula,data, nc=2,tRatio,MaxIt=200)
	{
		mycall = match.call();
		nx=ncol(data)-1; nobs = nrow(data)            # number of observations
			flag=0;
				n_keep=nobs*(1-tRatio)
				inds_in=sample(1:nobs,n_keep)
				llik=-Inf
				ccc=0
				while(flag==0){
					ccc=ccc+1
					llik_old=llik
					fres = flexmix_2(formula,data1=data[inds_in,],k=nc,mprior=0)            #
					nc_tmp=fres@k
					www = posterior(fres, newdata=data, unscale=TRUE)
					lll=rowSums(www)
					inds_in=sort(order(lll,decreasing=TRUE)[1:n_keep])
					llik=sum(lll[inds_in])
					print(llik)
					if(abs(llik_old-llik)/llik < 10^(-8) | ccc>MaxIt){
						flag=1
					}
					print(flag)
				}


		outliers=setdiff(1:nobs,inds_in)
                coffs_final=rbind(parameters(fres),apply(fres@posterior$scaled,2,sum)/length(inds_in))
                cates=vector("numeric",nobs)
                cates[inds_in]=fres@cluster
                cates[outliers]=-1
		pvals_final=sapply(summary(refit(fres))@components[[1]], function(x)(x[-1,4]))
                result = new("RobMixReg", inds_in=inds_in,indout=outliers,ctleclusters=cates,compcoef=coffs_final,comppvals=pvals_final)
                result@call <- mycall
                result
	}
)
