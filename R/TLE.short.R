
#' Method TLE.
#' @name TLE
#' @rdname TLE-methods
#' @exportMethod TLE
#' @param formula A symbolic description of the model to be fit.
#' @param data A data frame containing the variables in the model.
#' @param tRatio An tRatio parameter for TLE method.
#' @param MaxIt An MaxIt parameter for TLE method.
#' @param nc An optional number of clusters.
##########################################################################################
##### setGeneric function TLE
##########################################################################################
setGeneric("TLE",
	function(formula,data, nc=2,tRatio,MaxIt=200)
	standardGeneric("TLE"))

#' @rdname TLE-methods
#' @aliases CTLE,formula,ANY,numeric,numeric,numeric-method
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
					#print(llik)
					if(abs(llik_old-llik)/llik < 10^(-8) | ccc>MaxIt){
						flag=1
					}
					#print(flag)
				}


		outliers=setdiff(1:nobs,inds_in)
                coffs_final=rbind(parameters(fres),apply(fres@posterior$scaled,2,sum)/length(inds_in))
                cates=vector("numeric",nobs)
                cates[inds_in]=fres@cluster
                cates[outliers]=-1
		pvals_final=sapply(refit(fres)@components[[1]], function(x)(x[-1,4]))
                result = new("RobMixReg", inds_in=inds_in,indout=outliers,ctleclusters=cates,compcoef=coffs_final,comppvals=pvals_final)
                result@call <- mycall
                result
	}
)
