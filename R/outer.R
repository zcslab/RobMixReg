
#' Robust Mixture Regression using four methods.
#'
#' @param lr.method A robust mixture regression method to be used. Should be one of "flexmix","TLE","CTLE","biSquare".
#' @param formula A symbolic description of the model to be fit. The general form is y~x|g where y is the response, x the set of predictors and g an optional grouping factor for repeated measurements.
#' @param data A data frame containing the variables in the model.
#' @param nc An optional number of clusters.
#' @param nit An nit parameter for CTLE method.
#' @param tRatio An tRatio parameter for TLE method.
#' @param MaxIt An MaxIt parameter for TLE method.
#' @param numini An numini parameter for biSauqre method.
#' @return An S4 object about the regression result.
#' @examples
#' \dontrun{
#' library(RobMixReg)
#' #library(robust)
#' library(flexmix)
#' library(robustbase)
#' library(MASS)
#' library(gtools)
#' # gaussData
#' x=(gaussData$x);y=as.numeric(gaussData$y);
#' formula01=as.formula("y~x")
#' example_data01=data.frame(x,y)
#' res_rmr = rmr(lr.method='flexmix', formula=formula01, data=example_data01)
#' res_rmr = rmr(lr.method='CTLE', formula=formula01, data=example_data01)
#' }
rmr <- function(lr.method="flexmix", formula=NULL, data=NULL, nc=2, nit=20, tRatio=0.05, MaxIt=200, numini=20)
{
  if(is.null(formula)) stop('Please input formula!')
  if(is.null(data)) stop('data is null! Please input data')


  method.lib <- c('flexmix', 'TLE', 'CTLE', 'biSquare','mixLp')
  if(!lr.method %in% method.lib) stop('lr.method must be chosen from "flexmix","TLE","CTLE","biSquare","mixLp" !  ')

  if(lr.method=='flexmix')
  {
    res = flexmix(formula,data,k=nc)
  }

  if(lr.method=='CTLE')
  {
    res = CTLE(formula,data, nit,nc)
  }

  if(lr.method=='TLE')
  {
    res = TLE(formula,data, nc,tRatio,MaxIt)
  }

	if(lr.method=='biSquare')
	{
    res = mixlinrb_bi(formula,data,nc, numini)
	}

  if(lr.method=='mixLp')
  {
    res = mixLp(formula, data, nc, numini)
  }

  return(res)
}
