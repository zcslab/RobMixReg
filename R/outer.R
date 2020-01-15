
rmr <- function(lr.method=NULL, formula=NULL, data=NULL, nc=3, nit=20, tRatio=0.2, MaxIt=200, numini=20)
{
  if(is.null(formula)) stop('Please input formula!')
  if(is.null(data)) stop('data is null! Please input data')
  if(is.null(lr.method)) stop('Please input the regression method from ("flexmix","TLE","CTLE","biSquare"). ')

  method.lib <- c('flexmix', 'TLE', 'CTLE', 'biSquare')
  if(!lr.method %in% method.lib) stop('lr.method must be chosen from "flexmix","TLE","CTLE","biSquare"!  ')

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

  return(res)
}
