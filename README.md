# RMR
Robust Mixture Regression



# Installation
```
library("devtools")
devtools::install_github("changwn/RMR")
```

# Example
```
library("rmr")  #capital sensitive!
data_simu= data_simu_list[[1]]
                x=t(data_simu$x);y=as.numeric(data_simu$y);
                formula=as.formula(paste("y~",paste(colnames(x),collapse="+"),sep=""))
                data=data.frame((x),y)
                ##ctle
              
res_CRMR=CTLE(formula,data, nit=20,nc=3)
tRatio=0.05
res_TLE= TLE(formula,data, nc=3,tRatio,MaxIt=200)

```
# Questions & Problems

If you have any questions or problems, please feel free to open a new issue [here](https://github.com/changwn/RMR/issues). We will fix the new issue ASAP.  You can also email the maintainers and authors below.

- [Wennan Chang](https://zcslab.github.io/people/wennan/)
(wnchang@iu.edu)

PhD candidate at BDR group, Indiana University School of Medicine

- [Sha Cao](https://medicine.iu.edu/faculty/38873/cao-sha/)
(shacao@iu.edu)

Assistant Professor

Department of Biostatistics, Indiana University School of Medicine
