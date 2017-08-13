# CausalStump
This package fits a Gaussian process for treatment effect estimation. It is part of my MSc project in Computational Statistics and Machine Learning at UCL in London 2016/17.

This package can be installed using following command in R:
```R
install.packages("https://github.com/mazphilip/CausalStump/raw/master/tar/CausalStump_0.1.1.tar.gz", repos = NULL, type = "source")
```

As there is no documentation yet, here is an example of how to use the package to estimate treatment effects

```R
library(CausalStump)

#simualtion as in Hill (2011) and Hill & Su (2013)
set.seed(1231)
n=120
Z = rbinom(n, 1, 0.3)
X1 = rnorm(sum(Z), mean = 40, sd = 10)
X0 = rnorm(n-sum(Z), mean = 20, sd = 10)
X0[X0<0] = 0.01
X = data.frame(matrix(NaN,n,1))
X[Z==1,] = X1; X[Z==0,] = X0
mysort = sort(X[,1],index.return=TRUE)

y0_true = as.matrix(72 + 3 * sqrt(X))
y1_true = as.matrix(90 + exp(0.06 * X))
y0_true[is.nan(y0_true)] = 60
y0_true[y0_true<60] = 60; y1_true[y1_true<60] = 60
y0_true[y0_true>120] = 120; y1_true[y1_true>120] = 120

#add observation noise (this is different to Hill's example but preserves the additive normal noise assumption)
Y0 = rnorm(n, mean = y0_true, sd = 1)
Y1 = rnorm(n, mean = y1_true, sd = 1)
Y = Y0*(1-Z) + Y1*Z 

#parameter learning of Gaussian process (FALSE), student-t process (TRUE)
mycs = CausalStump(Y,X,Z,prior=FALSE)

#predict response surfaces
mypred0 = predict_surface(X,0,mycs)
mypred1 = predict_surface(X,1,mycs)

#plot surface fit
par(mfrow=c(1,1))
plot(X[,1],mypred1$map,ylim = c(70,120),ylab="Y",xlab="X")
points(X[,1],mypred0$map)
points(X[,1],Y,col=2,pch=20)
lines(mysort$x,y1_true[mysort$ix]); lines(mysort$x,y0_true[mysort$ix])
legend("topleft",c("true","estimates","osbervations"),pch=c(NA,1,20),lty=c(1,NA,NA),col=c(1,1,2) )

#predict treatment effect of training sample (GP: exact, TP: sampling)
mytreat = predict_treatment(X,mycs)

# the PEHE metric is used as in Hill (2011) (albeit with root):
PEHE = sqrt(mean( (y1_true-y0_true - mytreat$map)^2)); PEHE

```



