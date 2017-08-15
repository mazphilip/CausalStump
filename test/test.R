#simulation as in Hill (2011)
library(methods)
#source("R/GPSmain.R")
library(CausalStump)

set.seed(1231)

n=120
Z = rbinom(n, 1, 0.3)
Xt = rnorm(sum(Z), mean = 40, sd = 10)
Xc = rnorm(n-sum(Z), mean = 20, sd = 10)
Xc[Xc<0] = 0.01
X = data.frame(matrix(NaN,n,1))
X[Z==1,] = Xt
X[Z==0,] = Xc

y0_true = as.matrix(72 + 3 * sqrt(X))
y1_true = as.matrix(90 + exp(0.06 * X))
y0_true[is.nan(y0_true)] = 60
y0_true[y0_true<60] = 60
y1_true[y1_true<60] = 60
y0_true[y0_true>120] = 120
y1_true[y1_true>120] = 120

Y0 = rnorm(n, mean = y0_true, sd = 1)
Y1 = rnorm(n, mean = y1_true, sd = 1)
Y = Y0*(1-Z) + Y1*Z

CATE.true = mean(y1_true-y0_true)

mysort = sort(X[,1],index.return=TRUE)
X.sort = mysort$x

X2 = cbind(X)
X2 = data.frame(X2)

library(BayesTree)
myPSbart = BayesTree::bart(x.train = X,y.train = Z,x.test=X,binaryOffset=mean(Z),ntree=200)
pscore = pnorm(apply(myPSbart$yhat.test,2,mean))



CS_fit = CausalStump(Y,X,Z,pscore=pscore,maxiter=5000,tol=1e-4,learning_rate = 0.01,prior=FALSE,nu=100,myoptim = "Nadam")
#CS_fit = CausalStump(Y,X,Z,pscore=pscore,maxiter=5000,tol=1e-4,learning_rate = 0.01,prior=FALSE,myoptim = "GD")

CS_pred0 = predict(CS_fit,z=0)
CS_pred1 = predict(CS_fit,z=1)
par(mfrow=c(1,1))
plot(X[,1],CS_pred1$map,ylim = c(70,120),ylab="Y",xlab="X")
points(X[,1],CS_pred0$map)
points(X[,1],Y,col=2,pch=20)
lines(X.sort,y1_true[mysort$ix]); lines(X.sort,y0_true[mysort$ix])
legend("topleft",c("true","estimates","osbervations"),pch=c(NA,1,20),lty=c(1,NA,NA),col=c(1,1,2) )

lines(X.sort,CS_pred0$ci[mysort$ix,1])
lines(X.sort,CS_pred0$ci[mysort$ix,2])
lines(X.sort,CS_pred1$ci[mysort$ix,1])
lines(X.sort,CS_pred1$ci[mysort$ix,2])

mean(y1_true-y0_true)

CS_treat = treatment(CS_fit)
CS_treat$ate
CS_treat$ate_ci

par(mfrow=c(1,1))
plot(X[,1],CS_treat$map,ylim=c(0,30),pch=20)
lines(X.sort,CS_treat$ci[mysort$ix,1])
lines(X.sort,CS_treat$ci[mysort$ix,2])
lines(X.sort,(y1_true-y0_true)[mysort$ix],lty=2,lwd=2,col=2)

PEHE = sqrt(mean( (y1_true-y0_true - CS_treat$map)^2 )); PEHE




