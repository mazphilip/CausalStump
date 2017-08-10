#simulation as in Hill (2011)

source("R/GPSmain.R")

set.seed(1231)

n=120
Z = rbinom(n, 1, 0.3)
Xt = rnorm(sum(Z), mean = 30, sd = 10)
Xc = rnorm(n-sum(Z), mean = 30, sd = 20)
Xc[Xc<0] = 0.01
X = matrix(NaN,n,1)
X[Z==1,] = Xt
X[Z==0,] = Xc

y0_true = 72 + 3 * sqrt(X)
y1_true = 90 + exp(0.06 * X)
y0_true[is.nan(y0_true)] = 60
y0_true[y0_true<60] = 60
y1_true[y1_true<60] = 60
y0_true[y0_true>120] = 120
y1_true[y1_true>120] = 120

Y0 = rnorm(n, mean = y0_true, sd = 1)
Y1 = rnorm(n, mean = y1_true, sd = 1)
Y = Y0*(1-Z) + Y1*Z

CATE.true = mean(y1_true-y0_true)

mysort = sort(X,index.return=TRUE)
X.sort = mysort$x

X2 = cbind(X)
X2 = data.frame(X2)

library(BayesTree)
myPSbart = BayesTree::bart(x.train = X,y.train = Z,x.test=X,binaryOffset=mean(Z),ntree=200)
pscore = pnorm(apply(myPSbart$yhat.test,2,mean))

source("R/GPSmain.R")

CS_fit = CausalStump(Y,X,Z,pscore=pscore,maxiter=5000,tol=1e-4,learning_rate = 0.001,prior=TRUE,nu=2,myoptim = "GD")
CS_fit = CausalStump(Y,X,Z,pscore=pscore,maxiter=5000,tol=1e-4,learning_rate = 0.01,prior=FALSE,myoptim = "GD")

CS_pred0 = predict_surface(X2,0,CS_fit,pscore=pscore)
CS_pred1 = predict_surface(X2,1,CS_fit,pscore=pscore)
par(mfrow=c(1,1))
plot(X[,1],CS_pred1$map,ylim = c(70,120))
points(X[,1],CS_pred0$map)
lines(X.sort,y1_true[mysort$ix])
lines(X.sort,y0_true[mysort$ix])

CS_treat = predict_treatment(X2,CS_fit,pscore=pscore)
par(mfrow=c(1,1))
plot(X[,1],CS_treat$map)
lines(X.sort,CS_treat$ci[mysort$ix,1])
lines(X.sort,CS_treat$ci[mysort$ix,2])
lines(X.sort,(y1_true-y0_true)[mysort$ix],lty=2,lwd=2,col=2)

PEHE = sqrt(mean( (y1_true-y0_true - CS_treat$map)^2 )); PEHE

