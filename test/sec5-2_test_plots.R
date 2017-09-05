## Load Packages #####################################################################
library(methods)
#BART
library(BayesTree)

#GPS
library(CausalStump)

#CF and RF
library(grf)
library(ranger)

#propensity score:
library(gbm)

#parallelization
library(foreach)
library(parallel)

#exporting to TiKz:
library(latex2exp)
library(tikzDevice)


## Generate data #####################################################################
set.seed(1231)

n=120
Z = rbinom(n, 1, 0.3)
Xt = rnorm(sum(Z), mean = 40, sd = 10)
Xc = rnorm(n-sum(Z), mean = 20, sd = 10)
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

## Covariate matrices #####################################################################
X.train = data.frame(X)

X.test = -10:70
n.test = length(X.test)
y1_test = 90 + exp(0.06 * X.test)
y1_test[y1_test>120] = 120
y0_test =  72 + 3 * sqrt(X.test)
y0_test[is.nan(y0_test)] = 60
tau.test = y1_test - y0_test


#for BART (twice as long test file and with Z as covariate)
X.train.bart = data.frame(X = X.train,Z = Z)
X.test.bart  = data.frame( rbind(cbind(X.test,1),cbind(X.test,0)) )
colnames(X.test.bart) <- c("X","Z")

## Propensity score #####################################################################

#workaround as GBM can only use >1 input dimensions
tmpX = data.frame(X = X.train,Z = Z,V = rbinom(n,1,0.5))

myPSmodel = gbm(Z ~ X+V,data = tmpX,distribution="bernoulli",
                n.trees = 2000,
                shrinkage = .1,
                cv.folds = n,
                n.cores = 3,
                interaction.depth=1
)
best.fit = gbm.perf(myPSmodel,method="cv")
pscore0 = predict(myPSmodel,n.trees=best.fit,type="response",newdata = data.frame(X = X.test,V=0))
pscore1 = predict(myPSmodel,n.trees=best.fit,type="response",newdata = data.frame(X = X.test,V=1))
pscore.gbm = (pscore0 + pscore1) * 0.5

# BART
myPSbart = BayesTree::bart(x.train = X ,y.train = Z,x.test= rbind(X,matrix(X.test)), binaryOffset=mean(Z),ntree=100)
pscore.bart = apply(pnorm(myPSbart$yhat.test),2,mean)



## Append Pscore to matrices #####################################################################
pscore = pscore.bart[1:n]
pscore.test = pscore.bart[(n+1):(n+n.test)]

X.train.ps = data.frame(X = X.train,pscore = pscore)
X.test.ps = data.frame( X = X.test, pscore = pscore.test)

X.train.bart.ps = data.frame(X = X.train,pscore = pscore, Z = Z)
X.test.bart.ps = data.frame( rbind(cbind(X.test,pscore.test,1),cbind(X.test,pscore.test,0)) )
colnames(X.test.bart.ps) <- c("X","pscore","Z")


## Plots of sample #####################################################################
tikz(file = "ex1_surface.tex", width = 2.7, height = 2.3)
par(mfrow=c(1,1),mar=c(2.5,2.5,0.1,0.1),mgp=c(1.3,0.3,0),cex.lab=0.8,cex.axis=0.8,tcl=-0.3)
plot(X[Z==1],Y[Z==1],ylab=TeX("$Y$"),xlab=TeX("$X$"),pch=20,ylim=c(60,125))
points(X[Z==0],Y[Z==0],pch=21,cex=0.8)
lines(X.sort,y1_true[mysort$ix],lty=2)
lines(X.sort,y0_true[mysort$ix],lty=2)
legend("bottomright",c("Treated","Control","True"),lty=c(NA,NA,2),pch=c(20,21,NA),cex = 0.8,bg="white")
dev.off()

tikz(file = "ex1_treat.tex", width = 2.7, height = 2.3)
par(mfrow=c(1,1),mar=c(2.5,2.5,0.1,0.1),mgp=c(1.3,0.3,0),cex.lab=0.8,cex.axis=0.8,tcl=-0.3)
plot(X.sort,pscore[mysort$ix],ylim=c(0,1),xlim=c(min(X.test),max(X.test)),type = "n",lty=1,ylab="Propensity Score",xlab=TeX("$X$"))
grid(nx=0,ny=NULL)
abline(v=20,col="lightgrey",lty=3); abline(v=40,col="lightgrey",lty=3)
lines(X.test,pscore.gbm,lty=2)
lines(X.test,pscore.test,lty=3)
points(X[Z==1],Z[Z==1],col=1,pch=20)
points(X[Z==0],Z[Z==0],col=1,pch=21,cex=0.8)
legend("topleft",c("GBM","BART"),lty=c(2,3),pch=c(NA,NA),cex = 0.8,bg="white")
dev.off()

## Plot function #####################################################################
tau.plot <- function(tau,ci){
  par(mfrow=c(1,1),mar=c(2.5,2.5,0.1,0.1),mgp=c(1.3,0.3,0),cex.lab=0.8,cex.axis=0.8,tcl=-0.3)
  plot(NULL,NULL,ylim=c(0,30),xlim=c(min(X.test),max(X.test)),xlab=TeX("X"),ylab="Treatment effect")

  L = ci[,1]
  U = ci[,2]

  polygon(c(X.test,rev(X.test)),c(L,rev(U)),col = "grey75", border = FALSE)
  lines(X.test,tau.test,type="l",lty=2)
  lines(X.test,tau,lty=1)
  points(X[Z==1],rep(30,sum(Z)),pch=20,cex=0.5)
  points(X[Z==0],rep(0,sum(1-Z)),pch=21,cex=0.5)
}
plot_w = 3; plot_h = 1.3

## BART #####################################################################
## fit
BART_fit = BayesTree::bart(x.train = X.train.bart,y.train = Y,x.test = X.test.bart,
                           ntree=100,ndpost=3000,nskip=500,keepevery=1)

cat("Share of Z-splits for BART: ",apply(BART_fit$varcount/apply(BART_fit$varcount,1,sum),2,mean)[2])# 41.8% of splits are for Z



## obtain estimates
tau.hat.bart = BART_fit$yhat.test.mean[1:n.test] - BART_fit$yhat.test.mean[(n.test+1):(2*n.test)]
tau.ci.bart = {
  tmpsort = apply(BART_fit$yhat.test[,1:n.test]-BART_fit$yhat.test[,(n.test+1):(2*n.test)],2,sort)
  idx975 = round(nrow(BART_fit$yhat.test)*0.975,0)
  idx025 = round(nrow(BART_fit$yhat.test)*0.025,0)
  L = tmpsort[idx025,]
  U = tmpsort[idx975,]
  cbind(L,U)
}
## plot using test set
tikz(file = "ex1_bart.tex", width = plot_w, height = plot_h)
tau.plot(tau.hat.bart,ci = tau.ci.bart)
dev.off()

## BART with PS #####################################################################
## fit
BART_PS_fit = BayesTree::bart(x.train = X.train.bart.ps,y.train = Y,x.test = X.test.bart.ps,
                           ntree=100,ndpost=3000,nskip=500,keepevery=1)

apply(BART_PS_fit$varcount/apply(BART_PS_fit$varcount,1,sum),2,mean)# 30.9% of all splits are for Z

## obtain estimates
tau.hat.bart.ps = BART_PS_fit$yhat.test.mean[1:n.test] - BART_PS_fit$yhat.test.mean[(n.test+1):(2*n.test)]
tau.ci.bart.ps = {
  tmpsort = apply(BART_PS_fit$yhat.test[,1:n.test]-BART_PS_fit$yhat.test[,(n.test+1):(2*n.test)],2,sort)
  idx975 = round(nrow(BART_PS_fit$yhat.test)*0.975,0)
  idx025 = round(nrow(BART_PS_fit$yhat.test)*0.025,0)
  L = tmpsort[idx025,]
  U = tmpsort[idx975,]
  cbind(L,U)
}

## plot using test set
tikz(file = "ex1_bartps.tex", width = plot_w, height = plot_h)
tau.plot(tau.hat.bart.ps,ci = tau.ci.bart.ps)
dev.off()

## BCF / CF-BART via Lu et al (2017) #####################################################################
## no prior on variance of prior of the treatment effect (as in Hahn et al 2017)
#fit
BCF.fit0 = BayesTree::bart(x.train = X.train.ps[Z==0,], y.train = Y[Z==0],x.test = X.test.ps,
                           ntree=80,ndpost=3000,nskip=500,keepevery=1,power=2.0)
BCF.fit1 = BayesTree::bart(x.train = X.train.ps[Z==1,], y.train = Y[Z==1],x.test = X.test.ps,
                           ntree=20,ndpost=3000,nskip=500,keepevery=1,power=2.0)

## obtain estimates
tau.hat.bcf = BCF.fit1$yhat.test.mean- BCF.fit0$yhat.test.mean
tau.ci.bcf = {
  tmpsort = apply(abs(BCF.fit1$yhat.test - BCF.fit0$yhat.test),2,sort)
  idx975 = round(nrow(BCF.fit1$yhat.test)*0.975,0)
  idx025 = round(nrow(BCF.fit1$yhat.test)*0.025,0)
  L = tmpsort[idx025,]
  U = tmpsort[idx975,]
  #idx950 = round(nrow(BCF.fit1$yhat.test)*0.95,0)
  #U = tmpsort[idx950,]
  #L = tau.hat.bcf - (tmpsort[idx950,]-tau.hat.bcf)
  cbind(L,U)
}

## plot using test set
tikz(file = "ex1_bcf.tex", width = plot_w, height = plot_h)
tau.plot(tau.hat.bcf,ci = tau.ci.bcf)
dev.off()


## GP #####################################################################
## fit
GP_fit = CausalStump(Y,X.train,Z,maxiter=5000,tol=1e-4,learning_rate = 0.01,prior=FALSE,myoptim = "GD")

## predict
GS_treat = treatment(GP_fit,X=data.frame(X.test))

tau.hat.gp = -GS_treat$map
tau.ci.gp  = cbind(tau.hat.gp,tau.hat.gp) + (GS_treat$ci - cbind(GS_treat$map,GS_treat$map) )

tau.hat.gp = GS_treat$map
tau.ci.gp  = GS_treat$ci

GS_pred1 = predict(GP_fit,X=data.frame(X.test),z=1)
GS_pred0 = predict(GP_fit,X=data.frame(X.test),z=0)

tikz(file = "ex1_GP_uncertainty.tex", width = 4, height = 2.5)
par(mfrow=c(1,1),mar=c(2.5,2.5,0.1,0.1),mgp=c(1.3,0.3,0),cex.lab=0.8,cex.axis=0.8,tcl=-0.3)
plot(NULL,NULL,ylim=c(60,130),xlim=c(min(X.test),max(X.test)),ylab="Outcome",xlab=TeX("X"))
polygon(c(X.test,rev(X.test)),c(GS_pred0$ci[,1],rev(GS_pred0$ci[,2])),col = mycol60, border = FALSE)
polygon(c(X.test,rev(X.test)),c(GS_pred1$ci[,1],rev(GS_pred1$ci[,2])),col = mycol30, border = TRUE,density=20,angle=80)
lines(X.test,y1_test,lty=2)
lines(X.test,y0_test,lty=2)
lines(X.test,GS_pred1$map)
lines(X.test,GS_pred0$map)
#legendlines = legend("topleft",c("true","map"),lty=c(2,1),col=c(1,1))
#legend(x=0.5,y=legendlines$rect$top,c("95pc ci","95pc ci"),fill=c(mycol30,mycol60),density=c(20,NA),border=c(TRUE,FALSE))
#legend("topleft", legend = c("true", "map", "ci treat", "ci cntr"),
#       bty = "o", bg="white",
#       col = c(1,1,NA,NA),
#       pch=c(NA,NA,NA,NA),
#       lty = c(2,1,NA,NA),
#       density=c(NA,NA,20,NA),
#       fill = c(0,0,mycol30,mycol60),
#       angle = c(NA,NA,80,0),
#       border = c(NA,NA,mycol30,NA),
#       x.intersp=c(1,1,1,1)
#)
points(X[Z==1],rep(70+60,sum(Z)),col=1,pch=20,cex=0.5)
points(X[Z==0],rep(60,sum(1-Z)),col=1,pch=21,cex=0.5)
dev.off()

## Plot
tikz(file = "ex1_gp.tex", width = plot_w, height = plot_h)
tau.plot(tau.hat.gp,ci = tau.ci.gp)
dev.off()

## Extra plots for BCF #####################################################################
mycol60 <- rgb(191, 191, 191, max=255, alpha = (100-60)*255/100)
mycol30 <- rgb(191, 191, 191, max=255, alpha = (100-30)*255/100)

bcf1.ci = {
  tmpsort = apply((BCF.fit1$yhat.test),2,sort)
  idx975 = round(nrow(BCF.fit1$yhat.test)*0.975,0)
  idx025 = round(nrow(BCF.fit1$yhat.test)*0.025,0)
  L = tmpsort[idx025,]
  U = tmpsort[idx975,]
  cbind(L,U)
}

bcf0.ci = {
  tmpsort = apply((BCF.fit0$yhat.test),2,sort)
  idx975 = round(nrow(BCF.fit1$yhat.test)*0.975,0)
  idx025 = round(nrow(BCF.fit1$yhat.test)*0.025,0)
  L = tmpsort[idx025,]
  U = tmpsort[idx975,]
  cbind(L,U)
}

tikz(file = "ex1_BCF_uncertainty.tex", width = 4, height = 2.5)
  par(mfrow=c(1,1),mar=c(2.5,2.5,0.1,0.1),mgp=c(1.3,0.3,0),cex.lab=0.8,cex.axis=0.8,tcl=-0.3)
  plot(NULL,NULL,ylim=c(60,130),xlim=c(min(X.test),max(X.test)),ylab="Outcome",xlab=TeX("X"))
  polygon(c(X.test,rev(X.test)),c(bcf0.ci[,1],rev(bcf0.ci[,2])),col = mycol60, border = FALSE)
  polygon(c(X.test,rev(X.test)),c(bcf1.ci[,1],rev(bcf1.ci[,2])),col = mycol30, border = TRUE,density=20,angle=80)
  lines(X.test,y1_test,lty=2)
  lines(X.test,y0_test,lty=2)
  lines(X.test,BCF.fit0$yhat.test.mean)
  lines(X.test,BCF.fit1$yhat.test.mean)
  legend("topleft", legend = c("true", "map", "ci treat", "ci cntr"),
         bty = "o", bg="white",
         col = c(1,1,NA,NA),
         pch=c(NA,NA,NA,NA),
         lty = c(2,1,NA,NA),
         density=c(NA,NA,20,NA),
         fill = c(0,0,mycol30,mycol60),
         angle = c(NA,NA,80,0),
         border = c(NA,NA,mycol30,NA),
         x.intersp=c(1,1,1,1)
  )
  points(X[Z==1],rep(70+60,sum(Z)),col=1,pch=20)
  points(X[Z==0],rep(60,sum(1-Z)),col=1,pch=21)
dev.off()

## compare with BART
bart1.ci = {
  tmpsort = apply(BART_fit$yhat.test[,1:n.test],2,sort)
  idx975 = round(nrow(BART_fit$yhat.test)*0.975,0)
  idx025 = round(nrow(BART_fit$yhat.test)*0.025,0)
  L = tmpsort[idx025,]
  U = tmpsort[idx975,]
  cbind(L,U)
}

bart0.ci = {
  tmpsort = apply((BART_fit$yhat.test[,(n.test+1):(2*n.test)]),2,sort)
  idx975 = round(nrow(BART_fit$yhat.test)*0.975,0)
  idx025 = round(nrow(BART_fit$yhat.test)*0.025,0)
  L = tmpsort[idx025,]
  U = tmpsort[idx975,]
  cbind(L,U)
}

tikz(file = "ex1_BART_uncertainty.tex", width = 4, height = 2.5)
par(mfrow=c(1,1),mar=c(2.5,2.5,0.1,0.1),mgp=c(1.3,0.3,0),cex.lab=0.8,cex.axis=0.8,tcl=-0.3)
plot(NULL,NULL,ylim=c(60,130),xlim=c(min(X.test),max(X.test)),ylab="Outcome",xlab=TeX("X"))
polygon(c(X.test,rev(X.test)),c(bart0.ci[,1],rev(bart0.ci[,2])),col = mycol60, border = FALSE)
polygon(c(X.test,rev(X.test)),c(bart1.ci[,1],rev(bart1.ci[,2])),col = mycol30, border = TRUE,density=20,angle=80)
lines(X.test,y1_test,lty=2)
lines(X.test,y0_test,lty=2)
lines(X.test,BART_fit$yhat.test.mean[1:n.test])
lines(X.test,BART_fit$yhat.test.mean[(n.test+1):(2*n.test)])
#legend("topleft", legend = c("true", "map", "ci treat", "ci cntr"),
#       bty = "o", bg="white",
#       col = c(1,1,NA,NA),
#       pch=c(NA,NA,NA,NA),
#       lty = c(2,1,NA,NA),
#       density=c(NA,NA,20,NA),
#       fill = c(0,0,mycol30,mycol60),
#       angle = c(NA,NA,80,0),
#       border = c(NA,NA,mycol30,NA),
#       x.intersp=c(1,1,1,1)
#)
points(X[Z==1],rep(70+60,sum(Z)),col=1,pch=20)
points(X[Z==0],rep(60,sum(1-Z)),col=1,pch=21)
dev.off()

## compare with compare with CS-GP
GS_pred1 = predict(GP_fit,X=data.frame(X.test),z=1)
GS_pred0 = predict(GP_fit,X=data.frame(X.test),z=0)

tikz(file = "ex1_CSGP_uncertainty.tex", width = 4, height = 2.5)
par(mfrow=c(1,1),mar=c(2.5,2.5,0.1,0.1),mgp=c(1.3,0.3,0),cex.lab=0.8,cex.axis=0.8,tcl=-0.3)
plot(NULL,NULL,ylim=c(60,130),xlim=c(min(X.test),max(X.test)),ylab="Outcome",xlab=TeX("X"))
polygon(c(X.test,rev(X.test)),c(GS_pred0$ci[,1],rev(GS_pred0$ci[,2])),col = mycol60, border = FALSE)
polygon(c(X.test,rev(X.test)),c(GS_pred1$ci[,1],rev(GS_pred1$ci[,2])),col = mycol30, border = TRUE,density=20,angle=80)
lines(X.test,y1_test,lty=2)
lines(X.test,y0_test,lty=2)
lines(X.test,GS_pred1$map)
lines(X.test,GS_pred0$map)
#legendlines = legend("topleft",c("true","map"),lty=c(2,1),col=c(1,1))
#legend(x=0.5,y=legendlines$rect$top,c("95pc ci","95pc ci"),fill=c(mycol30,mycol60),density=c(20,NA),border=c(TRUE,FALSE))
legend("topleft", legend = c("true", "map", "ci treat", "ci cntr"),
       bty = "o", bg="white",
       col = c(1,1,NA,NA),
       pch=c(NA,NA,NA,NA),
       lty = c(2,1,NA,NA),
       density=c(NA,NA,20,NA),
       fill = c(0,0,mycol30,mycol60),
       angle = c(NA,NA,80,0),
       border = c(NA,NA,mycol30,NA),
       x.intersp=c(1,1,1,1)
)
points(X[Z==1],rep(70+60,sum(Z)),col=1,pch=20)
points(X[Z==0],rep(60,sum(1-Z)),col=1,pch=21)
dev.off()

## VT-GP / RR #####################################################################

## fit
VTGP_fit = CausalStump(Y,cbind(X.train,Z),rep(0,n),maxiter=5000,tol=1e-4,learning_rate = 0.01,prior=FALSE,myoptim = "GD")

#predict
VTGS_pred1 = predict(VTGP_fit,X=cbind(X.test,1),z=0)
VTGS_pred0 = predict(VTGP_fit,X=cbind(X.test,0),z=0)

tau.hat.vtgp = VTGS_pred1$map - VTGS_pred0$map

## CF-GP / RR #####################################################################
## fit
CFGP_fit1 = CausalStump(Y[Z==1],X=data.frame(X.train[Z==1,]),z=rep(0,sum(Z)),maxiter=5000,tol=1e-4,learning_rate = 0.01,prior=FALSE,myoptim = "GD")
CFGP_fit0 = CausalStump(Y[Z==0],X=data.frame(X.train[Z==0,]),z=rep(0,sum(1-Z)),maxiter=5000,tol=1e-4,learning_rate = 0.01,prior=FALSE,myoptim = "GD")

CFGS_pred1 = predict(CFGP_fit1,X=cbind(X.test),z=0)
CFGS_pred0 = predict(CFGP_fit0,X=cbind(X.test),z=0)

tau.hat.cfgp = CFGS_pred1$map - CFGS_pred0$map

## Print CS-GP v VT_RR v CF-RR ######

tikz(file = "ex1_GP_comparison.tex", width = 4, height = 2.5)
par(mfrow=c(1,1),mar=c(2.5,2.5,0.1,0.1),mgp=c(1.3,0.3,0),cex.lab=0.8,cex.axis=0.8,tcl=-0.3)
plot(NULL,NULL,ylim=c(0,30),xlim=c(min(X.test),max(X.test)),xlab=TeX("X"),ylab="Treatment effect")
#polygon(c(X.test,rev(X.test)),c(tau.ci.gp[,1],rev(tau.ci.gp[,2])),col = mycol60, border = FALSE)
lines(X.test,tau.test,lty=2,col=mycol30)
lines(X.test,tau.hat.gp,lty=1)
lines(X.test,tau.hat.cfgp,lty=3)
lines(X.test,tau.hat.vtgp,lty=4)
legend("topleft", legend = c("true", "CS-GP", "VT-GP", "CF-GP"),
       bty = "o", bg="white",
       col = c(mycol30,1,1,1),
       pch=c(NA,NA,NA,NA),
       lty = c(2,1,3,4)
)
points(X[Z==1],rep(30,sum(Z)),col=1,pch=20)
points(X[Z==0],rep(0,sum(1-Z)),col=1,pch=21)
dev.off()


## GP with PS #####################################################################
## fit
GP_PS_fit = CausalStump(Y,X.train,Z,pscore=pscore,maxiter=5000,tol=1e-4,learning_rate = 0.01,prior=FALSE,myoptim = "GD")

## predict
GS_PS_treat = treatment(GP_PS_fit, X = data.frame(X.test),pscore=pscore.test)

tau.hat.gp.ps = GS_PS_treat$map
tau.ci.gp.ps  = GS_PS_treat$ci

## plot
tikz(file = "ex1_gp.ps.tex", width = plot_w, height = plot_h)
tau.plot(tau.hat.gp.ps,tau.ci.gp.ps)
dev.off()

## Analysis GP w/wo PS #######

#par(mfrow=c(1,1),mar=c(2.5,2.5,0.1,0.1),mgp=c(1.3,0.3,0),cex.lab=0.8,cex.axis=0.8,tcl=-0.3)
#plot(NULL,NULL,ylim=c(0,30),xlim=c(0,60),xlab=TeX("X"),ylab="Treatment effect")
#polygon(c(X.test,rev(X.test)),c(tau.ci.gp.ps[,1],rev(tau.ci.gp.ps[,2])),col = "blue", border = TRUE)#,density=20,angle=80
#polygon(c(X.test,rev(X.test)),c(tau.ci.gp[,1],rev(tau.ci.gp[,2])),col = "red", border = TRUE)
#lines(X.test,tau.test,lty=2,col=mycol30)
#lines(X.test,tau.hat.gp,lty=1)
#lines(X.test,tau.hat.gp.ps,lty=3)
#legend("topleft", legend = c("true", "CS-GP", "VT-GP", "CF-GP"),
#       bty = "o", bg="white",
#       col = c(mycol30,1,1,1),
#       pch=c(NA,NA,NA,NA),
#       lty = c(2,1,3,4)
#)
#points(X[Z==1],rep(30,sum(Z)),col=1,pch=20)
#points(X[Z==0],rep(0,sum(1-Z)),col=1,pch=21)

## TP2 #####################################################################
## fit
TP2_fit = CausalStump(Y,X.train,Z,learning_rate = 0.0001,prior=TRUE,nu=2,myoptim = "Nesterov",momentum=0.1)

## predict
TP2_treat = treatment(TP2_fit,X=data.frame(X.test))

tau.hat.tp2 = TP2_treat$map
tau.ci.tp2  = TP2_treat$ci

## plot
tikz(file = "ex1_tp2.tex", width = plot_w, height = plot_h)
tau.plot(tau.hat.tp2,tau.ci.tp2)
dev.off()

### TP2 with PS #####################################################################
### fit
#TP2_PS_fit = CausalStump(Y,X.train,Z,pscore=pscore,maxiter=5000,tol=1e-4,learning_rate = 0.01,prior=TRUE,nu=2,myoptim = "Nadam")
#
### predict
#TP2_PS_pred = predict(TP2_PS_fit)
#TP2_PS_treat = treatment(TP2_PS_fit,X=data.frame(X.test),pscore=pscore.test)
#
#tau.hat.tp2.ps = TP2_PS_treat$map
#tau.ci.tp2.ps  = TP2_PS_treat$ci
#
### plot
#tikz(file = "ex1_tp2.ps.tex", width = plot_w, height = plot_h)
#tau.plot(tau.hat.tp2.ps,tau.ci.tp2.ps)
#dev.off()

## TP100 #####################################################################
## fit
TP100_fit = CausalStump(Y,X.train,Z,learning_rate = 0.0001,prior=TRUE,nu=100,myoptim = "Nesterov",momentum=0.1)

## predict
TP100_treat = treatment(TP100_fit,X=data.frame(X.test))
TP100_pred1 = predict(TP100_fit,X=data.frame(X.test),z=1)
TP100_pred0 = predict(TP100_fit,X=data.frame(X.test),z=0)

bcf1.ci = TP100_pred1$ci
bcf0.ci = TP100_pred0$ci

#par(mfrow=c(1,1),mar=c(2.5,2.5,0.1,0.1),mgp=c(1.3,0.3,0),cex.lab=0.8,cex.axis=0.8,tcl=-0.3)
#plot(NULL,NULL,ylim=c(60,130),xlim=c(min(X.test),max(X.test)),ylab="Outcome",xlab=TeX("X"))
#polygon(c(X.test,rev(X.test)),c(bcf0.ci[,1],rev(bcf0.ci[,2])),col = "blue", border = FALSE)
#polygon(c(X.test,rev(X.test)),c(bcf1.ci[,1],rev(bcf1.ci[,2])),col = "red", border = TRUE,density=20,angle=80)
#lines(X.test,y1_test,lty=2)
#lines(X.test,y0_test,lty=2)
#lines(X.test,TP100_pred0$map)
#lines(X.test,TP100_pred1$map)
#legend("topleft", legend = c("true", "map", "ci treat", "ci cntr"),
#       bty = "o", bg="white",
#       col = c(1,1,NA,NA),
#       pch=c(NA,NA,NA,NA),
#       lty = c(2,1,NA,NA),
#       density=c(NA,NA,20,NA),
#       fill = c(0,0,mycol30,mycol60),
#       angle = c(NA,NA,80,0),
#       border = c(NA,NA,mycol30,NA),
#       x.intersp=c(1,1,1,1)
#)
#points(X[Z==1],rep(70+60,sum(Z)),col=1,pch=20)
#points(X[Z==0],rep(60,sum(1-Z)),col=1,pch=21)

tau.hat.tp100 = TP100_treat$map
tau.ci.tp100  = TP100_treat$ci



tikz(file = "ex1_tp100.tex", width = plot_w, height = plot_h)
tau.plot(tau.hat.tp100,tau.ci.tp100)
dev.off()

## TP100 with PS #####################################################################
## fit
#TP100_PS_fit = CausalStump(Y,X.train,Z,pscore=pscore,maxiter=5000,tol=1e-4,learning_rate = 0.001,prior=TRUE,nu=100,myoptim = "Nadam")
#
#TP100_PS_fit$Kernel$parameters
#
#
### predict
#TP100_PS_pred = predict(TP100_PS_fit)
#TP100_PS_treat = treatment(TP100_PS_fit,X=X.test,pscore=pscore.test)
#
##
#
#y.hat.tp100.ps   = TP100_PS_pred$map
#tau.hat.tp100.ps = TP100_PS_treat$map
#tau.ci.tp100.ps  = TP100_PS_treat$ci
#ate.ci.tp100.ps  = TP100_PS_treat$ate_ci
#
#tikz(file = "ex1_tp100ps.tex", width = plot_w, height = plot_h)
#tau.plot(tau.hat.tp100.ps,tau.ci.tp100.ps)
#dev.off()


## CF #####################################################################
## fit
CF.forest = grf::causal_forest(X = as.matrix(X.train),Y =  Y, W = Z, num.trees = 5000,
                               precompute.nuisance = TRUE,num.threads=2,honesty=TRUE,min.node.size=10,mtry=1)

#split_frequencies(CF.forest)

## predict
CF.pred   = predict(CF.forest, newdata = X.test, estimate.variance = TRUE,num.threads=2)

tau.hat.cf = CF.pred$predictions
tau.ci.cf  = cbind(tau.hat.cf - 1.96 *sqrt(CF.pred$variance.estimates),
                   tau.hat.cf + 1.96 *sqrt(CF.pred$variance.estimates))

tikz(file = "ex1_cf.tex", width = plot_w, height = plot_h)
tau.plot(tau.hat.cf,ci = tau.ci.cf)
dev.off()

### CF.PS (obsolete -- too little data) #####################################################################
#
### fit
#CF.PS.forest = grf::causal_forest(X = as.matrix(X.train.ps),Y =  Y, W = Z, num.trees = 5000,precompute.nuisance = TRUE,num.threads=1,honesty = TRUE,min.node.size=10)
#
### predict
#CF.PS.pred   = predict(CF.PS.forest, newdata = as.matrix(X.test.ps), estimate.variance = TRUE,num.threads=2)
#
#
#tau.hat.cf.ps = CF.pred$predictions
#tau.ci.cf.ps  = cbind(tau.hat.cf.ps - 1.96 *sqrt(CF.PS.pred$variance.estimates),
#                   tau.hat.cf.ps + 1.96 *sqrt(CF.PS.pred$variance.estimates))
#
#
#tikz(file = "ex1_cfps.tex", width = plot_w, height = plot_h)
#tau.plot(tau.hat.cf.ps,tau.ci.cf.ps)
#dev.off()

## VT-RF #####################################################################

VTRF_fit = ranger::ranger(Y ~ ., data = data.frame(Y=Y,X=X.train,Z=Z),num.trees = 5000,min.node.size=5,mtry=2)
VTRF_pred1 = predict(VTRF_fit, data = data.frame(X=X.test,Z=1) , predict.all = TRUE)
VTRF_pred0 = predict(VTRF_fit, data = data.frame(X=X.test,Z=0) , predict.all = TRUE)

cat("The average VT-RF splitting rule over 1000 trees is")
splitvtrf = split_frequencies(grf::regression_forest(X = as.matrix(cbind(X.train,Z)),Y =  Y, num.trees = 1000,
                       num.threads=1,honesty=FALSE,min.node.size=5,mtry=1))
cat("Initial split: ",splitvtrf[1,2]/1000)
cat("Depth 2 share: ",splitvtrf[2,2]/1000)


tau.hat.vtrf = apply(VTRF_pred1$predictions,1,mean) - apply(VTRF_pred0$predictions,1,mean)
tau.ci.vtrf = {
  tmpsort = apply(VTRF_pred1$predictions - VTRF_pred0$predictions,1,sort)
  idx975 = round(ncol(VTRF_pred1$predictions)*0.975,0)
  idx025 = round(ncol(VTRF_pred1$predictions)*0.025,0)
  L = tmpsort[idx025,]
  U = tmpsort[idx975,]
  cbind(L,U)
}

tikz(file = "ex1_vtrf.tex", width = plot_w, height = plot_h)
tau.plot(tau.hat.vtrf, tau.ci.vtrf)
dev.off()

## VTRF response surface #############################
vtrf0_ci = {
  tmpsort = apply(VTRF_pred0$predictions,1,sort)
  idx975 = round(ncol(VTRF_pred1$predictions)*0.975,0)
  idx025 = round(ncol(VTRF_pred1$predictions)*0.025,0)
  L = tmpsort[idx025,]
  U = tmpsort[idx975,]
  cbind(L,U)
}

vtrf1_ci = {
  tmpsort = apply(VTRF_pred1$predictions,1,sort)
  idx975 = round(ncol(VTRF_pred1$predictions)*0.975,0)
  idx025 = round(ncol(VTRF_pred1$predictions)*0.025,0)
  L = tmpsort[idx025,]
  U = tmpsort[idx975,]
  cbind(L,U)
}

tikz(file = "ex1_VTRF_uncertainty.tex", width = 4, height = 2.5)
par(mfrow=c(1,1),mar=c(2.5,2.5,0.1,0.1),mgp=c(1.3,0.3,0),cex.lab=0.8,cex.axis=0.8,tcl=-0.3)
plot(NULL,NULL,ylim=c(60,130),xlim=c(min(X.test),max(X.test)),ylab="Outcome",xlab=TeX("X"))
polygon(c(X.test,rev(X.test)),c(vtrf0_ci[,1],rev(vtrf0_ci[,2])),col = mycol60, border = FALSE)
polygon(c(X.test,rev(X.test)),c(vtrf1_ci[,1],rev(vtrf1_ci[,2])),col = mycol30, border = TRUE,density=20,angle=80)
lines(X.test,y1_test,lty=2)
lines(X.test,y0_test,lty=2)
lines(X.test,apply(VTRF_pred0$predictions,1,mean))
lines(X.test,apply(VTRF_pred1$predictions,1,mean))
legend("topleft", legend = c("true", "map", "ci treat", "ci cntr"),
       bty = "o", bg="white",
       col = c(1,1,NA,NA),
       pch=c(NA,NA,NA,NA),
       lty = c(2,1,NA,NA),
       density=c(NA,NA,20,NA),
       fill = c(0,0,mycol30,mycol60),
       angle = c(NA,NA,80,0),
       border = c(NA,NA,mycol30,NA),
       x.intersp=c(1,1,1,1)
)
points(X[Z==1],rep(70+60,sum(Z)),col=1,pch=20)
points(X[Z==0],rep(60,sum(1-Z)),col=1,pch=21)
dev.off()

tikz(file = "ex1_VTRF_uncertainty_nolegend.tex", width = 4, height = 2.5)
par(mfrow=c(1,1),mar=c(2.5,2.5,0.1,0.1),mgp=c(1.3,0.3,0),cex.lab=0.8,cex.axis=0.8,tcl=-0.3)
plot(NULL,NULL,ylim=c(60,130),xlim=c(min(X.test),max(X.test)),ylab="Outcome",xlab=TeX("X"))
polygon(c(X.test,rev(X.test)),c(vtrf0_ci[,1],rev(vtrf0_ci[,2])),col = mycol60, border = FALSE)
polygon(c(X.test,rev(X.test)),c(vtrf1_ci[,1],rev(vtrf1_ci[,2])),col = mycol30, border = TRUE,density=20,angle=80)
lines(X.test,y1_test,lty=2)
lines(X.test,y0_test,lty=2)
lines(X.test,apply(VTRF_pred0$predictions,1,mean))
lines(X.test,apply(VTRF_pred1$predictions,1,mean))
points(X[Z==1],rep(70+60,sum(Z)),col=1,pch=20)
points(X[Z==0],rep(60,sum(1-Z)),col=1,pch=21)
dev.off()

## CF-RF #####################################################################
CFRF0_fit = ranger(Y ~ X, data = data.frame(Y=Y,X=X.train)[Z==0,],num.trees = 500,min.node.size=5,mtry=1)
CFRF1_fit = ranger(Y ~ X, data = data.frame(Y=Y,X=X.train)[Z==1,],num.trees = 500,min.node.size=5,mtry=1)

CFRF_pred0 = predict(CFRF0_fit, data = data.frame(X=X.test) , predict.all = TRUE)
CFRF_pred1 = predict(CFRF1_fit, data = data.frame(X=X.test) , predict.all = TRUE)

tau.hat.cfrf = apply(CFRF_pred1$predictions,1,mean) - apply(CFRF_pred0$predictions,1,mean)

tau.ci.cfrf = {
  tmpsort = apply(CFRF_pred1$predictions - CFRF_pred0$predictions,1,sort)
  idx975 = round(ncol(CFRF_pred0$predictions)*0.975,0)
  idx025 = round(ncol(CFRF_pred0$predictions)*0.025,0)
  L = tmpsort[idx025,]
  U = tmpsort[idx975,]
  cbind(L,U)
}

tikz(file = "ex1_cfrf.tex", width = plot_w, height = plot_h)
tau.plot(tau.hat.cfrf, tau.ci.cfrf)
dev.off()


## CFRF response surface #############################
cfrf0_ci = {
  tmpsort = apply(CFRF_pred0$predictions,1,sort)
  idx975 = round(ncol(CFRF_pred1$predictions)*0.975,0)
  idx025 = round(ncol(CFRF_pred1$predictions)*0.025,0)
  L = tmpsort[idx025,]
  U = tmpsort[idx975,]
  cbind(L,U)
}

cfrf1_ci = {
  tmpsort = apply(CFRF_pred1$predictions,1,sort)
  idx975 = round(ncol(CFRF_pred1$predictions)*0.975,0)
  idx025 = round(ncol(CFRF_pred1$predictions)*0.025,0)
  L = tmpsort[idx025,]
  U = tmpsort[idx975,]
  cbind(L,U)
}

tikz(file = "ex1_CFRF_uncertainty.tex", width = 4, height = 2.5)
par(mfrow=c(1,1),mar=c(2.5,2.5,0.1,0.1),mgp=c(1.3,0.3,0),cex.lab=0.8,cex.axis=0.8,tcl=-0.3)
plot(NULL,NULL,ylim=c(60,130),xlim=c(min(X.test),max(X.test)),ylab="Outcome",xlab=TeX("X"))
polygon(c(X.test,rev(X.test)),c(cfrf0_ci[,1],rev(cfrf0_ci[,2])),col = mycol60, border = FALSE)
polygon(c(X.test,rev(X.test)),c(cfrf1_ci[,1],rev(cfrf1_ci[,2])),col = mycol30, border = TRUE,density=20,angle=80)
lines(X.test,y1_test,lty=2)
lines(X.test,y0_test,lty=2)
lines(X.test,apply(CFRF_pred0$predictions,1,mean))
lines(X.test,apply(CFRF_pred1$predictions,1,mean))
#legend("topleft", legend = c("true", "map", "ci treat", "ci cntr"),
#       bty = "o", bg="white",
#       col = c(1,1,NA,NA),
#       pch=c(NA,NA,NA,NA),
#       lty = c(2,1,NA,NA),
#       density=c(NA,NA,20,NA),
#       fill = c(0,0,mycol30,mycol60),
#       angle = c(NA,NA,80,0),
#       border = c(NA,NA,mycol30,NA),
#       x.intersp=c(1,1,1,1)
#)
points(X[Z==1],rep(70+60,sum(Z)),col=1,pch=20)
points(X[Z==0],rep(60,sum(1-Z)),col=1,pch=21)
dev.off()
