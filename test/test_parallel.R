## Load Packages #####################################################################
#BART
library(BayesTree)

#GPS
library(Rcpp)
library(RcppArmadillo)
source("GPSbasics.R")

#CF
library(grf)

#for the propensity score:
library(rms)
library(gbm)

library(foreach)
library(parallel)

#exporting to TiKz:
#library(latex2exp)
#library(tikzDevice)

## Experiment length setup #####################################################################

max.iter=1000

## Predefine result files #####################################################################
#dimnames

mydimnames = list(c("PEHE.all","PEHE.cs","RMSE.all","RMSE.cs","E.ate","E.att","coverage.ite","coverage.ate"),c("BART","CF","GPS","BART.PS","CF.PS","GPS.PS"),1:max.iter)
results.mat = array(NaN, dim=c(8,6,max.iter),dimnames=mydimnames)

update.results <- function(iter,index,result.item,y.hat,tau.hat,tau.ci,ate.ci,y,z,y1,y0,pscore){
  threshold = 0.1

  result.item["PEHE.all",index] = sqrt(mean((y1 - y0 - tau.hat)^2))
  result.item["PEHE.cs",index]  = sqrt(mean((y1 - y0 - tau.hat)[pscore>threshold & pscore<(1-threshold)]^2))

  result.item["E.ate",index]  = (mean(y1-y0) - mean(tau.hat))
  result.item["E.att",index]  = (mean((y1-y0)[z==1]) - mean((tau.hat)[z==1]))

  result.item["RMSE.all",index] = sqrt(mean((y - y.hat)^2))
  result.item["RMSE.cs",index]  = sqrt(mean((y - y.hat)[pscore>threshold & pscore<(1-threshold)]^2))

  result.item["coverage.ite",index] = mean(((tau.ci[,1] < (y1-y0)) * (tau.ci[,2] > (y1-y0))))
  result.item["coverage.ate",index] = ((ate.ci[1]< mean(y1-y0)) * (ate.ci[2] >mean(y1-y0)))
  result.item
}

## Run experiment Loop #####################################################################
#for(iter in 1:max.iter){
run <- function(iter){



  ## Generate Data #####################################################################
  set.seed(1230+iter)

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

  ## Generate Covariate matrices #####################################################################
  X.train = X
  X.test  = X

  #for BART (twice as long test file and with Z as covariate)
  X.train.bart = cbind(X.train,Z)
  colnames(X.train.bart) = c("X","Z")
  X.test.bart = cbind(rbind(X,X),rbind(matrix(1,n,1),matrix(0,n,1)))
  colnames(X.test.bart) = c("X","Z")

  # ESTIMATE THE PROPENSITY SCORE --- BART
  myPSbart = BayesTree::bart(x.train = X,y.train = Z,x.test=X,binaryOffset=mean(Z),ntree=200)
  pscore = pnorm(apply(myPSbart$yhat.test,2,mean))

  #BART file with PS
  X.train.bart.ps = cbind(X.train.bart,pscore)
  colnames(X.train.bart.ps) = c("X","Z","pscore")
  X.test.bart.ps = cbind(X.test.bart,pscore)
  colnames(X.test.bart.ps) = c("X","Z","pscore")

  #data matrix for GPS and CF (Z is not a covariate in the matrix X)
  X.train.ps = cbind(X.train,pscore)
  colnames(X.train.ps) = c("X","pscore")
  X.test.ps  = cbind(X.test,pscore)
  colnames(X.test.ps) = c("X","pscore")

  result.item = matrix(NaN, 8,6,dimnames=list(c("PEHE.all","PEHE.cs","RMSE.all","RMSE.cs","E.ate","E.att","coverage.ite","coverage.ate"),c("BART","CF","GPS","BART.PS","CF.PS","GPS.PS")) )

  ## BART #####################################################################
  ## fit

  BART.fit = BayesTree::bart(x.train = X.train.bart,y.train = Y,x.test=X.test.bart,ntree=100,ndpost=3000,nskip=500,keepevery=1)

  ## obtain estimates
  y.hat.bart = apply( BART.fit$yhat.train,2,mean)
  tau.hat.bart = apply( (BART.fit$yhat.test[,1:n] - BART.fit$yhat.test[,(n+1):(2*n)]) ,2,mean)
  ate.ci.bart = {
    tmpvar = var(apply( BART.fit$yhat.train,1,mean))
    L = mean(tau.hat.bart) - 1.96 * sqrt(tmpvar)
    U = mean(tau.hat.bart) + 1.96 * sqrt(tmpvar)
    c(L,U)
  }
  tau.ci.bart = foreach(i = 1:n,.combine=rbind) %do% {
    tmpsort = sort(BART.fit$yhat.test[,i]-BART.fit$yhat.test[,n+i],index.return=TRUE)
    idx975 = round(nrow(BART.fit$yhat.test)*0.975,0)
    idx025 = round(nrow(BART.fit$yhat.test)*0.025,0)
    L = tmpsort$x[idx025]
    U = tmpsort$x[idx975]
    c(L,U)
  }

  ## evaluate and store results
  result.item = update.results(iter,"BART",result.item,
                               y.hat = y.hat.bart,
                               tau.hat = tau.hat.bart,
                               tau.ci = tau.ci.bart,
                               ate.ci = ate.ci.bart,
                               y = Y, z = Z, y1 = y1_true, y0 = y0_true,pscore = pscore)

  #rm(y.hat.bart,ate.ci.bart)
  ## BART EPS #####################################################################
  ## fit
  BART.PS.fit = BayesTree::bart(x.train = X.train.bart.ps,y.train = Y,x.test=X.test.bart.ps,ntree=100,ndpost=3000,nskip=500,keepevery=1)

  ## obtain estimates
  y.hat.bart.ps   = apply( BART.PS.fit$yhat.train,2,mean)
  tau.hat.bart.ps = apply( (BART.PS.fit$yhat.test[,1:n] - BART.PS.fit$yhat.test[,(n+1):(2*n)]) ,2,mean)
  ate.ci.bart.ps  = {
    tmpvar = var(apply( BART.PS.fit$yhat.train,1,mean))
    L = mean(tau.hat.bart.ps) - 1.96 * sqrt(tmpvar)
    U = mean(tau.hat.bart.ps) + 1.96 * sqrt(tmpvar)
    c(L,U)
  }
  tau.ci.bart.ps = foreach(i = 1:n,.combine=rbind) %do% {
    tmpsort = sort(BART.fit$yhat.test[,i]-BART.fit$yhat.test[,n+i],index.return=TRUE)
    idx975 = round(nrow(BART.fit$yhat.test)*0.975,0)
    idx025 = round(nrow(BART.fit$yhat.test)*0.025,0)
    L = tmpsort$x[idx025]
    U = tmpsort$x[idx975]
    c(L,U)
  }

  ## evaluate and store results
  result.item = update.results(iter,"BART.PS",result.item,
                               y.hat = y.hat.bart.ps,
                               tau.hat = tau.hat.bart.ps,
                               tau.ci = tau.ci.bart.ps,
                               ate.ci = ate.ci.bart.ps,
                               y = Y, z = Z, y1 = y1_true, y0 = y0_true,pscore = pscore)


  ## GPS #####################################################################
  ## need to load files as rcpp functions cannot be handed over to a parallel environment
  source("GPSbasics.R")
  source("GPStump_ARD.R")

  ## fit
  X.train.dummy = data.frame("X" = X.train,"dum" = rep(0,n))
  GPS.fit = GPSlearning(Y,X.train.dummy,Z,pscore=pscore,maxiter=5000,tol=1e-4,optim="Nadam",objfun = "noLOO",learnrate = 0.01,beta1=0.2)

  ## predict
  GPS.pred = GPStreatmenteffect(data.frame(X.train.ps),GPSobj = GPS.fit)

  y.hat.gps = GPS.fit$pred
  tau.hat.gps = GPS.pred$map
  tau.ci.gps  = GPS.pred$ci
  ate.ci.gps = c(GPS.pred$ATE - 1.96*sqrt(GPS.pred$ATEvar),
                 GPS.pred$ATE + 1.96*sqrt(GPS.pred$ATEvar))

  ## evaluate and store results
  result.item = update.results(iter,"GPS",result.item,
                               y.hat = y.hat.gps,
                               tau.hat = tau.hat.gps,
                               tau.ci = tau.ci.gps,
                               ate.ci = ate.ci.gps,
                               y = Y, z = Z, y1 = y1_true, y0 = y0_true,pscore = pscore)

  ## GPS.PS #####################################################################
  ## fit
  GPS.PS.fit = GPSlearning(Y,X.train.ps,Z,pscore=pscore,maxiter=5000,tol=1e-4,optim="Nadam",objfun = "noLOO",learnrate = 0.01,beta1=0.2)

  ## predict
  GPS.PS.pred = GPStreatmenteffect(data.frame(X.train.ps),GPSobj = GPS.PS.fit)

  y.hat.gps.ps   = GPS.PS.fit$pred
  tau.hat.gps.ps = GPS.PS.pred$map
  tau.ci.gps.ps  = GPS.PS.pred$ci
  ate.ci.gps.ps  = c(GPS.PS.pred$ATE - 1.96*sqrt(GPS.PS.pred$ATEvar),
                     GPS.PS.pred$ATE + 1.96*sqrt(GPS.PS.pred$ATEvar))
  #tau.plot("GPS",tau.hat.gps.ps,ci = tau.ci.gps.ps)
  ## evaluate and store results
  result.item = update.results(iter,"GPS.PS",result.item,
                               y.hat = y.hat.gps.ps,
                               tau.hat = tau.hat.gps.ps,
                               tau.ci = tau.ci.gps.ps,
                               ate.ci = ate.ci.gps.ps,
                               y = Y, z = Z, y1 = y1_true, y0 = y0_true,pscore = pscore)

  ## CF #####################################################################
  ## fit
  CF.forest = grf::causal_forest(X = X.train,Y =  Y, W = Z, num.trees = 5000,num.threads=NULL)

  ## predict
  CF.pred   = predict(CF.forest, X.test = X.test, estimate.variance = TRUE)

  y.hat.cf   = CF.forest$Y.hat
  tau.hat.cf = CF.pred$predictions
  tau.ci.cf  = cbind(tau.hat.cf - 1.96 *sqrt(CF.pred$variance.estimates),
                     tau.hat.cf + 1.96 *sqrt(CF.pred$variance.estimates))
  #cannot get confidence interval for CF

  ate.ci.cf  =  { tmp.ate = grf::estimate_average_effect(CF.forest, target.sample = "all");
  c(tmp.ate["estimate"] - 1.96 * tmp.ate["std.err"],
    tmp.ate["estimate"] + 1.96 * tmp.ate["std.err"])}

  ## evaluate and store results
  result.item = update.results(iter,"CF",result.item,
                               y.hat = y.hat.cf,
                               tau.hat = tau.hat.cf,
                               tau.ci = tau.ci.cf,
                               ate.ci = ate.ci.cf,
                               y = Y, z = Z, y1 = y1_true, y0 = y0_true,pscore = pscore)

  ## CF.PS #####################################################################

  ## fit
  CF.PS.forest = grf::causal_forest(X = X.train.ps,Y =  Y, W = Z, num.trees = 5000,num.threads=NULL)

  ## predict
  CF.PS.pred   = predict(CF.PS.forest, X.test = X.test, estimate.variance = TRUE)

  y.hat.cf.ps   = CF.PS.forest$Y.hat
  tau.hat.cf.ps = CF.PS.pred$predictions
  tau.ci.cf.ps  = cbind(tau.hat.cf - 1.96 *sqrt(CF.PS.pred$variance.estimates),
                        tau.hat.cf + 1.96 *sqrt(CF.PS.pred$variance.estimates))
  #cannot get confidence interval for CF

  ate.ci.cf.ps  =  { tmp.ate = grf::estimate_average_effect(CF.PS.forest, target.sample = "all");
  c(tmp.ate["estimate"] - 1.96 * tmp.ate["std.err"],
    tmp.ate["estimate"] + 1.96 * tmp.ate["std.err"])}

  ## evaluate and store results
  result.item = update.results(iter,"CF.PS",result.item,
                               y.hat = y.hat.cf.ps,
                               tau.hat = tau.hat.cf.ps,
                               tau.ci = tau.ci.cf.ps,
                               ate.ci = ate.ci.cf.ps,
                               y = Y, z = Z, y1 = y1_true, y0 = y0_true,pscore = pscore)


  cat("Iteration ",iter, " finished\n")
  #return results matrix
  result.item
  #list(results = result.item,ci =  tau.ci.gps.ps)
} # end of run


max.iter=1000
print(sprintf("Iterations: %d",max.iter))
cl <- makeCluster(detectCores()-1)
clusterExport(cl, list("update.results", "foreach", "%do%","sourceCpp"))
#Using load balancing as convergence with tolerance might lead to very different run times per node
results <- parSapply(cl,1:max.iter,run) #clusterApplyLB, parSapply
stopCluster(cl)

#results[[1]]$results

#cannot replicate?
# [21] 0.9083333 1.0000000 0.7000000 0.8416667 0.9583333 0.8250000 0.9666667        NA 0.7833333 0.3916667 [30]

#row means
print("means")
results.mean = matrix(apply(results,1,mean),8,6,dimnames = list(c("PEHE.all","PEHE.cs","RMSE.all","RMSE.cs","E.ate","E.att","coverage.ite","coverage.ate"),c("BART","CF","GPS","BART.PS","CF.PS","GPS.PS")))
results.mean
print("se")
results.se = matrix(sqrt(apply(results,1,var)),8,6,dimnames = list(c("PEHE.all","PEHE.cs","RMSE.all","RMSE.cs","E.ate","E.att","coverage.ite","coverage.ate"),c("BART","CF","GPS","BART.PS","CF.PS","GPS.PS")))
results.se
save.image(file = paste("ex1_",max.iter,"_evid.RData",sep=""))
#results[47,]

