#Hill (2011)'s IHDP simulation, surface B
#Sys.setenv(OMP_NUM_THREADS=2)
niters=1

## Libraries #######
library(CausalStump)
library(BayesTree)
library(grf)
library(ranger)
library(parallel)

## Set seed #######
set.seed(1234)

## Preprocessing Simulation data ####
load(file="sim.data")
obs <- imp1[!(imp1$treat==1 & imp1$momwhite==0),]

covs.cont.n=c("bw","b.head","preterm","birth.o","nnhealth","momage")
covs.cat.n=c("sex","twin","b.marr","mom.lths","mom.hs",	"mom.scoll","cig","first","booze","drugs","work.dur","prenatal","ark","ein","har","mia","pen","tex","was")
p=length(c(covs.cont.n,covs.cat.n))

### calculate pscores and weights for tot
Trt=obs$treat
form.qx=as.formula(obs[,c("treat",covs.cont.n,covs.cat.n)])
qx=glm(data=obs[,c("treat",covs.cont.n,covs.cat.n)],formula=form.qx,family=binomial)$fitted
wts=rep(1,nrow(obs))
# now *controls* get a weight of 1
wts[Trt==1]=(1-qx[Trt==1])/(qx[Trt==1])
# treated get weights equal to the probability of being untreated over
# probability of being treated (to weight them up to look like controls)

#### get data in right format for BART
xt=obs[,c(covs.cont.n,covs.cat.n,"treat")]
xt=as.matrix(xt)
xp1=xt[xt[,"treat"]==0,]
xp2=xp1
xp1[,ncol(xt)]=1
xp2[,ncol(xt)]=0
xp=rbind(xp1,xp2)

nc=sum(1-obs$treat)

##################### now simulate outcome data
##### covariate data, X
covs.ols = c(covs.cont.n,covs.cat.n)
X = obs[,covs.ols]
#X = na.omit(X)
# now standardize the continuous variables
X[,covs.cont.n]=as.data.frame(t((t(X[,covs.cont.n])-unlist(lapply(X[,covs.cont.n],mean)))/sqrt(unlist(lapply(X[covs.cont.n],var)))))


Hill11_surfaceB_limited <- function(X,z){
  N = nrow(X)
  dimx = ncol(X)
  Xmat = as.matrix(X)

  sigy = 1
  betaB = c(sample(c(.0,.1,.2,.3,.4),(dimx+1),replace=TRUE,prob=c(.6,.1,.1,.1,.1)))
  yb0hat = exp((cbind(rep(1, N), (Xmat+.5)) %*% betaB))
  yb1hat = cbind(rep(1, N), (Xmat+.5)) %*% betaB
  offset = c(mean(yb1hat[z==0] - yb0hat[z==0])) - 4
  yb1hat = cbind(rep(1, N), (Xmat+.5)) %*% betaB -offset
  YB0 = rnorm(N, yb0hat, sigy)
  YB1 = rnorm(N, yb1hat, sigy)

  YB = YB1; YB[z==0] = YB0[z==0]

  list(y = YB, X = Xmat ,z = z, tau = yb1hat - yb0hat, y1 = yb1hat, y0 = yb0hat )
}

## updating function ###########################

update.results <- function(index,result.item,testindex,
                           y_hat_all,
                           tau_hat_all,
                           tau_ci_all,
                           ate.ci.train,ate.ci.test,
                           y,z,y1,y0,pscore){
  threshold = 0.1
  y.hat.train = y_hat_all[-testindex]
  y.hat.test = y_hat_all[testindex]
  tau.hat.train = tau_hat_all[-testindex]
  tau.hat.test = tau_hat_all[testindex]

  z.train = z[-testindex]; z.test = z[testindex]
  y.train = y[-testindex]; y.test = y[testindex]

  tau_true = y1 - y0
  tau_true_train = tau_true[-testindex]; tau_true_test = tau_true[testindex]

  result.item["PEHE.all.train",index] = sqrt(mean((tau_true_train - tau.hat.train)^2))
  result.item["PEHE.tt.train",index]  = sqrt(mean((tau_true_train - tau.hat.train)[z.train=1]^2))
  result.item["PEHE.all.test",index] = sqrt(mean((tau_true_test - tau.hat.test)^2))
  result.item["PEHE.tt.test",index]  = sqrt(mean((tau_true_test - tau.hat.test)[z.test=1]^2))

  result.item["E.ate.train",index]  = (mean(tau_true_train) - mean(tau.hat.train))
  result.item["E.att.train",index]  = (mean((tau_true_train)[z.train==1]) - mean((tau.hat.train)[z.train==1]))
  result.item["E.ate.test",index]  = (mean(tau_true_test) - mean(tau.hat.test))
  result.item["E.att.test",index]  = (mean((tau_true_test)[z.test==1]) - mean((tau.hat.test)[z.test==1]))

  result.item["SE.ate.train",index]  = (mean(tau_true_train) - mean(tau.hat.train) )^2
  result.item["SE.att.train",index]  = (mean((tau_true_train)[z.train==1]) - mean((tau.hat.train)[z.train==1]))^2
  result.item["SE.ate.test",index]  = (mean(tau_true_test) - mean(tau.hat.test) )^2
  result.item["SE.att.test",index]  = (mean((tau_true_test)[z.test==1]) - mean((tau.hat.test)[z.test==1]))^2

  result.item["RMSE.all.train",index] = sqrt(mean((y.train - y.hat.train)^2))
  result.item["RMSE.tt.train",index]  = sqrt(mean((y.train - y.hat.train)[z.train==1]^2))
  result.item["RMSE.all.test",index] = sqrt(mean((y.test - y.hat.test)^2))
  result.item["RMSE.tt.test",index]  = sqrt(mean((y.test - y.hat.test)[z.test==1]^2))

  result.item["coverage.ite.train",index] = mean(((tau_ci_all[-testindex,1] < tau_true_train) * (tau_ci_all[-testindex,2] > tau_true_train)))
  result.item["coverage.ite.test",index] = mean(((tau_ci_all[testindex,1] < tau_true_test ) * (tau_ci_all[testindex,2] > tau_true_test)))

  result.item["coverage.ate.train",index] = ((ate.ci.train[1]< mean(tau_true_train)) * (ate.ci.train[2] >mean(tau_true_train)))
  result.item["coverage.ate.test",index] = ((ate.ci.test[1]< mean(tau_true_test)) * (ate.ci.test[2] >mean(tau_true_test)))

  result.item["range.ate.train",index] = (ate.ci.train[2] - ate.ci.train[1])
  result.item["range.ate.test",index] = (ate.ci.test[2] - ate.ci.test[1])

  result.item
}



## Simulation function ####
run <- function(i){
  #data generation within loop
  if(i<=500){set.seed(565 + i*5)}
  if(i>500){set.seed(7565 + i*5)}

  ## Set up result file ######
  mydimnames = list(c("PEHE.all.train","PEHE.tt.train","RMSE.all.train","RMSE.tt.train","E.ate.train","E.att.train","SE.ate.train","SE.att.train","coverage.ite.train","coverage.ate.train","range.ate.train",
                      "PEHE.all.test","PEHE.tt.test","RMSE.all.test","RMSE.tt.test","E.ate.test","E.att.test","SE.ate.test","SE.att.test","coverage.ite.test","coverage.ate.test","range.ate.test"),
                    c("BART","CF","VT-RF","CF-RF","GP","BCF","BART.PS","GP.PS","TP2.PS","TP2000.PS"))
  result.item = matrix(NaN, length(mydimnames[[1]]),length(mydimnames[[2]]),dimnames=list(mydimnames[[1]],mydimnames[[2]]) )
  ## Sampling data and train/test split ###########################
  mysample = Hill11_surfaceB_limited(X,Trt)
  Y = mysample$y; X = mysample$X; Z = mysample$z
  n = length(Y)
  idx.test = sample.int(n,ceiling(n*0.1))
  Y.train = Y[-idx.test]; Y.test = Y[idx.test]
  X.train = X[-idx.test,]; X.test = X[idx.test,]
  Z.train = Z[-idx.test]; Z.test = Z[idx.test]

  X.bart.train = cbind(X.train,Z=Z.train);
  X.bart.test = rbind(cbind(X,Z=1),cbind(X,Z=0));

  ## Propensity score estimation ###########################
  #Gradient Boosting Machines
  #myPSmodel = gbm::gbm(Z ~ .,data = Xdf,distribution="bernoulli",
  #                n.trees = 2000,
  #                shrinkage = .1,
  #                cv.folds = 20,
  #                n.cores = 1)
  #best.fit = gbm::gbm.perf(myPSmodel,method="cv")
  #pscore.gbm = gbm::predict(myPSmodel,n.trees=best.fit,type="response")

  #BART
  myPSbart = BayesTree::bart(x.train = X.train,y.train = Z.train ,x.test= X ,binaryOffset=mean(Z.train),ntree=200)
  pscore.bart = apply(pnorm(myPSbart$yhat.test),2,mean)

  pscore = pscore.bart
  pscore.train = pscore[-idx.test]; pscore.test = pscore[idx.test]

  X.bart.train.ps = cbind(X.bart.train,pscore = pscore.train);
  X.bart.test.ps = cbind(X.bart.test,pscore= rep(pscore,2));

  ## BART ################################
  bart.idx = (1:n)*Z + ((n+1):(2*n))*(1-Z)
  ## fit
  BART_fit = BayesTree::bart(x.train = X.bart.train, y.train = Y.train, x.test=X.bart.test,ntree=100,ndpost=3000,nskip=500,keepevery=1)

  ## obtain estimates
  y.hat.all = BART_fit$yhat.test.mean[bart.idx]
  tau.hat.all = BART_fit$yhat.test.mean[1:n] - BART_fit$yhat.test.mean[(n+1):(2*n)]

  ate.ci.train = {
    tmpset = BART_fit$yhat.test[,1:n] - BART_fit$yhat.test[,(n+1):(2*n)]
    tmpset = tmpset[,-idx.test]
    tmpvar = var(apply(tmpset,1,mean)) #variance over 1000 samples
    L = mean(tau.hat.all[-idx.test]) - 1.96 * sqrt(tmpvar)
    U = mean(tau.hat.all[-idx.test]) + 1.96 * sqrt(tmpvar)
    cbind(L,U)
  }

  ate.ci.test = {
    tmpset = BART_fit$yhat.test[,1:n] - BART_fit$yhat.test[,(n+1):(2*n)]
    tmpset = tmpset[,idx.test]
    tmpvar = var(apply(tmpset,1,mean)) #variance over 1000 samples
    L = mean(tau.hat.all[idx.test]) - 1.96 * sqrt(tmpvar)
    U = mean(tau.hat.all[idx.test]) + 1.96 * sqrt(tmpvar)
    cbind(L,U)
  }

  tau.ci.all = {
    tmpsort = apply(BART_fit$yhat.test[,1:n]-BART_fit$yhat.test[,(n+1):(2*n)],2,sort)
    idx975 = round(nrow(BART_fit$yhat.test)*0.975,0)
    idx025 = round(nrow(BART_fit$yhat.test)*0.025,0)
    L = tmpsort[idx025,]
    U = tmpsort[idx975,]
    cbind(L,U)
  }

  ## evaluate and store results
  result.item = update.results("BART",result.item,idx.test,
                               y_hat_all = y.hat.all,
                               tau_hat_all = tau.hat.all,
                               tau_ci_all = tau.ci.all,
                               ate.ci.train = ate.ci.train,
                               ate.ci.test = ate.ci.test,
                               y=Y,z=Z,y1 = mysample$y1, y0 = mysample$y0 ,pscore = pscore)

  rm(y.hat.all,tau.hat.all,ate.ci.train,ate.ci.test,tau.ci.all)
  ## BART with PS ####################################################
  ## fit
  BART_PS_fit = BayesTree::bart(x.train = X.bart.train.ps, y.train = Y.train, x.test=X.bart.test.ps,ntree=200,ndpost=3000,nskip=500,keepevery=1)

  ## obtain estimates
  y.hat.all = BART_PS_fit$yhat.test.mean[bart.idx]
  tau.hat.all = BART_PS_fit$yhat.test.mean[1:n] - BART_PS_fit$yhat.test.mean[(n+1):(2*n)]

  ate.ci.train = {
    tmpset = BART_PS_fit$yhat.test[,1:n] - BART_PS_fit$yhat.test[,(n+1):(2*n)]
    tmpset = tmpset[,-idx.test]
    tmpvar = var(apply(tmpset,1,mean)) #variance over 1000 samples
    L = mean(tau.hat.all[-idx.test]) - 1.96 * sqrt(tmpvar)
    U = mean(tau.hat.all[-idx.test]) + 1.96 * sqrt(tmpvar)
    cbind(L,U)
  }

  ate.ci.test = {
    tmpset = BART_PS_fit$yhat.test[,1:n] - BART_PS_fit$yhat.test[,(n+1):(2*n)]
    tmpset = tmpset[,idx.test]
    tmpvar = var(apply(tmpset,1,mean)) #variance over 1000 samples
    L = mean(tau.hat.all[idx.test]) - 1.96 * sqrt(tmpvar)
    U = mean(tau.hat.all[idx.test]) + 1.96 * sqrt(tmpvar)
    cbind(L,U)
  }

  tau.ci.all = {
    tmpsort = apply(BART_PS_fit$yhat.test[,1:n]-BART_PS_fit$yhat.test[,(n+1):(2*n)],2,sort)
    idx975 = round(nrow(BART_PS_fit$yhat.test)*0.975,0)
    idx025 = round(nrow(BART_PS_fit$yhat.test)*0.025,0)
    L = tmpsort[idx025,]
    U = tmpsort[idx975,]
    cbind(L,U)
  }

  ## evaluate and store results
  result.item = update.results("BART.PS",result.item,idx.test,
                               y_hat_all = y.hat.all,
                               tau_hat_all = tau.hat.all,
                               tau_ci_all = tau.ci.all,
                               ate.ci.train = ate.ci.train,
                               ate.ci.test = ate.ci.test,
                               y=Y,z=Z,y1 = mysample$y1, y0 = mysample$y0 ,pscore = pscore)

  rm(y.hat.all,tau.hat.all,ate.ci.train,ate.ci.test,tau.ci.all)
  ## BCF #####################################################################
  #no prior on variance of prior of the treatment effect (as in Hahn et al 2017)
  L=200
  L0 = ceiling(L * mean(Z.train))
  L1 = floor(L * mean(1-Z.train))
  BCF_fit0 = BayesTree::bart(x.train = cbind(X.train,pscore.train)[Z.train==0,],y.train = Y.train[Z.train==0],x.test=cbind(X,pscore),
                             ntree=L0,ndpost=3000,nskip=500,keepevery=1,power=3.0)
  BCF_fit1 = BayesTree::bart(x.train = cbind(X.train,pscore.train)[Z.train==1,],y.train = Y.train[Z.train==1],x.test=cbind(X,pscore),
                             ntree=L1,ndpost=3000,nskip=500,keepevery=1,power=3.0)


  ## obtain estimates
  y.hat.all = BCF_fit1$yhat.test.mean*Z + BCF_fit0$yhat.test.mean*(1-Z)
  tau.hat.all = BCF_fit1$yhat.test.mean - BCF_fit0$yhat.test.mean

  ate.ci.train = {
    tmpset = BCF_fit1$yhat.test - BCF_fit0$yhat.test
    tmpset = tmpset[,-idx.test]
    tmpvar = var(apply(tmpset,1,mean)) #variance over 1000 samples
    L = mean(tau.hat.all[-idx.test]) - 1.96 * sqrt(tmpvar)
    U = mean(tau.hat.all[-idx.test]) + 1.96 * sqrt(tmpvar)
    cbind(L,U)
  }

  ate.ci.test = {
    tmpset = BCF_fit1$yhat.test - BCF_fit0$yhat.test
    tmpset = tmpset[,idx.test]
    tmpvar = var(apply(tmpset,1,mean)) #variance over 1000 samples
    L = mean(tau.hat.all[idx.test]) - 1.96 * sqrt(tmpvar)
    U = mean(tau.hat.all[idx.test]) + 1.96 * sqrt(tmpvar)
    cbind(L,U)
  }

  tau.ci.all = {
    tmpsort = apply(BCF_fit1$yhat.test - BCF_fit0$yhat.test,2,sort)
    idx975 = round(nrow(BART_PS_fit$yhat.test)*0.975,0)
    idx025 = round(nrow(BART_PS_fit$yhat.test)*0.025,0)
    L = tmpsort[idx025,]
    U = tmpsort[idx975,]
    cbind(L,U)
  }

  ## evaluate and store results
  result.item = update.results("BCF",result.item,idx.test,
                               y_hat_all = y.hat.all,
                               tau_hat_all = tau.hat.all,
                               tau_ci_all = tau.ci.all,
                               ate.ci.train = ate.ci.train,
                               ate.ci.test = ate.ci.test,
                               y=Y,z=Z,y1 = mysample$y1, y0 = mysample$y0 ,pscore = pscore)

  rm(y.hat.all,tau.hat.all,ate.ci.train,ate.ci.test,tau.ci.all)
  ## GP #####################################################################
  GP_fit = CausalStump::CausalStump(y = Y.train,X = data.frame(X.train), z = Z.train,minimumRMSE=TRUE,
                                    maxiter=5000,tol=1e-2,learning_rate = 0.01,prior=FALSE,myoptim = "Nadam")

  ## predict
  GS_pred  =   predict(GP_fit, X = data.frame(X),z=Z)

  y.hat.all   = GS_pred$map
  tau.hat.all = rep(NaN,n)
  tau.ci.all  = matrix(NaN,n,2)

  GS_treat = treatment(GP_fit, X = data.frame(X.train))
  tau.hat.all[-idx.test] = GS_treat$map
  tau.ci.all[-idx.test,] = GS_treat$ci
  ate.ci.train= GS_treat$ate_ci

  GS_treat = treatment(GP_fit, X = data.frame(X.test))
  tau.hat.all[idx.test] = GS_treat$map
  tau.ci.all[idx.test,] = GS_treat$ci
  ate.ci.test = GS_treat$ate_ci

  ## evaluate and store results
  result.item = update.results("GP",result.item,idx.test,
                               y_hat_all = y.hat.all,
                               tau_hat_all = tau.hat.all,
                               tau_ci_all = tau.ci.all,
                               ate.ci.train = ate.ci.train,
                               ate.ci.test = ate.ci.test,
                               y=Y,z=Z,y1 = mysample$y1, y0 = mysample$y0 ,pscore = pscore)

  rm(y.hat.all,tau.hat.all,ate.ci.train,ate.ci.test,tau.ci.all)

  ## GP with PS #####################################################################
  GP_PS_fit = CausalStump::CausalStump(y = Y.train,X = data.frame(X.train), z = Z.train,pscore=pscore.train,minimumRMSE=TRUE,
                                       maxiter=5000,tol=1e-2,learning_rate = 0.01,prior=FALSE,myoptim = "Nadam")

  #exp(GP_PS_fit$Kernel$parameters$Lm[26])
  #exp(GP_PS_fit$Kernel$parameters$La[26])

  ## predict
  GS_PS_pred  = predict(GP_PS_fit, X = data.frame(X),z=Z,pscore=pscore)

  y.hat.all   = GS_PS_pred$map

  tau.hat.all = rep(NaN,n)
  tau.ci.all  = matrix(NaN,n,2)

  GS_PS_treat = treatment(GP_PS_fit, X = data.frame(X.train),pscore=pscore.train)
  tau.hat.all[-idx.test] = GS_PS_treat$map
  tau.ci.all[-idx.test,] = GS_PS_treat$ci
  ate.ci.train= GS_PS_treat$ate_ci

  GS_PS_treat = treatment(GP_PS_fit, X = data.frame(X.test),pscore=pscore.test)
  tau.hat.all[idx.test] = GS_PS_treat$map
  tau.ci.all[idx.test,] = GS_PS_treat$ci
  ate.ci.test = GS_PS_treat$ate_ci

  ## evaluate and store results
  result.item = update.results("GP.PS",result.item,idx.test,
                               y_hat_all = y.hat.all,
                               tau_hat_all = tau.hat.all,
                               tau_ci_all = tau.ci.all,
                               ate.ci.train = ate.ci.train,
                               ate.ci.test = ate.ci.test,
                               y=Y,z=Z,y1 = mysample$y1, y0 = mysample$y0 ,pscore = pscore)

  rm(y.hat.all,tau.hat.all,ate.ci.train,ate.ci.test,tau.ci.all)

  ## TP2 with PS #####################################################################
  TP2_PS_fit = CausalStump::CausalStump(y = Y.train,X = data.frame(X.train), z = Z.train,pscore=pscore.train,minimumRMSE=TRUE,
                                        maxiter=5000,tol=1e-2,learning_rate = 0.1,
                                        prior=TRUE,nu=2,
                                        myoptim = "Nadam")
  ## predict
  TP2_PS_pred  = predict(TP2_PS_fit, X = data.frame(X),z=Z,pscore=pscore)

  y.hat.all   = TP2_PS_pred$map

  tau.hat.all = rep(NaN,n)
  tau.ci.all  = matrix(NaN,n,2)

  TP2_PS_treat = treatment(TP2_PS_fit, X = data.frame(X.train),pscore=pscore.train)
  tau.hat.all[-idx.test] = TP2_PS_treat$map
  tau.ci.all[-idx.test,] = TP2_PS_treat$ci
  ate.ci.train= TP2_PS_treat$ate_ci

  TP2_PS_treat = treatment(TP2_PS_fit, X = data.frame(X.test),pscore=pscore.test)
  tau.hat.all[idx.test] = TP2_PS_treat$map
  tau.ci.all[idx.test,] = TP2_PS_treat$ci
  ate.ci.test = TP2_PS_treat$ate_ci

  ## evaluate and store results
  result.item = update.results("TP2.PS",result.item,idx.test,
                               y_hat_all = y.hat.all,
                               tau_hat_all = tau.hat.all,
                               tau_ci_all = tau.ci.all,
                               ate.ci.train = ate.ci.train,
                               ate.ci.test = ate.ci.test,
                               y=Y,z=Z,y1 = mysample$y1, y0 = mysample$y0 ,pscore = pscore)

  rm(y.hat.all,tau.hat.all,ate.ci.train,ate.ci.test,tau.ci.all)
  ## TP2000 with PS #####################################################################
  TP2000_PS_fit = CausalStump::CausalStump(y = Y.train,X = data.frame(X.train), z = Z.train,pscore=pscore.train,minimumRMSE=TRUE,
                                           maxiter=5000,tol=1e-2,learning_rate = 0.05,
                                           prior=TRUE,nu=2000,
                                           myoptim = "Nadam")

  ## predict
  TP2000_PS_pred  = predict(TP2000_PS_fit, X = data.frame(X),z=Z,pscore=pscore)

  y.hat.all   = TP2000_PS_pred$map

  tau.hat.all = rep(NaN,n)
  tau.ci.all  = matrix(NaN,n,2)

  TP2000_PS_treat = treatment(TP2000_PS_fit, X = data.frame(X.train),pscore=pscore.train)
  tau.hat.all[-idx.test] = TP2000_PS_treat$map
  tau.ci.all[-idx.test,] = TP2000_PS_treat$ci
  ate.ci.train= TP2000_PS_treat$ate_ci

  TP2000_PS_treat = treatment(TP2000_PS_fit, X = data.frame(X.test),pscore=pscore.test)
  tau.hat.all[idx.test] = TP2000_PS_treat$map
  tau.ci.all[idx.test,] = TP2000_PS_treat$ci
  ate.ci.test = TP2000_PS_treat$ate_ci

  ## evaluate and store results
  result.item = update.results("TP2000.PS",result.item,idx.test,
                               y_hat_all = y.hat.all,
                               tau_hat_all = tau.hat.all,
                               tau_ci_all = tau.ci.all,
                               ate.ci.train = ate.ci.train,
                               ate.ci.test = ate.ci.test,
                               y=Y,z=Z,y1 = mysample$y1, y0 = mysample$y0 ,pscore = pscore)

  rm(y.hat.all,tau.hat.all,ate.ci.train,ate.ci.test,tau.ci.all)

  ## CF #####################################################################
  ## fit
  CF_forest = grf::causal_forest(X = X.train,Y =  Y.train, W = Z.train, num.trees = 1000,precompute.nuisance = TRUE,num.threads=1,seed=(123+i))

  ## predict
  CF_pred   = predict(CF_forest, newdata = X, estimate.variance = TRUE,num.threads=1)

  y.hat.all = rep(NA,n)
  y.hat.all[-idx.test] = CF_forest$Y.hat

  tau.hat.all = CF_pred$predictions

  tau.ci.all  = cbind(tau.hat.all - 1.96 *sqrt(CF_pred$variance.estimates),
                      tau.hat.all + 1.96 *sqrt(CF_pred$variance.estimates))
  #cannot get confidence interval for CF

  ate.ci.train = { tmp.ate = grf::estimate_average_effect(CF_forest, target.sample = "all");
  c(tmp.ate["estimate"] - 1.96 * tmp.ate["std.err"],
    tmp.ate["estimate"] + 1.96 * tmp.ate["std.err"])}

  ate.ci.test = NA

  ## evaluate and store results
  result.item = update.results("CF",result.item,idx.test,
                               y_hat_all = y.hat.all,
                               tau_hat_all = tau.hat.all,
                               tau_ci_all = tau.ci.all,
                               ate.ci.train = ate.ci.train,
                               ate.ci.test = ate.ci.test,
                               y=Y,z=Z,y1 = mysample$y1, y0 = mysample$y0 ,pscore = pscore)

  rm(y.hat.all,tau.hat.all,ate.ci.train,ate.ci.test,tau.ci.all)

  ## VT-RF #####################################################################
  ## fit
  VTRF_fit = ranger::ranger(Y ~ ., data = data.frame(Y=Y.train,X=X.train,Z=Z.train),num.trees = 5000,num.threads=1,mtry=10)

  VTRF_pred1 = predict(VTRF_fit, data = data.frame(X=X,Z=1) , predict.all = TRUE,num.threads=1)
  VTRF_pred0 = predict(VTRF_fit, data = data.frame(X=X,Z=0) , predict.all = TRUE,num.threads=1)

  y.hat.all = apply(VTRF_pred1$predictions,1,mean) * Z + apply(VTRF_pred0$predictions,1,mean) * (1-Z)
  tau.hat.all = apply(VTRF_pred1$predictions,1,mean) - apply(VTRF_pred0$predictions,1,mean)

  ate.ci.train = {
    tmpset = VTRF_pred1$predictions[-idx.test,] - VTRF_pred0$predictions[-idx.test,]
    tmpvar = var(apply(tmpset,2,mean)) #variance over 1000 samples
    L = mean(tau.hat.all[-idx.test]) - 1.96 * sqrt(tmpvar)
    U = mean(tau.hat.all[-idx.test]) + 1.96 * sqrt(tmpvar)
    cbind(L,U)
  }

  ate.ci.test = {
    tmpset = VTRF_pred1$predictions[idx.test,] - VTRF_pred0$predictions[idx.test,]
    tmpvar = var(apply(tmpset,2,mean)) #variance over 1000 samples
    L = mean(tau.hat.all[idx.test]) - 1.96 * sqrt(tmpvar)
    U = mean(tau.hat.all[idx.test]) + 1.96 * sqrt(tmpvar)
    cbind(L,U)
  }

  tau.ci.all = {
    tmpmat = VTRF_pred1$predictions - VTRF_pred0$predictions
    tmpsort = apply(tmpmat,1,sort)
    idx975 = round(ncol(VTRF_pred1$predictions)*0.975,0)
    idx025 = round(ncol(VTRF_pred1$predictions)*0.025,0)
    L = tmpsort[idx025,]
    U = tmpsort[idx975,]
    cbind(L,U)
  }

  ## evaluate and store results
  result.item = update.results("VT-RF",result.item,idx.test,
                               y_hat_all = y.hat.all,
                               tau_hat_all = tau.hat.all,
                               tau_ci_all = tau.ci.all,
                               ate.ci.train = ate.ci.train,
                               ate.ci.test = ate.ci.test,
                               y=Y,z=Z,y1 = mysample$y1, y0 = mysample$y0 ,pscore = pscore)

  rm(y.hat.all,tau.hat.all,ate.ci.train,ate.ci.test,tau.ci.all)
  ## CF-RF #####################################################################
  ## fit
  CFRF_fit1 = ranger::ranger(Y ~ ., data = data.frame(Y=Y.train,X=X.train)[Z.train==1,],num.trees = 4000,num.threads=1,mtry=10)
  CFRF_fit0 = ranger::ranger(Y ~ ., data = data.frame(Y=Y.train,X=X.train)[Z.train==0,],num.trees = 4000,num.threads=1,mtry=10)

  CFRF_pred1 = predict(CFRF_fit1, data = data.frame(X=X) , predict.all = TRUE,num.threads=1)
  CFRF_pred0 = predict(CFRF_fit0, data = data.frame(X=X) , predict.all = TRUE,num.threads=1)

  y.hat.all = apply(CFRF_pred1$predictions,1,mean) * Z + apply(CFRF_pred0$predictions,1,mean) * (1-Z)
  tau.hat.all = apply(CFRF_pred1$predictions,1,mean) - apply(CFRF_pred0$predictions,1,mean)

  ate.ci.train = {
    tmpset = CFRF_pred1$predictions[-idx.test,] - CFRF_pred0$predictions[-idx.test,]
    tmpvar = var(apply(tmpset,2,mean)) #variance over 1000 samples
    L = mean(tau.hat.all[-idx.test]) - 1.96 * sqrt(tmpvar)
    U = mean(tau.hat.all[-idx.test]) + 1.96 * sqrt(tmpvar)
    cbind(L,U)
  }

  ate.ci.test = {
    tmpset = CFRF_pred1$predictions[idx.test,] - CFRF_pred0$predictions[idx.test,]
    tmpvar = var(apply(tmpset,2,mean)) #variance over 1000 samples
    L = mean(tau.hat.all[idx.test]) - 1.96 * sqrt(tmpvar)
    U = mean(tau.hat.all[idx.test]) + 1.96 * sqrt(tmpvar)
    cbind(L,U)
  }

  tau.ci.all = {
    tmpmat = CFRF_pred1$predictions - CFRF_pred0$predictions
    tmpsort = apply(tmpmat,1,sort)
    idx975 = round(ncol(CFRF_pred1$predictions)*0.975,0)
    idx025 = round(ncol(CFRF_pred1$predictions)*0.025,0)
    L = tmpsort[idx025,]
    U = tmpsort[idx975,]
    cbind(L,U)
  }

  ## evaluate and store results
  result.item = update.results("CF-RF",result.item,idx.test,
                               y_hat_all = y.hat.all,
                               tau_hat_all = tau.hat.all,
                               tau_ci_all = tau.ci.all,
                               ate.ci.train = ate.ci.train,
                               ate.ci.test = ate.ci.test,
                               y=Y,z=Z,y1 = mysample$y1, y0 = mysample$y0 ,pscore = pscore)

  rm(y.hat.all,tau.hat.all,ate.ci.train,ate.ci.test,tau.ci.all)
  result.item
}

## Parallel setup #####
#system.time(sapply(1,run))

max_iter=100
nr_cores = detectCores()
cat(sprintf("\n-----\nIterations: %d on %d cores\n",max_iter,nr_cores))
cl <- makeCluster(nr_cores,type="FORK")
system.time(results <- parSapplyLB(cl,1:max_iter,run)) #clusterApplyLB, parSapply
stopCluster(cl)


save.image(file = paste("H11_",max_iter,"_limited.RData",sep=""))
