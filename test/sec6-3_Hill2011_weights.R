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


## now create matrix of all interactions etc for third response surface C ####
#ytmp=rnorm(N)
#mod.bal <- glm(formula=ytmp~(bw+b.head+preterm+birth.o+nnhealth+momage+sex+twin+b.marr+mom.lths+mom.hs+mom.scoll+cig+first+booze+drugs+work.dur+prenatal+ark+ein+har+mia+pen+tex+was)^2 + I(bw^2) + I(b.head^2) + I(preterm^2) + I(birth.o^2) + I(nnhealth^2) + I(momage^2),x=T,data=cbind.data.frame(Xmat))
#coefs <- mod.bal$coef[-1]
#XX <- mod.bal$x[,-1]
#XX <- XX[,!is.na(coefs)]

#nouts=3
#os=c("YA","YB","YC")

# data with constant, only for response surface C
#XXXmat=cbind(rep(1,N),XX)
#rm(XX)

## updating function ###########################

update.results <- function(index,result.item,testindex,discard_idx,
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

  retain.train = discard_idx[-testindex]==0
  retain.test = discard_idx[testindex]==0

  tau_true_train = tau_true_train[retain.train]
  tau.hat.train = tau.hat.train[retain.train]
  z.train = z.train[retain.train]

  tau_true_test = tau_true_test[retain.test]
  tau.hat.test = tau.hat.test[retain.test]
  z.test = z.test[retain.test]

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

  y.train = y.train[retain.train]
  y.hat.train = y.hat.train[retain.train]

  y.test = y.test[retain.test]
  y.hat.test = y.hat.test[retain.test]

  result.item["RMSE.all.train",index] = sqrt(mean((y.train - y.hat.train)^2))
  result.item["RMSE.tt.train",index]  = sqrt(mean((y.train - y.hat.train)[z.train==1]^2))
  result.item["RMSE.all.test",index] = sqrt(mean((y.test - y.hat.test)^2))
  result.item["RMSE.tt.test",index]  = sqrt(mean((y.test - y.hat.test)[z.test==1]^2))

  result.item["coverage.ite.train",index] = mean((((tau_ci_all[-testindex,1])[retain.train] < tau_true_train) * ((tau_ci_all[-testindex,2])[retain.train] > tau_true_train)))
  result.item["coverage.ite.test",index] = mean((((tau_ci_all[testindex,1])[retain.test] < tau_true_test ) * ((tau_ci_all[testindex,2])[retain.test] > tau_true_test)))

  result.item["coverage.ate.train",index] = ((ate.ci.train[1]< mean(tau_true_train)) * (ate.ci.train[2] >mean(tau_true_train)))
  result.item["coverage.ate.test",index] = ((ate.ci.test[1]< mean(tau_true_test)) * (ate.ci.test[2] >mean(tau_true_test)))

  result.item["range.ate.train",index] = (ate.ci.train[2] - ate.ci.train[1])
  result.item["range.ate.test",index] = (ate.ci.test[2] - ate.ci.test[1])

  result.item["discard percent train",index] = mean(discard_idx[-testindex])
  result.item["discard percent test",index] = mean(discard_idx[testindex])

  result.item
}



## Simulation function ####
run <- function(i){
  #data generation within loop
  if(i<=500){set.seed(565 + i*5)}
  if(i>500){set.seed(7565 + i*5)}

  ## Set up result file ######
  mydimnames = list(c("PEHE.all.train","PEHE.tt.train","RMSE.all.train","RMSE.tt.train","E.ate.train","E.att.train","SE.ate.train","SE.att.train","coverage.ite.train","coverage.ate.train","range.ate.train","discard percent train",
                      "PEHE.all.test","PEHE.tt.test","RMSE.all.test","RMSE.tt.test","E.ate.test","E.att.test","SE.ate.test","SE.att.test","coverage.ite.test","coverage.ate.test","range.ate.test","discard percent test"),
                    c("BART w R1","BART w R2","BART w R3","GP w R1","GP w R2","GP w R3","GP weights","GP balancing","GP variance 1 weights","GP variance 2 weights","BART w GP-R2","BART w GP-R3","GP w BART R2"))
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


 #GP - variance 2 - IPW ######
 myweights = Z.train*pscore.train^2*(1-pscore.train)+(1-Z.train)*pscore.train*(1-pscore.train)^2
 GP_fit = CausalStump::CausalStump(y = Y.train,X = data.frame(X.train), z = Z.train,w = myweights,minimumRMSE=FALSE,
                                   maxiter=5000,tol=1e-2,learning_rate = 0.01,prior=FALSE,myoptim = "Nadam")

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

 DISCARD_NOTHING = rep(0,n)

 ## evaluate and store results
 result.item = update.results("GP variance 2 weights",result.item,idx.test,DISCARD_NOTHING,
                              y_hat_all = y.hat.all,
                              tau_hat_all = tau.hat.all,
                              tau_ci_all = tau.ci.all,
                              ate.ci.train = ate.ci.train,
                              ate.ci.test = ate.ci.test,
                              y=Y,z=Z,y1 = mysample$y1, y0 = mysample$y0 ,pscore = pscore)
 print(result.item)

  ## GP #####################################################################
  GP_fit <- CausalStump::CausalStump(y = Y.train,X = data.frame(X.train), z = Z.train,minimumRMSE=FALSE,
                                    maxiter=5000,tol=1e-2,learning_rate = 0.01,prior=FALSE,myoptim = "Nadam")

  ## predict
  GS_pred  =   predict(GP_fit, X = data.frame(X),z=Z)

  y.hat.all   = GS_pred$map
  tau.hat.all = rep(NaN,n)
  tau.ci.all  = matrix(NaN,n,2)

  GS_treat = treatment(GP_fit, X = data.frame(X.train))
  tau.hat.all[-idx.test] = GS_treat$map
  tau.ci.all[-idx.test,] = GS_treat$ci

  GS_treat = treatment(GP_fit, X = data.frame(X.test))
  tau.hat.all[idx.test] = GS_treat$map
  tau.ci.all[idx.test,] = GS_treat$ci

  gp_propsigma = tau.ci.all[,2]-tau.ci.all[,1]
  mgp = median(gp_propsigma)
  #mgp = quantile(gp_propsigma)[4]
  #mgp = mean(gp_propsigma)


  ## GP - Rule 1 - 10% ######
  GP_DISCARD_R1 = (gp_propsigma/mgp)^2 > 2.706
  #mean(GP_DISCARD_R1)

  GS_treat = treatment(GP_fit, X = data.frame(X.train[GP_DISCARD_R1[-idx.test]==0,]))
  ate.ci.train= GS_treat$ate_ci
  GS_treat = treatment(GP_fit, X = data.frame(X.test[GP_DISCARD_R1[idx.test]==0,]))
  ate.ci.test = GS_treat$ate_ci

  result.item = update.results("GP w R1",result.item,idx.test,GP_DISCARD_R1,
                               y_hat_all = y.hat.all,
                               tau_hat_all = tau.hat.all,
                               tau_ci_all = tau.ci.all,
                               ate.ci.train = ate.ci.train,
                               ate.ci.test = ate.ci.test,
                               y=Y,z=Z, y1 = mysample$y1, y0 = mysample$y0 ,pscore = pscore)

  ## GP - Rule 2 - 5% ######
  GP_DISCARD_R2 = (gp_propsigma/mgp)^2 > 3.841
  #mean(GP_DISCARD_R2)

  GS_treat = treatment(GP_fit, X = data.frame(X.train[GP_DISCARD_R2[-idx.test]==0,]))
  ate.ci.train= GS_treat$ate_ci
  GS_treat = treatment(GP_fit, X = data.frame(X.test[GP_DISCARD_R2[idx.test]==0,]))
  ate.ci.test = GS_treat$ate_ci

  result.item = update.results("GP w R2",result.item,idx.test,GP_DISCARD_R2,
                               y_hat_all = y.hat.all,
                               tau_hat_all = tau.hat.all,
                               tau_ci_all = tau.ci.all,
                               ate.ci.train = ate.ci.train,
                               ate.ci.test = ate.ci.test,
                               y=Y,z=Z, y1 = mysample$y1, y0 = mysample$y0 ,pscore = pscore)

  ## GP - Rule 3 - 1% ######
  GP_DISCARD_R3 = (gp_propsigma/mgp)^2 > 6.635
  #mean(GP_DISCARD_R2)

  GS_treat = treatment(GP_fit, X = data.frame(X.train[GP_DISCARD_R3[-idx.test]==0,]))
  ate.ci.train= GS_treat$ate_ci
  GS_treat = treatment(GP_fit, X = data.frame(X.test[GP_DISCARD_R3[idx.test]==0,]))
  ate.ci.test = GS_treat$ate_ci

  result.item = update.results("GP w R3",result.item,idx.test,GP_DISCARD_R3,
                               y_hat_all = y.hat.all,
                               tau_hat_all = tau.hat.all,
                               tau_ci_all = tau.ci.all,
                               ate.ci.train = ate.ci.train,
                               ate.ci.test = ate.ci.test,
                               y=Y,z=Z, y1 = mysample$y1, y0 = mysample$y0 ,pscore = pscore)


   ## BART with PS ####################################################
   ## fit
   bart.idx = (1:n)*Z + ((n+1):(2*n))*(1-Z)
   BART_PS_fit = BayesTree::bart(x.train = X.bart.train.ps, y.train = Y.train, x.test=X.bart.test.ps,ntree=200,ndpost=3000,nskip=500,keepevery=1)

   ## obtain estimates
   y.hat.all = BART_PS_fit$yhat.test.mean[bart.idx] # 0.5*(y.hat.all+ )
   tau.hat.all = BART_PS_fit$yhat.test.mean[1:n] - BART_PS_fit$yhat.test.mean[(n+1):(2*n)] #0.5*(tau.hat.all + )

   tau.ci.all = {
     tmpsort = apply(BART_PS_fit$yhat.test[,1:n]-BART_PS_fit$yhat.test[,(n+1):(2*n)],2,sort)
     idx975 = round(nrow(BART_PS_fit$yhat.test)*0.975,0)
     idx025 = round(nrow(BART_PS_fit$yhat.test)*0.025,0)
     L = tmpsort[idx025,]
     U = tmpsort[idx975,]
     cbind(L,U)
   }

  #c("BART w R1","BART w R2","BART w R3","GP w R1","GP w R2","GP w R3")

  ## BART - Rule 1 - 1SD ######
  bart.idx

  bart.s1 = sqrt(apply(BART_PS_fit$yhat.test[,1:n],2,var))
  bart.s0 = sqrt(apply(BART_PS_fit$yhat.test[,(n+1):(2*n)],2,var))
  sds1=sqrt(var(bart.s1))
  sds0=sqrt(var(bart.s0))
  m1 = max(bart.s1)
  m0 = max(bart.s0)

  tmp0 = bart.s0 > m1 + 2*sds1
  tmp1 = bart.s1 > m0 + 2*sds0

  BART_DISCARD_R1 = (Z*tmp0 + (1-Z)*tmp1)

  ate.ci.train = {
    tmpset = BART_PS_fit$yhat.test[,1:n] - BART_PS_fit$yhat.test[,(n+1):(2*n)]
    tmpset = tmpset[,-idx.test]
    tmpset = tmpset[,BART_DISCARD_R1[-idx.test]==0]
    tmpvar = var(apply(tmpset,1,mean)) #variance over 1000 samples
    L = mean((tau.hat.all[-idx.test])[BART_DISCARD_R1[-idx.test]==0]) - 1.96 * sqrt(tmpvar)
    U = mean((tau.hat.all[-idx.test])[BART_DISCARD_R1[-idx.test]==0]) + 1.96 * sqrt(tmpvar)
    cbind(L,U)
  }

  ate.ci.test = {
    tmpset = BART_PS_fit$yhat.test[,1:n] - BART_PS_fit$yhat.test[,(n+1):(2*n)]
    tmpset = tmpset[,idx.test]
    tmpset = tmpset[,BART_DISCARD_R1[idx.test]==0]
    tmpvar = var(apply(tmpset,1,mean)) #variance over 1000 samples
    L = mean((tau.hat.all[idx.test])[BART_DISCARD_R1[idx.test]==0]) - 1.96 * sqrt(tmpvar)
    U = mean((tau.hat.all[idx.test])[BART_DISCARD_R1[idx.test]==0]) + 1.96 * sqrt(tmpvar)
    cbind(L,U)
  }

  result.item = update.results("BART w R1",result.item,idx.test,BART_DISCARD_R1,
                               y_hat_all = y.hat.all,
                               tau_hat_all = tau.hat.all,
                               tau_ci_all = tau.ci.all,
                               ate.ci.train = ate.ci.train,
                               ate.ci.test = ate.ci.test,
                               y=Y,z=Z, y1 = mysample$y1, y0 = mysample$y0 ,pscore = pscore)

  ## BART - Rule 2 - 10% ######

  tmp0 = (bart.s1/bart.s0)^2 > 2.706
  tmp1 = (bart.s0/bart.s1)^2 > 2.706

  BART_DISCARD_R2 = Z*tmp1 + (1-Z)*tmp0

  ate.ci.train = {
    tmpset = BART_PS_fit$yhat.test[,1:n] - BART_PS_fit$yhat.test[,(n+1):(2*n)]
    tmpset = tmpset[,-idx.test]
    tmpset = tmpset[,BART_DISCARD_R2[-idx.test]==0]
    tmpvar = var(apply(tmpset,1,mean)) #variance over 1000 samples
    L = mean((tau.hat.all[-idx.test])[BART_DISCARD_R2[-idx.test]==0]) - 1.96 * sqrt(tmpvar)
    U = mean((tau.hat.all[-idx.test])[BART_DISCARD_R2[-idx.test]==0]) + 1.96 * sqrt(tmpvar)
    cbind(L,U)
  }

  ate.ci.test = {
    tmpset = BART_PS_fit$yhat.test[,1:n] - BART_PS_fit$yhat.test[,(n+1):(2*n)]
    tmpset = tmpset[,idx.test]
    tmpset = tmpset[,BART_DISCARD_R2[idx.test]==0]
    tmpvar = var(apply(tmpset,1,mean)) #variance over 1000 samples
    L = mean((tau.hat.all[idx.test])[BART_DISCARD_R2[idx.test]==0]) - 1.96 * sqrt(tmpvar)
    U = mean((tau.hat.all[idx.test])[BART_DISCARD_R2[idx.test]==0]) + 1.96 * sqrt(tmpvar)
    cbind(L,U)
  }

  result.item = update.results("BART w R2",result.item,idx.test,BART_DISCARD_R2,
                               y_hat_all = y.hat.all,
                               tau_hat_all = tau.hat.all,
                               tau_ci_all = tau.ci.all,
                               ate.ci.train = ate.ci.train,
                               ate.ci.test = ate.ci.test,
                               y=Y,z=Z, y1 = mysample$y1, y0 = mysample$y0 ,pscore = pscore)

  ## BART - Rule 3 - 5% ######

  tmp0 = (bart.s1/bart.s0)^2 > 3.841
  tmp1 = (bart.s0/bart.s1)^2 > 3.841

  BART_DISCARD_R3 = Z*tmp1 + (1-Z)*tmp0

  ate.ci.train = {
    tmpset = BART_PS_fit$yhat.test[,1:n] - BART_PS_fit$yhat.test[,(n+1):(2*n)]
    tmpset = tmpset[,-idx.test]
    tmpset = tmpset[,BART_DISCARD_R3[-idx.test]==0]
    tmpvar = var(apply(tmpset,1,mean)) #variance over 1000 samples
    L = mean((tau.hat.all[-idx.test])[BART_DISCARD_R3[-idx.test]==0]) - 1.96 * sqrt(tmpvar)
    U = mean((tau.hat.all[-idx.test])[BART_DISCARD_R3[-idx.test]==0]) + 1.96 * sqrt(tmpvar)
    cbind(L,U)
  }

  ate.ci.test = {
    tmpset = BART_PS_fit$yhat.test[,1:n] - BART_PS_fit$yhat.test[,(n+1):(2*n)]
    tmpset = tmpset[,idx.test]
    tmpset = tmpset[,BART_DISCARD_R3[idx.test]==0]
    tmpvar = var(apply(tmpset,1,mean)) #variance over 1000 samples
    L = mean((tau.hat.all[idx.test])[BART_DISCARD_R3[idx.test]==0]) - 1.96 * sqrt(tmpvar)
    U = mean((tau.hat.all[idx.test])[BART_DISCARD_R3[idx.test]==0]) + 1.96 * sqrt(tmpvar)
    cbind(L,U)
  }

  result.item = update.results("BART w R3",result.item,idx.test,BART_DISCARD_R3,
                               y_hat_all = y.hat.all,
                               tau_hat_all = tau.hat.all,
                               tau_ci_all = tau.ci.all,
                               ate.ci.train = ate.ci.train,
                               ate.ci.test = ate.ci.test,
                               y=Y,z=Z, y1 = mysample$y1, y0 = mysample$y0 ,pscore = pscore)

  ## BART - GP Rule 2 ######

  BART_DISCARD_R3 = GP_DISCARD_R2

  ate.ci.train = {
    tmpset = BART_PS_fit$yhat.test[,1:n] - BART_PS_fit$yhat.test[,(n+1):(2*n)]
    tmpset = tmpset[,-idx.test]
    tmpset = tmpset[,BART_DISCARD_R3[-idx.test]==0]
    tmpvar = var(apply(tmpset,1,mean)) #variance over 1000 samples
    L = mean((tau.hat.all[-idx.test])[BART_DISCARD_R3[-idx.test]==0]) - 1.96 * sqrt(tmpvar)
    U = mean((tau.hat.all[-idx.test])[BART_DISCARD_R3[-idx.test]==0]) + 1.96 * sqrt(tmpvar)
    cbind(L,U)
  }

  ate.ci.test = {
    tmpset = BART_PS_fit$yhat.test[,1:n] - BART_PS_fit$yhat.test[,(n+1):(2*n)]
    tmpset = tmpset[,idx.test]
    tmpset = tmpset[,BART_DISCARD_R3[idx.test]==0]
    tmpvar = var(apply(tmpset,1,mean)) #variance over 1000 samples
    L = mean((tau.hat.all[idx.test])[BART_DISCARD_R3[idx.test]==0]) - 1.96 * sqrt(tmpvar)
    U = mean((tau.hat.all[idx.test])[BART_DISCARD_R3[idx.test]==0]) + 1.96 * sqrt(tmpvar)
    cbind(L,U)
  }

  result.item = update.results("BART w GP-R2",result.item,idx.test,BART_DISCARD_R3,
                               y_hat_all = y.hat.all,
                               tau_hat_all = tau.hat.all,
                               tau_ci_all = tau.ci.all,
                               ate.ci.train = ate.ci.train,
                               ate.ci.test = ate.ci.test,
                               y=Y,z=Z, y1 = mysample$y1, y0 = mysample$y0 ,pscore = pscore)

  ## BART - GP Rule 3 ######

  BART_DISCARD_R3 = GP_DISCARD_R3

  ate.ci.train = {
    tmpset = BART_PS_fit$yhat.test[,1:n] - BART_PS_fit$yhat.test[,(n+1):(2*n)]
    tmpset = tmpset[,-idx.test]
    tmpset = tmpset[,BART_DISCARD_R3[-idx.test]==0]
    tmpvar = var(apply(tmpset,1,mean)) #variance over 1000 samples
    L = mean((tau.hat.all[-idx.test])[BART_DISCARD_R3[-idx.test]==0]) - 1.96 * sqrt(tmpvar)
    U = mean((tau.hat.all[-idx.test])[BART_DISCARD_R3[-idx.test]==0]) + 1.96 * sqrt(tmpvar)
    cbind(L,U)
  }

  ate.ci.test = {
    tmpset = BART_PS_fit$yhat.test[,1:n] - BART_PS_fit$yhat.test[,(n+1):(2*n)]
    tmpset = tmpset[,idx.test]
    tmpset = tmpset[,BART_DISCARD_R3[idx.test]==0]
    tmpvar = var(apply(tmpset,1,mean)) #variance over 1000 samples
    L = mean((tau.hat.all[idx.test])[BART_DISCARD_R3[idx.test]==0]) - 1.96 * sqrt(tmpvar)
    U = mean((tau.hat.all[idx.test])[BART_DISCARD_R3[idx.test]==0]) + 1.96 * sqrt(tmpvar)
    cbind(L,U)
  }

  result.item = update.results("BART w GP-R3",result.item,idx.test,BART_DISCARD_R3,
                               y_hat_all = y.hat.all,
                               tau_hat_all = tau.hat.all,
                               tau_ci_all = tau.ci.all,
                               ate.ci.train = ate.ci.train,
                               ate.ci.test = ate.ci.test,
                               y=Y,z=Z, y1 = mysample$y1, y0 = mysample$y0 ,pscore = pscore)

  rm(y.hat.all,tau.hat.all,ate.ci.train,ate.ci.test,tau.ci.all)

  ## GP - BART Rule 10% ######

  y.hat.all   = GS_pred$map
  tau.hat.all = rep(NaN,n)
  tau.ci.all  = matrix(NaN,n,2)

  GS_treat = treatment(GP_fit, X = data.frame(X.train))
  tau.hat.all[-idx.test] = GS_treat$map
  tau.ci.all[-idx.test,] = GS_treat$ci

  GS_treat = treatment(GP_fit, X = data.frame(X.test))
  tau.hat.all[idx.test] = GS_treat$map
  tau.ci.all[idx.test,] = GS_treat$ci

  GP_DISCARD_R3 = BART_DISCARD_R2

  GS_treat = treatment(GP_fit, X = data.frame(X.train[GP_DISCARD_R3[-idx.test]==0,]))
  ate.ci.train= GS_treat$ate_ci
  GS_treat = treatment(GP_fit, X = data.frame(X.test[GP_DISCARD_R3[idx.test]==0,]))
  ate.ci.test = GS_treat$ate_ci

  result.item = update.results("GP w BART R2",result.item,idx.test,GP_DISCARD_R3,
                               y_hat_all = y.hat.all,
                               tau_hat_all = tau.hat.all,
                               tau_ci_all = tau.ci.all,
                               ate.ci.train = ate.ci.train,
                               ate.ci.test = ate.ci.test,
                               y=Y,z=Z, y1 = mysample$y1, y0 = mysample$y0 ,pscore = pscore)

  ## GP - IPW ######
  myweights = pscore.train*Z.train+(1-pscore.train)*(1-Z.train)
  GP_fit = CausalStump::CausalStump(y = Y.train,X = data.frame(X.train), z = Z.train,w = myweights,minimumRMSE=FALSE,
                                    maxiter=5000,tol=1e-2,learning_rate = 0.01,prior=FALSE,myoptim = "Nadam")

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

  DISCARD_NOTHING = rep(0,n)

  ## evaluate and store results
  result.item = update.results("GP weights",result.item,idx.test,DISCARD_NOTHING,
                               y_hat_all = y.hat.all,
                               tau_hat_all = tau.hat.all,
                               tau_ci_all = tau.ci.all,
                               ate.ci.train = ate.ci.train,
                               ate.ci.test = ate.ci.test,
                               y=Y,z=Z,y1 = mysample$y1, y0 = mysample$y0 ,pscore = pscore)

  ## GP - subgroup balancing ######
  myweights = mean(Z.train)*Z.train+(1-mean(Z.train))*(1-Z.train)
  GP_fit = CausalStump::CausalStump(y = Y.train,X = data.frame(X.train), z = Z.train,w = myweights,minimumRMSE=FALSE,
                                    maxiter=5000,tol=1e-2,learning_rate = 0.01,prior=FALSE,myoptim = "Nadam")

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

  DISCARD_NOTHING = rep(0,n)

  ## evaluate and store results
  result.item = update.results("GP balancing",result.item,idx.test,DISCARD_NOTHING,
                               y_hat_all = y.hat.all,
                               tau_hat_all = tau.hat.all,
                               tau_ci_all = tau.ci.all,
                               ate.ci.train = ate.ci.train,
                               ate.ci.test = ate.ci.test,
                               y=Y,z=Z,y1 = mysample$y1, y0 = mysample$y0 ,pscore = pscore)

  ## GP - variance 1 ######
  myweights = pscore.train*(1-pscore.train)
  GP_fit = CausalStump::CausalStump(y = Y.train,X = data.frame(X.train), z = Z.train,w = myweights,minimumRMSE=FALSE,
                                    maxiter=5000,tol=1e-2,learning_rate = 0.01,prior=FALSE,myoptim = "Nadam")

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

  DISCARD_NOTHING = rep(0,n)

  ## evaluate and store results
  result.item = update.results("GP variance 1 weights",result.item,idx.test,DISCARD_NOTHING,
                               y_hat_all = y.hat.all,
                               tau_hat_all = tau.hat.all,
                               tau_ci_all = tau.ci.all,
                               ate.ci.train = ate.ci.train,
                               ate.ci.test = ate.ci.test,
                               y=Y,z=Z,y1 = mysample$y1, y0 = mysample$y0 ,pscore = pscore)

  result.item
}


max_iter=100
nr_cores = detectCores()
cat(sprintf("\n-----\nIterations: %d on %d cores\n",max_iter,nr_cores))
cl <- makeCluster(nr_cores,type="FORK")
system.time(results <- parSapplyLB(cl,1:max_iter,run)) #clusterApplyLB, parSapply
stopCluster(cl)

save.image(file = paste("H11_",max_iter,"_weights.RData",sep=""))


