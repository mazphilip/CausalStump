

niters=1

##### training
library(CausalStump)
library(BayesTree)
library(gbm)
set.seed(2659232)

## Preprocessing Simulation data ####
#From Hill 2011

load(file="test/sim.data")
obs <- imp1[!(imp1$treat==1 & imp1$momwhite==0),]

covs.cont.n=c("bw","b.head","preterm","birth.o","nnhealth","momage")
covs.cat.n=c("sex","twin","b.marr","mom.lths","mom.hs",	"mom.scoll","cig","first","booze","drugs","work.dur","prenatal","ark","ein","har","mia","pen","tex","was")
p=length(c(covs.cont.n,covs.cat.n))

### calculate pscores and weights for tot
Trt=obs$treat
form.qx=as.formula(obs[,c("treat",covs.cont.n,covs.cat.n)])
qx=glm(data=obs[,c("treat",covs.cont.n,covs.cat.n)],formula=form.qx,family=binomial)$fitted
wts=rep(1,nrow(obs))
wts[Trt==0]=qx[Trt==0]/(1-qx[Trt==0])

#### get data in right format for BART
xt=obs[,c(covs.cont.n,covs.cat.n,"treat")]
xt=as.matrix(xt)
xp1=xt[xt[,"treat"]==1,]
xp2=xp1
xp1[,ncol(xt)]=1
xp2[,ncol(xt)]=0
xp=rbind(xp1,xp2)

nt=sum(obs$treat)

##################### now simulate outcome data
##### covariate data, X
covs.ols = c(covs.cont.n,covs.cat.n)
X = obs[,covs.ols]
#X = na.omit(X)
# now standardize the continuous variables
X[,covs.cont.n]=as.data.frame(t((t(X[,covs.cont.n])-unlist(lapply(X[,covs.cont.n],mean)))/sqrt(unlist(lapply(X[covs.cont.n],var)))))
# record numbers of units and covariates
#N = nrow(X)
#dimx = ncol(X)
#Xmat = as.matrix(X)


Hill11_surfaceB <- function(X,z){
  N = nrow(X)
  dimx = ncol(X)
  Xmat = as.matrix(X)

  sigy = 1
  betaB = c(sample(c(.0,.1,.2,.3,.4),(dimx+1),replace=TRUE,prob=c(.6,.1,.1,.1,.1)))
  yb0hat = exp((cbind(rep(1, N), (Xmat+.5)) %*% betaB))
  yb1hat = cbind(rep(1, N), (Xmat+.5)) %*% betaB
  offset = c(mean(yb1hat[z==1] - yb0hat[z==1])) - 4
  yb1hat = cbind(rep(1, N), (Xmat+.5)) %*% betaB -offset
  YB0 = rnorm(N, yb0hat, sigy)
  YB1 = rnorm(N, yb1hat, sigy)

  YB = YB1; YB[z==0] = YB0[z==0]

  list(y = YB, X = Xmat ,z = z, tau = yb1hat - yb0hat, y1 = yb1hat, y0 = yb0hat )

}

### now create matrix of all interactions etc for third response surface
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

update.results <- function(iter,index,result.item,y.hat,tau.hat,tau.ci,ate.ci,y,z,y1,y0,pscore){
  threshold = 0.1

  result.item["PEHE.all",index] = sqrt(mean((y1 - y0 - tau.hat)^2))
  result.item["PEHE.cs",index]  = sqrt(mean((y1 - y0 - tau.hat)[pscore>threshold & pscore<(1-threshold)]^2))

  result.item["E.ate",index]  = (mean(y1-y0) - mean(tau.hat))
  result.item["E.att",index]  = (mean((y1-y0)[z==1]) - mean((tau.hat)[z==1]))

  result.item["MSE.ate",index]  = (mean(y1-y0) - mean(tau.hat) )^2
  result.item["MSE.att",index]  = (mean((y1-y0)[z==1]) - mean((tau.hat)[z==1]))^2

  result.item["RMSE.all",index] = sqrt(mean((y - y.hat)^2))
  result.item["RMSE.cs",index]  = sqrt(mean((y - y.hat)[pscore>threshold & pscore<(1-threshold)]^2))

  result.item["coverage.ite",index] = mean(((tau.ci[,1] < (y1-y0)) * (tau.ci[,2] > (y1-y0))))
  result.item["coverage.ate",index] = ((ate.ci[1]< mean(y1-y0)) * (ate.ci[2] >mean(y1-y0)))

  result.item["range.ate",index] = (ate.ci[2] - ate.ci[1])

  result.item
}

run <- function(i){
  #data generation within loop
  if(i<=500){set.seed(565 + i*5)}
  if(i>500){set.seed(7565 + i*5)}

  mydimnames = list(c("PEHE.all","PEHE.cs","RMSE.all","RMSE.cs","E.ate","E.att","MSE.ate","MSE.att","coverage.ite","coverage.ate","range.ate"),c("BART","CF","GP","TP2","TP100","BCF","BART.PS","CF.PS","GP.PS","TP2.PS","TP100.PS"))
  result.item = matrix(NaN, length(mydimnames[[1]]),length(mydimnames[[2]]),dimnames=list(mydimnames[[1]],mydimnames[[2]]) )

  ## sample data ####
  mysample = Hill11_surfaceB(X,Trt)
  Y = mysample$y; X = mysample$X; Z = mysample$z
  Xdf = data.frame(X)
  n = nrow(X)
  X.bart.train = cbind(X,Z);
  X.bart.test = rbind(cbind(X,rep(1,n)),cbind(X,rep(0,n)));


  ## Propensity score estimation ####
  #Gradient Boosting Machines
  myPSmodel = gbm(Z ~ .,data = Xdf,distribution="bernoulli",
                  n.trees = 2000,
                  shrinkage = .1,
                  cv.folds = 20,
                  n.cores = 1)
  best.fit = gbm.perf(myPSmodel,method="cv")
  pscore.gbm = predict(myPSmodel,n.trees=best.fit,type="response")

  #BART
  myPSbart = BayesTree::bart(x.train = X,y.train = Z,x.test=X,binaryOffset=mean(Z),ntree=200)
  pscore.bart = pnorm(apply(myPSbart$yhat.test,2,mean))

  pscore = pscore.bart

  ## BART ####
  BART_fit = BayesTree::bart(x.train = X.bart.train, y.train = Y, x.test=X.bart.test,ntree=100,ndpost=3000,nskip=500,keepevery=1)

  ## obtain estimates
  y.hat.bart = apply( BART_fit$yhat.train,2,mean)
  tau.hat.bart = apply( (BART_fit$yhat.test[,1:n] - BART_fit$yhat.test[,(n+1):(2*n)]) ,2,mean)
  ate.ci.bart = {
    tmpvar = var(apply( BART_fit$yhat.train,1,mean))
    L = mean(tau.hat.bart) - 1.96 * sqrt(tmpvar)
    U = mean(tau.hat.bart) + 1.96 * sqrt(tmpvar)
    cbind(L,U)
  }
  tau.ci.bart = {
    tmpsort = apply(BART_fit$yhat.test[,1:n]-BART_fit$yhat.test[,(n+1):(2*n)],2,sort)
    idx975 = round(nrow(BART_fit$yhat.test)*0.975,0)
    idx025 = round(nrow(BART_fit$yhat.test)*0.025,0)
    L = tmpsort[idx025,]
    U = tmpsort[idx975,]
    cbind(L,U)
  }

  ## evaluate and store results
  result.item = update.results(iter,"BART",result.item,
                               y.hat = y.hat.bart,
                               tau.hat = tau.hat.bart,
                               tau.ci = tau.ci.bart,
                               ate.ci = ate.ci.bart,
                               y = Y, z = Z, y1 = mysample$y1, y0 = mysample$y0 ,pscore = pscore)

  ## GP PS ####

  GP_PS_fit = CausalStump(Y,Xdf,Z,pscore=pscore,maxiter=1000,tol=1e-4,learning_rate = 0.1,prior=FALSE,myoptim = "Nadam")

  ## predict
  #GS_PS_pred = predict_surface(Xdf,Z,GP_PS_fit,pscore=pscore)
  GS_PS_treat = predict_treatment(Xdf,GP_PS_fit,pscore=pscore)
  plot(GS_PS_treat$map,mysample$tau); abline(0,1,col=3)
  PEHE = sqrt( mean( (GS_PS_treat$map - mysample$tau)^2 ) ); PEHE

  GS_PS_pred = predict_surface(X.CS,Z,GP_PS_fit)

  y.hat.gp.ps = GS_PS_pred$map
  tau.hat.gp.ps = GS_PS_treat$map
  tau.ci.gp.ps  = GS_PS_treat$ci
  ate.ci.gp.ps  = GS_PS_treat$ate_ci

  ## evaluate and store results
  result.item = update.results(iter,"GP.PS",result.item,
                               y.hat = y.hat.gps,
                               tau.hat = tau.hat.gps,
                               tau.ci = tau.ci.gps,
                               ate.ci = ate.ci.gps,
                               y = Y, z = Z, y1 = mysample$y1, y0 = mysample$y0 ,pscore = pscore)

}
