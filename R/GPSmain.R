
source("R/utilities.R")
source("R/kernel_classes.R")
source("R/optimizer_classes.R")

CausalStump <- function(y,X,z,w,pscore,kernelfun="SE",myoptim = "Nadam",maxiter=5000,tol=1e-4,learning_rate=0.01,beta1=0.2,beta2=0.999,momentum=0.0){
  #check dimensionality and class of X
  check_inputs(y,X,z);
  if(!missing(pscore)){ X = cbind(X, pscore); }

  p = ncol(X);
  #normalize variables
  norm_ret = norm_variables(y,X)
  moments = norm_ret$moments; y = norm_ret$y; X = norm_ret$X;

  #if not alternative weighting, use equal weights
  if(missing(w)){ w = rep(1,n) }

  if(kernelfun == "SE") {
    myKernel = KernelClass_SE$new(p = p,w = rep(1,n)) #object
  } else if(kernelfun == "Matern52") {
    myKernel = KernelClass_M52$new(p = p,w = rep(1,n)) #object
  } ## extend if necessary

  #initialize parameters
  myKernel$parainit(y);

  #select parameter
  if((myoptim=="Adam") || (myoptim=="Nadam")){
    if(myoptim=="Adam") {
      myOptimizer = optAdam$new(lr = learning_rate, beta1 = beta1, beta2 = beta2)
    } else {
      myOptimizer = optNadam$new(lr = learning_rate, beta1 = beta1, beta2 = beta2)
    }
  } else if(myoptim=="GD" || myoptim=="Nesterov"){
    if(myoptim=="GD"){ momentum=0.0 }
    myOptimizer = optNestorov$new(lr = learning_rate, momentum = momentum)
  }
  #set optimization variables
  myOptimizer$initOpt(myKernel);

  #store statistics for plotting
  stats = matrix(0,2,maxiter+2) #Evidence and RMSE

  for(iter in 1:maxiter){
    #print(myKernel$parameters)

    stats[,iter+1] = myKernel$para_update(iter,y,X,z,myOptimizer)

    #tolerance missing
    change = abs(stats[2,iter+1] - stats[2,iter])
    if((change < tol)){ print("Stopped due to tolerance"); break; }
  }
  if(iter == maxiter){ print("Stopped due to maximum iterations reached") }
  print("Final parameters")
  print(myKernel$parameters)
  #output changes

  #generate kernel once more and final stats
  stats[,iter+2] = myKernel$get_train_stats(y,X,z)



  #plot
  par(mfrow=c(1,2))
  plot(stats[2,3:(iter+2)],type="l",ylab="log Evidence",xlab="Iterations")
  plot(stats[1,3:(iter+2)],type="l",ylab="RMSE",xlab="Iterations")

  list(class="CSobject",Kernel = myKernel,moments=moments,train_data=list(y=y,X=X,z=z)) #generate S3 output class?
}

CS_fit = CausalStump(y,X2,z,maxiter=5000,learning_rate = 0.01,myoptim = "GD")





predict_surface <- function(y,X,z,CSobject){
  #this function returns the prediction for the fitted Gaussian process
  #only includes calculations with the new data

  #if(class(CSobject)!=""){ warning("CSobject: incorrect class", call. = FALSE) }
  n = length(y);
  #this makes it easier to plot both surfaces
  if(length(z)==1){ z = rep(z,n) }

  check_inputs(y,X,z);



  #normalize the non-binary variables
  norm_ret = norm_variables(y,X,CSobject$moments)
  y = norm_ret$y; X = norm_ret$X;

  #remaining kernel calculations
  K_xX = (CSobject$Kernel$kernel_mat(X,CSobject$train_data$X,z,CSobject$train_data$z));
  K_xX = K_xX$Kmat
  K_xx = (CSobject$Kernel$kernel_mat(X,X,z,z));
  K_xx = K_xx$Kmat

  #map
  map = moments$meanY + sqrt(moments$varY) * K_xX %*% CSobject$Kernel$invKmatn %*% (CSobject$train_data$y - CSobject$Kernel$parameters$mu)

  cov = sqrt(moments$varY) *(K_xx - K_xX %*% CSobject$Kernel$invKmatn %*% t(K_xX) )
  predvar = diag(cov);

  list(map = map,ci = cbind(map-1.96*sqrt(predvar),map+1.96*sqrt(predvar)))
}

#CS_pred0 = predict_surface(y,X2,0,CS_fit)
#CS_pred1 = predict_surface(y,X2,1,CS_fit)
#par(mfrow=c(1,1))
#plot(X[,1],CS_pred0$map,ylim=c(40,140))
#points(X[,1],CS_pred1$map)



predict_treatment <- function(){


}


#gradient only one dimensional!  ----- solved for Nadam and GD
#lambdam and lambdaanot converging, kernel?




X2 = cbind(X,rnorm(120))
X2 = data.frame(X2)

kernmat_SE_cpp(X2, X2, z, z, myKernel$parameters)


#sigma and sigma_z:works
#lambdam: rmse goes down but evid as well
#lambdaa: rmse goes down but evid as well
#Lm: evid up then down, RMSE up
#La: evid up, RMSE down then up
#mu: works

#kernel seems fine



# kernel not the problem
n2 = 3
myKernel = KernelClass_SE$new(p = 2,w = rep(1,n2))
myKernel$parainit(y)
myKernel$parameters$sigma=0
myKernel$parameters$sigma_z=0
myKernel$parameters$lambdam=log(2)
myKernel$parameters$lambdaa=log(2)
myKernel$parameters
X2 = data.frame(c(1,2,1.1))
