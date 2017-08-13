#require(RcppArmadillo)
#source("R/utilities.R")
#source("R/kernel_classes.R")
#source("R/optimizer_classes.R")

CausalStump <- function(y,X,z,w,pscore,kernelfun="SE",myoptim = "Nadam",prior=FALSE,nu=200,maxiter=5000,tol=1e-4,learning_rate=0.01,beta1=0.9,beta2=0.999,momentum=0.0){
  #check dimensionality and class of X
  check_inputs(y,X,z);
  if(!missing(pscore)){ X = cbind(X, pscore); }
  n = length(y);
  p = ncol(X);
  #normalize variables

  norm_ret = norm_variables(y,X)
  moments = norm_ret$moments; y = norm_ret$y; X = norm_ret$X;

  #if not alternative weighting, use equal weights
  if(missing(w)){ w = rep(1,n) }
  if(prior==TRUE){
    #if(kernelfun == "SE") {
    print("Fitting the Student-t Process")
      myKernel = KernelClass_TP_SE$new(p = p,w = rep(1,n)) #object
    #} else if(kernelfun == "Matern52") {
    #  myKernel = KernelClass_TP_M52$new(p = p,w = rep(1,n)) #object
    #} ## extend if necessary
      myKernel$parainit(y,nu);
  } else {
    #if(kernelfun == "SE") {
    print("Fitting the Gaussian process")
      myKernel = KernelClass_GP_SE$new(p = p,w = rep(1,n)) #object
    #} else if(kernelfun == "Matern52") {
    #  myKernel = KernelClass_GP_M52$new(p = p,w = rep(1,n)) #object
    #} ## extend if necessary
      myKernel$parainit(y);
  }

  #initialize optimizer
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
    #ubuntu initializes two instances for update!
    stats[,iter+1] = myKernel$para_update(iter,y,X,z,myOptimizer)

    #tolerance missing
    change = abs(stats[2,iter+1] - stats[2,iter])
    if((change < tol) && (iter > 3)){ print("Stopped due to tolerance"); break; }
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
  plot(stats[1,3:(iter+2)],type="l",ylab="RMSE",xlab="Iterations",ylim=c(0,5))

  Stump = list(Kernel = myKernel,moments=moments,train_data=list(y=y,X=X,z=z)) #generate S3 output class? and use overloaded predict etc. ?
  class(Stump) = "CausalStump"
  Stump
}

predict_surface <- function(X,z,CSobject,pscore){
  #this function returns the prediction for the fitted Gaussian process
  #only includes calculations with the new data
  if(!missing(pscore)){ X = cbind(X, pscore); }

  #if(class(CSobject)!=""){ warning("CSobject: incorrect class", call. = FALSE) }
  n = nrow(X);
  #this makes it easier to plot both surfaces
  if(length(z)==1){ z = rep(z,n) }

  #check_inputs(y,X,z);

  #normalize the non-binary variables
  X = norm_variables(X = X,moments = CSobject$moments)$X

  #remaining kernel calculations using the kernel class method
  pred_list = CSobject$Kernel$predict(CSobject$train_data$y,CSobject$train_data$X,CSobject$train_data$z,X,z)

  #add the moments
  map = CSobject$moments$meanY + sqrt(CSobject$moments$varY) * pred_list$map
  ci = CSobject$moments$varY * pred_list$ci+ cbind(map,map) #change to broadcasting

  list(map = map,ci = ci)
}

predict_treatment <- function(X,CSobject,pscore){
  #this function returns the prediction for the fitted Gaussian process
  #only includes calculations with the new data
  if(!missing(pscore)){ X = cbind(X, pscore); }

  #if(class(CSobject)!=""){ warning("CSobject: incorrect class", call. = FALSE) }
  n = nrow(X);

  #check_inputs(X,z);

  #normalize the non-binary variables
  X = norm_variables(X = X,moments = CSobject$moments)$X

  #remaining kernel calculations using the kernel class method
  pred_list = CSobject$Kernel$predict_treat(CSobject$train_data$y,CSobject$train_data$X,CSobject$train_data$z,X)


  #add the moments
  map = sqrt(CSobject$moments$varY) * pred_list$map
  ci = CSobject$moments$varY * pred_list$ci+ cbind(map,map) #change to broadcasting

  cate = sqrt(CSobject$moments$varY) *pred_list$ate_map;
  cate_ci = CSobject$moments$varY *pred_list$ate_ci + cbind(cate,cate)
  list(map = map,ci = ci, ate = cate , ate_ci = cate_ci )
}
