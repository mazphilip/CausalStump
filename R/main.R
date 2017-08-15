#' Fit a Gaussian process or Student-t process to the treatment and control surface.
#'
#' @param y A vector
#' @param X A data.frame
#' @param z A vector
#' @param w A vector (optional)
#' @param pscore A vector (optional)
#' @param kernelfun A string (default: SE)
#' @param myoptim A string (default: Gradient Descent -- GD)
#' @param maxiter (default: 5000)
#' @param tol (default: 1e-4)
#' @param prior A logic statement (default: FALSE)
#' @param nu A value (default: 200)
#' @param learning_rate (default: 0.01)
#' @param beta1 (default: 0.9)
#' @param beta2 (default: 0.999)
#' @param momentum (default: 0.0)
#' @return The function returns the fitted process as a CausalStump class object
#' @examples
#' #Generate data
#' n = 120
#' Z = rbinom(n, 1, 0.3)
# 'X1 = runif(sum(Z), min = 20, max = 40)
#' X0 = runif(n-sum(Z), min = 20, max = 40)
#' X = data.frame(matrix(NaN,n,1))
#' X[Z==1,] = X1; X[Z==0,] = X0
#' y0_true = as.matrix(72 + 3 * sqrt(X))
#' y1_true = as.matrix(90 + exp(0.06 * X))
#' Y0 = rnorm(n, mean = y0_true, sd = 1)
#' Y1 = rnorm(n, mean = y1_true, sd = 1)
#' Y = Y0*(1-Z) + Y1*Z
#' mystump <- CausalStump(Y,X,Z)

CausalStump <- function(y,X,z,w,pscore,kernelfun="SE",myoptim = "Nadam",maxiter=5000,tol=1e-4,prior=FALSE,nu=200,learning_rate=0.01,beta1=0.9,beta2=0.999,momentum=0.0){
  #check dimensionality and class of X
  check_inputs(y,X,z);
  if(!missing(pscore)){ X = cbind(X, pscore); }

  n = length(y); p = ncol(X);

  #normalize variables
  norm_ret = norm_variables(y,X)
  moments = norm_ret$moments; y = norm_ret$y; X = norm_ret$X;

  #if not alternative weighting, use equal weights
  if(missing(w)){ w = rep(1,n) }

  #set GP or TP
  if(prior==TRUE){
    cat(sprintf("\nFitting the Student-t Process with nu=%d:\n",nu))
    myKernel = KernelClass_TP_SE$new(p = p,w = rep(1,n)) #object
    myKernel$parainit(y,nu);
  } else {
    cat("\nFitting the Gaussian process:\n")
    myKernel = KernelClass_GP_SE$new(p = p,w = rep(1,n)) #object
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
    stats[,iter+1] = myKernel$para_update(iter,y,X,z,myOptimizer)

    #tolerance missing
    change = abs(stats[2,iter+1] - stats[2,iter])
    if((change < tol) && (iter > 3)){ cat( sprintf("Stopped: change smaller than tolerance after %d iterations",iter)); break; }
  }
  if(iter == maxiter){ cat("Stopped: maximum iterations reached") }

  #generate kernel once more and final stats
  stats[,iter+2] = myKernel$get_train_stats(y,X,z)

  #plot
  graphics::par(mfrow=c(1,2))
  graphics::plot(stats[2,3:(iter+2)],type="l",ylab="log Evidence",xlab="Iteration")
  graphics::plot(stats[1,3:(iter+2)],type="l",ylab="training RMSE",xlab="Iteration")

  mystump = structure(list(Kernel = myKernel,moments=moments,train_data=list(y=y,X=X,z=z)), class = "CausalStump")
}
