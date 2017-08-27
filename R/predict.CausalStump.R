#' Predict the response surfaces using a Gaussian process or Student-t process fit.
#'
#' @param cs_object A object from the CausalStump function
#' @param X A data.frame with new data. If not presented, using the training data.
#' @param z A vector or scalar with new treatment data. If it is a scalar, it predicts using the same value for all observations. If missing, it uses the training data.
#' @param pscore A vector with the propensity score. Throws an error if use is different to the CausalStump fit.
#' @return Returns the MAP and the 95 percent credible interval of the fitted process as a list.

predict.CausalStump <- function(cs_object,X,z,pscore){
  #this function returns the prediction for the fitted Gaussian process
  if(missing(X)){
    X = cs_object$train_data$X
    }  else {
    if(!missing(pscore)){ X = cbind(X, pscore); }
    if(ncol(X)!=ncol(cs_object$train_data$X)){ stop("Error: propensity score mismatch", call. = FALSE) }
    #normalize the non-binary variables
    X = norm_variables(X = X,moments = cs_object$moments)$X
    }
  n = nrow(X);
  if(missing(z)){ z = cs_object$train_data$z}
  else { if(length(z)==1){ z = rep(z,n) } }

  #remaining kernel calculations using the kernel class method
  pred_list = cs_object$Kernel$predict(cs_object$train_data$y,cs_object$train_data$X,cs_object$train_data$z,X,z)

  #get the output arguments
  map = cs_object$moments$meanY + sqrt(cs_object$moments$varY) * pred_list$map
  ci = cs_object$moments$varY * pred_list$ci+ cbind(map,map) #change to broadcasting

  list(map = map,ci = ci)
}
