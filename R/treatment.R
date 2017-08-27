treatment <- function(cs_object, ...)  UseMethod("treatment")

treatment.CausalStump <- function(cs_object,X,pscore){
  #this function returns the prediction for the fitted Gaussian process
  #only includes calculations with the new data

  if(missing(X)){
    X = cs_object$train_data$X
  }  else {
    if(!missing(pscore)){ X = cbind(X, pscore); }
    if(ncol(X)!=ncol(cs_object$train_data$X)){ stop("Error: propensity score mismatch", call. = FALSE) }
    #normalize the non-binary variables
    X = norm_variables(X = X,moments = cs_object$moments)$X
  }
  n = nrow(X);

  #remaining kernel calculations using the kernel class method
  pred_list = cs_object$Kernel$predict_treat(cs_object$train_data$y,cs_object$train_data$X,cs_object$train_data$z,X)

  #add the moments
  map = sqrt(cs_object$moments$varY) * pred_list$map
  ci = cs_object$moments$varY * pred_list$ci+ cbind(map,map) #change to broadcasting
  cate = sqrt(cs_object$moments$varY) *pred_list$ate_map;
  cate_ci = cs_object$moments$varY *pred_list$ate_ci + cbind(cate,cate)
  #return arguments as list
  list(map = map,ci = ci, ate = cate , ate_ci = cate_ci )
}
