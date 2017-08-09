norm_variables <- function(y,X,moments){
  p = ncol(X)

  #use function for both in and out-of-sample normalization
  if(missing(moments)){
    moments = list(meanX = rep(0,p),
                   varX  = rep(1,p),
                   meanY  = 0*mean(y), # we do not center the outcome surface
                   varY  = var(y))

    #only normalize non-binary variables
    nonbinary_varbs = unique(1:p * (apply(X,2,function(x) { all(na.omit(x) %in% 0:1) })==0))
    nonbinary_varbs = nonbinary_varbs[nonbinary_varbs!=0]

    tmp = data.frame(X[,nonbinary_varbs])

    moments$meanX[nonbinary_varbs] = apply(tmp,2,mean)
    moments$varX[nonbinary_varbs]  = apply(tmp,2,var)
  }

  y = (y - moments$meanY )/sqrt(moments$varY)
  mynorm <- function(i){ (X[,i]-moments$meanX[i]) / sqrt(moments$varX[i]) }

  X = data.frame(sapply(1:p,mynorm))

  list(moments = moments, y=y, X=X)
}

check_inputs <- function(y,X,z){
  if(length(y) != nrow(X) ){
    warning("y and X: nr of observations mismatch", call. = FALSE)
  }
  if(length(y) != length(z) ){
    warning("y and z: nr of observations mismatch", call. = FALSE)
  }
  if(class(X) != "data.frame"){
    warning("X is not of class data.frame", call. = FALSE)
  }
}
