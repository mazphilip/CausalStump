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
    within_unit_circle = unique(1:p * apply(Xp,2,function(x) { all(na.omit(x)<=1 & na.omit(x)>=0 ) }) )
    within_unit_circle = within_unit_circle[within_unit_circle!=0]

    normalize_varbs = setdiff(nonbinary_varbs,within_unit_circle)

    tmp = data.frame(X[,normalize_varbs])

    moments$meanX[normalize_varbs] = apply(tmp,2,mean)
    moments$varX[normalize_varbs]  = apply(abs(tmp),2,max)
  }
  mynorm <- function(i){ (X[,i]-moments$meanX[i]) / sqrt(moments$varX[i]) }
  X = data.frame(sapply(1:p,mynorm))
  if(missing(y)){
    out = list(X=X)
  } else {
    y = (y - moments$meanY ) / sqrt(moments$varY)
    out = list(moments=moments, y=y, X=X)
  }
  out
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
