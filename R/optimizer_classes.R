require(Rcpp)
require(RcppArmadillo)
sourceCpp("optimizer_cpp.cpp")

list2zero <- function(mylist){
  tmplist = within(stack(mylist), values <- 0.0)
  
  if(nrow(tmplist) == nlevels(tmplist[,2])){
    zerolist <- as.list(rep(0.0,nrow(tmplist)))
    zerolist = setNames(zerolist,tmplist[,"ind"])
  } else {
    zerolist = unstack(tmplist)  
  }
  zerolist
}

## define optimizer classes
optAdam <- setRefClass("AdamOpt",
                       fields = list(m = "list",
                                     v = "list",
                                     lr = "numeric",
                                     beta1 = "numeric",
                                     beta2 = "numeric"),
                       methods = list(
                         update = function(iter,parameters,gradients) { 
                           
                           mylist = Adam_cpp(iter,lr,beta1,beta2,1e-8,m,v,gradients,parameters)
                           
                           m <<- mylist[[1]]
                           v <<- mylist[[2]]
                           mylist[[3]]
                         },
                         initOpt = function(KernelObj){
                           
                           m <<- list2zero(KernelObj$parameters)
                           v <<- list2zero(KernelObj$parameters)
                         })
)

optNadam <- setRefClass("NadamOpt",
                        fields = list(m = "list",
                                      v = "list",
                                      lr = "numeric",
                                      beta1 = "numeric",
                                      beta2 = "numeric"),
                        methods = list(
                          update = function(iter,parameters,gradients) { 
                            mylist = Nadam_cpp(iter,lr,beta1,beta2,1e-8,m,v,gradients,parameters)
                            m <<- mylist[[1]];
                            v <<- mylist[[2]]
                            #print("Nadam: m")
                            #print(m$Lm)
                            #print(m$La)
                            #print("Nadam: v")                            
                            #print(v$Lm)
                            #print(v$La)
                            mylist[[3]]
                          },
                          initOpt = function(KernelObj){
                            m <<- list2zero(KernelObj$parameters)
                            v <<- list2zero(KernelObj$parameters)
                          })
)

optNestorov <- setRefClass("NesterovOpt",
                           fields = list(nu = "list",
                                         lr = "numeric",
                                         momentum = "numeric"),
                           methods = list(
                             update = function(iter,parameters,gradients) { 
                               mylist = Nesterov_cpp(lr,momentum,nu,gradients,parameters) 
                               nu <<- mylist[[1]]
                               mylist[[2]]
                             },
                             initOpt = function(KernelObj){
                               nu <<- list2zero(KernelObj$parameters)
                             })
)