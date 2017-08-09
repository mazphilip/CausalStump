require(Rcpp)
require(RcppArmadillo)
sourceCpp("src/kernel_SE_cpp.cpp")

KernelClass_SE <- setRefClass("SqExpKernel",
  fields = list(parameters = "list",
                invKmatn = "matrix",
                Kmat = "matrix",
                Km = "matrix",
                Ka = "matrix",
                w = "numeric",
                p = "numeric"),
  methods = list(
    parainit = function(y) {

      #for momentum and adam
      parameters <<- list(sigma=log(var(y)),
                        sigma_z=log(runif(1,min=1,max=2)),
                        lambdam=log(runif(1,min=0.2,max=1)),
                        lambdaa=log(runif(1,min=0.2,max=1)),
                        Lm=rep(-0.1,p),La=rep(0.1,p),
                        mu = mean(y)
                        )
                        #mu_z = mean(y[z==1]))

    },
    kernel_mat = function(X1,X2,z1,z2) {
      #intended use for the gardaient step and for prediction
      X1 = as.matrix(X1); #otherwise I need to rewrite the kernel function
      X2 = as.matrix(X2);

      Klist = kernmat_SE_cpp(X1,X2,z1,z2,parameters)

    },
    getinv_kernel = function(X,z) {
      #get matrices and return inverse for prediction, maybe return Ka in a later stage

      #get new kernel and invert with noise
      Klist = kernel_mat(X,X,z,z);

      Kmat <<- Klist[[1]]
      Km   <<- Klist[[2]]
      Ka   <<- Klist[[3]]

      invKmatn <<- invkernel_cpp(z,w,Kmat,parameters) #no error handling
      list(invKmatn,Kmat) #,Ka)?
    },
    para_update = function(iter,y,X,z,Optim) {
      #update Kmat and invKmat in the class environment
      getinv_kernel(X,z);

      gradlist = grad_SE_cpp(y,as.matrix(X),z,w,Kmat,Km,Ka,invKmatn,parameters)
      gradients = gradlist$gradients; stats = gradlist$stats
      #gradients$Lm = as.numeric(gradients$Lm)
      #gradients$La = as.numeric(gradients$La)
      parameters <<- Optim$update(iter,parameters,gradients)
      mean_solution(y,z)

      if(iter%%100 == 0){ cat(sprintf("%5d | log Evidence %9.4f | RMSE %9..4f \n", iter, stats[2], stats[1])) }
      get_train_stats(y,X,z)
      stats
    },
    get_train_stats = function(y,X,z){
      getinv_kernel(X,z);

      #gradlist = grad_SE_cpp(y,as.matrix(X),z,w,Kmat,Km,Ka,invKmatn,parameters)
      stats = stats_SE(y,Kmat, invKmatn, parameters)

    },
    mean_solution = function(y,z){
      #means using analytic solution
      mu = mu_solution_cpp(y, z,invKmatn,parameters)
      parameters$mu <<- mu #mu[[1]]; #parameters$mu_z <<- 0*mu[[2]]
      parameters <<- parameters
    }
    )
)

