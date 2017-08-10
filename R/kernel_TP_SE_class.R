require(Rcpp)
require(RcppArmadillo)
#require(invgamma)
sourceCpp("src/kernel_TP_SE_cpp.cpp")
#Using some prior-invariant functions from the GP cpp file:
sourceCpp("src/kernel_GP_SE_cpp.cpp")
require(mvnfast)


KernelClass_TP_SE <- setRefClass("SqExpKernel_TP",
                                 fields = list(parameters = "list",
                                               invSmatn = "matrix", #prior-scaled kernel matrices "lambda_0 * Sigma"
                                               Smat = "matrix",
                                               Sm = "matrix",
                                               Sa = "matrix",
                                               w = "numeric",
                                               p = "numeric"),
                                 methods = list(
                                   parainit = function(y,nu) {

                                     #for momentum and adam
                                     parameters <<- list(sigma=log(var(y)),
                                                         sigma_z=log(runif(1,min=1,max=2)),
                                                         lambdam=log(runif(1,min=0.2,max=1)),
                                                         nu=log(nu), #exp(nu)
                                                         lambda0=log(runif(1,min=0.1,max=0.7)), #exp( lambda0)
                                                         Lm=rep(-0.1,p),La=rep(0.1,p),
                                                         mu = mean(y)
                                                         )
                                   },
                                   kernel_mat = function(X1,X2,z1,z2) {
                                     #intended use for the gardaient step and for prediction
                                     X1 = as.matrix(X1); #otherwise I need to rewrite the kernel function
                                     X2 = as.matrix(X2);

                                     #unscaled kernel matrices "Sigma
                                     Slist = kernmat_TP_SE_cpp(X1,X2,z1,z2,parameters)

                                   },
                                   getinv_kernel = function(X,z) {
                                     #get matrices and return inverse for prediction, maybe return Ka in a later stage

                                     #get new kernel and invert with noise
                                     Klist = kernel_mat(X,X,z,z);

                                     Smat <<- Klist[[1]]
                                     Sm   <<- Klist[[2]]
                                     Sa   <<- Klist[[3]]

                                     invSmatn <<- invkernel_cpp(z,w,Smat,parameters) #no error handling
                                     list(invSmatn,Smat) #,Ka)?
                                   },
                                   para_update = function(iter,y,X,z,Optim) {
                                     #update Kmat and invKmat in the class environment
                                     getinv_kernel(X,z);

                                     gradlist = grad_TP_SE_cpp(y,as.matrix(X),z,w,Smat,Sm,Sa,invSmatn,parameters)
                                     gradients = gradlist$gradients; stats = gradlist$stats
                                     #gradients$Lm = as.numeric(gradients$Lm)
                                     #gradients$La = as.numeric(gradients$La)
                                     parameters <<- Optim$update(iter,parameters,gradients)
                                     mean_solution(y,z)

                                     if(iter%%100 == 0){ cat(sprintf("%5d | log Evidence %9.4f | RMSE %9.4f \n", iter, stats[2], stats[1])) }
                                     get_train_stats(y,X,z)
                                     stats
                                   },
                                   get_train_stats = function(y,X,z){
                                     getinv_kernel(X,z);

                                     stats = stats_TP_SE(y,Smat, invSmatn, parameters)

                                   },
                                   mean_solution = function(y,z){
                                     #means using analytic solution
                                     mu = mu_solution_cpp(y, z, invSmatn, parameters)
                                     parameters$mu <<- mu #mu[[1]]; #parameters$mu_z <<- 0*mu[[2]]
                                     parameters <<- parameters
                                   },
                                   predict = function(y,X,z,X2,z2){ #NEED CHANGE

                                     S_xX = kernel_mat(X2,X,z2,z)$Smat
                                     S_xx = kernel_mat(X2,X2,z2,z2)$Smat

                                     #map
                                     map = S_xX %*% invSmatn %*% (y - parameters$mu) + parameters$mu
                                     cov = S_xx - S_xX %*% invSmatn %*% t(S_xX) + diag(exp(parameters$sigma + parameters$sigma_z * z) )
                                     #in current specification, includes already lamnbda0
                                     #want:
                                     #ybar = y - parameters$mu;
                                     #using the mode of the sampled credible interval
                                     #mode = c( exp(parameters$nu + parameters$lambda0)) / (c(exp(parameters$nu)) + 2)
                                     #sqrt_var = sqrt( diag(cov) / mode );

                                     #uncentered_ci = cbind(-1.96 * sqrt_var,1.96 * sqrt_var )
                                     N_sample = 10000

                                     TP_samples = mvnfast::rmvt(N_sample, mu = map,  sigma = cov, df = c(exp(parameters$n)) )
                                     #symmetric,
                                     TP_samples = t(apply(TP_samples, 1, function(x) abs(x-map) ))

                                     U = apply(TP_samples,2,sort)[floor(N_sample*0.95),] #floor has a faster convergence as samples denser

                                     list(map=map,ci=c(-U,U))
                                   },
                                   predict_treat = function(y,X,z,X2){ #NEED CHANGE
                                     n = length(y);
                                     n2 = nrow(X2);
                                     z2 = rep(1,n)

                                     S_xX = kernel_mat(X2,X,z2,z)$Sa
                                     S_xx = kernel_mat(X2,X2,z2,z2)$Sa

                                     #map
                                     map = S_xX %*% invSmatn %*% (y - parameters$mu)
                                     ate = mean(map);

                                     #ci
                                     #median_lambdaa = qinvgamma(0.5, shape = c(exp(parameters$nu))/2, rate = c(exp(parameters$nu))/2 )
                                     #using the mode of the sampled credible interval

                                     #mode = c( exp(parameters$nu)) / (c(exp(parameters$nu)) + 2)
                                     cov = ( S_xx - S_xX %*% invSmatn %*% t(S_xX) + exp(parameters$sigma_z) * diag(n2);
                                     #ate_cov = sum(sum(cov))/(n^2);
                                     N_sample = 10000

                                     TP_samples = mvnfast::rmvt(N_sample, mu = map,  sigma = cov, df = c(exp(parameters$n)) )
                                     #symmetric,
                                     TP_samples = t(apply(TP_samples, 1, function(x) abs(x-map) ))
                                     U = apply(TP_samples,2,sort)[floor(N_sample*0.95),] #floor has a faster convergence as samples denser

                                     ate_U = sort(apply(TP_samples,1,mean))[floor(1000*0.95)]

                                     #sqrt_var = sqrt( diag(cov) );
                                     #uncentered_ci = cbind(-1.96 * sqrt_var, 1.96 * sqrt_var)
                                     #ate_uncentered_ci = cbind(-1.96 * sqrt(ate_cov), 1.96 * sqrt(ate_cov) )
                                     list(map = map,ci = c(-U,U),ate_map = ate , ate_ci = c(-ate_U,ate_U))
                                   }
                                 )
)
