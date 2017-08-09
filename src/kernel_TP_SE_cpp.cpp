// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <boost/math/special_functions/digamma.hpp>
#include <cmath>

using namespace arma;
using namespace Rcpp;

//same as for GP_SE
double evid_grad(arma::mat mymat, arma::mat dK) {
  //grad = - 0.5 * ( sum(diag( tmpK %*% dK ) ) )
  double grad = - 0.5 * trace( mymat * dK);
  return(grad);
}

//same as for GP_SE
Rcpp::List evid_scale_gradients(arma::mat X, arma::mat tmpK,arma::mat Km,arma::mat Ka, Rcpp::List parameters){
  // produce a (p x 2) matrix "grad" with the gradients of Lm and La in the first and second column, respectively
  unsigned int n = X.n_rows;
  unsigned int p = X.n_cols;

  //predefine matrices and temporary storage files
  arma::vec grad_m(p);
  arma::vec grad_a(p);
  arma::mat tmpX(n,n);
  arma::mat tmpM(n,n);
  arma::mat tmpA(n,n);

  for(unsigned int i=0; i<p; i++){
    //tmpX = abs(matrix(rep(X[,i],n),n,n) - t(matrix(rep(X[,i],n),n,n)))
    for(unsigned int r=0; r< n; r++){
      tmpX.col(r) = pow(X( r , i ) - X.col(i),2);
    }

    //Klist$Km * exp( 2*log(tmpX) - parameters$Lm[i] )
    tmpM = Km % tmpX * exp(- as<arma::vec>(parameters["Lm"])[i]);
    tmpA = Ka % tmpX * exp(- as<arma::vec>(parameters["La"])[i]);

    grad_m(i) = evid_grad(tmpK,tmpM);
    grad_a(i) = evid_grad(tmpK,tmpA);
  }

  return Rcpp::List::create(Named("m") = grad_m,Named("a") = grad_a );
}

// [[Rcpp::export]]
Rcpp::List kernmat_TP_SE_cpp(Rcpp::NumericMatrix X1,Rcpp::NumericMatrix X2,Rcpp::NumericVector Z1,Rcpp::NumericVector Z2, Rcpp::List para) {
  //input X1,X2,Z1,Z2,parameters

  //for(i in 1:p){
  //  tmpX = (matrix(rep(X1[,i],n2),n1,n2) - t(matrix(rep(X2[,i],n1),n2,n1)))^2
  //  Km        = Km        + tmpX        * exp(-parameters$Lm[i])
  //  Ka[Zmask] = Ka[Zmask] + tmpX[Zmask] * exp(-parameters$La[i])
  //}

  unsigned int n1 = X1.nrow();
  unsigned int n2 = X2.nrow();
  unsigned int p = X2.ncol();

  arma::mat tmpKm(n1,n2); tmpKm.zeros();
  arma::mat tmpKa(n1,n2); tmpKa.zeros();

  arma::rowvec tmprow(n2);

  NumericVector Lm = as<NumericVector>(para["Lm"]);
  NumericVector La = as<NumericVector>(para["La"]);

  for(unsigned int i=0;i < p;i++){
    //for every output row
    for(unsigned int r=0;r < n1;r++){
      tmprow = pow(X1( r , i ) - X2( _, i),2);

      tmpKm.row(r) += tmprow * exp(- Lm(i) );
      tmpKa.row(r) += tmprow * exp(- La(i) );

    }
  }

  double lambda_m = as<double>(para["lambdam"]);
  //double lambda_a = as<double>(para["lambdaa"]);

  Rcpp::NumericMatrix tmpK(n1, n2);

  for(unsigned int r = 0; r < n1; r++){
    for(unsigned int c = 0; c < n2; c++){
      tmpKm(r,c) = exp(lambda_m - tmpKm(r,c));
      tmpKa(r,c) = exp( - tmpKa(r,c)) * Z1(r) * Z2(c); //"lambdaa" = 0
      tmpK(r,c) = tmpKm(r,c) + tmpKa(r,c) ;
    }
  }

  //Put matrices in a list
  Rcpp::List out; out["Smat"] = tmpK; out["Sm"] = tmpKm; out["Sa"] = tmpKa;
  return(out);
  //return 0;
}

// [[Rcpp::export]]
Rcpp::List grad_TP_SE_cpp(arma::colvec y, arma::mat X, arma::colvec z,arma::colvec w, arma::mat Smat, arma::mat Sm, arma::mat Sa, arma::mat invSmatn, Rcpp::List parameters) {
  Rcpp::List gradients = clone(parameters); //for the same list structure

  unsigned int n = X.n_rows;
  //unsigned int p = X.n_cols;

  //preallocate memory
  arma::colvec stats(2);
  arma::colvec ybar(n);
  arma::colvec alpha(n);
  arma::mat tmpS(n,n);
  arma::mat dS(n,n);

  /*
  parameters <<- list(sigma=log(var(y)),
                      sigma_z=log(runif(1,min=1,max=2)),
                      lambdam=log(runif(1,min=0.2,max=1)),
                      nu=log(runif(1,min=0.2,max=1)), #exp(nu)
                        lambda0=log(runif(1,min=0.2,max=1)), #exp( lambda0)
                        Lm=rep(-0.1,p),La=rep(0.1,p),
                          mu = mean(y)
  )*/

  double M_norm;
  double TP_adjust;
  double lambda0 = exp( as<double>(parameters["lambda0"]));
  double nu = exp( as<double>(parameters["nu"]));

  ybar = y - as<double>(parameters["mu"]); // - parameters["mu_z"]
  alpha = invSmatn * ybar;
  M_norm = arma::dot(ybar, alpha) / lambda0;
  TP_adjust = (nu + n) / ((nu + M_norm) * lambda0 );
  tmpS = invSmatn - TP_adjust * alpha * alpha.t();

  //mu - gradient approach
  gradients["mu"] = TP_adjust * arma::sum( invSmatn * ybar );

  //nu
  gradients["nu"] = 0 * 0.5*boost::math::digamma( 0.5*(as<double>(parameters["nu"])+n) ) - 0.5*boost::math::digamma( 0.5*as<double>(parameters["nu"]) ) - 0.5*( ((n - M_norm)/(nu + M_norm))  + log(nu + M_norm) - log(nu) );

  //lambda0
  gradients["lambda0"] = - 0.5* lambda0 * (nu * (n - M_norm) )/(lambda0 * (nu + M_norm)) ;

  //Same as for GP
  //sigma
  dS = arma::diagmat( exp( as<double>(parameters["sigma"]) + as<double>(parameters["sigma_z"]) * z ) % w );
  gradients["sigma"] = evid_grad(tmpS, dS);

  //sigma_z
  dS.diag() = dS.diag() % z;
  gradients["sigma_z"] = evid_grad(tmpS, dS);

  //lambdam
  gradients["lambdam"] = evid_grad(tmpS, Sm);

  //Lm and La
  Rcpp::List tmp = evid_scale_gradients(X, tmpS, Sm, Sa, parameters);
  gradients["Lm"] = as<std::vector<double> >(tmp["m"]);
  gradients["La"] = as<std::vector<double> >(tmp["a"]);

  //RMSE
  stats(0) = pow(arma::norm(y - (Smat * alpha + as<double>(parameters["mu"]))),2);

  //Evidence
  double val;
  double sign;
  arma::log_det(val,sign,invSmatn);  //logdet of INVERSE -> -val
  stats(1) = lgamma( 0.5*(nu+n) ) - lgamma( 0.5 * nu) - 0.5 * (n * log(nu  * arma::datum::pi ) - val + n * log(lambda0) + (n + nu) * ( 1 + M_norm/nu ) );

  //output
  return Rcpp::List::create(Named("gradients") = gradients,Named("stats") = stats );
}

// [[Rcpp::export]]
arma::rowvec stats_TP_SE(arma::colvec y, arma::mat Kmat, arma::mat invKmatn, Rcpp::List parameters) {
  Rcpp::List gradients = clone(parameters); //for the same list structure

  unsigned int n = y.size();
  //preallocate memory
  arma::rowvec stats(2);
  arma::colvec ybar(n);
  arma::colvec alpha(n);
  double lambda0 = exp( as<double>(parameters["lambda0"]));
  double nu = exp( as<double>(parameters["nu"]));
  double M_norm = arma::dot(ybar, alpha) / lambda0;

  ybar = y - as<double>(parameters["mu"]); // - parameters["mu_z"]
  alpha = invKmatn * ybar;

  //RMSE
  stats(0) = pow(arma::norm(y - (Kmat * alpha + as<double>(parameters["mu"]))),2);

  //Evidence
  double val;
  double sign;
  arma::log_det(val,sign,invKmatn);  //logdet of INVERSE -> -val
  stats(1) = lgamma( 0.5*(nu+n) ) - lgamma( 0.5 * nu) - 0.5 * (n * log(nu  * arma::datum::pi ) - val + n * log(lambda0) + (n + nu) * ( 1 + M_norm/nu ) );

  return stats;
}
