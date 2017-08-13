// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
//#include <boost/math/special_functions/digamma.hpp>
//#include <boost/math/distributions/inverse_gamma.hpp>
//#include <cmath>

using namespace arma;
using namespace Rcpp;

//same as for GP_SE
double evid_grad_TP(arma::mat mymat, arma::mat dK) {
  //grad = - 0.5 * ( sum(diag( tmpK %*% dK ) ) )
  double grad = - 0.5 * trace( mymat * dK);
  return(grad);
}

//same as for GP_SE
Rcpp::List evid_scale_gradients_TP(arma::mat X, arma::mat tmpK,arma::mat Km,arma::mat Ka, Rcpp::List parameters){
  // produce a (p x 2) matrix "grad" with the gradients of Lm and La in the first and second column, respectively
  unsigned int n = X.n_rows;
  unsigned int p = X.n_cols;

  //predefine matrices and temporary storage files
  arma::vec grad_m(p);
  arma::vec grad_a(p);
  arma::mat tmpX(n,n);
  arma::mat tmpM(n,n);
  arma::mat tmpA(n,n);
  //double lambda0 = exp( as<double>(parameters["lambda0"]) );

  for(unsigned int i=0; i<p; i++){
    //tmpX = abs(matrix(rep(X[,i],n),n,n) - t(matrix(rep(X[,i],n),n,n)))
    for(unsigned int r=0; r< n; r++){
      tmpX.col(r) = pow(X( r , i ) - X.col(i),2);
    }

    //Klist$Km * exp( 2*log(tmpX) - parameters$Lm[i] )
    tmpM = Km % tmpX * exp(- as<arma::vec>(parameters["Lm"])[i]);
    tmpA = Ka % tmpX * exp(- as<arma::vec>(parameters["La"])[i]);

    grad_m(i) = evid_grad_TP(tmpK,tmpM);
    grad_a(i) = evid_grad_TP(tmpK,tmpA);
  }

  return Rcpp::List::create(Named("m") = grad_m,Named("a") = grad_a );
}

// lambda0 instead of lambdaa but otherwise identical to GP
// [[Rcpp::export]]
Rcpp::List kernmat_TP_SE_cpp(Rcpp::NumericMatrix X1,Rcpp::NumericMatrix X2,Rcpp::NumericVector Z1,Rcpp::NumericVector Z2, Rcpp::List para) {
  //input X1,X2,Z1,Z2,parameters

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
  double lambda0 = as<double>(para["lambda0"]);

  Rcpp::NumericMatrix tmpK(n1, n2);

  for(unsigned int r = 0; r < n1; r++){
    for(unsigned int c = 0; c < n2; c++){

      tmpKm(r,c) = exp(lambda_m - tmpKm(r,c)); //- lambda0
      tmpKa(r,c) = exp(lambda0 - tmpKa(r,c)) * Z1(r) * Z2(c); //"lambdaa" = 0
      tmpK(r,c)  = tmpKm(r,c) + tmpKa(r,c) ;

    }
  }

  //Put matrices in a list
  Rcpp::List out; out["Smat"] = tmpK; out["Sm"] = tmpKm; out["Sa"] = tmpKa;
  return(out);
}

//double log_multivariate_t_density(arma::colvec y,double mu, double nu, arma::mat Omega){
//  unsigned int n = y.size();
//
//  logp = lgamma((nu+n)/2 ) - lgamma( nu/2 ) - 0.5 * (n * log(nu*arma::datum::pi) + );
//
//  return(logp);
//}

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

  double M_norm;
  double TP_adjust;
  double lambda0 = exp( as<double>(parameters["lambda0"]) );
  double nu = exp( as<double>(parameters["nu"]) );

  ybar = y - as<double>(parameters["mu"]); // - parameters["mu_z"]
  alpha = invSmatn * ybar;
  M_norm = arma::dot(ybar,alpha);
  TP_adjust = ( nu + n ) / ( (nu + M_norm) * lambda0 ) ;
  //Rcpp::Rcout << M_norm << std::endl;
  tmpS = (invSmatn - TP_adjust * alpha * alpha.t());

  //mu - gradient approach
  gradients["mu"] =  TP_adjust * arma::sum( invSmatn * ybar ) ;

  //nu
  gradients["nu"] =  0.0; // 0.5 * nu * (2*boost::math::digamma( 0.5*(nu+n) ) - 2*boost::math::digamma( 0.5*nu ) - ( ((n - M_norm)/(nu + M_norm))  + log(1 + M_norm/nu) ));

  //lambda0 - similar as lambdaa for GP, when considering that hatlambdam cancels out
  gradients["lambda0"] = evid_grad_TP(tmpS, Sa);

  //Same as for GP
  //sigma
  dS = arma::diagmat( exp(as<double>(parameters["sigma"]) + as<double>(parameters["sigma_z"]) * z ) );
  gradients["sigma"] = evid_grad_TP(tmpS, dS);

  //sigma_z
  dS.diag() = dS.diag() % z;
  gradients["sigma_z"] = evid_grad_TP(tmpS, dS);

  //lambdam
  gradients["lambdam"] = evid_grad_TP(tmpS, Sm);

  //Lm and La
  Rcpp::List tmp = evid_scale_gradients_TP(X, tmpS, Sm, Sa, parameters);
  gradients["Lm"] = as<std::vector<double> >(tmp["m"]);
  gradients["La"] = as<std::vector<double> >(tmp["a"]);

  //RMSE
  stats(0) = pow(arma::norm(y - (Smat * alpha + as<double>(parameters["mu"]))),2);

  //Evidence
  double val;
  double sign;
  arma::log_det(val,sign,invSmatn);  //logdet of INVERSE -> -val
  stats(1) = lgamma( 0.5 * (nu+n) ) - lgamma( 0.5 * nu ) - 0.5*(n*log(nu * arma::datum::pi) - val + (nu+n) * log( 1 + M_norm/nu ) );

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

  ybar = y - as<double>(parameters["mu"]); // - parameters["mu_z"]
  alpha = invKmatn * ybar;

  //double lambda0 = exp( as<double>(parameters["lambda0"]));
  double nu = exp( as<double>(parameters["nu"]));
  double M_norm = arma::dot(ybar,alpha);

  //RMSE
  stats(0) = pow(arma::norm(y - (Kmat * alpha + as<double>(parameters["mu"]))),2);

  //Evidence
  double val;
  double sign;
  arma::log_det(val,sign,invKmatn);  //logdet of INVERSE -> -val
  stats(1) = lgamma( 0.5 * (nu+n) ) - lgamma( 0.5 * nu ) - 0.5*(n*log(nu * arma::datum::pi) - val + (nu+n) * log( 1 + M_norm/nu ) );

  return stats;
}
