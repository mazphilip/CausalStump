// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <cmath>

using namespace arma;
using namespace Rcpp;

double evid_grad(arma::mat mymat, arma::mat dK) {
  //grad = - 0.5 * ( sum(diag( tmpK %*% dK ) ) )

  double grad = - 0.5 * trace( mymat * dK);

  return(grad);
}

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
    //Rcpp::Rcout << grad_m(i) << std::endl;
    grad_a(i) = evid_grad(tmpK,tmpA);
  }

  return Rcpp::List::create(Named("m") = grad_m,Named("a") = grad_a );
}

// [[Rcpp::export]]
arma::mat invkernel_cpp(arma::colvec z,arma::colvec w, arma::mat mymat, Rcpp::List parameters){
  //unsigned int n = z.n_cols;
  //add to diagonal, for loop to minimize unnecessary memory access
  /*for(unsigned int i=0; i < n; i++){
    mat(i,i) = mat(i,i) + exp(as<double>(parameters["sigma"]) + as<double>(parameters["sigma_z"]) * z(i) );
  }*/
  //neater:

  mymat.diag() += exp(as<double>(parameters["sigma"]) + as<double>(parameters["sigma_z"]) * z ) % w;
  mymat = inv_sympd(mymat);

  return(mymat);
}

// [[Rcpp::export]]
double mu_solution_cpp(arma::colvec y, arma::colvec z, arma::mat invKmat, Rcpp::List parameters) {
  //calculates the exact solutions to the maximization problem
  double mu = 0.5 * arma::sum(invKmat * y ) / arma::sum(arma::sum(invKmat));

  //mu(0) = 0.5 * arma::sum(invKmat * ( y - as<double>(parameters["mu_z"]) * z )) / arma::sum(arma::sum(invKmat));
  //mu(1) = 0.5 * arma::sum( (z.t() * invKmat * ( y - as<double>(parameters["mu"]) ) ) / ( z.t() * invKmat * z ) );

  return(mu);
}

// [[Rcpp::export]]
Rcpp::List kernmat_GP_SE_cpp(Rcpp::NumericMatrix X1,Rcpp::NumericMatrix X2,Rcpp::NumericVector Z1,Rcpp::NumericVector Z2, Rcpp::List para) {
  //input X1,X2,Z1,Z2,parameters

  //for(i in 1:p){
  //  tmpX = (matrix(rep(X1[,i],n2),n1,n2) - t(matrix(rep(X2[,i],n1),n2,n1)))^2
  //  Km        = Km        + tmpX        * exp(-parameters$Lm[i])
  //  Ka[Zmask] = Ka[Zmask] + tmpX[Zmask] * exp(-parameters$La[i])
  //}

  unsigned int n1 = X1.nrow();
  unsigned int n2 = X2.nrow();
  unsigned int p = X2.ncol();
  /*
  NumericMatrix tmpKm(n1, n2);
  NumericMatrix tmpKa(n1, n2);
  std::fill(tmpKm.begin(), tmpKm.end(), 0);
  std::fill(tmpKa.begin(), tmpKa.end(), 0);
  */
  arma::mat tmpKm(n1,n2); tmpKm.zeros();
  arma::mat tmpKa(n1,n2); tmpKa.zeros();

  arma::rowvec tmprow(n2);

  NumericVector Lm = as<NumericVector>(para["Lm"]);
  NumericVector La = as<NumericVector>(para["La"]);

  //Rcpp::Rcout << p <<std::endl;
  //Rcpp::Rcout << Lm.size() << std::endl;
  //Rcpp::Rcout << La.size() << std::endl;
  //double tmp;
  //arma::rowvec tmp2(n2);

  for(unsigned int i=0;i < p;i++){
    //for every output row
    for(unsigned int r=0;r < n1;r++){
      tmprow = pow(X1( r , i ) - X2( _, i),2);
      //Rcpp::Rcout << tmprow << std::endl;

      //Rcpp::Rcout << std::endl << tmp << std::endl;
      //Rcpp::Rcout << std::endl << tmp2 << std::endl;
      tmpKm.row(r) += tmprow * exp(- Lm(i) );
      tmpKa.row(r) += tmprow * exp(- La(i) );
      //tmpKm(r,_) = tmpKm(r,_) + tmp;
      //tmpKa(r,_) = tmpKa(r,_) + tmp2;

    }
  }

  double lambda_m = as<double>(para["lambdam"]);
  double lambda_a = as<double>(para["lambdaa"]);
  //Rcpp::Rcout << "lambdas" <<std::endl;
  //Rcpp::Rcout << lambda_m << "  " << lambda_a <<std::endl;

  Rcpp::NumericMatrix tmpK(n1, n2);

  for(unsigned int r = 0; r < n1; r++){
    for(unsigned int c = 0; c < n2; c++){
      tmpKm(r,c) = exp(lambda_m - tmpKm(r,c));
      tmpKa(r,c) = exp(lambda_a - tmpKa(r,c)) * Z1(r) * Z2(c);
      tmpK(r,c) = tmpKm(r,c) + tmpKa(r,c) ;
    }
  }

  //Put matrices in a list
  Rcpp::List out; out["Kmat"] = tmpK; out["Km"] = tmpKm; out["Ka"] = tmpKa;
  return(out);
  //return 0;
}


// [[Rcpp::export]]
Rcpp::List grad_GP_SE_cpp(arma::colvec y, arma::mat X, arma::colvec z,arma::colvec w, arma::mat Kmat, arma::mat Km, arma::mat Ka, arma::mat invKmatn, Rcpp::List parameters) {
  Rcpp::List gradients = clone(parameters); //for the same list structure

  unsigned int n = X.n_rows;
  //unsigned int p = X.n_cols;

  //preallocate memory
  arma::colvec stats(2);
  arma::colvec ybar(n);
  arma::colvec alpha(n);
  arma::mat tmpK(n,n);
  arma::mat dK(n,n);

  ybar = y - as<double>(parameters["mu"]); // - parameters["mu_z"]
  alpha = invKmatn * ybar;
  tmpK = invKmatn - alpha * alpha.t();

  //sigma
  dK = arma::diagmat( exp( as<double>(parameters["sigma"]) + as<double>(parameters["sigma_z"]) * z ) % w );
  gradients["sigma"] = evid_grad(tmpK, dK);

  //sigma_z
  dK.diag() = dK.diag() % z;
  gradients["sigma_z"] = evid_grad(tmpK, dK);

  //lambdam
  gradients["lambdam"] = evid_grad(tmpK, Km);

  //lambdaa
  gradients["lambdaa"] = evid_grad(tmpK, Ka);

  //Lm and La
  Rcpp::List tmp = evid_scale_gradients(X, tmpK, Km, Ka, parameters);
  //arma::vec tmp(2,fill::zeros);

  gradients["Lm"] = as<std::vector<double> >(tmp["m"]);
  gradients["La"] = as<std::vector<double> >(tmp["a"]);

  //mu - gradient approach

  gradients["mu"] = arma::sum( invKmatn * ybar );

  //mu_z?
  //gradients["mu"] = z.t() * invKmatn * ybar ;

  //Might remove stats but it is more efficient to calculate it here

  //RMSE
  stats(0) = pow(arma::norm(y - (Kmat * alpha + as<double>(parameters["mu"]))),2);

  //Evidence
  double val;
  double sign;
  arma::log_det(val,sign,invKmatn);  //logdet of INVERSE -> -val
  stats(1) = - 0.5 * (n * log( 2.0 * arma::datum::pi ) - val + arma::dot( ybar, alpha ) ) ;

  //output
  return Rcpp::List::create(Named("gradients") = gradients,Named("stats") = stats );
}

// [[Rcpp::export]]
arma::rowvec stats_GP_SE(arma::colvec y, arma::mat Kmat, arma::mat invKmatn, Rcpp::List parameters) {
  Rcpp::List gradients = clone(parameters); //for the same list structure

  unsigned int n = y.size();
  //preallocate memory
  arma::rowvec stats(2);
  arma::colvec ybar(n);
  arma::colvec alpha(n);

  ybar = y - as<double>(parameters["mu"]); // - parameters["mu_z"]
  alpha = invKmatn * ybar;

  //RMSE
  stats(0) = pow(arma::norm(y - (Kmat * alpha + as<double>(parameters["mu"]))),2);

  //Evidence
  double val;
  double sign;
  arma::log_det(val,sign,invKmatn);  //logdet of INVERSE -> -val
  stats(1) = - 0.5 * (n * log( 2.0 * arma::datum::pi ) - val + arma::dot( ybar, alpha ) ) ;

  return stats;
}





