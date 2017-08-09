// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;
//#include <Rcpp.h>
//using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List Nesterov_cpp(double learn_rate,double momentum, Rcpp::List nu, Rcpp::List grad, Rcpp::List para){
  
  unsigned int no_L = nu.size(); //excludes means
  Rcpp::List output_nu = clone(nu);
  Rcpp::List output_para = clone(para);
  //Rcpp::NumericVector tmp;
    
  //for each vector in the list
  for(unsigned int i=0; i < no_L; i++){
    output_nu[i] = momentum * as<NumericVector>(nu[i]) + learn_rate * as<NumericVector>(grad[i]);
    output_para[i] = as<NumericVector>(para[i]) + as<NumericVector>(nu[i]);
    
    //tmp = as<NumericVector>(nu[i]) + learn_rate * as<NumericVector>(grad[i]);
      
    //Rcpp::Rcout << "in " << as<Rcpp::NumericVector>(para[i]) << std::endl;  
    //Rcpp::Rcout << "nu " << as<Rcpp::NumericVector>(nu[i]) << std::endl;  
    ////Rcpp::Rcout << "nu + grad " << tmp << std::endl;  
    //Rcpp::Rcout << "outnu " << as<Rcpp::NumericVector>(output_nu[i]) << std::endl;     
    //Rcpp::Rcout << "grad " << as<Rcpp::NumericVector>(grad[i]) << std::endl; 
    //Rcpp::Rcout << "out " << as<Rcpp::NumericVector>(output_para[i]) << std::endl;  
  }
  
  //Output
  Rcpp::List out; out["nu"]=output_nu; out["parameters"] = output_para;
  return(out);
}

// [[Rcpp::export]]
Rcpp::List Nadam_cpp(double iter,double learn_rate,double beta1, double beta2, double eps, Rcpp::List m, Rcpp::List v, Rcpp::List grad, Rcpp::List para){
  //number of list elements
  unsigned int no_L = m.size() ; //excludes means
  //Rcpp::List grad = clone(gradin);
  Rcpp::List output_m = clone(m);
  Rcpp::List output_v = clone(v);
  Rcpp::List output_para = clone(para);
  
  //for(unsigned int i=0; i < no_L; i++){
  //  Rcpp::Rcout << "grad --  " << as<Rcpp::NumericVector>(grad[i]) << std::endl; 
  //}
    
  //for each vector in the list
  for(unsigned int i=0; i < no_L; i++){
    output_m[i] = (beta1 * as<Rcpp::NumericVector>(m[i]) + (1-beta1) * as<Rcpp::NumericVector>(grad[i]));
    output_v[i] = (beta2 * as<Rcpp::NumericVector>(v[i]) + (1-beta2) * Rcpp::pow(as<Rcpp::NumericVector>(grad[i]),2));     
    output_para[i] = as<Rcpp::NumericVector>(para[i]) + learn_rate * ( ( beta1 * (as<Rcpp::NumericVector>(output_m[i]))  + (1-beta1)*(as<Rcpp::NumericVector>(grad[i])) ) / (1-pow(beta1,iter)) ) / ( sqrt(as<Rcpp::NumericVector>(output_v[i]) / (1-pow(beta2,iter)) ) + eps);
    
    //Rcpp::Rcout << "in " << as<Rcpp::NumericVector>(para[i]) << std::endl;  
    //Rcpp::Rcout << "m " << as<Rcpp::NumericVector>(m[i]) << std::endl;  
    //Rcpp::Rcout << "outm " << as<Rcpp::NumericVector>(output_m[i]) << std::endl;     
    //Rcpp::Rcout << "v " << as<Rcpp::NumericVector>(v[i]) << std::endl;  
    //Rcpp::Rcout << "grad " << as<Rcpp::NumericVector>(grad[i]) << std::endl; 
    //Rcpp::Rcout << "out " << as<Rcpp::NumericVector>(output_para[i]) << std::endl;  
  }
  return Rcpp::List::create(Named("m") = output_m, Named("v") = output_v,Named("parameters") = output_para );
}

// [[Rcpp::export]]
Rcpp::List Adam_cpp(double iter,double learn_rate,double beta1, double beta2, double eps, Rcpp::List m, Rcpp::List v, Rcpp::List grad, Rcpp::List para){
  //number of list elements
  unsigned int no_L = m.size() ; //excludes means
  //Rcpp::List grad = clone(gradin);
  Rcpp::List output_m = clone(m);
  Rcpp::List output_v = clone(v);
  Rcpp::List output_para = clone(para);
  
  //for(unsigned int i=0; i < no_L; i++){
  //  Rcpp::Rcout << "grad --  " << as<Rcpp::NumericVector>(grad[i]) << std::endl; 
  //}
  
  //for each vector in the list
  for(unsigned int i=0; i < no_L; i++){
    output_m[i] = (beta1 * as<Rcpp::NumericVector>(m[i]) + (1-beta1) * as<Rcpp::NumericVector>(grad[i]));
    output_v[i] = (beta2 * as<Rcpp::NumericVector>(v[i]) + (1-beta2) * Rcpp::pow(as<Rcpp::NumericVector>(grad[i]),2));     
    output_para[i] = as<Rcpp::NumericVector>(para[i]) + learn_rate * ( (as<Rcpp::NumericVector>(output_m[i]))  / (1-pow(beta1,iter)) ) / ( sqrt(as<Rcpp::NumericVector>(output_v[i]) / (1-pow(beta2,iter)) ) + eps);
    
  }
  return Rcpp::List::create(Named("m") = output_m, Named("v") = output_v,Named("parameters") = output_para );
}
