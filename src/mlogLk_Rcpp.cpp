// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <cmath>

//using namespace std;
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double mlogLk_Rcpp(arma::mat allprobs, arma::mat egamma, arma::rowvec foo, int n){

  //if(allprobs.has_nan()) stop("Probs must be a number");
  //if(allprobs.has_inf()) stop("Probs must be a finite number");

  //if(egamma.has_nan()) stop("Transition Probs must be a number");
  //if(egamma.has_inf()) stop("Transition Probs must be a finite number");

  //if(foo.has_nan()) stop("Initial Probs must be a number");
  //if(foo.has_inf()) stop("Initial Probs must be a finite number");

  double lscale=0;
  double sumfoo;
  rowvec pivot;

  for(int i=1; i<n; i++){
    pivot = foo*egamma;
    foo = pivot%allprobs.row(i);

    sumfoo=sum(foo);
    if(sumfoo!=0){
      foo=foo/sumfoo;
    }else{
      stop("sumfoo must be positive!");
    }

    lscale=lscale+log(sumfoo);
  }

  return lscale;
}
