// [[Rcpp::depends(RcppArmadillo, RcppGSL, RcppParallel)]]

#include <RcppArmadillo.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_errno.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;

struct NormPdfWorker : public Worker {
  const RVector<double> y;
  RVector<double> ndens;

  // Construtor
  NormPdfWorker(const NumericVector& y, NumericVector& ndens)
    : y(y), ndens(ndens) {}

  void operator()(std::size_t begin, std::size_t end) {
    // A função gsl_ran_gaussian_pdf calcula:
    // (1 / (sigma * sqrt(2*pi))) * exp(-y^2 / (2 * sigma^2))

    for (std::size_t i = begin; i < end; i++) {
      ndens[i] = gsl_ran_gaussian_pdf(y[i], 1.0);
    }
  }
};

// [[Rcpp::export]]
NumericVector pdf_n(NumericVector y) {
  // 1. Segurança: Desativa handler de erro da GSL
  gsl_set_error_handler_off();

  int n = y.size();
  NumericVector ndens(n);

  // 2. Execução Paralela
  NormPdfWorker worker(y, ndens);
  parallelFor(0, n, worker);

  return ndens;
}
