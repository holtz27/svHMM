// [[Rcpp::depends(RcppArmadillo, RcppGSL, RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_errno.h>

using namespace Rcpp;
using namespace RcppParallel;

struct TDistPdfWorker : public Worker {
  const RVector<double> y;
  const double df;
  RVector<double> result;

  TDistPdfWorker(const NumericVector& y, double df, NumericVector& result)
    : y(y), df(df), result(result) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      // gsl_ran_tdist_pdf é thread-safe se o error handler estiver off
      result[i] = gsl_ran_tdist_pdf(y[i], df);
    }
  }
};

// [[Rcpp::export]]
NumericVector pdf_t(NumericVector y, double df) {
  // 1. Essencial: Desativa o manipulador de erro global da GSL
  gsl_set_error_handler_off();

  NumericVector result(y.size());
  TDistPdfWorker worker(y, df, result);
  parallelFor(0, y.size(), worker);

  return result;
}
