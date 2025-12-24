// [[Rcpp::depends(RcppArmadillo, RcppGSL, RcppParallel)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_errno.h>
#include <Rmath.h>

using namespace Rcpp;
using namespace RcppParallel;

struct PdfS2Worker : public Worker {
  const RVector<double> y;
  const double nu;
  const double log_a; // Pré-calculado
  const double s;     // Pré-calculado
  RVector<double> result;

  PdfS2Worker(const NumericVector& y, double nu, double log_a, double s, NumericVector& result)
    : y(y), nu(nu), log_a(log_a), s(s), result(result) {}

  void operator()(std::size_t begin, std::size_t end) {
    double tol = 1e-8;

    for (std::size_t i = begin; i < end; i++) {
      double x = 0.5 * y[i] * y[i];
      double z;

      if (x < tol) {
        // Limite analítico quando y -> 0
        z = std::exp(log_a - std::log(s));
      } else {
        /* Sua fórmula original: a * (1/x)^s * P(s, x) * Gamma(s)
         Para evitar overflow, calculamos no espaço logarítmico:
         log(z) = log(a) - s*log(x) + log(P(s, x)) + log(Gamma(s))
         */
        double log_P = std::log(gsl_sf_gamma_inc_P(s, x));
        double log_gamma_s = gsl_sf_lngamma(s);

        z = std::exp(log_a - s * std::log(x) + log_P + log_gamma_s);
      }

      result[i] = z;
    }
  }
};

// [[Rcpp::export]]
NumericVector pdf_s(NumericVector y, double nu) {
  // 1. Desativa o manipulador de erros global da GSL (ESSENCIAL para paralelo)
  gsl_set_error_handler_off();

  // 2. Pré-cálculo de constantes fora do loop paralelo
  double log_a = std::log(nu) - 0.5 * std::log(2.0 * M_PI);
  double s = nu + 0.5;

  NumericVector result(y.size());

  // 3. Execução
  PdfS2Worker worker(y, nu, log_a, s, result);
  parallelFor(0, y.size(), worker);

  return result;
}
