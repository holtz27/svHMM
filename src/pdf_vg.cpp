#include <RcppArmadillo.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <RcppParallel.h>

// [[Rcpp::depends(RcppArmadillo, RcppGSL, RcppParallel)]]

using namespace Rcpp;
using namespace RcppParallel;

struct PdfVG2Worker : public Worker {
  const RVector<double> y;
  const double nu;
  const double log_const_fix; // Log da parte constante da densidade

  // Output
  RVector<double> vgdens;

  // Construtor: Inicializa os dados e pré-calcula o que for constante
  PdfVG2Worker(const NumericVector& y, double nu, NumericVector& vgdens, double log_const_fix)
    : y(y), nu(nu), log_const_fix(log_const_fix), vgdens(vgdens) {}

  void operator()(std::size_t begin, std::size_t end) {
    double tol = 1e-8;
    double sqrt_nu = std::sqrt(nu);
    double nu_minus_1_half = 0.5 * (nu - 1.0);

    for(std::size_t i = begin; i < end; i++) {
      double abs_yi = std::abs(y[i]);

      if(abs_yi < tol) {
        // Caso limite para y = 0
        // f(0) = (sqrt(nu)/sqrt(pi)) * (Gamma((nu-1)/2) / (2 * Gamma(nu/2)))
        double log_val = std::log(0.5) + 0.5 * std::log(nu/M_PI) +
          gsl_sf_lngamma(nu_minus_1_half) - gsl_sf_lngamma(0.5 * nu);
        vgdens[i] = std::exp(log_val);
      } else {
        // Cálculo via Log-Space para evitar Underflow/Overflow
        // O log da densidade VG (simplificada) envolve:
        // log(const) + log((|y|*sqrt(nu)/2)^((nu-1)/2)) + log(K_{ (nu-1)/2 }( |y|*sqrt(nu) ))

        double bessel_ln_val = gsl_sf_bessel_lnKnu(nu_minus_1_half, abs_yi * sqrt_nu);
        double term_pow = nu_minus_1_half * std::log(0.5 * abs_yi * sqrt_nu);

        vgdens[i] = std::exp(log_const_fix + term_pow + bessel_ln_val);
      }
    }
  }
};

// [[Rcpp::export]]
NumericVector pdf_vg(NumericVector y, double nu) {
  // 1. Desativa o manipulador de erros da GSL (Crucial para multi-thread)
  gsl_set_error_handler_off();

  int n = y.size();
  NumericVector vgdens(n);

  // 2. Pré-calculo da constante logarítmica: log( sqrt(nu/pi) / Gamma(nu/2) )
  double log_const_fix = 0.5 * std::log(nu / M_PI) - gsl_sf_lngamma(0.5 * nu);

  // 3. Criar o Worker e executar em paralelo
  PdfVG2Worker worker(y, nu, vgdens, log_const_fix);
  parallelFor(0, n, worker);

  return vgdens;
}
