# svHMM

**svHMM** is an R package for fast approximate Bayesian inference in **Stochastic Volatility in Mean models with heavy-tailed distributions**, using **Hidden Markov Models (HMMs)**.

The package is based on the methodology proposed in:

> *Stochastic Volatility in Mean Models with Heavy Tails: A Fast Approximate Bayesian Inference Using Hidden Markov Models*

---

## Installation

You can install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("https://github.com/holtz27/svHMM.git")
```

## Example
```r
library(svHMM)

# Simulate data
log.ret=svmt.sim(mu=1.0, 
                 phi=0.98, 
                 sigma=0.2, 
                 nu=10, 
                 beta=c(0.2, -0.01, -0.2),
                 y0=0.2,
                 g_dim=1e3)
plot(log.ret$y, type='l')

# Initial parameter values
theta0=list(mu0=log(var(log.ret$y)), 
            phi0=0.98, 
            sigma0=0.15, 
            nu0=8, 
            beta0=c(mean(log.ret$y),-.01,-0.2))

# Fit SVM-t HMM model
fit=svmtHMM(y=log.ret$y, m=200, gmax=5, theta_init=theta0, y0=0.2)

summary(fit)

lv=logvol(fit, plot=TRUE)
lines(exp(0.5*log.ret$h), col='red', lwd=2)
