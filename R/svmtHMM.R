
svmt.pn2pw = function(beta, mu, phi, sigma, nu){
  lbeta1 = beta[1]
  lbeta2 =  atanh(beta[2]) #log((1+beta[2])/(1-beta[2]))
  lbeta3 = beta[3]
  lmu = mu
  lphi = atanh(phi) #log((1+phi)/(1-phi))
  lsigma = log(sigma)
  # 2<nu<40
  #lnu = log(nu-2)-log(40-nu)
  alpha=0.1
  lnu=(2/alpha)*atanh( (2*nu-40-2)/(40-2))
  parvect = c(lbeta1,lbeta2,lbeta3,lmu,lphi,lsigma,lnu)
  return(parvect)
}
svmt.pw2pn <- function(parvect){
  beta=array(0,dim=3)
  beta[1]= parvect[1]
  beta[2]=tanh(parvect[2]) #(exp(parvect[2])-1)/(exp(parvect[2])+1)
  beta[3]=parvect[3]
  mu=parvect[4]
  phi=tanh(parvect[5]) #(exp(parvect[5])-1)/(exp(parvect[5])+1)
  sigma=exp(parvect[6])
  #nu = exp(parvect[7]) + 2
  #nu = (40*exp(parvect[7])+2)/(1+exp(parvect[7]))
  alpha=0.1
  nu=0.5*((40-2)*tanh(0.5*alpha*parvect[7])+(40+2))
  #return(list(beta=beta,mu=mu,phi=phi,sigma=sigma,nu=nu))
  return(c(beta, mu, phi, sigma, nu))
}
## function that will be used to compute 'allprobs' in mllk below
fillallprobs_t <- function(x,beg,beta,nu,y){
  z = (x-beta[1]-beta[2]*y-beta[3]*beg^2)/beg
  #w = pdf_vg(y=z, mu=0, sigma=1, nu=nu)
  w = pdf_t(y=z, df=nu)
  return(w/beg)
}
svmt.mllk <-function(parvect,y,y0,m,gmax){
  ny = length(y)
  p = svmt.pw2pn(parvect)
  K = m+1
  b=seq(-gmax,gmax,length=K)
  bs=(b[-1]+b[-K])*0.5
  E=p[4]+p[5]*(bs-p[4])
  intlen <- b[2]-b[1]
  sey = exp(bs/2)
  Gamma=matrix(0,m,m) #06
  for (i in 1:m){
    Gamma[i,]=dnorm(bs, mean=E[i], sd=p[6])
    sg = sum(Gamma[i,])
    if(sg == 0){
      stop('Built Gamma error')
    }else{
      Gamma[i,] = Gamma[i,]/sg
    }
  }
  Gamma = intlen*Gamma
  #Gamma = Gamma/apply(Gamma,1,sum)

  xx<-y
  yy<-c(y0,y[1:(ny-1)])
  allprobs = outer(xx,sey,"fillallprobs_t", beta=p[1:3], nu=p[7], yy)
  delta=dnorm(bs,p[4],p[6]/sqrt(1-p[5]^2))*intlen
  foo = delta*allprobs[1,]
  lscale = mlogLk_Rcpp(allprobs,Gamma,foo,ny) #Rcpp function
  return(-lscale)
}
svmt.prior = function(parvect){

  # b0
  lprior = dnorm(parvect[1], 0, sqrt(10), log=TRUE )
  # b1
  x=0.5*(tanh(parvect[2])+1)
  j=abs(0.5/cosh(parvect[2])^2)
  lprior=lprior+dbeta(x, shape1=5, shape2=1.5, log=TRUE)+log(j)
  # b2
  lprior=lprior+dnorm(parvect[3], 0, sqrt(10), log=TRUE)
  # mu
  lprior=lprior+dnorm(parvect[4], 0, sqrt(10), log=TRUE)
  # phi
  x=0.5*(tanh(parvect[5])+1)
  j=abs(0.5/cosh(parvect[5])^2)
  lprior=lprior+dbeta(x, shape1=20, shape2=1.5, log=TRUE)+log(j)
  # sigma
  x=exp(2*parvect[6])
  j=2*x
  lprior=lprior+invgamma::dinvgamma(x, shape=2.5, rate=0.025, log=TRUE)+log(j)
  # nu
  lprior=lprior+dnorm(parvect[7], -10, 10, log=TRUE)


  #+ log(dnorm(parvect[2], 0.5, 10))
  #+ log(dnorm(parvect[3], 0, 10))
  #+ log(dnorm(parvect[4], 0, 10))
  #+ log(dnorm(parvect[5], 4.5, 10))
  #+ log(dnorm(parvect[6], -1.5, 10))
  #+ log(dnorm(parvect[7], -10, 10))

  return(-lprior)
}
svmt.posterior <-function(parvect,y,y0,m,gmax){
  return(svmt.mllk(parvect,y,y0,m,gmax)+svmt.prior(parvect))
}
svmt.map <- function(y, m, parvect0, y0, gmax){

  #parvect <- svm.pn2pw(beta=beta0,mu=mu0,phi=phi0,sigma=sigma0,nu=nu0)
  mod <- optim(parvect0,svmt.posterior,y=y,y0=y0,m=m,gmax=gmax,hessian=TRUE)
  mode <- mod$par
  conv = mod$convergence
  return(list(mode=mode,
              lpostsvvg=-mod$value,
              hessian=mod$hessian,
              conv = conv))
}

#' @export
svmt.sim = function(mu,phi,sigma,nu,beta,y0,g_dim){

  y=array(0,dim=g_dim)
  h=array(0,dim=g_dim)
  l=rgamma(g_dim,shape=0.5*nu,rate=0.5*nu)

  b0=beta[1]
  b1=beta[2]
  b2=beta[3]

  h[1]=rnorm(1,mu,sigma/sqrt(1-phi^2))
  y[1]= b0+b1*y0+b2*exp(h[1])+(exp(0.5*h[1])/sqrt(l[1]))*rnorm(1)
  for (j in 2:g_dim){
    h[j]=mu+phi*(h[j-1]-mu)+rnorm(1,0,sigma)
    y[j]=b0+b1*y[j-1]+b2*exp(h[j])+(exp(0.5*h[j])/sqrt(l[j]))*rnorm(1)
  }

  return (list(y=y,h=h,l=l))
}

#' Fit SV-t HMM model
#'
#' @export
svmtHMM <- function(y, y0=0.2, theta_init, m, gmax, nIS=1e3){

  time <- Sys.time()

  mu0    <- theta_init$mu0
  phi0   <- theta_init$phi0
  sigma0 <- theta_init$sigma0
  beta0  <- theta_init$beta0
  nu0    <- theta_init$nu0

  ##########################################################################
  # Optimization
  eigenvalues <- rep(0, 7)
  cont <- 0
  parvect0 <- svmt.pn2pw(beta0, mu0, phi0, sigma0, nu0)

  cat("Optimization step...\n")
  optim_time <- Sys.time()
  while (any(eigenvalues <= 0) && cont < 5) {
    optim_res <- svmt.map(y=y, m=m, parvect0=parvect0, y0=y0, gmax=gmax)
    H1 <- signif(solve(optim_res$hessian), 6)
    eigenvalues <- eigen(H1, only.values=TRUE)$values
    parvect0 <- optim_res$mode
    cont <- cont + 1
  }
  optim_time <- as.numeric(Sys.time() - optim_time, units = "mins")
  cat("Done. Elapsed time:", optim_time, "min\n")

  if(any(eigenvalues <= 0)){
    stop("Hessian is not positive definite")
  }

  ##########################################################################
  # Importance Sampling
  cat("\nImportance Sampling step!\n")
  is_time <- Sys.time()

  k   <- -optim_res$lpostsvvg
  map <- optim_res$mode

  X <- mvtnorm::rmvnorm(nIS, map, H1)
  #df=10
  #X=mvtnorm::rmvt(nIS, delta=map, sigma=H1, df=df, type='shifted')
  Weigth <- numeric(nIS)

  largura=40
  ultimo_print=-5 # garante que 0% seja impresso
  for (j in seq_len(nIS)) {
    #############################
    progresso=j/nIS
    pct=floor(progresso*100)
    if (pct %% 5 == 0 && pct != ultimo_print || j == nIS){
      preenchido = round(largura * progresso)
      barra = paste0( "[", paste(rep("=", preenchido), collapse = ""),
                      paste(rep(" ", largura - preenchido), collapse = ""), "]" )
      cat(sprintf("\r%s %3.0f%%", barra, pct))
      ultimo_print=pct
    }
    Sys.sleep(0.02)
    #############################
    Weigth[j] <- exp(
      k -
      svmt.posterior(X[j,], y, y0, m, gmax) -
      #mvtnorm::dmvt(X[j,],delta=map,sigma=H1,type='shifted',log=TRUE, df=df)
      mvtnorm::dmvnorm(X[j,], map, H1, log=TRUE)
    )
  }

  s=sum(Weigth)
  if((s != 0) && !is.nan(s)){
    Weigth=Weigth/s
  }else{
    stop('Error normalize constante weigths!')
  }
  # ESS
  #ess=1/sum(Weigth^2)
  #cat('\n', 'Raw ess: ',ess, '\n')

  ### Pareto Smooth
  #psis=loo::psis(log(Weigth), r_eff=ess)
  #Weigth=exp(psis$log_weights)
  #Weigth=Weigth/sum(Weigth)
  #cat('\n', 'Pareto k: ', psis$diagnostics$pareto_k, '\n')

  # Resample
  indx = sample(seq_len(nIS), prob=Weigth, replace=TRUE)
  X = X[indx, ]
  Weigth = rep(1 / nIS, nIS)

  is_time <- as.numeric(Sys.time() - is_time, units = "mins")
  cat('\n', "Done. Elapsed time:", is_time, "min\n")

  times <- as.numeric(Sys.time() - time, units = "mins")

  structure(
    list(
      Weigth=Weigth,
      X=X,
      y=y,
      y0=y0,
      m=m,
      gmax=gmax,
      times=times,
      call=match.call()
    ),
    class = "svmtHMM"
  )
}

#' @export
summary.svmtHMM <- function(object, ...) {
  ISdiag(
    Weigth = object$Weigth,
    X = object$X,
    normal = FALSE,
    p = svmt.pw2pn
  )
}

#' @export
info_crit <- function(object, ...) {
  UseMethod("info_crit")
}

#' @export
info_crit.svmtHMM <- function(object, ...) {

  theta_hat <- colMeans(object$X)

  D <- 2 * svmt.mllk(
    parvect = theta_hat,
    y = object$y,
    y0 = object$y0,
    m = object$m,
    gmax = object$gmax
  )

  Dbar <- 0
  for (j in seq_along(object$Weigth)) {
    Dbar <- Dbar +
      object$Weigth[j] *
      svmt.mllk(object$X[j, ], object$y, object$y0, object$m, object$gmax)
  }

  Dbar <- 2 * Dbar
  pd <- Dbar - D

  list(
    DIC = D + 2 * pd,
    LPS = 0.5 * D / length(object$y)
  )
}

#' @export
logvol <- function(object, ...) {
  UseMethod("logvol")
}

#' @export
logvol.svmtHMM <- function(object, plot=FALSE, ...) {

  X      <- object$X
  Weigth <- object$Weigth

  # Transform parameters from working space to natural space
  Thetas <- t(apply(X = X, MARGIN = 1, FUN = svmt.pw2pn))

  # Weighted posterior mean
  wThetas <- apply(X = Thetas, MARGIN = 2, FUN = "*", Weigth)
  theta_hat <- apply(X = wThetas, MARGIN = 2, sum)

  lv=svm.viterbi(
    y = object$y,
    y0 = object$y0,
    theta_hat = theta_hat,
    m = object$m,
    gmax = object$gmax,
    model = "fillallprobs_t"
  )

  ## Optional plot
  if (isTRUE(plot)) {

    plot(
      abs(object$y),
      type = "l",
      col  = "gray",
      ylab = "Magnitude",
      xlab = "Time",
      main = "Observed returns vs fitted volatility",
      ...
    )

    lines(
      exp(0.5*lv),
      lwd=2
    )
  }

  return(lv)
}


