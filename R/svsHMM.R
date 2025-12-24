

svs.pn2pw <- function(mu, phi, sigma, nu){

  lmu = mu
  lphi = atanh(phi) #log((1+phi)/(1-phi))
  lsigma = log(sigma)
  lnu=log(nu)
  parvect = c(lmu,lphi,lsigma,lnu)
  return(parvect)
}
svs.pw2pn <- function(parvect){

  mu=parvect[1]
  phi=tanh(parvect[2]) #(exp(parvect[5])-1)/(exp(parvect[5])+1)
  sigma=exp(parvect[3])
  nu=exp(parvect[4])
  #return(list(beta=beta,mu=mu,phi=phi,sigma=sigma,nu=nu))
  return(c(mu, phi, sigma, nu))
}
## function that will be used to compute 'allprobs' in mllk below
fillallprobs_s2 <- function(x,beg,nu){
  z = x/beg
  w = pdf_s(y=z, nu=nu)
  #w = pdf_s(z,0.0,1.0,nu)
  return(w/beg)
}
svs.mllk <-function(parvect,y,m,gmax){
  ny = length(y)
  p = svs.pw2pn(parvect)
  K = m+1
  b=seq(-gmax,gmax,length=K)
  bs=(b[-1]+b[-K])*0.5
  E=p[1]+p[2]*(bs-p[1])
  intlen <- b[2]-b[1]
  sey = exp(bs/2)
  Gamma=matrix(0,m,m) #06
  for (i in 1:m){
    Gamma[i,]=dnorm(bs, mean=E[i], sd=p[3])
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
  #yy<-c(y0,y[1:(ny-1)])
  allprobs = outer(xx,sey,"fillallprobs_s2",nu=p[4])
  delta=dnorm(bs,p[1],p[3]/sqrt(1-p[2]^2))*intlen
  foo = delta*allprobs[1,]
  lscale = mlogLk_Rcpp(allprobs,Gamma,foo,ny) #Rcpp function
  return(-lscale)
}
svs.prior <-function(parvect){

  # mu
  lprior=dnorm(parvect[1], 0, sqrt(10), log=TRUE)
  # phi
  x=0.5*(tanh(parvect[2])+1)
  j=abs(0.5/cosh(parvect[2])^2)
  lprior=lprior+dbeta(x, shape1=20, shape2=1.5, log=TRUE)+log(j)
  # sigma
  x=exp(2*parvect[3])
  j=2*x
  lprior=lprior+invgamma::dinvgamma(x, shape=2.5, rate=0.025, log=TRUE)+log(j)
  # nu
  v=exp(parvect[4])
  j=abs(v)
  lprior=lprior+dgamma(v, shape=0.08, rate=0.04) + log(j)

  #+ log(dnorm(parvect[2], 0.5, 10))
  #+ log(dnorm(parvect[3], 0, 10))
  #+ log(dnorm(parvect[4], 0, 10))
  #+ log(dnorm(parvect[5], 4.5, 10))
  #+ log(dnorm(parvect[6], -1.5, 10))
  #+ log(dnorm(parvect[7], -10, 10))
  return(-lprior)
}
svs.posterior <-function(parvect,y,m,gmax){
  return(svs.mllk(parvect,y,m,gmax)+svs.prior(parvect))
}
svs.map <- function(y, m, parvect0, gmax){

  mod <- optim(parvect0, svs.posterior,y=y,m=m,gmax=gmax,hessian=T)
  mode <- mod$par
  conv = mod$convergence
  return(list(mode=mode,
              lpostsvs=-mod$value,
              hessian=mod$hessian,
              conv = conv))
}

#' @export
svs.sim <-function(mu,phi,sigma,nu,g_dim){

  y=array(0,dim=g_dim)
  h=array(0,dim=g_dim)
  l=array(0,dim=g_dim)
  l=rbeta(g_dim,nu,1,ncp=0)

  h[1]=rnorm(1,mu,sigma/sqrt(1-phi^2))
  y[1]= (exp(0.5*h[1])/sqrt(l[1]))*rnorm(1)
  for (j in 2:g_dim){
    h[j]=mu+phi*(h[j-1]-mu)+rnorm(1,0,sigma)
    y[j]=(exp(0.5*h[j])/sqrt(l[j]))*rnorm(1)
  }
  return (list(y=y,h=h,l=l))
}

#' @export
svsHMM=function(y,  theta_init, m, gmax, nIS=1e3){

  mu0=theta_init$mu0
  phi0=theta_init$phi0
  sigma0=theta_init$sigma0
  nu0=theta_init$nu0

  time = Sys.time()
  ############################################################################
  # Optimize
  eigenvalues = rep(0, 4)
  cont = 0
  parvect0 = svs.pn2pw(mu0, phi0, sigma0, nu0)
  cat("Optimization step...\n")
  optim_time <- Sys.time()
  while( any(eigenvalues <= 0) && cont<5 ){

    optim_res = svs.map(y=y,m=m, parvect0=parvect0, gmax=gmax)
    H1 = signif(solve(optim_res$hessian), 6)
    eigenvalues = eigen(H1, only.values=TRUE)$values
    parvect0 = optim_res$mode
    cont = cont + 1

  }
  optim_time <- as.numeric(Sys.time() - optim_time, units = "mins")
  cat("Done. Elapsed time:", optim_time, "min\n")
  ############################################################################
  # Test if H1 is positive definite
  if(any(eigenvalues <= 0)){
    stop('Hessian is not positive definite')
  }
  ############################################################################
  cat("\nImportance Sampling step!\n")
  is_time <- Sys.time()

  k=-optim_res$lpostsvs
  map=optim_res$mode
  ##########################################################################
  # Weigths Evaluation
  #n=1e3
  X=mvtnorm::rmvnorm(nIS, map, H1)
  Weigth <- array(0,dim=c(nIS,1))
  largura=40
  ultimo_print=-5 # garante que 0% seja impresso
  for(j in 1:nIS){
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
    Weigth[j,1]=exp(k
                    -svs.posterior(parvect=X[j,],y=y,m=m,gmax=gmax)
                    -mvtnorm::dmvnorm(X[j,],mean=map,sigma=H1,log=TRUE)
                    )
  }
  s = sum(Weigth)
  if((s != 0) && !is.nan(s)){
    Weigth=Weigth/s
  }else{
    stop('Error normalize constante weigths!')
  }
  #times=as.numeric(Sys.time()-time, units='mins')
  ### Resample
  indx=sample(1:nIS, prob=Weigth, replace=TRUE)
  X=X[indx,]
  Weigth=rep(1/nIS, nIS)
  #Results=ISdiag(Weigth=Weigth, X=X, normal=FALSE, p=svs.pw2pn, svm=FALSE)

  is_time <- as.numeric(Sys.time() - is_time, units = "mins")
  cat('\n', "Done. Elapsed time:", is_time, "min\n")

  times <- as.numeric(Sys.time() - time, units = "mins")

  ### -log p(y|\theta)
  # DIC
  #theta_hat=apply(X, 2, mean)
  #D=2*svs.mllk(parvect=theta_hat,y=y,m=m,gmax=gmax)

  # Evaluating \bar{D(\theta)}
  #Dbar=0
  #for(j in 1:nIS){
  #  pv=X[j,]
  #  Dbar=Dbar+Weigth[j]*svs.mllk(parvect=pv,y=y,m=m,gmax=gmax)
  #}
  #Dbar=2*Dbar
  #pd=Dbar-D
  #DIC=D+2*pd

  # Log Predictive Score
  #LPS=0.5*D/length(y)

  structure(
    list(
      Weigth=Weigth,
      X=X,
      y=y,
      m=m,
      gmax=gmax,
      times=times,
      call=match.call()
    ),
    class = "svsHMM"
  )
}

#' @export
summary.svsHMM <- function(object, ...) {
  ISdiag(
    Weigth=object$Weigth,
    X=object$X,
    normal=FALSE,
    p=svs.pw2pn,
    svm=FALSE
  )
}

#' @export
info_crit <- function(object, ...) {
  UseMethod("info_crit")
}

#' @export
info_crit.svsHMM <- function(object, ...) {

  theta_hat <- colMeans(object$X)

  D <- 2 * svs.mllk(
    parvect = theta_hat,
    y = object$y,
    m = object$m,
    gmax = object$gmax
  )

  Dbar <- 0
  for (j in seq_along(object$Weigth)) {
    Dbar <- Dbar +
      object$Weigth[j] *
      svs.mllk(object$X[j, ], object$y, object$m, object$gmax)
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
logvol.svsHMM <- function(object, plot=FALSE, ...) {

  X      <- object$X
  Weigth <- object$Weigth

  # Transform parameters from working space to natural space
  Thetas <- t(apply(X = X, MARGIN = 1, FUN = svs.pw2pn))

  # Weighted posterior mean
  wThetas <- apply(X = Thetas, MARGIN = 2, FUN = "*", Weigth)
  theta_hat <- apply(X = wThetas, MARGIN = 2, sum)

  lv=sv.viterbi(
    y = object$y,
    theta_hat = theta_hat,
    m = object$m,
    gmax = object$gmax,
    model = "fillallprobs_s2"
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
