svmn.pn2pw = function(beta,mu,phi,sigma){
  lbeta1<- beta[1]
  lbeta2<-atanh(beta[2]) #log((1+beta[2])/(1-beta[2]))
  lbeta3<-beta[3]
  lmu<-mu
  lphi <- atanh(phi) #log((1+phi)/(1-phi))
  lsigma <- log(sigma)
  parvect <- c(lbeta1,lbeta2,lbeta3,lmu,lphi,lsigma)
  return(parvect)
}
svmn.pw2pn = function(parvect){
  beta=array(0,dim=3)
  beta[1]= parvect[1]
  beta[2]= tanh(parvect[2]) #(exp(parvect[2])-1)/(exp(parvect[2])+1)
  beta[3]=parvect[3]
  mu=parvect[4]
  phi = tanh(parvect[5]) #(exp(parvect[5])-1)/(exp(parvect[5])+1)
  sigma = exp(parvect[6])
  return(c(beta, mu, phi, sigma))
}
## function that will be used to compute 'allprobs' in mllk below
fillallprobs_n = function(x,beg,beta,y){
  z = (x-beta[1]-beta[2]*y-beta[3]*beg^2)/beg
  #w = pdf_vg(y=z, mu=0, sigma=1, nu=nu)
  w = pdf_n(y=z)
  return(w/beg)
}
svmn.mllk =function(parvect,y,y0,m,gmax){
  ny <- length(y)
  p <- svmn.pw2pn(parvect)
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
  allprobs <- outer(xx,sey,"fillallprobs_n",beta=p[1:3],yy)
  delta <-dnorm(bs,p[4],p[6]/sqrt(1-p[5]^2))*intlen
  foo <- delta*allprobs[1,]
  lscale = mlogLk_Rcpp(allprobs,Gamma,foo,ny) #Rcpp function
  return(-lscale)
}
svmn.prior =function(parvect){

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

  #+ log(dnorm(parvect[2], 0.5, 10))
  #+ log(dnorm(parvect[3], 0, 10))
  #+ log(dnorm(parvect[4], 0, 10))
  #+ log(dnorm(parvect[5], 4.5, 10))
  #+ log(dnorm(parvect[6], -1.5, 10))
  #+ log(dnorm(parvect[7], -10, 10))
  return(-lprior)
}
svmn.posterior =function(parvect,y,y0,m,gmax){
  return(svmn.mllk(parvect,y,y0,m,gmax)+svmn.prior(parvect))
}
svmn.map = function(y, m, parvect0, y0, gmax){

  #parvect <- svm.pn2pw(beta=beta0,mu=mu0,phi=phi0,sigma=sigma0,nu=nu0)
  mod = optim(parvect0,svmn.posterior,y=y,y0=y0,m=m,gmax=gmax,hessian=TRUE)
  mode = mod$par
  conv = mod$convergence
  return(list(mode=mode,
              lpost=-mod$value,
              hessian=mod$hessian,
              conv = conv))
}

#' @export
svmn.sim =function(y0,beta,mu,phi,sigma,g_dim){
  y=array(0,dim=g_dim)
  h=array(0,dim=g_dim)
  h[1]=rnorm(1,mu,sigma/sqrt(1-phi^2))
  y[1]=beta[1]+beta[2]*y0+beta[3]*exp(h[1])+exp(h[1]/2)*rnorm(1,0,1)
  for (j in 2:g_dim){
    h[j]=mu+phi*(h[j-1]-mu)+rnorm(1,0,sigma)
    y[j]=beta[1]+beta[2]*y[j-1]+beta[3]*exp(h[j])+exp(h[j]/2)*rnorm(1,0,1)
  }
  return (list(y=y,h=h,y0=y0))
}

#' @export
svmnHMM = function(y, y0=0.2, m, gmax, theta_init=NULL, nIS=1e3){

  time_start = Sys.time()

  if(is.null(theta_init)){
    mu0     = log(var(y))
    phi0    = 0.98
    sigma0  = 0.15
    beta0   = c(mean(y), 0, -0.05)
  }else{
    # Extracting initial parameters
    mu0     = theta_init$mu0
    phi0    = theta_init$phi0
    sigma0  = theta_init$sigma0
    beta0   = theta_init$beta0
  }

  # ----------------------------------------------------------------------------
  # STEP 1: Numerical Optimization (MAP Estimation)
  # ----------------------------------------------------------------------------
  eigenvalues = rep(0, 6)
  cont = 0
  parvect0 = svmn.pn2pw(beta0, mu0, phi0, sigma0)

  cat("\n[1/2] Starting numerical optimization (MAP)...\n")
  optim_time_start <- Sys.time()

  while( any(eigenvalues <= 0) && cont < 5 ){
    optim_res = svmn.map(y=y, m=m, parvect0=parvect0, y0=y0, gmax=gmax)

    # Inverting the Hessian matrix
    H1 = signif(solve(optim_res$hessian), 6)
    eigenvalues = eigen(H1, only.values=TRUE)$values
    parvect0 = optim_res$mode
    cont = cont + 1
  }

  if(any(eigenvalues <= 0)){
    cat("\n")
    stop('Error: Hessian matrix is not positive definite after 5 attempts. Check initial values.')
  }

  optim_duration <- as.numeric(Sys.time() - optim_time_start, units = "mins")
  cat(sprintf("      Optimization completed in %.2f min (%d iterations).\n", optim_duration, cont))

  # ----------------------------------------------------------------------------
  # STEP 2: Importance Sampling
  # ----------------------------------------------------------------------------
  cat(sprintf("\n[2/2] Starting Importance Sampling (n = %d)...\n", nIS))
  is_time_start <- Sys.time()

  k   = -optim_res$lpost
  map = optim_res$mode

  # Proposal Sampling (Multivariate Normal)
  X = mvtnorm::rmvnorm(nIS, map, H1)
  Weight = numeric(nIS)

  largura = 40
  ultimo_print = -1

  for(j in 1:nIS){
    # Weight evaluation
    Weight[j] = exp(k -
                    svmn.posterior(parvect=X[j,], y=y, y0=y0, m=m, gmax=gmax) -
                    mvtnorm::dmvnorm(X[j,], mean=map, sigma=H1, log=TRUE))

    # Progress bar logic
    progresso = j / nIS
    pct = floor(progresso * 100)
    if (pct %% 5 == 0 && pct != ultimo_print || j == nIS){
      preenchido = round(largura * progresso)
      barra = paste0("[", paste(rep("=", preenchido), collapse = ""),
                     paste(rep(" ", largura - preenchido), collapse = ""), "]")
      cat(sprintf("\r      Progress: %s %3.0f%%", barra, pct))
      ultimo_print = pct
    }
  }

  # Normalize weights
  s = sum(Weight)
  if((s != 0) && !is.nan(s)){
    Weight = Weight / s
  } else {
    cat("\n")
    stop('Critical Error: Failed to normalize weights (Sum is zero or NaN).')
  }

  # Resampling
  indx = sample(1:nIS, prob=Weight, replace=TRUE)
  X = X[indx, ]
  Weight = rep(1/nIS, nIS)

  is_duration <- as.numeric(Sys.time() - is_time_start, units = "mins")
  cat(sprintf("\n      Sampling completed in %.2f min.\n", is_duration))

  # ----------------------------------------------------------------------------
  # Finalization
  # ----------------------------------------------------------------------------
  total_time <- as.numeric(Sys.time() - time_start, units = "mins")
  cat(sprintf("\nProcess completed successfully! Total elapsed time: %.2f min.\n\n", total_time))

  structure(
    list(
      Weigth=Weight,
      X=X,
      y=y,
      y0=y0,
      m=m,
      gmax=gmax,
      times=total_time,
      call=match.call()
    ),
    class = "svmnHMM"
    )


}

#' @export
summary.svmnHMM <- function(object, ...) {
  ISdiag(
    Weigth = object$Weigth,
    X = object$X,
    normal=TRUE,
    p = svmn.pw2pn
  )
}

#' @export
info_crit <- function(object, ...) {
  UseMethod("info_crit")
}

#' @export
info_crit.svmnHMM <- function(object, ...) {

  theta_hat <- colMeans(object$X)

  D <- 2 * svmn.mllk(
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
      svmn.mllk(object$X[j, ], object$y, object$y0, object$m, object$gmax)
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
logvol.svmnHMM <- function(object, plot=FALSE, ...) {

  X      <- object$X
  Weigth <- object$Weigth

  # Transform parameters from working space to natural space
  Thetas <- t(apply(X = X, MARGIN = 1, FUN = svmn.pw2pn))

  # Weighted posterior mean
  wThetas <- apply(X = Thetas, MARGIN = 2, FUN = "*", Weigth)
  theta_hat <- apply(X = wThetas, MARGIN = 2, sum)

  lv=svm.viterbi(
    y = object$y,
    y0 = object$y0,
    theta_hat = theta_hat,
    m = object$m,
    gmax = object$gmax,
    svmn=TRUE,
    model = "fillallprobs_n"
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

lforward.svmnHMM=function(x, m, gbmax, p, y0){

  nx <- length(x)
  yy = c(y0, x[1:(nx-1)])
  #p<-svmn.pw2pn(pvect)
  lalpha <- matrix(NA,m,nx)
  K <- m+1
  gb <- seq(-gbmax,gbmax,length=K)
  g <- (gb[-1]+gb[-K])*0.5
  beg <-exp(g/2)
  gamma <- matrix(0,m,m)
  E=p[4]+p[5]*(g-p[4])
  intlen <- gb[2]-gb[1]
  for (i in 1:m){
    goo <- dnorm(g,E[i],p[6])*intlen
    gamma[i,] <- goo/sum(goo)
  }

  lscale=0
  delta <- dnorm(g,p[4],p[6]/sqrt(1-p[5]^2))*intlen

  for (i in 1:nx){

    if (i==1) foo <- delta*1/beg*dnorm((x[1]-p[1]-p[2]*yy[1]-p[3]*beg^2)/beg)

    if ( i>1)
      foo <- foo%*%gamma*1/beg*dnorm((x[i]-p[1]-p[2]*yy[i]-p[3]*beg^2)/beg)
    lscale <- lscale+log(sum(foo));
    foo <- foo/sum(foo)
    lalpha[,i] <- log(foo)+lscale
  }

  return(lalpha)

}

#' @export
state.forecast.SVMnHMM <- function(object, h=1) {
  # 1. Obter as probabilidades forward (alfa) até o tempo final T
  # Usamos a sua função lalpha já definida
  x=object$y
  m=object$m
  gbmax=object$gmax
  y0=object$y0
  nx <- length(x)
  p <- summary(object)[,1]

  lalpha <- lforward.svmnHMM(x=x, m=m, gbmax=gbmax, p=p, y0=y0)

  # 2. Reconstruir a matriz de transição (gamma) discretizada
  K <- m + 1
  gb <- seq(-gbmax, gbmax, length=K)
  g <- (gb[-1] + gb[-K]) * 0.5
  intlen <- gb[2] - gb[1]

  gamma <- matrix(0, m, m)
  E <- p[4] + p[5] * (g - p[4])
  for (i in 1:m) {
    goo <- dnorm(g, E[i], p[6]) * intlen
    gamma[i,] <- goo / sum(goo)
  }

  # 3. Preparar a matriz de saída
  statepreds <- matrix(NA, nrow = m, ncol = h)

  # 4. Probabilidade atual dos estados no tempo T (filtrada)
  # Para evitar sobrefluxo (underflow), normalizamos no log
  c <- max(lalpha[, nx])
  foo <- exp(lalpha[, nx] - c)
  foo <- foo / sum(foo) # Distribuição de probabilidade no tempo T

  # 5. Projeção para h passos à frente
  for (i in 1:h) {
    foo <- foo %*% gamma
    statepreds[, i] <- foo
  }

  # Retorna a distribuição de probabilidade dos bins para cada horizonte h
  # Você pode usar colSums(statepreds * g) para obter a volatilidade esperada
  return(list(probs = statepreds, grids = g))
}

svmn.HMM.quantiles=function(x, res.time, m, gbmax, p, lalpha, y0){

  nx = length(x)
  yy = c(y0, x[1:(nx-1)])
  K = m+1
  #p=svm.pw2pn(pvect)
  gb = seq(-gbmax,gbmax,length=K)
  g = (gb[-1]+gb[-K])*0.5
  beg = exp(g/2)
  gamma = matrix(0,m,m)
  E=p[4]+p[5]*(g-p[4])
  intlen = gb[2]-gb[1]
  for(i in 1:m){
    goo = dnorm(g,E[i],p[6])*intlen
    gamma[i,] = goo/sum(goo)
  }
  c=max(lalpha[,res.time-1])   # to avoid numerical underflow
  a=exp(lalpha[,res.time-1]-c) # scalar factor cancels in Pro computation below
  npsr=t(a)%*%(gamma/sum(a))%*%pnorm((x[res.time]-p[1]-p[2]*yy[res.time]-p[3]*beg^2)/beg)

  return(npsr)
}

#' @export
svmnHMM.pseudo_residuals=function(object, ytest){

  ytrain=object$y
  m=object$m
  gbmax=object$gmax
  y0=object$y0
  p=summary(object)[,1]

  lalpha=lforward.svmnHMM(x=c(ytrain, ytest), m=m, gbmax=gbmax, p=p, y0=y0)

  psressvn=rep(NA,length(ytest))
  for(k in 1:length(ytest)){
    # note that the function 'SV.HMM.quantiles' could alternatively be written
    # such that no loop would be required here, but I didn't bother since it's
    # fast enough as it stands
    psressvn[k]=svmn.HMM.quantiles(x=c(ytrain, ytest), res.time=length(ytrain)+k,
                                   m=m, gbmax=gbmax, p=p, lalpha=lalpha, y0=y0)
  }

  return(psressvn)
}
