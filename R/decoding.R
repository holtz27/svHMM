#' @export
sv.viterbi = function(y, theta_hat, m, gmax, svn=FALSE, model){
  n = length(y)
  p = theta_hat
  K = m+1
  b = seq(-gmax,gmax,length=K)
  bs = (b[-1]+b[-K])*0.5
  E = p[1] + p[2]*(bs-p[1])
  intlen = b[2]-b[1]
  sey = exp(bs/2)
  Gamma = matrix(0,m,m)
  for(i in 1:m){
    Gamma[i,] = dnorm(bs,E[i],p[3])*intlen
  }
  Gamma = Gamma/apply(Gamma,1,sum)
  #xx = y
  #yy = c(y0,y[1:(n-1)])
  if(svn){
    allprobs = outer(y, sey, model)
  }else{
    allprobs = outer(y, sey, model, nu=p[4])
  }
  delta = dnorm(bs, p[1], p[3]/sqrt(1-p[2]^2))*intlen

  # defining xi
  xi = matrix(0, n, m)
  foo = delta*allprobs[1,]
  xi[1,] = foo/sum(foo)
  for(i in 2:n){
    foo = apply(xi[i-1,]*Gamma, 2, max) * allprobs[i,]
    xi[i,] = foo/sum(foo)
  }

  # evaluating i
  iv = numeric(n)
  iv[n] = which.max(xi[n,])
  for(i in (n-1):1){
    iv[i] = which.max(Gamma[,iv[i+1]]*xi[i,])
  }

  return(bs[iv])
}


#' @export
svm.viterbi = function(y, y0, theta_hat, m, gmax, svmn=FALSE, model){
  n = length(y)
  p = theta_hat
  K = m+1
  b = seq(-gmax,gmax,length=K)
  bs = (b[-1]+b[-K])*0.5
  E = p[4] + p[5]*(bs-p[4])
  intlen = b[2]-b[1]
  sey = exp(bs/2)
  Gamma = matrix(0,m,m)
  for(i in 1:m){
    Gamma[i,] = dnorm(bs,E[i],p[6])*intlen
  }
  Gamma = Gamma/apply(Gamma,1,sum)
  xx = y
  yy = c(y0,y[1:(n-1)])
  if(svmn){
    allprobs = outer(xx, sey, model, beta=p[1:3], yy)
  }else{
    allprobs = outer(xx, sey, model, beta=p[1:3], nu=p[7], yy)
  }
  delta = dnorm(bs, p[4], p[6]/sqrt(1-p[5]^2))*intlen

  # defining xi
  xi = matrix(0, n, m)
  foo = delta*allprobs[1,]
  xi[1,] = foo/sum(foo)
  for(i in 2:n){
    foo = apply(xi[i-1,]*Gamma, 2, max) * allprobs[i,]
    xi[i,] = foo/sum(foo)
  }

  # evaluating i
  iv = numeric(n)
  iv[n] = which.max(xi[n,])
  for(i in (n-1):1){
    iv[i] = which.max(Gamma[,iv[i+1]]*xi[i,])
  }

  return(bs[iv])
}
