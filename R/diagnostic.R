quant = function(x, weights, probs=c(0.025, 0.5, 0.975)){

  ix = order(x)
  x_sorted = x[ix]
  weights_sorted = weights[ix]
  weights_sorted = weights_sorted / sum(weights_sorted)
  cumulative_weights = cumsum(weights_sorted)
  quants = numeric(length(probs))
  for(i in seq_along(probs)){
    index = which(cumulative_weights >= probs[i])[1]
    quants[i] = x_sorted[index]
  }
  return(quants)
}
ISdiag = function(Weigth, X, normal=FALSE, svm=TRUE, p){

  #p=svmn.pw2pn
  ### mean #####################################################################
  Thetas=t(apply(X=X, MARGIN=1, FUN=p))
  wThetas = apply(X=Thetas, MARGIN = 2, FUN='*', Weigth)
  theta_hat = apply(X=wThetas, MARGIN = 2, sum)
  ### var ######################################################################
  Vars = t(apply(X=Thetas, 1, '-', theta_hat))
  wVars = apply(X=Vars, MARGIN = 2, FUN='*', Weigth)
  wVars = wVars^2
  Vars_hat = apply(X=wVars, MARGIN = 2, sum)
  ### quantiles ################################################################
  Quants = t(apply(Thetas, 2, quant, Weigth))
  Results = cbind(theta_hat, Vars_hat, Quants)
  colnames(Results) = c('mean', 'var', '2.5%', '50%', '97.5%')

  if(normal){
    if(svm){
      row.names(Results) = c('b0','b1','b2','mu','phi','sigma')
    }else{
      row.names(Results) = c('mu','phi','sigma')
    }
  }else{
    if(svm){
      row.names(Results) = c('b0','b1','b2','mu','phi','sigma','nu')
    }else{
      row.names(Results) = c('mu','phi','sigma','nu')
    }
  }
  Results = round(Results, 4)#signif(Results, 3)

  return(Results=Results)
}
