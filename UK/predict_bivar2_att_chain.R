
#predict_bivar2 <- function(y, m0, C0, h, V, Ft, Gt, delta, ages, lower, upper){
  predict_bivar2_chain <- function(y, m0, C0, h, V, Ft, Gt, delta, ages){
  n = dim(V)[1]
  N = nrow(y)
  q = ncol(y)
  p = length(m0)
  if(length(delta) > 1){ delta_pred = delta[N] } else{ delta_pred = delta }
  sim <- array(NA, dim = c(n, h, q))
  
  for(i in 1:n){
    aux = ff(m0 = m0, C0 = C0, y = as.matrix(y),
             V = V[i,,], Ft = Ft, Gt = Gt, delta = delta)
    
    Wt = aux$C[N,,] * (1 - delta_pred) / delta_pred
    at = Gt %*% aux$m[N,] 
    Rt = Gt %*% aux$C[N,,] %*% t(Gt) + Wt
    ft = Ft %*% at
    Qt = Ft %*% Rt %*% t(Ft) + V[i,,]
    At = Rt %*% t(Ft) %*% solve(Qt)
    Ct = Rt - At %*% Ft %*% Rt   ## second moment
    
   
    sim[i,1,] <- MASS::mvrnorm(1, ft, Qt)
    
    #### verificar entradas de lower e upper
    
    if(h > 1) for(k in 2:h){
      Wt = Ct #* (1 - delta_pred) / delta_pred ## Confirmar essa linha
      at = Gt %*% at 
      Rt = Gt %*% Rt %*% t(Gt) + Wt
      ft = Ft %*% at
      Qt = Ft %*% Rt %*% t(Ft) + V[i]
      At = Rt %*% t(Ft) %*% solve(Qt)
      Ct = Rt - At %*% Ft %*% Rt   ## second moment
      
      sim[i,k,] <- MASS::mvrnorm(1, ft, Qt)
      
    }
  }
  
  qx_sim = exp(sim)
  qx_sim[qx_sim < 0] = 0
  #qx_sim[qx_sim > 1] = 1
  #qx_fitted = apply(qx_sim, 2:3, median, na.rm = T)
  #qx_lim = apply(qx_sim, 2:3, quantile, probs = c(0.025, 0.975), na.rm = T)
  #qx_fitted = data.frame(Ages = (ages[N]+1):(ages[N]+h), qx_fitted = qx_fitted)
  #return(list(qx_pred = qx_fitted, qx_int_pred = qx_lim))
  #return(list(qxpred.1 = cbind(qx_fitted[,1:2], t(qx_lim[,,1])),
  #qxpred.2 = cbind(qx_fitted[,c(1,3)], t(qx_lim[,,2]))))
  return(qx_sim)
}


