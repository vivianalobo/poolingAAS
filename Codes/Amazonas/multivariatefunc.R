#######################################################################
###### funcoes auxiliares para rodar o modelo para dados da amazonas
#######################################################################


#' Filtragem com V fixo e fator de desconto para W
#'
#' @param m0 vetor de dimensao p com os valores iniciais de mt da filtragem
#' @param C0 matriz de dimensao pxp com os valores iniciais de Ct da filtragem
#' @param y vetor de tamanho N
#' @param V variância observacional - escalar
#' @param Ft0 matriz de dimensao Nxp
#' @param Gt matriz de dimensao pxp
#' @param delta fator de desconto - escalar
#'
#' @return Uma lista contendo mt, Ct, at, Rt e Wt
#'
#' @examples
#' n = length(Nile)
#' mld = ff(m0 = 0, C0 = 100000, y = Nile, V = 15099.8,
#'          Ft0 = matrix(rep(1,n),n,1), Gt = 1, delta = 0.85)
#' xg = 1:length(Nile)
#' q = qnorm(1 - 0.05/2)
#' qfinf = mld$m - q*sqrt(as.vector(mld$C))
#' qfsup = mld$m + q*sqrt(as.vector(mld$C))
#' ts.plot(as.vector(Nile))
#' polygon(x = c(xg, rev(xg)), y = c(qfinf, rev(qfsup)),
#'         col = adjustcolor('lightcoral', alpha.f = 0.3), border = NA)
#' lines(as.vector(mld$m), col = 2, lwd = 2)
#'
# ff = function(m0, C0, y, V, Ft, Gt, delta){
# 
#   # if (!is.matrix(C0)) stop("C0 deve ser uma matriz")
#   # if (dim(C0)[1] != dim(C0)[2]) stop("C0 deve ser uma matriz quadrada")
#   # if (!is.vector(m0)) stop("m0 deve ser um vetor")
#   # if (length(m0) != dim(C0)[1]) stop("Dimensões de m0 e C0 incompativeis")
#   # if (!is.vector(y)) stop("y deve ser um vetor")
#   # if (!is.matrix(Ft0)) stop("Ft0 deve ser uma matriz")
#   # if (dim(Ft0)[1] != length(y)) stop("Dimensões de Ft0 e y imcompativeis")
#   # if (!is.matrix(Gt)) stop("Gt deve ser uma matriz")
#   # if (delta < 0 | delta > 1) stop("delta deve estar entre 0 e 1")
# 
#   N = nrow(y)
#   p = length(m0)
#   resultado.m = matrix(NA, N, p)
#   resultado.C = array(NA,c(N,p,p))
#   resultado.W = array(NA,c(N,p,p))
#   resultado.a = matrix(NA, N, p)
#   resultado.R = array(NA,c(N,p,p))
#   if(length(delta) == 1){ delta = rep(delta, N) }
# 
#   ## Filtro de Kalman
#   ### passo 1 (inicializacao)
# 
#   Wt = C0 * (1 - delta[1]) / delta[1]
# 
#   #  if(is.matrix(Ft0) == TRUE){ Ft = Ft0[1,] }
#   at = Gt %*% m0
#   Rt = Gt %*% C0 %*% t(Gt) + Wt
#   ft = Ft %*% at
#   Qt = Ft %*% Rt %*% t(Ft) + V
#   et = y[1,] - ft
#   At = Rt %*% t(Ft) %*% solve(Qt)
#   mt = at + At %*% et   ### first moment
#   Ct = Rt - At %*% Ft %*% Rt   ## second moment
# 
#   resultado.m[1,] = mt
#   resultado.C[1,,] = Ct
#   resultado.W[1,,] = Ct * (1 - delta[1]) / delta[1]
#   resultado.a[1,] = at
#   resultado.R[1,,] = Rt
# 
#   ### passo 2 (atualizacao)
#   for (j in 2:N) {
#     Wt = Ct * (1 - delta[j]) / delta[j]
#     # if(is.matrix(Ft0) == TRUE){ Ft = Ft0[j,] }
#     at = Gt %*% mt
#     Rt = Gt %*% Ct %*% t(Gt) + Wt
#     ft = Ft %*% at
#     Qt = Ft %*% Rt %*% t(Ft) + V
#     et = y[j,] - ft
#     At = Rt %*% t(Ft) %*% solve(Qt)
#     mt = at + At %*% et  ### mean
#     Ct = Rt - At %*% Ft %*% Rt ### variance
# 
#     resultado.m[j,] = mt
#     resultado.C[j,,] = Ct
#     resultado.W[j,,] = Wt
#     resultado.a[j,] = at
#     resultado.R[j,,] = Rt
#   }
# 
#   return(list(m = resultado.m, C = resultado.C, a = resultado.a, R = resultado.R, W = resultado.W))
# }

# ff = function(m0, C0, y, V, Ft, Gt, delta){
#   
#   N = nrow(y)
#   p = length(m0)
#   J = (p - 1) / 2  # número de populações
#   
#   resultado.m = matrix(NA, N, p)
#   resultado.C = array(NA, c(N, p, p))
#   resultado.W = array(NA, c(N, p, p))
#   resultado.a = matrix(NA, N, p)
#   resultado.R = array(NA, c(N, p, p))
#   
#   # Se delta vier como escalar - replica
#   if(is.null(dim(delta))){
#     delta = matrix(delta, nrow = N, ncol = J + 1)
#   }
#   
#   make_Wt <- function(Ct, delta_row){
#     
#     deltas_j <- delta_row[1:J]      # população
#     delta_g  <- delta_row[J + 1]    # termo comum
#     
#     vec_delta <- numeric(p)  # vetor completo respeitando ordem do estado
#     vec_delta[1:J] <- (1 - deltas_j) / deltas_j ### mu
#     vec_delta[(J+1):(2*J)] <- (1 - deltas_j) / deltas_j ### beta
#     
#     vec_delta[2*J + 1] <- (1 - delta_g) / delta_g
#     D <- diag(sqrt(vec_delta))
#     Wt <- D %*% Ct %*% D
#     return(Wt)
#   }
#   
#   Wt = make_Wt(C0, delta[1,])
#   
#   at = Gt %*% m0
#   Rt = Gt %*% C0 %*% t(Gt) + Wt
#   
#   ft = Ft %*% at
#   Qt = Ft %*% Rt %*% t(Ft) + V
#   et = y[1,] - ft
#   
#   At = Rt %*% t(Ft) %*% solve(Qt)
#   mt = at + At %*% et
#   Ct = Rt - At %*% Ft %*% Rt
#   
#   resultado.m[1,] = mt
#   resultado.C[1,,] = Ct
#   resultado.W[1,,] = Wt
#   resultado.a[1,] = at
#   resultado.R[1,,] = Rt
#   
#   for (t in 2:N) {
#     
#     Wt = make_Wt(Ct, delta[t,])
#     
#     at = Gt %*% mt
#     Rt = Gt %*% Ct %*% t(Gt) + Wt
#     
#     ft = Ft %*% at
#     Qt = Ft %*% Rt %*% t(Ft) + V
#     et = y[t,] - ft
#     
#     At = Rt %*% t(Ft) %*% solve(Qt)
#     mt = at + At %*% et
#     Ct = Rt - At %*% Ft %*% Rt
#     
#     resultado.m[t,] = mt
#     resultado.C[t,,] = Ct
#     resultado.W[t,,] = Wt
#     resultado.a[t,] = at
#     resultado.R[t,,] = Rt
#   }
#   
#   return(list(m = resultado.m, C = resultado.C, a = resultado.a, R = resultado.R, W = resultado.W))
#   
# }

ff = function(m0, C0, y, V, Ft, Gt, delta){
  
  N = nrow(y)
  p = length(m0)
  
  # Detecta se existe termo comum
  common = (p %% 2 != 0)
  
  if(common){
    J = (p - 1) / 2
  } else {
    J = p / 2
  }
  
  resultado.m = matrix(NA, N, p)
  resultado.C = array(NA, c(N, p, p))
  resultado.W = array(NA, c(N, p, p))
  resultado.a = matrix(NA, N, p)
  resultado.R = array(NA, c(N, p, p))
  
  # Ajusta dimensão de delta
  if(is.null(dim(delta))){
    
    if(common){
      delta = matrix(delta, nrow = N, ncol = J + 1)
    } else {
      delta = matrix(delta, nrow = N, ncol = J)
    }
    
  }
  
  make_Wt <- function(Ct, delta_row){
    
    vec_delta <- numeric(p)
    
    deltas_j <- delta_row[1:J]
    
    # mu
    vec_delta[1:J] <- (1 - deltas_j) / deltas_j
    
    # beta
    vec_delta[(J+1):(2*J)] <- (1 - deltas_j) / deltas_j
    
    # termo comum (se existir)
    if(common){
      
      delta_g <- delta_row[J + 1]
      
      vec_delta[2*J + 1] <-
        (1 - delta_g) / delta_g
    }
    
    D <- diag(sqrt(vec_delta))
    
    Wt <- D %*% Ct %*% D
    
    return(Wt)
  }
  
  Wt = make_Wt(C0, delta[1,])
  
  at = Gt %*% m0
  Rt = Gt %*% C0 %*% t(Gt) + Wt
  
  ft = Ft %*% at
  Qt = Ft %*% Rt %*% t(Ft) + V
  et = y[1,] - ft
  
  At = Rt %*% t(Ft) %*% solve(Qt)
  
  mt = at + At %*% et
  Ct = Rt - At %*% Ft %*% Rt
  
  resultado.m[1,] = mt
  resultado.C[1,,] = Ct
  resultado.W[1,,] = Wt
  resultado.a[1,] = at
  resultado.R[1,,] = Rt
  
  for (t in 2:N) {
    
    Wt = make_Wt(Ct, delta[t,])
    
    at = Gt %*% mt
    Rt = Gt %*% Ct %*% t(Gt) + Wt
    
    ft = Ft %*% at
    Qt = Ft %*% Rt %*% t(Ft) + V
    
    et = y[t,] - ft
    
    At = Rt %*% t(Ft) %*% solve(Qt)
    
    mt = at + At %*% et
    Ct = Rt - At %*% Ft %*% Rt
    
    resultado.m[t,] = mt
    resultado.C[t,,] = Ct
    resultado.W[t,,] = Wt
    resultado.a[t,] = at
    resultado.R[t,,] = Rt
  }
  
  return(list(
    m = resultado.m,
    C = resultado.C,
    a = resultado.a,
    R = resultado.R,
    W = resultado.W
  ))
}




#' Backward Sampling
#'
#' Usada internamente na funcao ffbs()
#'
#' @param m vetor de dimensao p com o resultado da filtragem
#' @param C matriz de dimensao pxp com o resultado da filtragem
#' @param a vetor de dimensao p com o resultado da filtragem
#' @param R matriz de dimensao pxp com o resultado da filtragem
#' @param Gt matriz de dimensao pxp
#'
#' @return Uma lista contendo theta, mt, Ct, as e Rs
#'
#' @importFrom MASS mvrnorm
#'
bs = function(m,C,a,R,Gt){
  
  N = nrow(m)
  p = ncol(m)
  
  as = matrix(NA, N, p)
  Rs = array(NA,c(N,p,p))
  theta <- matrix(NA,N,p)
  
  as[N,] = m[N,]
  Rs[N,,] = C[N,,]
  
  ### draw theta_T - page 162 petris petroni
  theta[N,] <- MASS::mvrnorm(1, as[N,], Rs[N,,])
  ### step 3 - algorithm 4.1 Backward Sampling
  for (t in (N - 1):1) {
    
    Bt = C[t,,] %*% t(Gt) %*% solve(R[t + 1,,])
    
    # Rs[t,,] = C[t,,] + Bt %*% (Rs[t + 1,,] - R[t + 1,,]) %*% t(Bt)
    
    # as[t,] = m[t,] + Bt %*% (as[t + 1,] - a[t + 1,])
    
    ht <- m[t,] + Bt %*% (theta[t + 1,] - a[t + 1,])
    Ht = C[t,,] - Bt %*% R[t + 1,,] %*% t(Bt)
    ### draw theta_t
    theta[t,]  = MASS::mvrnorm(1,ht, Ht)
    #### theta[t,]  = MASS::mvrnorm(1,as[t,], Rs[t,,])
  }
  return(list(d = theta, m = m, C = C))
  #m.m = as, C.C = Rs))
  
}

#' Filtragem e suavização com FFBS e fator de desconto para W
#'
#' Usada internamente na funcao gibbsSigma2
#'
#' @param m0 vetor de dimensao p com os valores iniciais de mt da filtragem
#' @param C0 matriz de dimensao pxp com os valores iniciais de Ct da filtragem
#' @param y vetor de tamanho N
#' @param V variância observacional - escalar
#' @param Ft0 matriz de dimensao Nxp
#' @param Gt matriz de dimensao pxp
#' @param delta fator de desconto - escalar
#'
#' @return Uma lista contendo theta, mt, Ct, as, Rs e W
#'
ffbs <- function(m0, C0, y, V, Ft, Gt, delta){
  
  aux.f = ff(m0, C0, y, V, Ft, Gt, delta)
  
  res = bs(aux.f$m,aux.f$C,aux.f$a,aux.f$R,Gt)
  res$W = aux.f$W
  
  return(res)
  
}


#' Estimação de V constante via Gibss
#'
#' Utlizada ffbs e um passo de Gibbs para estimar uma variancia constante no tmepo.
#'
#' @param m0 vetor de dimensao p com os valores iniciais de mt da filtragem
#' @param C0 matriz de dimensao pxp com os valores iniciais de Ct da filtragem
#' @param y vetor de tamanho N
#' @param Ft0 matriz de dimensao Nxp
#' @param Gt matriz de dimensao pxp
#' @param delta fator de desconto - escalar
#' @param V valor inicial para V
#' @param v0 priori inversa Wishart para matriz de covariancia (ou precisao, confirmar)
#' @param s0 priori inversa Wishart para matriz de  covariancia (ou precisao, confirmar)
#' @param nit numero de iteracoes
#' @param shiny habilita as funcoes para exibir progresso no shiny e cancelar o calculo
#' @param status_file só é usado quando shiny = T.
#'
#' @return Uma lista contendo uma cadeia da posteriori de mu, beta e sigma2.
#'
#' @examples
#'
#' data_hp$qx <- data_hp$dx/data_hp$nx
#' data_hp$qx <- 1 - exp(-data_hp$qx)
#' y <- log(data_hp$qx)
#' plot(y, t='l')
#' n <- length(y)
#'
#' Gt <- matrix(c(1,1,0,1),  ncol=2, byrow=TRUE)
#' Ft0 <- matrix(c(1,0), n , ncol=2, byrow=TRUE)
#' m0 <-  rep(0, 2)
#' C0 <- diag(100,2)
#' delta <- 0.65
#' V <- 0.002674713
#'
#' res <- gibbsSigma2(m0,C0,y,Ft0,Gt,delta,0.002674713,0.01,0.01,500)
#'
#' plot(res$sig2, t='l')
#'
#' mu <- apply(res$mu, 2, median)
#' Ft <- t(c(1,0))
#' media <- mu
#' med = exp(media)[2:76]
#' plot(data_hp$x, med, type = "l", log = "y", ylim=c(8e-05, 1e-01))
#' points(data_hp$x,data_hp$qx,pch = 19, cex = 0.5, col = "red")
#'
#' @import progress
#'
dlm.multivariate <- function(m0, C0, y, Ft0, Gt, delta, V, v0, s0, nit,
                             bn,
                             bn.sample = NULL, thin,shiny = F, status_file = NULL){
  n <- nrow(y)
  q <- ncol(y)
  p = length(m0)
  Ft <- Ft0
  V.post <- array(NA, dim = c(nit, q, q))
  mu.post <- array(NA, dim = c(nit, n, q))
  # theta.post <- matrix(NA, nit, ncol = n)
  theta.post <- array(dim = c(nit, n, ncol(Ft)))
  Wt = array(NA,dim = c(nit, n, p, p))
  
  pb  = progress::progress_bar$new(format = "Simulating [:bar] :percent in :elapsed",
                                   total = nit, clear = FALSE, width = 60)
  
  ## first parameter from posterior of V
  alpha.star <- (v0 + 1 + n)/2
  V0 <- (v0-2)*s0
  
  for (k in 1:nit) {
    
    pb$tick()
    
    ## FFBS for thetas (each age x)
    mld <- ffbs(m0 , C0 , y, V = V, Ft0 , Gt , delta)
    dt <- mld$d
    Wt[k,,,] = mld$W
    
    mu.post[k,,] <- t(Ft%*%t(dt)) ;   theta.post[k,,] <- dt ;
    muk <- mu.post[k,,]
    
    
    ####### gibbs for V
    SSy = t(y-muk)%*%(y-muk)
    beta.star <- 0.5*(SSy + V0)
    V.post[k,,] <- solve(rWishart(1, alpha.star, solve(beta.star))[,,1])
    V <- V.post[k,,]
  }
 # return(list(mu = mu.post, theta = theta.post, V = V.post, Wt = Wt))
   return(list(mu = mu.post[seq(bn+1, nit, by = thin),,],
              theta = theta.post[seq(bn+1, nit, by = thin),,],
              V = V.post[seq(bn+1, nit, by = thin),,], 
              Wt = Wt[seq(bn+1, nit, by = thin),,,]))
  
}



gibbs_missing_data <- function(m0, C0, y, Ft0, Gt, delta, V, v0, s0, nit, lower = NULL, upper = NULL, bn,
                               bn.sample = NULL, thin, y_missing = NULL,
                               ...){
  n <- nrow(y)
  q <- ncol(y)
  p = length(m0)
  ##### aqui mostra em qual linha (age) e coluna (sex) temos missing
  ind_missing = which(is.na(y), arr.ind = T) 
  
  Ft <- Ft0
  V.post <- array(NA, dim = c(nit, q, q))
  mu.post <- array(NA, dim = c(nit, n, q))
  # theta.post <- matrix(NA, nit, ncol = n)
  theta.post <- array(dim = c(nit, n, ncol(Ft)))
  Wt = array(NA, dim = c(nit, n, p, p))
  
  ### parte do missing data
  input.post = matrix(NA, nrow = nit, ncol = nrow(ind_missing))
  
  ## chutes iniciais
  if(length(ind_missing) > 0 & is.null(y_missing)){
    y_missing = rowMeans(y[ind_missing[,1],], na.rm = T)
  }
  y[ind_missing] = y_missing
  t_missing = unique(ind_missing[,1])
  
  pb  = progress::progress_bar$new(format = "Simulating [:bar] :percent in :elapsed",
                                   total = nit, clear = FALSE, width = 60)
  
  ## first parameter from posterior of V
  alpha.star <- (v0 + 1 + n)/2
  V0 <- (v0-2)*s0
  
  for (k in 1:nit) {
    
    pb$tick()
    
    # FFBS for thetas (each age x)
    mld <- ffbs(m0, C0, y, V = V, Ft0, Gt, delta, ...)
    dt <- mld$d
    Wt[k,,,] = mld$W
    
    mu.post[k,,] <- t(Ft%*%t(dt)) ;   theta.post[k,,] <- dt ;
    muk <- mu.post[k,,]
    
    ####### gibbs for V (pagina 174 do petris distribuicao a posteriori)
    SSy = t(y-muk)%*%(y-muk)
    beta.star <- 0.5*(SSy + V0)  ### S0= 0.5*V0     
    V.post[k,,] <- solve(rWishart(1, alpha.star, solve(beta.star))[,,1])
    V <- V.post[k,,]
    
    ## condicional completa y
    input = c()
    ## indice com 1 se refere ao dado faltante e indice 2 ao dado observado
    for(i in 1:length(t_missing)){
      
      # índices
      cols_miss = ind_missing[ind_missing[,1] == t_missing[i], 2]
      cols_obs  = setdiff(1:q, cols_miss)
      # dados observados (força vetor coluna)
      aux2 = matrix(y[t_missing[i], cols_obs], ncol = 1)
      
      # médias
      mu1 = matrix(muk[t_missing[i], cols_miss], ncol = 1)
      mu2 = matrix(muk[t_missing[i], cols_obs], ncol = 1)
      
      # covariâncias (FORÇANDO matriz sempre)
      V11 = V[cols_miss, cols_miss, drop = FALSE]
      V12 = V[cols_miss, cols_obs, drop = FALSE]
      V22 = V[cols_obs, cols_obs, drop = FALSE]
      inv_V22 = solve(V22)
      V21 = t(V12)
      
      ## Parametros da condicional
      mu_condicional = mu1 + V12 %*% inv_V22 %*% (aux2 - mu2)
      sigma_condicional = V11 - V12 %*% inv_V22 %*% V21
      
      y[t_missing[i], cols_miss] =
        MASS::mvrnorm(1, mu = as.vector(mu_condicional),Sigma = sigma_condicional)
      input = append(input, y[t_missing[i], cols_miss])
    }
    input.post[k,] = input
    
    # for(i in 1:length(t_missing)){
    #   ## valores dos y observados da i-esima idade com dados faltantes
    #   aux2 = y[t_missing[i], -ind_missing[ind_missing[,1] == t_missing[i], 2]]
    #   ## Medias
    #   mu1 = muk[t_missing[i], ind_missing[ind_missing[,1] == t_missing[i], 2]]
    #   mu2 = muk[t_missing[i], -ind_missing[ind_missing[,1] == t_missing[i], 2]]
    #   ## Variancias
    #   V11 = V[ind_missing[ind_missing[,1] == t_missing[i], 2], ind_missing[ind_missing[,1] == t_missing[i], 2]]
    #   V12 = V[ind_missing[ind_missing[,1] == t_missing[i], 2], -ind_missing[ind_missing[,1] == t_missing[i], 2]]
    #   V22 = V[-ind_missing[ind_missing[,1] == t_missing[i], 2], -ind_missing[ind_missing[,1] == t_missing[i], 2]]
    #   inv_V22 = solve(V22)
    #   V21 = t(V12)
    #   
    #   ## Parametros da condicional
    #   mu_condicional = mu1 + V12%*%inv_V22%*%(aux2-mu2) ## aux 2 eh yobs e mu2 é Fxthetaxobs
    #   sigma_condicional = V11 - V12%*%inv_V22%*%V21
    #   
    #   y[t_missing[i], ind_missing[ind_missing[,1] == t_missing[i], 2]] = MASS::mvrnorm(1, mu = mu_condicional, Sigma = sigma_condicional)
    #   input = append(input, y[t_missing[i], ind_missing[ind_missing[,1] == t_missing[i], 2]])
    # }
    # input.post[k,] = input
    
  }
  return(list(mu = mu.post[seq(bn+1, nit, by = thin),,],
              theta = theta.post[seq(bn+1, nit, by = thin),,],
              V = V.post[seq(bn+1, nit, by = thin),,], 
              Wt = Wt[seq(bn+1, nit, by = thin),,,], 
              input = input.post[seq(bn+1, nit, by = thin),]))
}



#### calculo do qx

qx_fitted <- function(fit){
  # fit$mu: array [idades x períodos x populações] com médias
  # fit$V: array [idades x populações x populações] com covariâncias para cada idade
  
  n_i <- dim(fit$mu)[1]   # número de idades
  n_t <- dim(fit$mu)[2]   # número de períodos
  n_j <- dim(fit$mu)[3]   # número de populações
  
  samples <- fit$mu
  V.samples <- fit$V
  
  # array para guardar os draws simulados
  fitted <- array(NA, dim = dim(samples))
  
  # Simulação multivariada para cada idade e período
  for(i in 1:n_i){
    for(t in 1:n_t){
      fitted[i,t,] <- exp(MASS::mvrnorm(1, mu = samples[i,t,], Sigma = V.samples[i,,]))
    }
  }
  
  # Calcula os quantis (mediana e intervalos 95%)
  qx_fitted_array <- apply(fitted, 2:3, quantile, probs = c(0.5, 0.025, 0.975))
  
  # Lista para armazenar resultados de cada população
  result_list <- vector("list", n_j)
  
  for(j in 1:n_j){
    mat <- t(qx_fitted_array[,,j])
    colnames(mat) <- c("qx.fitted", "qx.lower", "qx.upper")
    result_list[[j]] <- mat
  }
  
  # Retorna a lista, cada elemento corresponde a uma população
  return(result_list)
}



### valores iniciais de V
make_pd <- function(S, eps = 1e-6){
  S <- (S + t(S)) / 2
  eig <- eigen(S)
  eig$values[eig$values < eps] <- eps
  S_pd <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
  return(S_pd)
}
