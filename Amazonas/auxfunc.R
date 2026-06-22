######################################################
#### authors: figueiredo, lobo, alves and fonseca
#### output: auxiliary functions
####  modelling multiple populations
### data: amazonas, 2022-2024
######################################################

### funcoes uteis 
buildFtGt <- function(J, alpha = TRUE){
  if(alpha){
    # dimensão total:
    # mu_1,...,mu_J,
    # beta_1,...,beta_J,
    # alpha
    m <- 2*J + 1
    
    Gt <- matrix(0, nrow = m, ncol = m)
    
    # nível
    for(j in 1:J){
      Gt[j, j] <- 1
      Gt[j, J + j] <- 1
    }
    
    # tendência
    for(j in 1:J){
      Gt[J + j, J + j] <- 1
    }
    
    # termo comum
    Gt[m, m] <- 1
    
    Ft <- matrix(0, nrow = J, ncol = m)
    
    for(j in 1:J){
      Ft[j, j] <- 1
      Ft[j, m] <- 1
    }
  } else {
    
    # sem termo comum
    m <- 2*J
    Gt <- matrix(0, nrow = m, ncol = m)
    # nível
    for(j in 1:J){
      Gt[j, j] <- 1
      Gt[j, J + j] <- 1
    }
    # tendência
    for(j in 1:J){
      Gt[J + j, J + j] <- 1
    }
    Ft <- matrix(0, nrow = J, ncol = m)
    for(j in 1:J){
      Ft[j, j] <- 1
    }
  }
  return(list(Ft = Ft, Gt = Gt))
}

## convergencia entre curvas 
predict_conv_chain <- function(obj_pred, age, h, y_base, delta.aux,
                               V = 0.01, final_age ){
  
  step_age <- diff(age)[1]
  k <- (final_age - max(age)) / step_age
  
  cw <- seq(0, 1, length.out = k)
  
  
  ### eh o que faz o blend de uma tabua media - o delta nao eh o mesma matrix delta do multivariado
  fit_base <- predict2.DLM(
    dlm(y_base, delta = delta.aux, ages = age,
        M = nrow(obj_pred)),
    h = h, V = V
  )
  
  # peso ao longo das idades...
  peso <- matrix(c(cw, rep(1, h - k)),
                 ncol = h,
                 nrow = nrow(obj_pred),
                 byrow = TRUE)
  # médias preditas dos dois modelos
  M_pred <- obj_pred[,,1]
  F_pred <- obj_pred[,,2]
  # blending
  M_blend <- M_pred * (1 - peso) + fit_base * peso
  F_blend <- F_pred * (1 - peso) + fit_base * peso
  
  center <- (M_blend + F_blend) / 2 ### nivel medio da mortalodade
  diff   <- M_blend - F_blend # gap entre as populacoes
  # intensidade do shrink pode variar na idade
  shrink_w <- peso
  diff_new <- (1 - shrink_w) * diff ### shrinkage da diferenca, neste aso
  # reduz o quanto duas curvas são diferentes entre si, sem mexer diretamente no nível médio delas.
  M_final <- center + diff_new / 2
  F_final <- center - diff_new / 2
  
  mx.m <- t(apply(M_final, 2, quantile, c(0.5, 0.025, 0.975)))
  mx.f <- t(apply(F_final, 2, quantile, c(0.5, 0.025, 0.975)))
  
  return(list(
    mx.m = mx.m,
    mx.f = mx.f
  ))
}
