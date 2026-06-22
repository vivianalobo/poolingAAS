######################################################
#### authors: figueiredo, lobo, alves and fonseca
#### output: metrics 
####  modelling multiple populations
### data: amazonas, 2022-2024
######################################################


## measure functions
### log posterior predictive density
library(mvtnorm)
log_pred_density <- function(y, mu_post, V_post){
  
  nsim <- dim(mu_post)[1]
  N    <- dim(mu_post)[2]
  J    <- dim(mu_post)[3]
  
  logscore <- numeric(N)
  
  for(t in 1:N){
    dens <- numeric(nsim)
    yt <- y[t, ]
    for(m in 1:nsim){
      mut <- mu_post[m, t, ]
      Vt <- V_post[m, , ]
      dens[m] <- dmvnorm(yt,mean = mut, sigma = Vt,log = FALSE)
    }
    logscore[t] <- log(mean(dens))
  }
  
  return(sum(logscore))
}

extract_var <- function(data, var,regions,sex = NULL,log.transform = FALSE){
  
  out <- lapply(regions, function(r){
    
    df <- data %>%
      filter(Microrregião.IBGE == r)
    
    if(!is.null(sex)){
      df <- df %>% filter(sex == !!sex)
    }
    
    x <- df[[var]]
    if(log.transform)
      x <- log(x)
    x
  })
  names(out) <- regions
  return(out)
}

### entra com os dados observados, ajuste e exposicao nas idades para cada populacao j
weighted_metrics <- function(obs, fit, expo){
  
  ### population j
  wlmse.pop <- sapply(seq_len(ncol(obs)), function(j){
    
    sum(
      expo[,j] * (obs[,j] - fit[,j])^2,
      na.rm = TRUE
    ) /
      sum(expo[,j], na.rm = TRUE)
    
  })
  
  
  wmae.pop <- sapply(seq_len(ncol(obs)), function(j){
    sum(
      expo[,j] * abs(obs[,j] - fit[,j]),
      na.rm = TRUE
    ) /
      sum(expo[,j], na.rm = TRUE)
  })
  
  ### medidas globais
  wlmse.total <- sum(
    as.matrix(expo) *
      (as.matrix(obs) - as.matrix(fit))^2,
    na.rm = TRUE
  ) /
    sum(as.matrix(expo), na.rm = TRUE)
  
  
  wmae.total <- sum(
    as.matrix(expo) *
      abs(as.matrix(obs) - as.matrix(fit)),
    na.rm = TRUE
  ) /
    sum(as.matrix(expo), na.rm = TRUE)
  
  
  out <- data.frame(
    Population = c(colnames(obs), "Total"),
    WLMSE = c(wlmse.pop, wlmse.total),
    WMAE  = c(wmae.pop,  wmae.total),
    row.names = NULL
  )
  
  return(out)
}
