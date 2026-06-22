

source("ffbs_-_bivar.R")
source('gibbs_missing_data.R')


### no missing data - multivariate model without common term
Ft = matrix(c(1, 0, 0, 0,
              0, 1, 0, 0), nrow = 2, ncol = 4, byrow = T)
Gt = matrix(c(1,0,1,0,
              0,1,0,1,
              0,0,1,0,
              0,0,0,1), nrow = 4, ncol = 4, byrow = T)

m0 = rep(0, nrow(Gt))
C0 = diag(100, nrow(Gt))
s0 <- diag(2)*0.01  ###priori do petris vaga
v0 <- 3

d12 <- rep(c(0.99, 0.85, 0.99), c(5, 85, 14)) ## Ate 104
fitB<- gibbsV_corrigido(y = yt, delta = d12, m0 = m0, C0 = C0,
                        Ft0 = Ft, Gt = Gt, nit = 5000, v0 = v0, s0 = s0, V = V)

## no missing data - multivariate model without common term
Gt = matrix(c(1,0,1,0,0,
              0,1,0,1,0,
              0,0,1,0,0,
              0,0,0,1,0,
              1,1,0,0,1), nrow = 5, ncol = 5, byrow = T)
Ft = matrix(c(1, 0, 0, 0, 1,
              0, 1, 0, 0, 1), nrow = 2, ncol = 5, byrow = T)
m0 = rep(0, nrow(Gt))
C0 = diag(100, nrow(Gt))
fitB.TC<- gibbsV_corrigido(y = yt, delta = d12, m0 = m0, C0 = C0,
                           Ft0 = Ft, Gt = Gt, nit = 5000, v0 = v0, s0 = s0, V = V)




## fitted no missing
qxB.m = data.frame(age=1:(w), rbind(qx_fitted(fitB)[[1]]))
qxB.f = data.frame(age=1:(w), rbind(qx_fitted(fitB)[[2]]))
qxBTC.m = data.frame(age=1:(w), rbind(qx_fitted(fitB.TC)[[1]]))
qxBTC.f = data.frame(age=1:(w), rbind(qx_fitted(fitB.TC)[[2]]))


outAll<-bind_rows(qxB.all,qxBTC.all,
                  .id="id1")
outAll = outAll %>% 
  mutate(sex = ifelse(id == 1, 'male', 'female')) 

outAll. = outAll %>% 
  mutate(sex = ifelse(id == 1, 'male', 'female')) %>%
  filter(age <=35)

