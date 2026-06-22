
source("ffbs_-_bivar.R")
source('gibbs_missing_data.R')


### missing data - multivariate model without common term
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
fitB.miss<- gibbs_missing_data(y = y, Ft = Ft, Gt = Gt, nit = 5000, bn = 3000, thin = 1,
                               v0 = v0, s0 = s0, m0 = m0, C0 = C0, delta = d12, V = V)

## missing data - multivariate model without common term
Gt = matrix(c(1,0,1,0,0,
              0,1,0,1,0,
              0,0,1,0,0,
              0,0,0,1,0,
              1,1,0,0,1), nrow = 5, ncol = 5, byrow = T)
Ft = matrix(c(1, 0, 0, 0, 1,
              0, 1, 0, 0, 1), nrow = 2, ncol = 5, byrow = T)
m0 = rep(0, nrow(Gt))
C0 = diag(100, nrow(Gt))
fitB.TCmiss<- gibbs_missing_data(y = y, Ft = Ft, Gt = Gt, nit = 5000, bn = 3000, thin = 1,
                                 v0 = v0, s0 = s0, m0 = m0, C0 = C0, delta = d12, V = V)


### missing
qxBmiss.m = data.frame(age=1:(w), rbind(qx_fitted(fitB.miss)[[1]]))
qxBmiss.f = data.frame(age=1:(w), rbind(qx_fitted(fitB.miss)[[2]]))
qxBTCmiss.m = data.frame(age=1:(w), rbind(qx_fitted(fitB.TCmiss)[[1]]))
qxBTCmiss.f = data.frame(age=1:(w), rbind(qx_fitted(fitB.TCmiss)[[2]]))


### missing 
qxBmiss.all<- bind_rows(qxBmiss.m,qxBmiss.f, .id="id") 
qxBTCmiss.all<- bind_rows(qxBTCmiss.m,qxBTCmiss.f, .id="id") 


out<-bind_rows(qxBmiss.all,qxBTCmiss.all,
               .id="id1")
out = out %>% 
  mutate(sex = ifelse(id == 1, 'male', 'female')) 

out. = out %>% 
  mutate(sex = ifelse(id == 1, 'male', 'female')) %>%
  filter(age <=35)
