######################################
##### output univariate modelling

#---------------------------------------------------------------
### cenario 1 

## delta fixo
fitU.2010m =  fitU.2010f= list()
fitU.2011m =  fitU.2011f= list()
fitU.2012m =  fitU.2012f= list()


for(i in 1:4){
  fitU.2010m[[i]] <- dlm(y.2010m$log.mx, delta = delta[i], prior = list(m0 = c(0,0), C0 = diag(1, 2)), ages = 1:104)
  fitU.2010f[[i]] <- dlm(y.2010f$log.mx, delta = delta[i], prior = list(m0 = c(0,0), C0 = diag(1, 2)), ages = 1:104)
  
  fitU.2011m[[i]] <- dlm(y.2011m$log.mx, delta = delta[i], prior = list(m0 = c(0,0), C0 = diag(1, 2)), ages = 1:104)
  fitU.2011f[[i]] <- dlm(y.2011f$log.mx, delta = delta[i], prior = list(m0 = c(0,0), C0 = diag(1, 2)), ages = 1:104)
  
  fitU.2012m[[i]] <- dlm(y.2012m$log.mx, delta = delta[i], prior = list(m0 = c(0,0), C0 = diag(1, 2)), ages = 1:104)
  fitU.2012f[[i]] <- dlm(y.2012f$log.mx, delta = delta[i], prior = list(m0 = c(0,0), C0 = diag(1, 2)), ages = 1:104)
}


### Montando as estimativas do modelo com as extrapolacoes
qx.m = lapply(fitU.2010m, fitted)
qx.f = lapply(fitU.2010f, fitted)
## Aplicacao do predict 
qx.pred.m = lapply(fitU.2010m, predict.DLM, h = 16)
qx.pred.f = lapply(fitU.2010f, predict.DLM, h = 16)

qx.m2010 = qx.f2010 =  list()
for(i in 1:4){
  qx.m2010[[i]] = rbind(qx.m[[i]], qx.pred.m[[i]])
  qx.f2010[[i]] = rbind(qx.f[[i]], qx.pred.f[[i]])
}

qx.m = lapply(fitU.2011m, fitted)
qx.f = lapply(fitU.2011f, fitted)
## Aplicacao do predict 
qx.pred.m = lapply(fitU.2011m, predict.DLM, h = 16)
qx.pred.f = lapply(fitU.2011f, predict.DLM, h = 16)

qx.m2011 = qx.f2011 =  list()
for(i in 1:4){
  qx.m2011[[i]] = rbind(qx.m[[i]], qx.pred.m[[i]])
  qx.f2011[[i]] = rbind(qx.f[[i]], qx.pred.f[[i]])
}

qx.m = lapply(fitU.2012m, fitted)
qx.f = lapply(fitU.2012f, fitted)
## Aplicacao do predict 
qx.pred.m = lapply(fitU.2012m, predict.DLM, h = 16)
qx.pred.f = lapply(fitU.2012f, predict.DLM, h = 16)

qx.m2012 = qx.f2012 =  list()
for(i in 1:4){
  qx.m2012[[i]] = rbind(qx.m[[i]], qx.pred.m[[i]])
  qx.f2012[[i]] = rbind(qx.f[[i]], qx.pred.f[[i]])
}

## Cenario 1: delta variando 
fitU.2010m2 <- dlm(y.2010m$log.mx, delta = d12, prior = list(m0 = c(0,0), C0 = diag(1, 2)), ages = 1:104)
fitU.2010f2 <- dlm(y.2010f$log.mx, delta = d12, prior = list(m0 = c(0,0), C0 = diag(1, 2)), ages = 1:104)

fitU.2011m2 <- dlm(y.2011m$log.mx, delta = d12, prior = list(m0 = c(0,0), C0 = diag(1, 2)), ages = 1:104)
fitU.2011f2 <- dlm(y.2011f$log.mx, delta = d12, prior = list(m0 = c(0,0), C0 = diag(1, 2)), ages = 1:104)

fitU.2012m2 <- dlm(y.2012m$log.mx, delta = d12, prior = list(m0 = c(0,0), C0 = diag(1, 2)), ages = 1:104)
fitU.2012f2 <- dlm(y.2012f$log.mx, delta = d12, prior = list(m0 = c(0,0), C0 = diag(1, 2)), ages = 1:104)

qx.m2 = fitted(fitU.2010m2)
qx.f2 = fitted(fitU.2010f2)
## Aplicacao do predict 
qx.pred.m2 = predict.DLM(fitU.2010m2,  h = 16)
qx.pred.f2 = predict.DLM(fitU.2010f2,  h = 16)
qx.m2010.d12 = rbind(qx.m2, qx.pred.m2)
qx.f2010.d12 = rbind(qx.f2, qx.pred.f2)

qx.m2 = fitted(fitU.2011m2)
qx.f2 = fitted(fitU.2011f2)
## Aplicacao do predict 
qx.pred.m2 = predict.DLM(fitU.2011m2,  h = 16)
qx.pred.f2 = predict.DLM(fitU.2011f2,  h = 16)
qx.m2011.d12 = rbind(qx.m2, qx.pred.m2)
qx.f2011.d12 = rbind(qx.f2, qx.pred.f2)

qx.m2 = fitted(fitU.2012m2)
qx.f2 = fitted(fitU.2012f2)
## Aplicacao do predict 
qx.pred.m2 = predict.DLM(fitU.2012m2,  h = 16)
qx.pred.f2 = predict.DLM(fitU.2012f2,  h = 16)
qx.m2012.d12 = rbind(qx.m2, qx.pred.m2)
qx.f2012.d12 = rbind(qx.f2, qx.pred.f2)

#---------------------------------------------------
#---------------------------------------------------------------
### cenario 2 

## delta fixo
fitU2.2010m =  fitU2.2010f= list()
fitU2.2011m =  fitU2.2011f= list()
fitU2.2012m =  fitU2.2012f= list()


for(i in 1:4){
  fitU2.2010m[[i]] <- dlm(y.2010m$log.mx[1:100], delta = delta[i], prior = list(m0 = c(0,0), C0 = diag(1, 2)), ages = 1:100)
  fitU2.2010f[[i]] <- dlm(y.2010f$log.mx[1:100], delta = delta[i], prior = list(m0 = c(0,0), C0 = diag(1, 2)), ages = 1:100)
  
  fitU2.2011m[[i]] <- dlm(y.2011m$log.mx[1:100], delta = delta[i], prior = list(m0 = c(0,0), C0 = diag(1, 2)), ages = 1:100)
  fitU2.2011f[[i]] <- dlm(y.2011f$log.mx[1:100], delta = delta[i], prior = list(m0 = c(0,0), C0 = diag(1, 2)), ages = 1:100)
  
  fitU2.2012m[[i]] <- dlm(y.2012m$log.mx[1:100], delta = delta[i], prior = list(m0 = c(0,0), C0 = diag(1, 2)), ages = 1:100)
  fitU2.2012f[[i]] <- dlm(y.2012f$log.mx[1:100], delta = delta[i], prior = list(m0 = c(0,0), C0 = diag(1, 2)), ages = 1:100)
}


### Montando as estimativas do modelo com as extrapolacoes
qx2.m = lapply(fitU2.2010m, fitted)
qx2.f = lapply(fitU2.2010f, fitted)
## Aplicacao do predict 
qx2.pred.m = lapply(fitU2.2010m, predict.DLM, h = 20)
qx2.pred.f = lapply(fitU2.2010f, predict.DLM, h = 20)

qx2.m2010 = qx2.f2010 =  list()
for(i in 1:4){
  qx2.m2010[[i]] = rbind(qx2.m[[i]], qx2.pred.m[[i]])
  qx2.f2010[[i]] = rbind(qx2.f[[i]], qx2.pred.f[[i]])
}

qx2.m = lapply(fitU2.2011m, fitted)
qx2.f = lapply(fitU2.2011f, fitted)
## Aplicacao do predict 
qx2.pred.m = lapply(fitU2.2011m, predict.DLM, h = 20)
qx2.pred.f = lapply(fitU2.2011f, predict.DLM, h = 20)

qx2.m2011 = qx2.f2011 =  list()
for(i in 1:4){
  qx2.m2011[[i]] = rbind(qx2.m[[i]], qx2.pred.m[[i]])
  qx2.f2011[[i]] = rbind(qx2.f[[i]], qx2.pred.f[[i]])
}

qx2.m = lapply(fitU2.2012m, fitted)
qx2.f = lapply(fitU2.2012f, fitted)
## Aplicacao do predict 
qx2.pred.m = lapply(fitU2.2012m, predict.DLM, h = 20)
qx2.pred.f = lapply(fitU2.2012f, predict.DLM, h = 20)

qx2.m2012 = qx2.f2012 =  list()
for(i in 1:4){
  qx2.m2012[[i]] = rbind(qx2.m[[i]], qx2.pred.m[[i]])
  qx2.f2012[[i]] = rbind(qx2.f[[i]], qx2.pred.f[[i]])
}

## Cenario 2: delta variando 
fitU2.2010m2 <- dlm(y.2010m$log.mx[1:100], delta = d22, prior = list(m0 = c(0,0), C0 = diag(1, 2)), ages = 1:100)
fitU2.2010f2 <- dlm(y.2010f$log.mx[1:100], delta = d22, prior = list(m0 = c(0,0), C0 = diag(1, 2)), ages = 1:100)

fitU2.2011m2 <- dlm(y.2011m$log.mx[1:100], delta = d22, prior = list(m0 = c(0,0), C0 = diag(1, 2)), ages = 1:100)
fitU2.2011f2 <- dlm(y.2011f$log.mx[1:100], delta = d22, prior = list(m0 = c(0,0), C0 = diag(1, 2)), ages = 1:100)

fitU2.2012m2 <- dlm(y.2012m$log.mx[1:100], delta = d22, prior = list(m0 = c(0,0), C0 = diag(1, 2)), ages = 1:100)
fitU2.2012f2 <- dlm(y.2012f$log.mx[1:100], delta = d22, prior = list(m0 = c(0,0), C0 = diag(1, 2)), ages = 1:100)

qx2.m2 = fitted(fitU2.2010m2)
qx2.f2 = fitted(fitU2.2010f2)
## Aplicacao do predict 
qx2.pred.m2 = predict.DLM(fitU2.2010m2,  h = 20)
qx2.pred.f2 = predict.DLM(fitU2.2010f2,  h = 20)
qx2.m2010.d22 = rbind(qx2.m2, qx2.pred.m2)
qx2.f2010.d22 = rbind(qx2.f2, qx2.pred.f2)

qx2.m2 = fitted(fitU2.2011m2)
qx2.f2 = fitted(fitU2.2011f2)
## Aplicacao do predict 
qx2.pred.m2 = predict.DLM(fitU2.2011m2,  h = 20)
qx2.pred.f2 = predict.DLM(fitU2.2011f2,  h = 20)
qx2.m2011.d22 = rbind(qx2.m2, qx2.pred.m2)
qx2.f2011.d22 = rbind(qx2.f2, qx2.pred.f2)

qx2.m2 = fitted(fitU2.2012m2)
qx2.f2 = fitted(fitU2.2012f2)
## Aplicacao do predict 
qx2.pred.m2 = predict.DLM(fitU2.2012m2,  h = 20)
qx2.pred.f2 = predict.DLM(fitU2.2012f2,  h = 20)
qx2.m2012.d22 = rbind(qx2.m2, qx2.pred.m2)
qx2.f2012.d22 = rbind(qx2.f2, qx2.pred.f2)

#---------------------------------------------------
