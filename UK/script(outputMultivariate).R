library(dplyr)
library(ggplot2)
library(BayesMortalityPlus)
library(patchwork)
library(tidyverse)
require(stringr)
library(prodlim)
library(purrr)
library(tidyr)
library(scales)
require(gridExtra)
library(tmvtnorm)
library(matrixcalc)
library(corpcor)


grad <- seq(0,110,20)
grad[(grad/5)%%2 != 0] <- ""
point <- format_format(big.mark = " ", decimal.mark = ".", scientific = FALSE)


#####==========#####
##### dataset  #####
#####==========#####
dx <- read.csv("study2/deaths_MYB2.csv")
Ex <- read.csv("study2/MYB2_mais_OLDAGE.csv")

dx.all = dx %>%
  dplyr::select(!X) %>% 
  pivot_longer(cols = starts_with("X"), values_to = "dx", names_to = "year") %>%
  mutate(year = substring(year, 2)) %>%
  filter(Age != '105+' & Age != '0') %>%
  mutate(age = as.numeric(Age)) %>%
  dplyr::select(age, sex, year, dx)

Ex.all = Ex %>%
  dplyr::select(!X) %>% 
  pivot_longer(cols = starts_with("population_"), values_to = "Ex", names_to = "year") %>%
  mutate(year = str_remove(year, "population_")) %>%
  filter(age != 105 & age !=0 )

mx = dx.all %>% 
  left_join(Ex.all, by = c("sex", "age", "year")) %>%
  mutate(mx = dx/Ex) %>%
  mutate(sex = as.character(sex))

#---------------------------------------------------
qx_fitted = function(fit){
  n = dim(fit$mu)[1]
  # samples = fit$mu[seq(bn+1, n, by = thin),,]
  # V.samples = fit$V[seq(bn+1, n, by = thin),,]
  samples = fit$mu
  V.samples = fit$V
  fitted <- array(NA, dim = dim(samples))
  
  for(i in 1:dim(samples)[1]){
    for(t in 1:dim(samples)[2])
      fitted[i,t,] = exp(MASS::mvrnorm(1, mu = samples[i,t,], Sigma = V.samples[i,,]))
  }
  
  qx_fitted = apply(fitted, 2:3, quantile, probs = c(0.5, 0.025, 0.975))
  qx_fitted_m = t(qx_fitted[,,1])
  qx_fitted_f = t(qx_fitted[,,2])
  
  colnames(qx_fitted_m) = colnames(qx_fitted_f) = c("qx.fitted", "qx.lower", "qx.upper")
  
  return(list(qx_fitted_m, qx_fitted_f))
  
}

####### model
Ft = matrix(c(1, 1, 0, 0,
              1, 0, 0, 0), nrow = 2, ncol = 4, byrow = T)

Gt = matrix(c(1, 0, 1, 0,
              0, 1, 0, 1,
              0, 0, 1, 0,
              0, 0, 0, 1), nrow = 4, ncol = 4, byrow = T)

#d12 = 0.85
d12 <- rep(c(0.99, 0.8, 0.99), c(5, 85, 14)) ## Ate 104 
w<-104

m0 = rep(0, nrow(Gt))
C0 = diag(100, nrow(Gt))
s0 = diag(2)*0.01
v0 = 3

bn = 3000; thin = 2;  bn.sample=4000 ; nit=5000
aux.names <- c("qx.fitted", "qx.lower", "qx.upper")


### USUAL
ytime=2010
df.mx<- mx %>% mutate(log.mx=log(mx))
y.m<- df.mx %>%
  filter(year==!!ytime, sex=='1') %>%
  select(log.mx) 
y.f<- df.mx %>% 
  filter(year==!!ytime, sex=='2') %>% 
  select(log.mx) 
y_2010<- as.matrix(cbind(y.m,y.f)) 

ytime=2011
df.mx<- mx %>% mutate(log.mx=log(mx))
y.m<- df.mx %>%
  filter(year==!!ytime, sex=='1') %>%
  select(log.mx) 
y.f<- df.mx %>% 
  filter(year==!!ytime, sex=='2') %>% 
  select(log.mx) 
y_2011<- as.matrix(cbind(y.m,y.f))

ytime=2012
df.mx<- mx %>% mutate(log.mx=log(mx))
y.m<- df.mx %>%
  filter(year==!!ytime, sex=='1') %>%
  select(log.mx) 
y.f<- df.mx %>% 
  filter(year==!!ytime, sex=='2') %>% 
  select(log.mx) 
y_2012<- as.matrix(cbind(y.m,y.f))


source("study2/predict_dlm_at.R")

predict_conv_chain <- function(obj_pred, age, h, y_base, delta.aux,
                               V = 0.01, final_age ){
  k <- (final_age - max(age))
  cw = seq(0, 1, length.out = k)
  fit_base <- predict2.DLM(dlm(y_base, delta = delta.aux, ages = age,
                               M = nrow(obj_pred)), h = h, V = V) ##returns qx
  
 
  peso <- matrix(c(cw, rep(1, h - k)),
                 ncol = h,
                 nrow = nrow(obj_pred),
                 byrow = TRUE)
  M_pred <- obj_pred[,,1]
  F_pred <- obj_pred[,,2]
 
   M_blend <- M_pred * (1 - peso) + fit_base * peso
  F_blend <- F_pred * (1 - peso) + fit_base * peso
  
  center <- (M_blend + F_blend) / 2 
  diff   <- M_blend - F_blend 
  shrink_w <- peso
  diff_new <- (1 - shrink_w) * diff
  M_final <- center + diff_new / 2
  F_final <- center - diff_new / 2
  
  mx.m <- t(apply(M_final, 2, quantile, c(0.5, 0.025, 0.975)))
  mx.f <- t(apply(F_final, 2, quantile, c(0.5, 0.025, 0.975)))
  
  return(list(
    mx.m = mx.m,
    mx.f = mx.f
  ))
  
}


mx.mar = mx %>% select(age, year, dx, Ex) %>%
  group_by(age, year) %>%
  summarise(dx = sum(dx), Ex = sum(Ex)) %>%
  mutate(logmx = log(dx/Ex))


source("study2/ffbs_-_bivar.R")
source("study2/predict_bivar2_att_chain.R")
### ajuste para os anos 2010, 2011 e 2012
fitB_10 <- readRDS("study2/fitB2010_nontrunc.RDS")
fitB_11 <- readRDS("study2/fitB2011_nontrunc.RDS")
fitB_12 <- readRDS("study2/fitB2012_nontrunc.RDS")


h = 16
fitB.pred.10 <- predict_bivar2_chain(y_2010, m0, C0, h, V = fitB_10$V,
                               Ft, Gt, delta=d12, ages = 1:w)

fitB.pred.11 <- predict_bivar2_chain(y_2011, m0, C0, h, V = fitB_11$V,
                               Ft, Gt, delta=d12, ages = 1:w)

fitB.pred.12 <- predict_bivar2_chain(y_2012, m0, C0, h, V = fitB_12$V,
                               Ft, Gt, delta=d12, ages = 1:w)


pred.10 <- predict_conv_chain(obj_pred = fitB.pred.10,
                              age = 1:95,
                              h = 16,
                              y_base = mx.mar %>% filter(year == 2010) %>% pull(logmx),
                              delta = d12,
                              V = 0.01,
                              final_age = 95)







qxB.m1 = data.frame(age=1:(w+h), rbind(qx_fitted(fitB_10)[[1]],
                                       setNames(pred.10[[1]], aux.names)) )
qxB.f1 = data.frame(age=1:(w+h), rbind(qx_fitted(fitB_10)[[2]],
                                       setNames(pred.10[[2]], aux.names)) )



pred.11 <- predict_conv_chain(obj_pred = fitB.pred.11,
                              age = 1:104,
                              h = 16,
                              y_base = mx.mar %>% filter(year == 2011) %>% pull(logmx),
                              delta = d12,
                              V = 0.01,
                              final_age = 115)


qxB.m2 = data.frame(age=1:(w+h), rbind(qx_fitted(fitB_11)[[1]],
                                       setNames(pred.11[[1]], aux.names)) )
qxB.f2 = data.frame(age=1:(w+h), rbind(qx_fitted(fitB_11)[[2]],
                                       setNames(pred.11[[2]], aux.names)) )


pred.12 <- predict_conv_chain(obj_pred = fitB.pred.12,
                              age = 1:104,
                              h = 16,
                              y_base = mx.mar %>% filter(year == 2012) %>% pull(logmx),
                              delta = d12,
                              V = 0.01,
                              final_age = 115)


qxB.m3 = data.frame(age=1:(w+h), rbind(qx_fitted(fitB_12)[[1]],
                                       setNames(pred.12[[1]], aux.names)) )
qxB.f3 = data.frame(age=1:(w+h), rbind(qx_fitted(fitB_12)[[2]],
                                       setNames(pred.12[[2]], aux.names)) )

qxB.all<- bind_rows(qxB.m1,qxB.f1,
                    qxB.m2,qxB.f2,
                    qxB.m3,qxB.f3, .id="id")  

out<-bind_rows(qxB.all,
               .id="id1")
out. = out %>% 
  mutate(sex = ifelse(id %in% c(1,3,5), 'male', 'female'))
out.[241:480,1] <- rep(2, length(241:480)) 
out.[481:720,1] <- rep(3, length(481:720)) 

## crude mortality rate 
mx.<- mx %>% 
  mutate(id1 = rep(1:3, length.out = 624)) %>%
  select(age,mx,sex, id1) 
mx. = mx. %>% 
  mutate(sex = ifelse(sex == 1, 'male', 'female')) 

### ZOOM
mx..<- mx. %>%
  filter(age >=80)
out..<- out. %>%
  filter(age>=80)

ggplot() +
  geom_point(data = mx.., aes(x = age, y = mx, color=sex)) + #, shape=sex), col="gray") +
  #geom_line(data=out.., aes(x=age, y = -log(1-qx.fitted), color = sex)) + 
  geom_line(data=out.., aes(x=age, y = qx.fitted, color = sex)) +
  geom_hline(yintercept=1, linetype='dashed', col = 'black')+
  geom_ribbon(data = out.., aes(x = age, ymin = qx.lower, ymax = qx.upper, fill=sex), alpha = 0.25) + 
  theme_classic(base_size = 20) + scale_y_continuous("qx", limits = 10^-c(NA,NA), 
                                                     trans = 'log10',
                                                     labels = scales::comma) + 
  scale_x_continuous("Age", breaks = seq(0, 200, by = 10)) + 
  scale_color_manual(values=c( "tomato","steelblue"), labels=c("female", "male"))+
  scale_fill_manual(values=c("tomato","steelblue"), labels=c("female", "male")) + 
  theme(legend.position = c(0.93,0.15), strip.background=element_rect(colour="black", fill="gray87"), 
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1), 
        legend.title=element_blank(), legend.key.height = unit(.6, "cm") , 
        legend.text=element_text(color="black", size=16), 
        axis.title = element_text(color = "black", size = 14), 
        axis.text = element_text(color="black",size=14))+ 
  facet_wrap(~id1, ncol=3, labeller = labeller(id1 = c("1" = "2010",
                                                       "2" = "2011",
                                                       "3" = "2012"))) -> pt
plotly::ggplotly(pt)
graphics.off()


##### PLOTS
## crude mortality rate 
pdf("Fig15a.pdf", width=12, height=5)
ggplot() +
  geom_point(data = mx., aes(x = age, y = mx, color=sex)) + #, shape=sex), col="gray") +
  #geom_line(data=out.., aes(x=age, y = -log(1-qx.fitted), color = sex)) + 
  geom_line(data=out., aes(x=age, y = qx.fitted, color = sex)) +
  #geom_hline(yintercept=1, linetype='dashed', col = 'black')+
  geom_ribbon(data = out., aes(x = age, ymin = qx.lower, ymax = qx.upper, fill=sex), alpha = 0.25) + 
  theme_classic(base_size = 20) + scale_y_continuous(expression(m[x]), limits = 10^-c(NA,NA), 
                                                     trans = 'log10', labels = scales::comma) + 
  scale_x_continuous("Age", breaks = seq(0, 200, by = 20)) + 
  scale_color_manual(values=c( "tomato","steelblue"), labels=c("female", "male"))+
  scale_fill_manual(values=c("tomato","steelblue"), labels=c("female", "male")) + 
  theme(legend.position = c(0.93,0.15), strip.background=element_rect(colour="black", fill="gray87"), 
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1), 
        legend.title=element_blank(), legend.key.height = unit(.6, "cm") , 
        legend.text=element_text(color="black", size=16), 
        axis.title = element_text(color = "black", size = 14), 
        axis.text = element_text(color="black",size=14))+ 
  facet_wrap(~id1, ncol=3, labeller = labeller(id1 = c("1" = "2010",
                                                       "2" = "2011",
                                                       "3" = "2012"))) 
graphics.off()

### ZOOM
pdf("Fig15b.pdf", width=12, height=5)
ggplot() +
  geom_point(data = mx.., aes(x = age, y = mx, color=sex)) + #, shape=sex), col="gray") +
  #geom_line(data=out.., aes(x=age, y = -log(1-qx.fitted), color = sex)) + 
  geom_line(data=out.., aes(x=age, y = qx.fitted, color = sex)) +
  #geom_hline(yintercept=1, linetype='dashed', col = 'black')+
  geom_ribbon(data = out.., aes(x = age, ymin = qx.lower, ymax = qx.upper, fill=sex), alpha = 0.25) + 
  theme_classic(base_size = 20) + scale_y_continuous(expression(m[x]), limits = 10^-c(NA,NA), 
                                                     trans = 'log10', labels = scales::comma) + 
  scale_x_continuous("Age", breaks = seq(0, 200, by = 10)) + 
  scale_color_manual(values=c( "tomato","steelblue"), labels=c("female", "male"))+
  scale_fill_manual(values=c("tomato","steelblue"), labels=c("female", "male")) + 
  theme(legend.position = c(0.93,0.15), strip.background=element_rect(colour="black", fill="gray87"), 
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1), 
        legend.title=element_blank(), legend.key.height = unit(.6, "cm") , 
        legend.text=element_text(color="black", size=16), 
        axis.title = element_text(color = "black", size = 14), 
        axis.text = element_text(color="black",size=14))+ 
  facet_wrap(~id1, ncol=3, labeller = labeller(id1 = c("1" = "2010",
                                                       "2" = "2011",
                                                       "3" = "2012"))) 
graphics.off()


mx..<- mx. %>%
  filter(age<=30)
out..<- out. %>%
  filter(age<=30)

pdf("Fig7c.pdf", width=12, height=5)
ggplot() +
  geom_point(data=mx.., aes(x=age, y =  mx, col=sex)) +
  geom_line(data=out.., aes(x=age, y = qx.fitted, color = sex)) +
  theme_classic(base_size = 20) +
  geom_hline(yintercept=1, linetype='dashed', col = 'black')+
  geom_ribbon(data = out.., aes(x = age, ymin = qx.lower, ymax = qx.upper, fill=sex), alpha = 0.25) +
  scale_y_continuous(expression(m[x]), limits = c(0.00003, 0.0030),
                     trans = 'log10', labels = scales::comma) +
  scale_x_continuous("Age", breaks = seq(0, 40, by = 10))+
  scale_color_manual(values=c( "tomato","steelblue"), labels=c("female", "male"))+
  scale_fill_manual(values=c("tomato","steelblue"), labels=c("female", "male")) +
  theme(legend.position = c(0.92, 0.15),
        strip.background=element_rect(colour="black",
                                      fill="gray87"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.title=element_blank(),
        legend.key.height = unit(.6, "cm") ,
        legend.text=element_text(color="black", size=16),
        axis.title = element_text(color = "black", size = 14),
        axis.text = element_text(color="black",size=14))+
  facet_wrap(~id1, ncol=3, labeller = labeller(id1 = c("1" = "2010",
                                                       "2" = "2011",
                                                       "3" = "2012"))) 
graphics.off()
