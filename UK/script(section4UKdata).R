######################################################
#### authors: figueiredo, lobo, alves and fonseca
#### output: modelling the UK data
### data: England and Wales, 2010-2012
######################################################


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

grad <- seq(0,110,20)
grad[(grad/5)%%2 != 0] <- ""
point <- format_format(big.mark = " ", decimal.mark = ".", scientific = FALSE)


### 4.2. England and Wales mortality data

#The England and Wales dataset \citep{OFS2019} contains death counts and mid-year population estimates for males and females
#by single year of age. Since \citet{forster22} analyse each year separately, we present only the 2010 results for illustration
#ses. Figure~\ref{fig:fig1} displays the crude mortality rates together with a moving-average smoother for both populations. 
#The smoothed curves reveal the well-known mortality differential between males and females throughout most ages, with male
#mortality generally exceeding female mortality. However, at the youngest and oldest ages, where exposure is lower and 
#mortality rates are more variable, the crude observations exhibit substantial fluctuations and occasional local crossings 
#between the two populations.

#-------------------------------------------------------------------------------------------------------------
### read mortality dataset
dx <- read.csv("deaths_MYB2.csv")
Ex <- read.csv("MYB2_mais_OLDAGE.csv")

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


mx.2010<- mx %>%filter(year == 2010)

mx2010.<- mx.2010 %>%
  filter(age <=30)
mx2010..<- mx.2010 %>%
  filter(age >=80)

mx.2010$group<- "age 1-104"
mx2010.$group<- "age < 30"
mx2010..$group<- "age  80+"


##### moving average
library(zoo)
mx.smooth2010 <- mx.2010 %>%
  group_by(sex) %>%
  arrange(age) %>%
  group_modify(~{
    
    dados <- .x
    
    d1 <- dados %>%
      filter(age <= 30)
    
    d1$mx.fit <- rollapply(
      d1$mx,
      width = 2.5,
      FUN = mean,
      partial = TRUE,
      align = "center"
    )
    
    d2 <- dados %>%
      filter(age > 30 & age <= 80)
    
    d2$mx.fit <- rollapply(
      d2$mx,
      width = 5,
      FUN = mean,
      partial = TRUE,
      align = "center"
    )
    
    d3 <- dados %>%
      filter(age > 80)
    
    d3$mx.fit <- rollapply(
      d3$mx,
      width = 3,
      FUN = mean,
      partial = TRUE,
      align = "center"
    )
    
    bind_rows(d1, d2, d3)
    
  })


mx2010smooth.<- mx.smooth2010 %>%
  filter(age <=30)
mx2010smooth..<- mx.smooth2010 %>%
  filter(age >=80)

mx.smooth2010$group<- "age 1-104"
mx2010smooth.$group<- "age < 30"
mx2010smooth..$group<- "age  80+"

plot.UKdata<- function(data, smooth){
  
  ymin <- 10^floor(log10(min(data$mx[data$mx > 0], na.rm = TRUE)))
  ymax <- 10^ceiling(log10(max(data$mx, na.rm = TRUE)))
  
  exp_min <- floor(log10(ymin))
  exp_max <- ceiling(log10(ymax))
  
  if((exp_max - exp_min) < 2){
    exp_min <- exp_min - 1
    exp_max <- exp_max + 1
  }
  
  breaks_y <- 10^(exp_min:exp_max)
  
  g1<- ggplot(data=data, aes(x=age, y =  mx, col = sex)) + 
    geom_point(size = 1.8) + 
    theme_classic(base_size = 20) + 
    geom_line(data=smooth,aes(x = age,y = mx.fit, col = sex),linewidth = 1.2) +
    scale_y_continuous(expression("log rate"), limits = c(min(breaks_y), max(breaks_y)),
                       breaks = breaks_y) +
    scale_x_continuous("Age", breaks = seq(0, 110, by = 20), labels = grad) +
    scale_color_manual(name="sex",
                       labels=c("male","female"),
                       values=c("steelblue","tomato"))+
    theme(legend.position = "none", strip.background=element_rect(colour="black",
                                                                  fill="gray87"), panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          legend.title=element_blank(),
          legend.text=element_text(size=rel(1)),
          axis.title = element_text(color = "black"),#, size = 14),
          axis.text = element_text(color="black"))+ #,size=14))+
    facet_wrap(~group)
  return(g1)
}



plot.UKdata <- function(data, smooth){
  
  ymin <- floor(min(log(data$mx), na.rm = TRUE))
  ymax <- ceiling(max(log(data$mx), na.rm = TRUE))
  
  if((ymax - ymin) < 4){
    ymin <- ymin - 1
    ymax <- ymax + 1
  }
  
  breaks_y <- seq(ymin, ymax, by = 2)
  
  g1 <- ggplot(data = data,
               aes(x = age,
                   y = log(mx),
                   col = sex)) +
    
    geom_point(size = 1.8) +
    
    geom_line(data = smooth,
              aes(x = age,
                  y = log(mx.fit),
                  col = sex),
              linewidth = 1.2) +
    
    scale_y_continuous(
      expression("log rate"),
      limits = c(ymin, ymax),
      breaks = breaks_y
    ) +
    
    scale_x_continuous(
      "Age",
      breaks = seq(0,110,by=20),
      labels = grad
    ) +
    
    scale_color_manual(
      name = "sex",
      labels = c("male","female"),
      values = c("steelblue","tomato")
    ) +
    
    theme_classic(base_size = 20) +
    
    theme(
      legend.position = "none",
      strip.background = element_rect(
        colour = "black",
        fill = "gray87"
      ),
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 1
      ),
      legend.title = element_blank(),
      legend.text = element_text(size = rel(1)),
      axis.title = element_text(color = "black"),
      axis.text = element_text(color = "black")
    ) +
    
    facet_wrap(~group)
  
  return(g1)
}

# Figure 9: Mortality rates for England and Wales (E+W) for ages 1-104 and year 2010. 
#The second and third columns zoom in on ages 40- and 80+, respectively. 
#Blue and red dots represent crude mortality rates for male and female populations. 
#Solid lines represent a moving average smoothing.
plot.UKdata(mx.2010, mx.smooth2010)
plot.UKdata(mx2010., mx2010smooth.)
plot.UKdata(mx2010..,mx2010smooth..)

#====================================================================================
### 4.2.1 Mortality modelling and extrapolation at advanced ages

#In this study, the univariate and multivariate state-space models are applied to UK male and female mortality data, 
#with emphasis on extrapolation at advanced ages and the impact of information sharing between populations. 
#Following \citet{forster22}, we consider both the complete dataset (ages 1--104) and a restricted dataset (ages 1--100), 
#illustrating how arbitrary termination ages may affect old-age mortality extrapolation under separate univariate fitting. 



#Following \citet{forster22}, we consider both the complete dataset (ages 1--104) and a restricted dataset (ages 1--100), 
#illustrating how arbitrary termination ages may affect old-age mortality extrapolation under separate univariate fitting.
#Figure 10 illustrates the fits resulting from different discounting strategies for the restricted data set (1-100) and 
#complete dataset (0-104), with unplausible crossing at older ages for male and female mortality curves. 
delta = c(0.8, 0.85, 0.9, 0.95)
## configuration discount factor - scenario 1: dataset (ages 1--104)
d12 = rep(c(0.99, 0.8, 0.85, 0.99), c(5, 30, 50, 19))
#### configuration discount factor - scenario 2:  dataset (ages 1--100)
d22 = rep(c(0.99, 0.8, 0.85, 0.99), c(5, 30, 50, 15))


source("fitUnivariate.R")
source("script(outputUnivariate).R")
source("predict_dlm.R")

#--------------------------------
## Figure 10: The univariate SSM fits extrapolated to age 120 with 95\% predictive intervals (shaded area) for
#fixed discount factors $\delta=(0.80,0.90,0.95)$ and age-varying factors
#($\delta_{1,1:5}=0.99$, $\delta_{2,6:35}=0.80$, $\delta_{3,36:85}=0.85$, $\delta_{4,86+}=0.99$),
#for the restricted dataset (0-100) and complete dataset (0-104) for male (blue) and female (red) populations.

mx.m2010<- mx %>%
  filter(year=='2010', sex=='1') %>%
  dplyr::  select(age,mx)

mx.f2010<- mx %>%
  filter(year=='2010', sex=='2') %>%
  dplyr::select(age,mx)


### data output - scenario 1 - male 2010
d1<-qx.m2010[[1]]; d2<-qx.m2010[[2]]; d3<-qx.m2010[[3]]; d4<-qx.m2010[[4]];
d5<- qx.m2010.d12
qx1.m.all<- bind_rows(d1,d3,d4,d5, .id="id")
### data output - scenario 2 - male 2010
d1<-qx2.m2010[[1]]; d2<-qx2.m2010[[2]]; d3<-qx2.m2010[[3]]; d4<-qx2.m2010[[4]];
d5<- qx2.m2010.d22
qx2.m.all<- bind_rows(d1,d3,d4,d5, .id="id")


### data output - scenario 1 - female 2010
d1<-qx.f2010[[1]]; d2<-qx.f2010[[2]]; d3<-qx.f2010[[3]]; d4<-qx.f2010[[4]];
d5<- qx.f2010.d12
qx1.f.all<- bind_rows(d1,d3,d4,d5, .id="id")
### data output - scenario 2 - female 2010
d1<-qx2.f2010[[1]]; d2<-qx2.f2010[[2]]; d3<-qx2.f2010[[3]]; d4<-qx2.f2010[[4]];
d5<- qx2.f2010.d22
qx2.f.all<- bind_rows(d1,d3,d4,d5, .id="id")

#> 80+
###scenario 1: 1-104 ano 2010
ggplot() +
  geom_point(data = mx.f2010 %>%
               filter(age>=80), aes(x = age, y = log(mx)), col = "tomato") +
  geom_point(data = mx.m2010 %>%
               filter(age>=80), aes(x = age, y = log(mx)), col = "steelblue", ) +
  geom_line(data=qx1.f.all %>%
              filter(age>=80), aes(x=age, y =  log(-log(1 - qx.fitted)), col = id), col="tomato") +
  geom_line(data=qx1.m.all %>%
              filter(age>=80), aes(x=age, y =  log(-log(1 - qx.fitted)), col = id), col="steelblue") +
  geom_ribbon(data = qx1.f.all %>%
                filter(age>=80), aes(x = age,  ymin = log(-log(1 - qx.lower)),
                                     ymax = log(-log(1 - qx.upper)), fill="tomato"), alpha = 0.25) +
  geom_ribbon(data = qx1.m.all %>%
                filter(age>=80), aes(x = age,  ymin = log(-log(1 - qx.lower)),
                                     ymax = log(-log(1 - qx.upper)), fill= "steelblue"), alpha = 0.25) +
  scale_color_manual(values=c("steelblue", "tomato"), labels=c("male", "female"))+
  scale_fill_manual(values=c("steelblue","tomato"), labels=c("male", "female")) +
  theme_classic(base_size = 20) +
  scale_y_continuous(expression("log rate")) +
  scale_x_continuous("Age", breaks = seq(80, 120, by = 20)) +
  theme(legend.position = c(0.95, 0.20),
        strip.background=element_rect(colour="black",
                                      fill="gray87"), panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.title=element_blank(),
        legend.key.height = unit(.6, "cm") ,
        legend.text=element_text(color="black", size=16),
        axis.title = element_text(color = "black", size = 16),
        axis.text = element_text(color="black",size=16))+
  facet_wrap(~id, ncol=5, labeller = labeller(id =
                                                c("1" = "0.80",
                                                  
                                                  "2" = "0.90",
                                                  "3" = "0.95",
                                                  "4" = "different per age")))


#> 80+
### scenario 2: 1-100 ano 2010
ggplot() + 

  geom_point(data = filter(mx.f2010 %>%
                             filter(age>=80), age %in% 1:100), aes(x = age, y = log(mx)), col = "tomato") +
  geom_point(data = filter(mx.m2010 %>%
                             filter(age>=80), age %in% 1:100), aes(x = age, y = log(mx)), col = "steelblue") +
  geom_point(data = filter(mx.f2010 %>%
                             filter(age>=80), age %in% 101:104), aes(x = age, y = log(mx)), col = "tomato", shape=1) +
  geom_point(data = filter(mx.m2010 %>%
                             filter(age>=80), age %in% 101:104), aes(x = age, y = log(mx)), col = "steelblue", shape=1) +
  geom_line(data=qx2.f.all %>%
              filter(age>=80), aes(x=age, y =  log(-log(1 - qx.fitted)), col = id), col="tomato") + 
  geom_line(data=qx2.m.all %>%
              filter(age>=80), aes(x=age, y =  log(-log(1 - qx.fitted)), col = id), col="steelblue") + 
  geom_ribbon(data = qx2.f.all %>%
                filter(age>=80), aes(x = age, ymin = log(-log(1 - qx.lower)),
                                     ymax = log(-log(1 - qx.upper)), fill="tomato"), alpha = 0.25) + 
  geom_ribbon(data = qx2.m.all %>%
                filter(age>=80), aes(x = age, ymin = log(-log(1 - qx.lower)),
                                     ymax = log(-log(1 - qx.upper)), fill= "steelblue"), alpha = 0.25) + 
  scale_color_manual(values=c("steelblue", "tomato"), labels=c("male", "female"))+
  scale_fill_manual(values=c("steelblue","tomato"), labels=c("male", "female")) + 
  theme_classic(base_size = 20) +
  scale_y_continuous(expression("log rate")) +
  scale_x_continuous("Age", breaks = seq(80, 120, by = 20)) + 
  theme(legend.position = c(0.95, 0.20),
        strip.background=element_rect(colour="black",
                                      fill="gray87"), panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.title=element_blank(),
        legend.key.height = unit(.6, "cm") ,
        legend.text=element_text(color="black", size=16),
        axis.title = element_text(color = "black", size = 16),
        axis.text = element_text(color="black",size=16))+
  facet_wrap(~id, ncol=5, labeller = labeller(id = 
                                                c("1" = "0.80",
                                                  
                                                  "2" = "0.90",
                                                  "3" = "0.95",
                                                  "4" = "different per age")))

#----------------------------------------

### Figure 11: The multivariate SSM fit extrapolated with 95\% predictive interval (shaded area) up to  the age of 120, 
#applying age-varying $\delta_x$ and  considering convergence between male and female populations. 
#The second and third columns zoom in on ages 30- and 80+, respectively.

source("ffbs_-_bivar.R")
source("predict_bivar2_att_chain.R")

qx_fitted = function(fit){
  n = dim(fit$mu)[1]
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

d12 <- rep(c(0.99, 0.8, 0.99), c(5, 85, 14))
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


### running multivariate 2010
 fitB_10<- gibbsV_corrigido(y = y_2010, delta = d12, m0 = m0, C0 = C0,
                            Ft0 = Ft, Gt = Gt, nit = 5000, v0 = v0, s0 = s0, V = cov(y_2010)) 

 ### predict up to age 120
 h = 16
 fitB.pred.10 <- predict_bivar2_chain(y_2010, m0, C0, h, V = fitB_10$V,
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
 
 
 
 qxB.all<- bind_rows(qxB.m1,qxB.f1,
.id="id")  
 
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
 
### 
## 1-104 - year 2010 - multivariate
ggplot() +
  geom_point(data = mx.%>% filter(id1==1), aes(x = age, y = log(mx), color=sex)) + #, shape=sex), col="gray") +
  geom_line(data=out.%>% filter(id1==1), aes(x=age, y = log(qx.fitted), color = sex)) +
  geom_ribbon(data = out.%>% filter(id1==1), aes(x = age, ymin = log(qx.lower),
                                                 ymax = log(qx.upper), fill=sex), alpha = 0.25) + 
  theme_classic(base_size = 20) + 
  scale_y_continuous(expression("log rate")) +
  scale_x_continuous("Age", breaks = seq(0, 200, by = 20)) + 
  scale_color_manual(values=c( "tomato","steelblue"), labels=c("female", "male"))+
  scale_fill_manual(values=c("tomato","steelblue"), labels=c("female", "male")) + 
  theme(legend.position = "none", strip.background=element_rect(colour="black", fill="gray87"), 
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1), 
        legend.title=element_blank(), legend.key.height = unit(.6, "cm") , 
        legend.text=element_text(color="black", size=16), 
        axis.title = element_text(color = "black", size = 14), 
        axis.text = element_text(color="black",size=14))+ 
  facet_wrap(~id1, ncol=1, labeller = labeller(id1 = c("1" = "age 1-104"
  ))) 


mx..<- mx. %>%
  filter(age<=30)
out..<- out. %>%
  filter(age<=30)

##### < 30
ggplot() +
  geom_point(data=mx..%>% filter(id1==1), aes(x=age, y =  log(mx), col=sex)) +
  geom_line(data=out..%>% filter(id1==1), aes(x=age, y = log(qx.fitted), color = sex)) +
  theme_classic(base_size = 20) +
  geom_hline(yintercept=1, linetype='dashed', col = 'black')+
  geom_ribbon(data = out..%>% filter(id1==1), aes(x = age, ymin = log(qx.lower), ymax = log(qx.upper), fill=sex), alpha = 0.25) +
  scale_y_continuous(expression("log rate"), limits = c(-10, -6)) +
  scale_x_continuous("Age", breaks = seq(0, 30, by = 10))+
  scale_color_manual(values=c( "tomato","steelblue"), labels=c("female", "male"))+
  scale_fill_manual(values=c("tomato","steelblue"), labels=c("female", "male")) +
  theme(legend.position = "none",
        strip.background=element_rect(colour="black",
                                      fill="gray87"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.title=element_blank(),
        legend.key.height = unit(.6, "cm") ,
        legend.text=element_text(color="black", size=16),
        axis.title = element_text(color = "black", size = 16),
        axis.text = element_text(color="black",size=16))+
  facet_wrap(~id1, ncol=1, labeller = labeller(id1 = c("1" = "age < 30"))) 

#### > 80+
mx..<- mx. %>%
  filter(age >=80)
out..<- out. %>%
  filter(age>=80)
ggplot() + 
  geom_point(data = mx..%>% filter(id1==1), aes(x = age, y = log(mx), color=sex)) + #, shape=sex), col="gray") +
  geom_line(data=out..%>% filter(id1==1), aes(x=age, y = log(qx.fitted), color = sex)) +
  geom_ribbon(data = out..%>% filter(id1==1), aes(x = age, ymin = log(qx.lower), ymax = log(qx.upper), fill=sex), alpha = 0.25) + 
  theme_classic(base_size = 20) + 
  scale_y_continuous(expression("log rate")) +
  scale_x_continuous("Age", breaks = seq(0, 200, by = 10)) + 
  scale_color_manual(values=c( "tomato","steelblue"), labels=c("female", "male"))+
  scale_fill_manual(values=c("tomato","steelblue"), labels=c("female", "male")) + 
  theme(legend.position = c(0.8,0.15), strip.background=element_rect(colour="black", fill="gray87"), 
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1), 
        legend.title=element_blank(), legend.key.height = unit(.6, "cm") , 
        legend.text=element_text(color="black", size=16), 
        axis.title = element_text(color = "black", size = 14), 
        axis.text = element_text(color="black",size=14))+ 
  facet_wrap(~id1, ncol=3, labeller = labeller(id1 = c("1" = "age 80+"))) 
graphics.off()



##------------------------------------------------------------
### 4.2.2. Joint modelling under missing data

#To evaluate the robustness of the proposed framework under incomplete information, 
#we consider a set of missing-data scenarios based on the England and Wales female mortality data for 2010 (ages 1--104).
#Missing observations are progressively introduced at younger ages, as summarised in Table~\ref{st3:tab1}. 
#The objective is to assess the ability of the joint model to recover plausible mortality patterns by borrowing information
#from related populations.

ytime=2010
df.mx<- mx %>% mutate(log.mx=log(mx))
y.m<- df.mx %>%
  filter(year==!!ytime, sex=='1') %>%
  select(log.mx) 
y.f<- df.mx %>% 
  filter(year==!!ytime, sex=='2') %>% 
  select(log.mx) 
y<- as.matrix(cbind(y.m,y.f)) 
yt= y ### no missing


## run models
source("fitallmodelsMissing.R") ### multivariate with and without common term and missing data
source("fitallmodelsNoMissing.R") ### multivariate with and without common term and no missing data


# For each scenario, the matrix y changes. Therefore, it is necessary to run 
# source("fitallmodelsMissing.R") separately for each scenario.
### recommendation to save the output for each scenario

## scenario a: 5% missing data for females 2010 
y[c(4:8),2] = NA
V = cov(y[-c(4:8),])

## scenario b: 10% missing data for females 2010 
y[c(4:10,15:17),2] = NA
V = cov(y[-c(4:10,15:17),])

## scenario c: 15% missing data for females 2010 
y[c(3:16),2]= NA
V = cov(y[-c(3:16),])

## scenario d: 25% missing data for females 2010 
y[c(1:25),2] = NA
V = cov(y[-c(1:25),])

## scenario e: 33% missing data female for females 2010 
y[c(1:16,23:41),2] = NA
V = cov(y[-c(1:16,23:41),])

## scenario f: 43% missing data for females 2010 
y[c(1:45),2] = NA
V = cov(y[-c(1:45),])

### Figure 12: Missing data in the female population under representative scenarios of increasing missingness. 
#SSM fits with age-varying discount factors $\delta_x$ and 95\% predictive intervals. 
#Empty dots indicate missing observations. Black dashed lines represent the fit obtained from the complete dataset.

## crude mortality rate 
mx.<- mx %>% 
  filter(year==!!ytime) %>%
  dplyr::select(age,mx,sex) 
mx.<- mx. %>%
  mutate(sex = ifelse(sex == 1, 'male', 'female'))
mx.. = mx. %>% 
  filter(age <= 30)

mx..f<- mx. %>%
  filter(sex=="female", age<=30)
mx..m<- mx. %>%
  filter(sex=="male", age<=30)

### fitted missing 
out<-bind_rows(qxBmiss.all,qxBTCmiss.all,
               .id="id1")
out = out %>% 
  mutate(sex = ifelse(id == 1, 'male', 'female')) 
out. = out %>% 
  mutate(sex = ifelse(id == 1, 'male', 'female')) %>%
  filter(age <=35)
out..<- out. %>%
  filter(sex=="female")

## fitted no missing 
outAll<-bind_rows(qxB.all,qxBTC.all,
                  .id="id1")
outAll = outAll %>% 
  mutate(sex = ifelse(id == 1, 'male', 'female')) 
outAll.. = outAll %>% 
  mutate(sex = ifelse(id == 1, 'male', 'female')) %>%
  filter(sex=='female', age <=35)

ggplot(NULL, aes(x = 1:30)) +
  geom_point(data = filter(mx..f, (age %in% 4:8)), aes(x = age, y = log(mx), color="female"), size=1.7, shape=1) +
  geom_point(data = filter(mx..f, !(age %in% 4:8)), aes(x = age, y = log(mx), color="female")) +
  geom_point(data = mx..m, aes(x = age, y = log(mx), color = "male"))+ #, linetype = "male")) +
  geom_line(data = out., aes(x = age, y = log(qx.fitted), color = sex)) +#, linetype = sex)) +
  geom_line(data = outAll.., aes(x = age, y = log(qx.fitted), color = "no missing"), linetype = "dashed",linewidth = 0.7) +
  geom_ribbon(data = out.., aes(x = age, ymin = log(qx.lower), ymax = log(qx.upper)), fill = "tomato", alpha = 0.25, inherit.aes = FALSE) +
  
  theme_classic(base_size = 20) +
  scale_y_continuous(expression("log rate"), labels = scales::comma) +
  scale_x_continuous("Age", breaks = seq(0, 40, by = 10)) +
  scale_color_manual(values = c("female" = "tomato", "male" = "steelblue", "no missing" = "black")) +
  scale_linetype_manual(values = c("female" = "solid", "male" = "solid", "no missing" = "dashed"), guide="none") +
  
  theme( legend.position = c(0.88, 0.2),
         strip.background = element_rect(colour = "black", fill = "gray87"),
         panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
         legend.title = element_blank(),
         legend.key.height = unit(.6, "cm"),
         legend.text = element_text(color = "black", size = 16),
         axis.title = element_text(color = "black", size = 16),
         axis.text = element_text(color = "black", size = 16)) +
  facet_wrap(~id1, ncol = 4, labeller = labeller(id1 = c(
    "1" = "Usual",
    "2" = "Common term")))


### Table 5: Model comparison based on Mean Square Prediction Error (MSPE), Mean Absolute Prediction Error (MAPE)
#and Width of Credibility Interval (WCI) for scenarios (a)-(f) a varying percentage of missing data for both competing models.

### for each scenario run 
out<- out %>%
  mutate(sex = ifelse(id == 1, 'male', 'female')) 

qx.missF.fit <- out %>%
  filter(sex== "female") %>%
  select(age, id1, qx.fitted) %>%
  spread(key = id1, value = qx.fitted)


qx.missF.fit025 <- out %>%
  filter(sex== "female") %>%
  select(age, id1, qx.lower) %>%
  spread(key = id1, value = qx.lower)

qx.missF.fit975 <- out %>%
  filter(sex== "female") %>%
  select(age, id1, qx.upper) %>%
  spread(key = id1, value = qx.upper)


df.mx. <- df.mx %>%
  mutate(qx= 1-exp(-mx))
ytime=2010
qxF<-  df.mx. %>% 
  filter(year==!!ytime, sex=='2') %>% 
  select(qx) 
qxM<-  df.mx. %>% 
  filter(year==!!ytime, sex=='1') %>% 
  select(qx) 

coverage.func <- function(y_real, lower, upper) {
  mean(y_real >= lower & y_real <= upper)
}

mae.func <- function(y_true, y_pred) {
  mean(abs(y_true - y_pred))
}



### female
mse<-NULL ; coverage<- NULL ; rmse<- NULL ; wci<- NULL ; mae<- NULL
for(j in 2:5){
  mse[j]= sum((qxF$qx - qx.missF.fit[,j])^2)/w
  mae[j]= mae.func(qxF$qx, qx.missF.fit[,j])
  rmse[j] = sqrt(mse[j])
  coverage[j]= coverage.func(qxF$qx, qx.missF.fit025[,j],qx.missF.fit975[,j])
  wci[j] = mean(qx.missF.fit975[,j] - qx.missF.fit025[,j])
  
}


### comparison
dfF.comp <- data.frame(
  MSE = round(mse[2:5], 5),
  MAE = round(mae[2:5], 5),
  WCI = round(wci[2:5], 5))#,
# IS = round(is, 5))
dfF.comp



### Figure 13: The multivariate SSM with common term considering missing data. Estimated female mortality in 2010 under 
#scenarios (c) and (f),
#using pooling with one auxiliary population (male, 2010) and with two
#auxiliary populations (male, 2010; female, 2012).
#SSM fits with age-varying discount factors $\delta_x$ and $\delta_{1,1:5}=0.9995$.
#Black dashed lines: fit using the complete dataset.
#Empty red dots: missing observations for females (2010).
#Shaded areas: 95\% predictive intervals.


### the new dataset
df.mx<- mx %>% mutate(log.mx=log(mx))
y.m2010<- df.mx %>%
  filter(year==!!2010, sex=='1') %>%
  select(log.mx) 
y.f2010<- df.mx %>% 
  filter(year==!!2010, sex=='2') %>% 
  select(log.mx) 
y.f2012<- df.mx %>% 
  filter(year==!!2012, sex=='2') %>% 
  select(log.mx) 
y<- as.matrix(cbind(y.m2010,y.f2010, y.f2012)) 
yt= y # no missing data


## for scenario c
y[c(3:16),2] = NA
V = cov(y[-c(3:16),])
V
## and again run source("fitallmodelsMissing.R") and source("fitallmodelsNoMissing.R")
outAll.f <- outAll %>%
  filter(sex=="female", age <=35)

ggplot(NULL, aes(x = 1:35)) + 
  geom_point(data = filter(mx..f, (age %in% c(3:16))), aes(x = age, y = log(mx), color="female"), size=1.7, shape=1) +
  geom_point(data = filter(mx..f, !(age %in% c(3:16))), aes(x = age, y = log(mx), color="female")) +#, color="tomato") +
  geom_point(data = mx..m, aes(x = age, y = log(mx), color="male")) + #, color="steelblue") +
  geom_point(data = mx..f2, aes(x = age, y = log(mx), color= "female2")) + #, color="orange") +
  geom_line(data=out.missAll., aes(x=age, y = log(qx.fitted), color = sex)) +
  geom_line(data = outAll.f, aes(x = age, y = log(qx.fitted),color = "no missing" ),linetype = "dashed",linewidth = 0.7) +
  geom_ribbon(data = out..missf, aes(x = age, ymin = log(qx.lower), ymax = log(qx.upper), fill=sex), alpha = 0.25) +
  theme_classic(base_size = 20) + 
  scale_y_continuous(expression("log rate")) + 
  scale_x_continuous("Age", breaks = seq(0, 40, by = 10)) + 
  guides( fill = "none") +
  scale_color_manual(values=c( "tomato", "orange","steelblue", "black"), labels=c("female (2010)", "female (2012)","male (2010)", "no missing"))+
  scale_fill_manual(values=c("tomato","orange","steelblue", "black"), labels=c("female (2010)","female (2012)", "male (2010)", "no missing")) +
  scale_linetype_manual(values = c("female (2010)" = "solid", "female (2012)" = "solid", "male (2010)" = "solid", "no missing" = "dashed"), guide="none") +
  
  theme(legend.position = c(0.88,0.18), strip.background=element_rect(colour="black", fill="gray87"), 
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1), 
        legend.title=element_blank(), legend.key.height = unit(.6, "cm") , 
        legend.text=element_text(color="black", size=16), 
        axis.title = element_text(color = "black", size = 16), 
        axis.text = element_text(color="black",size=16))+ 
  facet_wrap(~id1, ncol=4, labeller = labeller(id1 = c("1" = "Usual",
                                                       "2" = "Common term")))


## for scenario f
y[c(1:45),2] = NA
V = cov(y[-c(1:45),])
## and again run source("fitallmodelsMissing.R") and source("fitallmodelsNoMissing.R")
out.missAll. <- out.missAll %>%
  filter(age <=60)

outAll. <- outAll %>%
  filter(age <=60)

mx.. = mxtot. %>% 
  filter(age <= 60)

mx..f<- mxtot. %>%
  filter(sex=="female", age<=60)
mx..m<- mxtot. %>%
  filter(sex=="male", age<=60)
mx..f2<- mxtot. %>%
  filter(sex=="female2", age<=60)

out..missf<- out.missAll %>%
  filter(sex=="female", age <=60)

outAll.f <- outAll %>%
  filter(sex=="female", age <=60)


ggplot(NULL, aes(x = 1:60)) + 
  geom_point(data = filter(mx..f, (age %in% c(1:45))), aes(x = age, y = log(mx), color="female"), size=1.7, shape=1) +
  geom_point(data = filter(mx..f, !(age %in% c(1:45))), aes(x = age, y = log(mx), color="female")) +
  geom_point(data = mx..m, aes(x = age, y = log(mx), color="male")) + #, color="steelblue") +
  geom_point(data = mx..f2, aes(x = age, y = log(mx), color="female2"))+ #, color="orange") +
  geom_line(data=out.missAll., aes(x=age, y = log(qx.fitted), color = sex)) +
  geom_line(data=outAll.f, aes(x=age, y = log(qx.fitted), color = "no missing"), linetype="dashed", linewidth=0.7) +
  geom_ribbon(data = out..missf, aes(x = age, ymin = log(qx.lower), ymax = log(qx.upper), fill=sex), alpha = 0.25) +
  theme_classic(base_size = 20) + 
  scale_y_continuous(expression("log rate")) + 
  scale_x_continuous("Age", breaks = seq(0, 60, by = 10)) + 
  guides( fill = "none") +
  scale_color_manual(values=c( "tomato", "orange","steelblue", "black"), labels=c("female (2010)", "female (2012)","male (2010)", "no missing"))+
  scale_fill_manual(values=c("tomato","orange","steelblue", "black"), labels=c("female (2010)","female (2012)", "male (2010)", "no missing")) + 
  theme(legend.position = c(0.88,0.18), strip.background=element_rect(colour="black", fill="gray87"), 
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1), 
        legend.title=element_blank(), legend.key.height = unit(.6, "cm") , 
        legend.text=element_text(color="black", size=16), 
        axis.title = element_text(color = "black", size = 16), 
        axis.text = element_text(color="black",size=16))+ 
  facet_wrap(~id1, ncol=4, labeller = labeller(id1 = c("1" = "Usual",
                                                       "2" = "Common term")))


