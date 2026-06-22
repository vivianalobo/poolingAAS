######################################################
#### output: univariate modelling 
### data: England + Wales, male and female, 2010-2012
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



#####==========#####
##### dataset  #####
#####==========#####
dx <- read.csv("Results/deaths_MYB2.csv")
Ex <- read.csv("Results/MYB2_mais_OLDAGE.csv")

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


#####===============================================#####
##### modelagem univariada dos sexos por ano a ano  #####
#####===============================================#####


### Cenario 1: Age 1-104
### Cenario 2: Age 1-100


### Dados do cenario 1
df.mx<- mx %>%
     mutate(log.mx=log(mx))

y.2010m<- df.mx %>%
     filter(year=='2010', sex=='1') %>%
  dplyr:: select(age,log.mx)
y.2010f<- df.mx %>%
  filter(year=='2010', sex=='2') %>%
  dplyr:: select(age, log.mx)
y.2011m<- df.mx %>%
  filter(year=='2011', sex=='1') %>%
  dplyr:: select(age,log.mx)
y.2011f<- df.mx %>%
  filter(year=='2011', sex=='2') %>%
  dplyr::  select(age, log.mx)
y.2012m<- df.mx %>%
  filter(year=='2012', sex=='1') %>%
  dplyr::  select(age,log.mx)
y.2012f<- df.mx %>%
  filter(year=='2012', sex=='2') %>%
  dplyr::  select(age, log.mx)

## Cenario 1 e delta variando
delta = c(0.8, 0.85, 0.9, 0.95)
## 0.99 de 1-5 anos; 0.8 de 6-90 anos; 0.99 de 91+
d12 = rep(c(0.99, 0.8, 0.85, 0.99), c(5, 30, 50, 19))
#### Cenario 2 e delta variando
d22 = rep(c(0.99, 0.8, 0.85, 0.99), c(5, 30, 50, 15)) ## Ate 100


########################################################################
### output do ajuste univariado - ajuste com extrapolacao até 120 anos

### nao precisa rodar, ja tem um .RData abaixo salvo os ajustes
# source("predict_dlm.R")
# source("fitUnivariate.R")
#save.image("outUnivariate.RData")
### resultados do ajuste
load("~/Dropbox/Semestre 2024.1 UFRJ/Projeto Multivariado (viviana)/outUnivariate.RData")

#### Fig2a ; Fig3a ; Fig4a - male 
#### Fig2b; Fig3b; Fig4b - female 
#### Fig5a; Fig5b ; Fig5c - comparison male x female scenario 1
#### Fig 6a; Fig6b; Fig6c - comparion male x female scenario 2

qx.m2010 ; qx.f2010 ; qx.m2011; qx.f2011; qx.m2012;qx.f2012
qx.m2010.d12 ; qx.f2010.d12 ; qx.m2011.d12; qx.f2011.d12; qx.m2012.d12;qx.f2012.d12

qx2.m2010 ; qx2.f2010 ; qx2.m2011; qx2.f2011; qx2.m2012;qx2.f2012
qx2.m2010.d22 ; qx2.f2010.d22 ; qx2.m2011.d22; qx2.f2011.d22; qx2.m2012.d22;qx2.f2012.d22


#################################################################
##### Rodar os resultados a partir daqui para gerar graficos
### data output cenario 1 - masculino 2010
d1<-qx.m2010[[1]]; d2<-qx.m2010[[2]]; d3<-qx.m2010[[3]]; d4<-qx.m2010[[4]];
d5<- qx.m2010.d12
qx1.m.all<- bind_rows(d1,d2,d3,d4,d5, .id="id")
mx.m2010<- mx %>%
  filter(year=='2010', sex=='1') %>%
  dplyr::  select(age,mx)
### data output cenario 2 - masculino 2010
d1<-qx2.m2010[[1]]; d2<-qx2.m2010[[2]]; d3<-qx2.m2010[[3]]; d4<-qx2.m2010[[4]];
d5<- qx2.m2010.d22
qx2.m.all<- bind_rows(d1,d2,d3,d4,d5, .id="id")

### data output cenario 1 - feminino 2010
d1<-qx.f2010[[1]]; d2<-qx.f2010[[2]]; d3<-qx.f2010[[3]]; d4<-qx.f2010[[4]];
d5<- qx.f2010.d12
qx1.f.all<- bind_rows(d1,d2,d3,d4,d5, .id="id")
mx.f2010<- mx %>%
  filter(year=='2010', sex=='2') %>%
  dplyr::select(age,mx)
### data output cenario 2 - feminino 2010
d1<-qx2.f2010[[1]]; d2<-qx2.f2010[[2]]; d3<-qx2.f2010[[3]]; d4<-qx2.f2010[[4]];
d5<- qx2.f2010.d22
qx2.f.all<- bind_rows(d1,d2,d3,d4,d5, .id="id")

### cenario 1 versus cenario 2 - 2010 masculino 
pdf("Fig2a.pdf", width=16, height=5)
ggplot() + 
  geom_point(data = filter(mx.m2010, age %in% 1:100), aes(x = age, y = mx), col = "steelblue") +
  geom_point(data = filter(mx.m2010, age %in% 101:104), aes(x = age, y = mx), col = "steelblue", shape=1) +
  geom_line(data=qx1.m.all, aes(x=age, y =  -log(1-qx.fitted), color = id, lty = "1-104"),col = "steelblue") + 
  geom_line(data=qx2.m.all, aes(x=age, y =  -log(1-qx.fitted), color = id, lty = "1-100"),col = "steelblue") + 
  
 # geom_hline(yintercept=1, linetype='dashed', col = 'black')+
 
   geom_ribbon(data = qx1.m.all, aes(x = age, ymin = -log(1-qx.lower), ymax = -log(1-qx.upper), fill="steelblue"), alpha = 0.25) + 
  geom_ribbon(data = qx2.m.all, aes(x = age, ymin = -log(1- qx.lower), ymax = -log(1- qx.upper), fill="steelblue"), alpha = 0.25) + 

    scale_fill_manual(values=c("steelblue", "steelblue"), guide="none")+
  theme_classic(base_size = 20) + 
   scale_y_continuous(expression(m[x]), limits = 10^-c(NA,NA), 
                       trans = 'log10', labels = scales::comma) +
    scale_x_continuous("Age", breaks = seq(0, 200, by = 20)) +
  scale_color_manual(id, values=c(rep("steelblue",5)), guide="none")+
  theme(legend.position = c(0.95, 0.2),
        strip.background=element_rect(colour="black",
        fill="gray87"), panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.title=element_blank(),
        legend.key.height = unit(.6, "cm") ,
        legend.text=element_text(color="black", size=16),
        axis.title = element_text(color = "black", size = 14),
        axis.text = element_text(color="black",size=14))+
        facet_wrap(~id, ncol=5, labeller = labeller(id = 
                                                c("1" = "0.80",
                                                  "2" = "0.85",
                                                  "3" = "0.90",
                                                  "4" = "0.95",
                                                  "5" = "different per age")))
graphics.off()

### colocando o grafico com zoom a partir de 80 anos

mx.f2010.<- mx.f2010 %>%
  filter(age >=80)
mx.m2010.<- mx.m2010 %>%
  filter(age >=80)
qx1.m.all. <- qx1.m.all %>%
  filter(age >=80)
qx1.f.all. <- qx1.f.all %>%
  filter(age >=80)
qx2.m.all. <- qx2.m.all %>%
  filter(age >=80)
qx2.f.all. <- qx2.f.all %>%
  filter(age >=80)

pdf("Fig2aOlder.pdf", width=16, height=5)
ggplot() + 
  # geom_point(data = mx.m2010, aes(x = age, y = mx), col = "gray") +
  geom_point(data = filter(mx.m2010., age %in% 80:100), aes(x = age, y = mx), col = "steelblue") +
  geom_point(data = filter(mx.m2010., age %in% 80:104), aes(x = age, y = mx), col = "steelblue", shape=1) +
  geom_line(data=qx1.m.all., aes(x=age, y =  -log(1-qx.fitted), color = id, lty = "1-104"),col = "steelblue") + 
  geom_line(data=qx2.m.all., aes(x=age, y =  -log(1-qx.fitted), color = id, lty = "1-100"),col = "steelblue") + 
  
#  geom_hline(yintercept=1, linetype='dashed', col = 'black')+
  
  geom_ribbon(data = qx1.m.all., aes(x = age, ymin = -log(1-qx.lower), ymax = -log(1-qx.upper), fill="steelblue"), alpha = 0.25) + 
  geom_ribbon(data = qx2.m.all., aes(x = age, ymin = -log(1- qx.lower), ymax = -log(1- qx.upper), fill="steelblue"), alpha = 0.25) + 
  
  scale_fill_manual(values=c("steelblue", "steelblue"), guide="none")+
  theme_classic(base_size = 20) + 
  # scale_y_continuous(expression(q[x]), limits = 10^-c(NA,NA), 
  #                    trans = 'log10', labels = scales::comma) +
  # scale_x_continuous("Age", breaks = seq(0, 200, by = 20)) +
  scale_y_continuous(expression(m[x]), limits = 10^-c(NA,NA), 
                     trans = 'log10', labels = scales::comma) +
  scale_x_continuous("Age", breaks = seq(80, 200, by = 10))+
  scale_color_manual(id, values=c(rep("steelblue",5)), guide="none")+
  theme(legend.position = c(0.95, 0.2),
        strip.background=element_rect(colour="black",
                                      fill="gray87"), panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.title=element_blank(),
        legend.key.height = unit(.6, "cm") ,
        legend.text=element_text(color="black", size=16),
        axis.title = element_text(color = "black", size = 14),
        axis.text = element_text(color="black",size=14))+
  facet_wrap(~id, ncol=5, labeller = labeller(id = 
                                                c("1" = "0.80",
                                                  "2" = "0.85",
                                                  "3" = "0.90",
                                                  "4" = "0.95",
                                                  "5" = "different per age")))
graphics.off()



pdf("Fig2b.pdf", width=16, height=5)
ggplot() + 
 # geom_point(data = mx.f2010, aes(x = age, y = mx), col = "gray") +
  geom_point(data = filter(mx.f2010, age %in% 1:100), aes(x = age, y = mx), col = "tomato") +
  geom_point(data = filter(mx.f2010, age %in% 101:104), aes(x = age, y = mx), col = "tomato", shape=1) +
  geom_line(data=qx1.f.all, aes(x=age, y =  -log(1-qx.fitted), col = id, lty = "1-104"), col="tomato") + 
  geom_line(data=qx2.f.all, aes(x=age, y =  -log(1-qx.fitted), col = id, lty = "1-100"), col="tomato") + 
#  geom_hline(yintercept=1, linetype='dashed', col = 'black')+
  geom_ribbon(data = qx1.f.all, aes(x = age, ymin = -log(1-qx.lower), ymax = -log(1-qx.upper), fill="tomato"), alpha = 0.25) + 
  geom_ribbon(data = qx2.f.all, aes(x = age, ymin = -log(1- qx.lower), ymax = -log(1- qx.upper), fill="tomato"), alpha = 0.25) + 
  scale_fill_manual(values=c("tomato", "tomato"), guide="none")+
   theme_classic(base_size = 20) + 
  scale_y_continuous(expression(m[x]), limits = 10^-c(NA,NA), trans = 'log10', labels = scales::comma) +
  scale_x_continuous("Age", breaks = seq(0, 200, by = 20)) +
  theme(legend.position = c(0.95, 0.20),
        strip.background=element_rect(colour="black",
                                      fill="gray87"), panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.title=element_blank(),
        legend.key.height = unit(.6, "cm") ,
        legend.text=element_text(color="black", size=16),
        axis.title = element_text(color = "black", size = 14),
        axis.text = element_text(color="black",size=14))+
   facet_wrap(~id, ncol=5, labeller = labeller(id = 
                                                c("1" = "0.80",
                                                  "2" = "0.85",
                                                  "3" = "0.90",
                                                  "4" = "0.95",
                                                  "5" = "different per age")))
graphics.off()


pdf("Fig2bOlder.pdf", width=16, height=5)
ggplot() + 
  # geom_point(data = mx.m2010, aes(x = age, y = mx), col = "gray") +
  geom_point(data = filter(mx.f2010., age %in% 80:100), aes(x = age, y = mx), col = "tomato") +
  geom_point(data = filter(mx.f2010., age %in% 80:104), aes(x = age, y = mx), col = "tomato", shape=1) +
  geom_line(data=qx1.f.all., aes(x=age, y =  -log(1-qx.fitted), color = id, lty = "1-104"),col = "tomato") + 
  geom_line(data=qx2.f.all., aes(x=age, y =  -log(1-qx.fitted), color = id, lty = "1-100"),col = "tomato") + 
  
  #geom_hline(yintercept=1, linetype='dashed', col = 'black')+
  
  geom_ribbon(data = qx1.f.all., aes(x = age, ymin = -log(1-qx.lower), ymax = -log(1-qx.upper), fill="tomato"), alpha = 0.25) + 
  geom_ribbon(data = qx2.f.all., aes(x = age, ymin = -log(1- qx.lower), ymax = -log(1- qx.upper), fill="tomato"), alpha = 0.25) + 
  
  scale_fill_manual(values=c("tomato", "tomato"), guide="none")+
  theme_classic(base_size = 20) + 
  # scale_y_continuous(expression(q[x]), limits = 10^-c(NA,NA), 
  #                    trans = 'log10', labels = scales::comma) +
  # scale_x_continuous("Age", breaks = seq(0, 200, by = 20)) +
  scale_y_continuous(expression(m[x]), limits = 10^-c(NA,NA), 
                     trans = 'log10', labels = scales::comma) +
  scale_x_continuous("Age", breaks = seq(80, 200, by = 10))+
  scale_color_manual(id, values=c(rep("tomato",5)), guide="none")+
  theme(legend.position = c(0.95, 0.2),
        strip.background=element_rect(colour="black",
                                      fill="gray87"), panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.title=element_blank(),
        legend.key.height = unit(.6, "cm") ,
        legend.text=element_text(color="black", size=16),
        axis.title = element_text(color = "black", size = 14),
        axis.text = element_text(color="black",size=14))+
  facet_wrap(~id, ncol=5, labeller = labeller(id = 
                                                c("1" = "0.80",
                                                  "2" = "0.85",
                                                  "3" = "0.90",
                                                  "4" = "0.95",
                                                  "5" = "different per age")))
graphics.off()


#### zoomed
mx.m2010.<- mx %>%
  filter(year=='2010', sex=='1') %>%
  dplyr::select(age,mx) %>%
  filter(age>=80)

mx.f2010.<- mx %>%
  
  filter(year=='2010', sex=='2') %>%
  dplyr::select(age,mx) %>%
  filter(age>=80)
qx1.f.all.<- qx1.f.all %>%
  filter(age>=80)

qx2.f.all.<- qx2.f.all %>%
  filter(age>=80)

qx1.m.all.<- qx1.m.all %>%
  filter(age>=80)

qx2.m.all.<- qx2.m.all %>%
  filter(age>=80)

### 1-104 ano 2010
pdf("Fig5a.pdf", width=16, height=5)
ggplot() + 
  geom_point(data = mx.f2010., aes(x = age, y = mx), col = "tomato") +
  geom_point(data = mx.m2010., aes(x = age, y = mx), col = "steelblue", ) +
  geom_line(data=qx1.f.all., aes(x=age, y =  -log(1-qx.fitted), col = id), col="tomato") + 
  geom_line(data=qx1.m.all., aes(x=age, y =  -log(1-qx.fitted), col = id), col="steelblue") + 
  #geom_hline(yintercept=1, linetype='dashed', col = 'black')+
  geom_ribbon(data = qx1.f.all., aes(x = age, ymin = -log(1-qx.lower), ymax = -log(1-qx.upper), fill="tomato"), alpha = 0.25) + 
  geom_ribbon(data = qx1.m.all., aes(x = age, ymin = -log(1-qx.lower), ymax = -log(1-qx.upper), fill= "steelblue"), alpha = 0.25) + 
#  scale_fill_manual(values=c("steelblue", "tomato"), guide="none")+
  scale_color_manual(values=c("steelblue", "tomato"), labels=c("male", "female"))+
  scale_fill_manual(values=c("steelblue","tomato"), labels=c("male", "female")) + 
  theme_classic(base_size = 20) + 
  scale_y_continuous(expression(m[x]), limits = 10^-c(NA,NA), 
                                                     trans = 'log10', labels = scales::comma) + 
  scale_x_continuous("Age", breaks = seq(80, 200, by = 20)) + 
  theme(legend.position = c(0.95, 0.20),
        strip.background=element_rect(colour="black",
                                      fill="gray87"), panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.title=element_blank(),
        legend.key.height = unit(.6, "cm") ,
        legend.text=element_text(color="black", size=16),
        axis.title = element_text(color = "black", size = 14),
        axis.text = element_text(color="black",size=14))+
  facet_wrap(~id, ncol=5, labeller = labeller(id = 
                                                c("1" = "0.80",
                                                  "2" = "0.85",
                                                  "3" = "0.90",
                                                  "4" = "0.95",
                                                  "5" = "different per age")))


graphics.off()


#### zoomed
mx.f2010.<- mx.f2010 %>%
  filter(age <=40)
mx.m2010.<- mx.m2010 %>%
  filter(age <=40)
qx1.m.all. <- qx1.m.all %>%
  filter(age <=40)
qx1.f.all. <- qx1.f.all %>%
  filter(age <=40)
qx2.m.all. <- qx2.m.all %>%
  filter(age<=40)
qx2.f.all. <- qx2.f.all %>%
  filter(age <=40)

pdf("Fig2bYoung.pdf", width=16, height=5)
ggplot() + 
  # geom_point(data = mx.m2010, aes(x = age, y = mx), col = "gray") +
  geom_point(data = mx.f2010., aes(x = age, y = mx), col = "tomato") +
#  geom_point(data = filter(mx.f2010., age %in% 80:104), aes(x = age, y = mx), col = "tomato", shape=1) +
  geom_line(data=qx1.f.all., aes(x=age, y =  -log(1-qx.fitted), color = id, lty = "1-104"),col = "tomato") + 
  geom_line(data=qx2.f.all., aes(x=age, y =  -log(1-qx.fitted), color = id, lty = "1-100"),col = "tomato") + 
 # geom_hline(yintercept=1, linetype='dashed', col = 'black')+
  geom_ribbon(data = qx1.f.all., aes(x = age, ymin = -log(1-qx.lower), ymax = -log(1-qx.upper), fill="tomato"), alpha = 0.25) + 
  geom_ribbon(data = qx2.f.all., aes(x = age, ymin = -log(1- qx.lower), ymax = -log(1- qx.upper), fill="tomato"), alpha = 0.25) + 
  
  scale_fill_manual(values=c("tomato", "tomato"), guide="none")+
  theme_classic(base_size = 20) + 
  scale_y_continuous(expression(m[x]), limits = 10^-c(NA,NA), 
                     trans = 'log10', labels = scales::comma) +
  scale_x_continuous("Age", breaks = seq(0, 40, by = 10))+
  scale_color_manual(id, values=c(rep("tomato",5)), guide="none")+
  theme(legend.position = c(0.95, 0.2),
        strip.background=element_rect(colour="black",
                                      fill="gray87"), panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.title=element_blank(),
        legend.key.height = unit(.6, "cm") ,
        legend.text=element_text(color="black", size=16),
        axis.title = element_text(color = "black", size = 14),
        axis.text = element_text(color="black",size=14))+
  facet_wrap(~id, ncol=5, labeller = labeller(id = 
                                                c("1" = "0.80",
                                                  "2" = "0.85",
                                                  "3" = "0.90",
                                                  "4" = "0.95",
                                                  "5" = "different per age")))
graphics.off()

pdf("Fig2aYoung.pdf", width=16, height=5)
ggplot() + 
  # geom_point(data = mx.m2010, aes(x = age, y = mx), col = "gray") +
  geom_point(data = mx.m2010., aes(x = age, y = mx), col = "steelblue") +
  #  geom_point(data = filter(mx.f2010., age %in% 80:104), aes(x = age, y = mx), col = "tomato", shape=1) +
  geom_line(data=qx1.m.all., aes(x=age, y =  -log(1-qx.fitted), color = id, lty = "1-104"),col = "steelblue") + 
  geom_line(data=qx2.m.all., aes(x=age, y =  -log(1-qx.fitted), color = id, lty = "1-100"),col = "steelblue") + 
  # geom_hline(yintercept=1, linetype='dashed', col = 'black')+
  geom_ribbon(data = qx1.m.all., aes(x = age, ymin = -log(1-qx.lower), ymax = -log(1-qx.upper), fill="steelblue"), alpha = 0.25) + 
  geom_ribbon(data = qx2.m.all., aes(x = age, ymin = -log(1- qx.lower), ymax = -log(1- qx.upper), fill="steelblue"), alpha = 0.25) + 
  
  scale_fill_manual(values=c("steelblue", "steelblue"), guide="none")+
  theme_classic(base_size = 20) + 
  scale_y_continuous(expression(m[x]), limits = 10^-c(NA,NA), 
                     trans = 'log10', labels = scales::comma) +
  scale_x_continuous("Age", breaks = seq(0, 40, by = 10))+
  scale_color_manual(id, values=c(rep("steelblue",5)), guide="none")+
  theme(legend.position = c(0.95, 0.2),
        strip.background=element_rect(colour="black",
                                      fill="gray87"), panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.title=element_blank(),
        legend.key.height = unit(.6, "cm") ,
        legend.text=element_text(color="black", size=16),
        axis.title = element_text(color = "black", size = 14),
        axis.text = element_text(color="black",size=14))+
  facet_wrap(~id, ncol=5, labeller = labeller(id = 
                                                c("1" = "0.80",
                                                  "2" = "0.85",
                                                  "3" = "0.90",
                                                  "4" = "0.95",
                                                  "5" = "different per age")))
graphics.off()

ggplot() +
  geom_point(data = mx.f2010., aes(x = age, y = mx), col = "tomato") +
  geom_point(data = mx.m2010., aes(x = age, y = mx), col = "steelblue", ) +
  geom_line(data=qx1.f.all., aes(x=age, y =  -log(1-qx.fitted), col = id), col="tomato") +
  geom_line(data=qx1.m.all., aes(x=age, y =  -log(1-qx.fitted), col = id), col="steelblue") +
  #geom_hline(yintercept=1, linetype='dashed', col = 'black')+
  geom_ribbon(data = qx1.f.all., aes(x = age, ymin = -log(1-qx.lower), ymax = -log(1-qx.upper), fill="tomato"), alpha = 0.25) +
  geom_ribbon(data = qx1.m.all., aes(x = age, ymin = -log(1-qx.lower), ymax = -log(1-qx.upper), fill= "steelblue"), alpha = 0.25) +
  #  scale_fill_manual(values=c("steelblue", "tomato"), guide="none")+
  scale_color_manual(values=c("steelblue", "tomato"), labels=c("male", "female"))+
  scale_fill_manual(values=c("steelblue","tomato"), labels=c("male", "female")) +
  theme_classic(base_size = 20) +
  scale_y_continuous(expression(m[x]), limits = 10^-c(NA,NA),
                     trans = 'log10', labels = scales::comma) +
  scale_x_continuous("Age", breaks = seq(0, 40, by = 20)) +
  theme(legend.position = c(0.95, 0.20),
        strip.background=element_rect(colour="black",
                                      fill="gray87"), panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.title=element_blank(),
        legend.key.height = unit(.6, "cm") ,
        legend.text=element_text(color="black", size=16),
        axis.title = element_text(color = "black", size = 14),
        axis.text = element_text(color="black",size=14))+
  facet_wrap(~id, ncol=5, labeller = labeller(id =
                                                c("1" = "0.80",
                                                  "2" = "0.85",
                                                  "3" = "0.90",
                                                  "4" = "0.95",
                                                  "5" = "different per age")))
# 


### 1-100 ano 2010
pdf("Fig6a.pdf", width=16, height=5)
ggplot() + 
 # geom_point(data = mx.f2010., aes(x = age, y = mx), col = "tomato") +
#  geom_point(data = mx.m2010., aes(x = age, y = mx), col = "steelblue", ) +
  geom_point(data = filter(mx.f2010., age %in% 1:100), aes(x = age, y = mx), col = "tomato") +
  geom_point(data = filter(mx.m2010., age %in% 1:100), aes(x = age, y = mx), col = "steelblue") +
  geom_point(data = filter(mx.f2010., age %in% 101:104), aes(x = age, y = mx), col = "tomato", shape=1) +
  geom_point(data = filter(mx.m2010., age %in% 101:104), aes(x = age, y = mx), col = "steelblue", shape=1) +
  geom_line(data=qx2.f.all., aes(x=age, y =  -log(1-qx.fitted), col = id), col="tomato") + 
  geom_line(data=qx2.m.all., aes(x=age, y =  -log(1-qx.fitted), col = id), col="steelblue") + 
 
 # geom_hline(yintercept=1, linetype='dashed', col = 'black')+
  geom_ribbon(data = qx2.f.all., aes(x = age, ymin = -log(1-qx.lower), ymax = -log(1-qx.upper), fill="tomato"), alpha = 0.25) + 
  geom_ribbon(data = qx2.m.all., aes(x = age, ymin = -log(1-qx.lower), ymax = -log(1-qx.upper), fill= "steelblue"), alpha = 0.25) + 
  #scale_fill_manual(values=c("steelblue", "tomato"), guide="none")+
  scale_color_manual(values=c("steelblue", "tomato"), labels=c("male", "female"))+
  scale_fill_manual(values=c("steelblue","tomato"), labels=c("male", "female")) + 
  
  theme_classic(base_size = 20) +
  scale_y_continuous(expression(m[x]), limits = 10^-c(NA,NA),                                                    trans = 'log10', labels = scales::comma) + 
  scale_x_continuous("Age", breaks = seq(80, 200, by = 20)) + 
  theme(legend.position = c(0.95, 0.20),
        strip.background=element_rect(colour="black",
                                      fill="gray87"), panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.title=element_blank(),
        legend.key.height = unit(.6, "cm") ,
        legend.text=element_text(color="black", size=16),
        axis.title = element_text(color = "black", size = 14),
        axis.text = element_text(color="black",size=14))+
  facet_wrap(~id, ncol=5, labeller = labeller(id = 
                                                c("1" = "0.80",
                                                  "2" = "0.85",
                                                  "3" = "0.90",
                                                  "4" = "0.95",
                                                  "5" = "different per age")))


graphics.off()

#-------------

### data output cenario 1 - masculino 2011
d1<-qx.m2011[[1]]; d2<-qx.m2011[[2]]; d3<-qx.m2011[[3]]; d4<-qx.m2011[[4]];
d5<- qx.m2011.d12
qx1.m.all<- bind_rows(d1,d2,d3,d4,d5, .id="id")
mx.m2011<- mx %>%
  filter(year=='2011', sex=='1') %>%
  dplyr::  select(age,mx)
### data output cenario 2 - masculino 2011
d1<-qx2.m2011[[1]]; d2<-qx2.m2011[[2]]; d3<-qx2.m2011[[3]]; d4<-qx2.m2011[[4]];
d5<- qx2.m2011.d22
qx2.m.all<- bind_rows(d1,d2,d3,d4,d5, .id="id")


### data output cenario 1 - feminino 2011
d1<-qx.f2011[[1]]; d2<-qx.f2011[[2]]; d3<-qx.f2011[[3]]; d4<-qx.f2011[[4]];
d5<- qx.f2011.d12
qx1.f.all<- bind_rows(d1,d2,d3,d4,d5, .id="id")
mx.f2011<- mx %>%
  filter(year=='2011', sex=='2') %>%
  dplyr::  select(age,mx)
### data output cenario 2 - feminino 2011
d1<-qx2.f2011[[1]]; d2<-qx2.f2011[[2]]; d3<-qx2.f2011[[3]]; d4<-qx2.f2011[[4]];
d5<- qx2.f2011.d22
qx2.f.all<- bind_rows(d1,d2,d3,d4,d5, .id="id")


### cenario 1 versus cenario 2 - 2011 masculino 
pdf("Fig3a.pdf", width=16, height=5)
ggplot() + 
 # geom_point(data = mx.m2011, aes(x = age, y = mx), col = "gray") +
  geom_point(data = filter(mx.m2011, age %in% 1:100), aes(x = age, y = mx), col = "steelblue") +
  geom_point(data = filter(mx.m2011, age %in% 101:104), aes(x = age, y = mx), col = "steelblue", shape=1) +
  
  geom_line(data=qx1.m.all, aes(x=age, y =  -log(1-qx.fitted), color = id, lty = "1-104"),col = "steelblue") + 
  geom_line(data=qx2.m.all, aes(x=age, y =  -log(1-qx.fitted), color = id, lty = "1-100"),col = "steelblue") + 
 # geom_hline(yintercept=1, linetype='dashed', col = 'black')+
  geom_ribbon(data = qx1.m.all, aes(x = age, ymin = -log(1-qx.lower), ymax = -log(1-qx.upper), fill= "steelblue"), alpha = 0.25) + 
  geom_ribbon(data = qx2.m.all, aes(x = age, ymin = -log(1- qx.lower), ymax = -log(1- qx.upper), fill= "steelblue"), alpha = 0.25) + 
  scale_fill_manual(values=c("steelblue", "steelblue"), guide="none")+
  theme_classic(base_size = 20) + 
  scale_y_continuous(expression(m[x]), limits = 10^-c(NA,NA), trans = 'log10', labels = scales::comma) +
  scale_x_continuous("Age", breaks = seq(0, 200, by = 20)) +
  scale_color_manual(id, values=c(rep("steelblue",5)), guide="none")+
  theme(legend.position = c(0.95, 0.20),
        strip.background=element_rect(colour="black",
                                      fill="gray87"), panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.title=element_blank(),
        legend.key.height = unit(.6, "cm") ,
        legend.text=element_text(color="black", size=16),
        axis.title = element_text(color = "black", size = 16),
        axis.text = element_text(color="black",size=14))+
  facet_wrap(~id, ncol=5, labeller = labeller(id = 
                                                c("1" = "0.80",
                                                  "2" = "0.85",
                                                  "3" = "0.90",
                                                  "4" = "0.95",
                                                  "5" = "different per age")))
graphics.off()



pdf("Fig3b.pdf", width=16, height=5)
ggplot() + 
#  geom_point(data = mx.f2011, aes(x = age, y = mx), col = "gray") +
  geom_point(data = filter(mx.f2011, age %in% 1:100), aes(x = age, y = mx), col = "tomato") +
  geom_point(data = filter(mx.f2011, age %in% 101:104), aes(x = age, y = mx), col = "tomato", shape=1) +
   geom_line(data=qx1.f.all, aes(x=age, y =  -log(1-qx.fitted), col = id, lty = "1-104"), col="tomato") + 
  geom_line(data=qx2.f.all, aes(x=age, y =  -log(1-qx.fitted), col = id, lty = "1-100"), col="tomato") + 
#  geom_hline(yintercept=1, linetype='dashed', col = 'black')+
  geom_ribbon(data = qx1.f.all, aes(x = age, ymin = -log(1-qx.lower), ymax = -log(1-qx.upper), fill="tomato"), alpha = 0.25) + 
  geom_ribbon(data = qx2.f.all, aes(x = age, ymin = -log(1- qx.lower), ymax = -log(1- qx.upper), fill="tomato"), alpha = 0.25) + 
  scale_fill_manual(values=c("tomato", "tomato"), guide="none")+
  theme_classic(base_size = 20) + 
  scale_y_continuous(expression(m[x]), limits = 10^-c(NA,NA), trans = 'log10', labels = scales::comma) +
  scale_x_continuous("Age", breaks = seq(0, 200, by = 20)) +
  theme(legend.position = c(0.95, 0.20),
        strip.background=element_rect(colour="black",
                                      fill="gray87"), panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.title=element_blank(),
        legend.key.height = unit(.6, "cm") ,
        legend.text=element_text(color="black", size=16),
        axis.title = element_text(color = "black", size = 16),
        axis.text = element_text(color="black",size=14))+
  facet_wrap(~id, ncol=5, labeller = labeller(id = 
                                                c("1" = "0.80",
                                                  "2" = "0.85",
                                                  "3" = "0.90",
                                                  "4" = "0.95",
                                                  "5" = "different per age")))
graphics.off()





#### zoomed
mx.m2011.<- mx %>%
  filter(year=='2011', sex=='1') %>%
  dplyr::  select(age,mx) %>%
  filter(age>=80)

mx.f2011.<- mx %>%
  
  filter(year=='2011', sex=='2') %>%
  dplyr::  select(age,mx) %>%
  filter(age>=80)
qx1.f.all.<- qx1.f.all %>%
  filter(age>=80)

qx2.f.all.<- qx2.f.all %>%
  filter(age>=80)

qx1.m.all.<- qx1.m.all %>%
  filter(age>=80)

qx2.m.all.<- qx2.m.all %>%
  filter(age>=80)

### 1-104 ano 2011
pdf("Fig5b.pdf", width=16, height=5)
ggplot() + 
  geom_point(data = mx.f2011., aes(x = age, y = mx), col = "tomato") +
  geom_point(data = mx.m2011., aes(x = age, y = mx), col = "steelblue", ) +
   geom_line(data=qx1.f.all., aes(x=age, y =  -log(1-qx.fitted), col = id), col="tomato") + 
  geom_line(data=qx1.m.all., aes(x=age, y =  -log(1-qx.fitted), col = id), col="steelblue") + 
  #geom_hline(yintercept=1, linetype='dashed', col = 'black')+
   geom_ribbon(data = qx1.f.all., aes(x = age, ymin = -log(1-qx.lower), ymax = -log(1-qx.upper), fill="tomato"), alpha = 0.25) + 
  geom_ribbon(data = qx1.m.all., aes(x = age, ymin = -log(1-qx.lower), ymax = -log(1-qx.upper), fill= "steelblue"), alpha = 0.25) + 
  #scale_fill_manual(values=c("steelblue", "tomato"), guide="none")+
  scale_color_manual(values=c("steelblue", "tomato"), labels=c("male", "female"))+
  scale_fill_manual(values=c("steelblue","tomato"), labels=c("male", "female")) + 
  theme_classic(base_size = 20) + 
  scale_y_continuous(expression(m[x]), limits = 10^-c(NA,NA), 
                                                     trans = 'log10', labels = scales::comma) + 
  scale_x_continuous("Age", breaks = seq(80, 200, by = 20)) + 
  theme(legend.position = c(0.95, 0.20),
        strip.background=element_rect(colour="black",
                                      fill="gray87"), panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.title=element_blank(),
        legend.key.height = unit(.6, "cm") ,
        legend.text=element_text(color="black", size=16),
        axis.title = element_text(color = "black", size = 14),
        axis.text = element_text(color="black",size=14))+
  facet_wrap(~id, ncol=5, labeller = labeller(id = 
                                                c("1" = "0.80",
                                                  "2" = "0.85",
                                                  "3" = "0.90",
                                                  "4" = "0.95",
                                                  "5" = "different per age")))


graphics.off()



#### zoomed
mx.m2011.<- mx %>%
  filter(year=='2011', sex=='1') %>%
  dplyr::select(age,mx) %>%
  filter(age<=35)

mx.f2011.<- mx %>%
  
  filter(year=='2011', sex=='2') %>%
  dplyr::select(age,mx) %>%
  filter(age<=35)
qx1.f.all.<- qx1.f.all %>%
  filter(age<=35)

qx1.m.all.<- qx1.m.all %>%
  filter(age<=35)

ggplot() +
  geom_point(data = mx.f2011., aes(x = age, y = mx), col = "tomato") +
  geom_point(data = mx.m2011., aes(x = age, y = mx), col = "steelblue", ) +
  geom_line(data=qx1.f.all., aes(x=age, y =  -log(1-qx.fitted), col = id), col="tomato") +
  geom_line(data=qx1.m.all., aes(x=age, y =  -log(1-qx.fitted), col = id), col="steelblue") +
  #geom_hline(yintercept=1, linetype='dashed', col = 'black')+
  geom_ribbon(data = qx1.f.all., aes(x = age, ymin = -log(1-qx.lower), ymax = -log(1-qx.upper), fill="tomato"), alpha = 0.25) +
  geom_ribbon(data = qx1.m.all., aes(x = age, ymin = -log(1-qx.lower), ymax = -log(1-qx.upper), fill= "steelblue"), alpha = 0.25) +
  #  scale_fill_manual(values=c("steelblue", "tomato"), guide="none")+
  scale_color_manual(values=c("steelblue", "tomato"), labels=c("male", "female"))+
  scale_fill_manual(values=c("steelblue","tomato"), labels=c("male", "female")) +
  theme_classic(base_size = 20) +
  scale_y_continuous(expression(m[x]), limits = 10^-c(NA,NA),
                     trans = 'log10', labels = scales::comma) +
  scale_x_continuous("Age", breaks = seq(0, 40, by = 20)) +
  theme(legend.position = c(0.95, 0.20),
        strip.background=element_rect(colour="black",
                                      fill="gray87"), panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.title=element_blank(),
        legend.key.height = unit(.6, "cm") ,
        legend.text=element_text(color="black", size=16),
        axis.title = element_text(color = "black", size = 14),
        axis.text = element_text(color="black",size=14))+
  facet_wrap(~id, ncol=5, labeller = labeller(id =
                                                c("1" = "0.80",
                                                  "2" = "0.85",
                                                  "3" = "0.90",
                                                  "4" = "0.95",
                                                  "5" = "different per age")))
# 



### 1-100 ano 2011
pdf("Fig6b.pdf", width=16, height=5)
ggplot() + 
  geom_point(data = filter(mx.f2011., age %in% 1:100), aes(x = age, y = mx), col = "tomato") +
  geom_point(data = filter(mx.m2011., age %in% 1:100), aes(x = age, y = mx), col = "steelblue") +
  geom_point(data = filter(mx.f2011., age %in% 101:104), aes(x = age, y = mx), col = "tomato", shape=1) +
  geom_point(data = filter(mx.m2011., age %in% 101:104), aes(x = age, y = mx), col = "steelblue", shape=1) +
  geom_line(data=qx2.f.all., aes(x=age, y =  -log(1-qx.fitted), col = id), col="tomato") + 
  geom_line(data=qx2.m.all., aes(x=age, y =  -log(1-qx.fitted), col = id), col="steelblue") + 
 # geom_hline(yintercept=1, linetype='dashed', col = 'black')+
   geom_ribbon(data = qx2.f.all., aes(x = age, ymin = -log(1-qx.lower), ymax = -log(1-qx.upper), fill="tomato"), alpha = 0.25) + 
  geom_ribbon(data = qx2.m.all., aes(x = age, ymin = -log(1-qx.lower), ymax = -log(1-qx.upper), fill= "steelblue"), alpha = 0.25) + 
#  scale_fill_manual(values=c("steelblue", "tomato"), guide="none")+
  scale_color_manual(values=c("steelblue", "tomato"), labels=c("male", "female"))+
  scale_fill_manual(values=c("steelblue","tomato"), labels=c("male", "female")) + 
  
    theme_classic(base_size = 20) + scale_y_continuous(expression(m[x]), limits = 10^-c(NA,NA), 
                                                     trans = 'log10', labels = scales::comma) + 
   scale_x_continuous("Age", breaks = seq(80, 200, by = 20)) + 
  theme(legend.position = c(0.95, 0.20),
        strip.background=element_rect(colour="black",
                                      fill="gray87"), panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.title=element_blank(),
        legend.key.height = unit(.6, "cm") ,
        legend.text=element_text(color="black", size=16),
        axis.title = element_text(color = "black", size = 14),
        axis.text = element_text(color="black",size=14))+
  facet_wrap(~id, ncol=5, labeller = labeller(id = 
                                                c("1" = "0.80",
                                                  "2" = "0.85",
                                                  "3" = "0.90",
                                                  "4" = "0.95",
                                                  "5" = "different per age")))


graphics.off()

#-------------

### data output cenario 1 - masculino 2012
d1<-qx.m2012[[1]]; d2<-qx.m2012[[2]]; d3<-qx.m2012[[3]]; d4<-qx.m2012[[4]];
d5<- qx.m2012.d12
qx1.m.all<- bind_rows(d1,d2,d3,d4,d5, .id="id")
mx.m2012<- mx %>%
  filter(year=='2012', sex=='1') %>%
  dplyr::  select(age,mx)
### data output cenario 2 - masculino 2012
d1<-qx2.m2012[[1]]; d2<-qx2.m2012[[2]]; d3<-qx2.m2012[[3]]; d4<-qx2.m2012[[4]];
d5<- qx2.m2012.d22
qx2.m.all<- bind_rows(d1,d2,d3,d4,d5, .id="id")


### data output cenario 1 - feminino 2012
d1<-qx.f2012[[1]]; d2<-qx.f2012[[2]]; d3<-qx.f2012[[3]]; d4<-qx.f2012[[4]];
d5<- qx.f2012.d12
qx1.f.all<- bind_rows(d1,d2,d3,d4,d5, .id="id")
mx.f2012<- mx %>%
  filter(year=='2012', sex=='2') %>%
  dplyr::  select(age,mx)
### data output cenario 2 - feminino 2012
d1<-qx2.f2012[[1]]; d2<-qx2.f2012[[2]]; d3<-qx2.f2012[[3]]; d4<-qx2.f2012[[4]];
d5<- qx2.f2012.d22
qx2.f.all<- bind_rows(d1,d2,d3,d4,d5, .id="id")


### cenario 1 versus cenario 2 - 2012 masculino 
pdf("Fig4a.pdf", width=16, height=5)
ggplot() + 
#  geom_point(data = mx.m2012, aes(x = age, y = mx), col = "gray") +
  geom_point(data = filter(mx.m2012, age %in% 1:100), aes(x = age, y = mx), col = "steelblue") +
  geom_point(data = filter(mx.m2012, age %in% 101:104), aes(x = age, y = mx), col = "steelblue", shape=1) +
  
    geom_line(data=qx1.m.all, aes(x=age, y =  -log(1-qx.fitted), color = id, lty = "1-104"),col = "steelblue") + 
  geom_line(data=qx2.m.all, aes(x=age, y =  -log(1-qx.fitted), color = id, lty = "1-100"),col = "steelblue") + 
  #geom_hline(yintercept=1, linetype='dashed', col = 'black')+
   geom_ribbon(data = qx1.m.all, aes(x = age, ymin = -log(1-qx.lower), ymax = -log(1-qx.upper),fill= "steelblue"), alpha = 0.25) + 
  geom_ribbon(data = qx2.m.all, aes(x = age, ymin = -log(1- qx.lower), ymax = -log(1- qx.upper),fill= "steelblue"), alpha = 0.25) + 
  scale_fill_manual(values=c("steelblue", "steelblue"), guide="none")+
  theme_classic(base_size = 20) + 
  scale_y_continuous(expression(m[x]), limits = 10^-c(NA,NA), trans = 'log10', labels = scales::comma) +
  scale_x_continuous("Age", breaks = seq(0, 200, by = 20)) +
  scale_color_manual(id, values=c(rep("steelblue",5)), guide="none")+
  theme(legend.position = c(0.95, 0.20),
        strip.background=element_rect(colour="black",
                                      fill="gray87"), panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.title=element_blank(),
        legend.key.height = unit(.6, "cm") ,
        legend.text=element_text(color="black", size=16),
        axis.title = element_text(color = "black", size = 16),
        axis.text = element_text(color="black",size=14))+
  facet_wrap(~id, ncol=5, labeller = labeller(id = 
                                                c("1" = "0.80",
                                                  "2" = "0.85",
                                                  "3" = "0.90",
                                                  "4" = "0.95",
                                                  "5" = "different per age")))
graphics.off()



pdf("Fig4b.pdf", width=16, height=5)
ggplot() + 
 # geom_point(data = mx.f2012, aes(x = age, y = mx), col = "gray") +
  geom_point(data = filter(mx.f2012, age %in% 1:100), aes(x = age, y = mx), col = "tomato") +
  geom_point(data = filter(mx.f2012, age %in% 101:104), aes(x = age, y = mx), col = "tomato", shape=1) +
  
  geom_line(data=qx1.f.all, aes(x=age, y =  -log(1-qx.fitted), col = id, lty = "1-104"), col="tomato") + 
  geom_line(data=qx2.f.all, aes(x=age, y =  -log(1-qx.fitted), col = id, lty = "1-100"), col="tomato") + 
 # geom_hline(yintercept=1, linetype='dashed', col = 'black')+
  geom_ribbon(data = qx1.f.all, aes(x = age, ymin = -log(1-qx.lower), ymax = -log(1-qx.upper), fill="tomato"), alpha = 0.25) + 
  geom_ribbon(data = qx2.f.all, aes(x = age, ymin = -log(1- qx.lower), ymax = -log(1- qx.upper), fill="tomato"), alpha = 0.25) + 
  scale_fill_manual(values=c("tomato", "tomato"), guide="none")+
   theme_classic(base_size = 20) + 
  scale_y_continuous(expression(m[x]), limits = 10^-c(NA,NA), trans = 'log10', labels = scales::comma) +
  scale_x_continuous("Age", breaks = seq(0, 200, by = 20)) +
  theme(legend.position = c(0.95, 0.20),
        strip.background=element_rect(colour="black",
                                      fill="gray87"), panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.title=element_blank(),
        legend.key.height = unit(.6, "cm") ,
        legend.text=element_text(color="black", size=16),
        axis.title = element_text(color = "black", size = 16),
        axis.text = element_text(color="black",size=14))+
  facet_wrap(~id, ncol=5, labeller = labeller(id = 
                                                c("1" = "0.80",
                                                  "2" = "0.85",
                                                  "3" = "0.90",
                                                  "4" = "0.95",
                                                  "5" = "different per age")))
graphics.off()




#### zoomed
mx.m2012.<- mx %>%
  filter(year=='2012', sex=='1') %>%
  dplyr::  select(age,mx) %>%
  filter(age>=80)

mx.f2012.<- mx %>%
  
  filter(year=='2012', sex=='2') %>%
  dplyr::  select(age,mx) %>%
  filter(age>=80)
qx1.f.all.<- qx1.f.all %>%
  filter(age>=80)

qx2.f.all.<- qx2.f.all %>%
  filter(age>=80)

qx1.m.all.<- qx1.m.all %>%
  filter(age>=80)

qx2.m.all.<- qx2.m.all %>%
  filter(age>=80)

### 1-104 ano 2012
pdf("Fig5c.pdf", width=16, height=5)
ggplot() + 
  geom_point(data = mx.f2012., aes(x = age, y = mx), col = "tomato") +
  geom_point(data = mx.m2012., aes(x = age, y = mx), col = "steelblue", ) +
  geom_line(data=qx1.f.all., aes(x=age, y =  -log(1-qx.fitted), col = id), col="tomato") + 
  geom_line(data=qx1.m.all., aes(x=age, y =  -log(1-qx.fitted), col = id), col="steelblue") + 
 # geom_hline(yintercept=1, linetype='dashed', col = 'black')+
  geom_ribbon(data = qx1.f.all., aes(x = age, ymin = -log(1-qx.lower), ymax = -log(1-qx.upper), fill="tomato"), alpha = 0.25) + 
  geom_ribbon(data = qx1.m.all., aes(x = age, ymin = -log(1-qx.lower), ymax = -log(1-qx.upper), fill= "steelblue"), alpha = 0.25) + 
 # scale_fill_manual(values=c("steelblue", "tomato"), guide="none")+
  scale_color_manual(values=c("steelblue", "tomato"), labels=c("male", "female"))+
  scale_fill_manual(values=c("steelblue","tomato"), labels=c("male", "female")) + 
  
    theme_classic(base_size = 20) + scale_y_continuous(expression(m[x]), limits = 10^-c(NA,NA), 
                                                    trans = 'log10', labels = scales::comma) + 
  scale_x_continuous("Age", breaks = seq(80, 200, by = 20)) + 
  theme(legend.position = c(0.95, 0.20),
        strip.background=element_rect(colour="black",
                                      fill="gray87"), panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.title=element_blank(),
        legend.key.height = unit(.6, "cm") ,
        legend.text=element_text(color="black", size=16),
        axis.title = element_text(color = "black", size = 14),
        axis.text = element_text(color="black",size=14))+
  facet_wrap(~id, ncol=5, labeller = labeller(id = 
                                                c("1" = "0.80",
                                                  "2" = "0.85",
                                                  "3" = "0.90",
                                                  "4" = "0.95",
                                                  "5" = "different per age")))


graphics.off()



#### zoomed
mx.m2012.<- mx %>%
  filter(year=='2012', sex=='1') %>%
  dplyr::select(age,mx) %>%
  filter(age<=35)

mx.f2012.<- mx %>%
  
  filter(year=='2012', sex=='2') %>%
  dplyr::select(age,mx) %>%
  filter(age<=35)
qx1.f.all.<- qx1.f.all %>%
  filter(age<=35)

qx1.m.all.<- qx1.m.all %>%
  filter(age<=35)

ggplot() +
  geom_point(data = mx.f2012., aes(x = age, y = mx), col = "tomato") +
  geom_point(data = mx.m2012., aes(x = age, y = mx), col = "steelblue", ) +
  geom_line(data=qx1.f.all., aes(x=age, y =  -log(1-qx.fitted), col = id), col="tomato") +
  geom_line(data=qx1.m.all., aes(x=age, y =  -log(1-qx.fitted), col = id), col="steelblue") +
  #geom_hline(yintercept=1, linetype='dashed', col = 'black')+
  geom_ribbon(data = qx1.f.all., aes(x = age, ymin = -log(1-qx.lower), ymax = -log(1-qx.upper), fill="tomato"), alpha = 0.25) +
  geom_ribbon(data = qx1.m.all., aes(x = age, ymin = -log(1-qx.lower), ymax = -log(1-qx.upper), fill= "steelblue"), alpha = 0.25) +
  #  scale_fill_manual(values=c("steelblue", "tomato"), guide="none")+
  scale_color_manual(values=c("steelblue", "tomato"), labels=c("male", "female"))+
  scale_fill_manual(values=c("steelblue","tomato"), labels=c("male", "female")) +
  theme_classic(base_size = 20) +
  scale_y_continuous(expression(m[x]), limits = 10^-c(NA,NA),
                     trans = 'log10', labels = scales::comma) +
  scale_x_continuous("Age", breaks = seq(0, 40, by = 20)) +
  theme(legend.position = c(0.95, 0.20),
        strip.background=element_rect(colour="black",
                                      fill="gray87"), panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.title=element_blank(),
        legend.key.height = unit(.6, "cm") ,
        legend.text=element_text(color="black", size=16),
        axis.title = element_text(color = "black", size = 14),
        axis.text = element_text(color="black",size=14))+
  facet_wrap(~id, ncol=5, labeller = labeller(id =
                                                c("1" = "0.80",
                                                  "2" = "0.85",
                                                  "3" = "0.90",
                                                  "4" = "0.95",
                                                  "5" = "different per age")))
# 

### 1-100 ano 2012
pdf("Fig6c.pdf", width=16, height=5)
ggplot() + 
  geom_point(data = filter(mx.f2012., age %in% 1:100), aes(x = age, y = mx), col = "tomato") +
  geom_point(data = filter(mx.m2012., age %in% 1:100), aes(x = age, y = mx), col = "steelblue") +
  geom_point(data = filter(mx.f2012., age %in% 101:104), aes(x = age, y = mx), col = "tomato", shape=1) +
  geom_point(data = filter(mx.m2012., age %in% 101:104), aes(x = age, y = mx), col = "steelblue", shape=1) +
  geom_line(data=qx2.f.all., aes(x=age, y =  -log(1-qx.fitted), col = id), col="tomato") + 
  geom_line(data=qx2.m.all., aes(x=age, y =  -log(1-qx.fitted), col = id), col="steelblue") + 
 # geom_hline(yintercept=1, linetype='dashed', col = 'black')+
   geom_ribbon(data = qx2.f.all., aes(x = age, ymin = -log(1-qx.lower), ymax = -log(1-qx.upper), fill="tomato"), alpha = 0.25) + 
  geom_ribbon(data = qx2.m.all., aes(x = age, ymin = -log(1-qx.lower), ymax = -log(1-qx.upper), fill= "steelblue"), alpha = 0.25) + 
 # scale_fill_manual(values=c("steelblue", "tomato"), guide="none")+
  scale_color_manual(values=c("steelblue", "tomato"), labels=c("male", "female"))+
  scale_fill_manual(values=c("steelblue","tomato"), labels=c("male", "female")) + 
  
    theme_classic(base_size = 20) + scale_y_continuous(expression(m[x]), limits = 10^-c(NA,NA), 
                                                     trans = 'log10', labels = scales::comma) + 
  scale_x_continuous("Age", breaks = seq(80, 200, by = 20)) + 
  theme(legend.position = c(0.95, 0.20),
        strip.background=element_rect(colour="black",
                                      fill="gray87"), panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.title=element_blank(),
        legend.key.height = unit(.6, "cm") ,
        legend.text=element_text(color="black", size=16),
        axis.title = element_text(color = "black", size = 14),
        axis.text = element_text(color="black",size=14))+
  facet_wrap(~id, ncol=5, labeller = labeller(id = 
                                                c("1" = "0.80",
                                                  "2" = "0.85",
                                                  "3" = "0.90",
                                                  "4" = "0.95",
                                                  "5" = "different per age")))


graphics.off()
