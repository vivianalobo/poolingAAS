######################################################
#### authors: figueiredo, lobo, alves and fonseca
#### output: modelling amazonas data
### data: amazonas, 2022-2024
######################################################

#### section4: Amazonas mortality data

### run amazon dataset and auxiliary functions 
source("amazoniadata.R")
source("multivariatefunc.R")
source("predict_dlm_at.R")
source("predict_bivar2_att_chain.R")
source("funcmetrics.R")
source("auxfunc.R")


#--------------------------------------------------------------------------------------------------------------------------------------
###### Section 4: Mortality data analysis

### 4.1 Amazonas mortality data
# We next focus on an empirical application based on mortality data from micro-regions in the Brazilian Amazon.
#The ages range from 0 to 80+ years, and are available in five-year intervals. 
#We use this setting to examine three complementary aspects of the proposed joint modelling framework.
#Together, these analyses illustrate how the multivariate state-space structure can accommodate extrapolation and mortality closure at advanced ages,
#improve estimation in data-sparse settings, and exploit dependence across related populations while preserving local demographic features.  

########==================================================================================
###4.1.1 Joint modelling of male and female mortality in Manaus
########==================================================================================
### dataset male and female - mortality rates
df.m<- df.agregado %>%
  filter(sex == "Male")
df.f<- df.agregado %>%
  filter(sex == "Female")

df.m<- df.m %>%
  dplyr::select(Microrregião.IBGE, age_start, log_mx) %>%  
  pivot_wider(
    names_from = Microrregião.IBGE,  
    values_from = log_mx          
  ) %>%
  arrange(age_start)  

df.f<- df.f %>%
  dplyr::select(Microrregião.IBGE, age_start, log_mx) %>% 
  pivot_wider(
    names_from = Microrregião.IBGE,  
    values_from = log_mx            
  ) %>%
  arrange(age_start)  

age<- seq(0,80, by=5)


yM.manaus<- as.vector(df.m[,8]) ; yF.manaus<- as.vector(df.f[,8]) ;
y.manaus<-  data.frame(yF.manaus, yM.manaus)
y.manaus<- as.matrix(y.manaus)

V4= cov(y.manaus)
V4
J=2

mx.obs <- df.agregado %>%
  filter( Microrregião.IBGE == "MANAUS") 
mx.obs <-  mx.obs[, c("sex","age_start", "mx")]
colnames(mx.obs)<- c("regiao", "age", "mx")


age<- seq(0,80, by=5)


### priors, Ft and Gt
res <- buildFtGt(J, alpha= TRUE)
Ft <- res$Ft
Gt <- res$Gt
m0 = rep(0, nrow(Gt))
C0 = diag(100, nrow(Gt))
s0 <- diag(J)*0.01  ###priori do petris vaga
v0 <- 5


#--------------------------------------------------------------------
####  modelling with common term and different discount factor matrix 

## discount configuration 1
d_male   <- c(rep(0.99, 4), rep(0.95,7), rep(0.999, 6))
d_female <- c(rep(0.99, 3), rep(0.95,8), rep(0.999, 6))
d_alpha  <- c(rep(0.95, 5), rep(0.999, 12))
deltaJ1<- cbind(d_female, d_male, d_alpha)
fit.t1<- dlm.multivariate(y = y.manaus, Ft = Ft, Gt = Gt, nit = 30000, bn = 10000, thin = 1,
                               v0 = v0, s0 = s0, m0 = m0, C0 = C0, delta = deltaJ1, V = V4)


## discount configuration 2
d_male   <- c(rep(0.99, 4), rep(0.95,7), rep(0.999, 6))
d_female <- c(rep(0.99, 3), rep(0.95,8), rep(0.999, 6))
d_alpha  <- c(rep(0.999, 5), rep(0.999, 12))
deltaJ2<- cbind(d_female, d_male, d_alpha)

fit.t2<- dlm.multivariate(y = y.manaus, Ft = Ft, Gt = Gt, nit = 30000, bn = 10000, thin = 1,
                               v0 = v0, s0 = s0, m0 = m0, C0 = C0, delta = deltaJ2, V = V4)


## discount configuration 3
### best configuration following log predictive posterior density
d_male   <- c(rep(0.90, 4), rep(0.95,7), rep(0.999, 6))
d_female <- c(rep(0.90, 3), rep(0.95,8), rep(0.999, 6))
d_alpha  <- c(rep(0.95, 5), rep(0.999, 12))
deltaJ3<- cbind(d_female, d_male, d_alpha)
fit.t3<- dlm.multivariate(y = y.manaus, Ft = Ft, Gt = Gt, nit = 30000, bn = 10000, thin = 1,
                               v0 = v0, s0 = s0, m0 = m0, C0 = C0, delta = deltaJ3, V = V4)


## discount configuration 4
d_male   <- c(rep(0.95, 4), rep(0.95,7), rep(0.999, 6))
d_female <- c(rep(0.95, 3), rep(0.95,8), rep(0.999, 6))
d_alpha  <- c(rep(0.95, 5), rep(0.999, 12))
deltaJ4<- cbind(d_female, d_male, d_alpha)
fit.t4<- dlm.multivariate(y = y.manaus, Ft = Ft, Gt = Gt, nit = 30000, bn = 10000, thin = 1,
                          v0 = v0, s0 = s0, m0 = m0, C0 = C0, delta = deltaJ4, V = V4)


# Table 1: Comparison of candidate discount matrix specifications based on the log posterior predictive density (LPPD). 
#The best-performing model is indicated in boldface.
ls1 <- log_pred_density(y.manaus, fit.t1$mu, fit.t1$V)
ls2 <- log_pred_density(y.manaus, fit.t2$mu, fit.t2$V)
ls3 <- log_pred_density(y.manaus, fit.t3$mu, fit.t3$V)
ls4 <- log_pred_density(y.manaus, fit.t4$mu, fit.t4$V)

 data.frame(
  model = c("proposta 1",
            "proposta 2",
            "proposta 3",
            "proposta 4"),
  logscore = c(ls1, ls2, ls3,ls4)
)

 #---------------------------------------------
#### modelling without common term

 ### priors, Ft and Gt
 res <- buildFtGt(J=2, alpha= FALSE)
 Ft2 <- res$Ft
 Gt2 <- res$Gt
 
 m02 = rep(0, nrow(Gt2))
 C02 = diag(100, nrow(Gt2))
 s02 <- diag(J)*0.01  ###priori do petris vaga
 v02 <- 5

####  modelling with common term and different discount factor matrix 
## discount configuration 1
 
 d_male   <- c(rep(0.99, 4), rep(0.95,7), rep(0.999, 6))
 d_female <- c(rep(0.99, 3), rep(0.95,8), rep(0.999, 6))
 deltaJ1wtct<- cbind(d_female, d_male)
 
 fit.t1wtct<- dlm.multivariate(y = y.manaus, Ft = Ft2, Gt = Gt2, nit = 30000, bn = 10000, thin = 1,
                           v0 = v02, s0 = s02, m0 = m02, C0 = C02, delta = deltaJ1wtct, V = V4)

 ## discount configuration 2
 d_male   <- c(rep(0.99, 4), rep(0.95,7), rep(0.999, 6))
 d_female <- c(rep(0.99, 3), rep(0.95,8), rep(0.999, 6))
 deltaJ2wtct<- cbind(d_female, d_male)
 
 fit.t2wtct<- dlm.multivariate(y = y.manaus, Ft = Ft2, Gt = Gt2, nit = 30000, bn = 10000, thin = 1,
                           v0 = v02, s0 = s02, m0 = m02, C0 = C02, delta = deltaJ2wtct, V = V4)
 
 
 ## discount configuration 3
 ### best configuration following log predictive posterior density
 d_male   <- c(rep(0.90, 4), rep(0.95,7), rep(0.999, 6))
 d_female <- c(rep(0.90, 3), rep(0.95,8), rep(0.999, 6))
 deltaJ3wtct<- cbind(d_female, d_male)
 fit.t3wtct<- dlm.multivariate(y = y.manaus, Ft = Ft2, Gt = Gt2, nit = 30000, bn = 10000, thin = 1,
                           v0 = v02, s0 = s02, m0 = m02, C0 = C02, delta = deltaJ3wtct, V = V4)
 ## discount configuration 4
 d_male   <- c(rep(0.95, 4), rep(0.95,7), rep(0.999, 6))
 d_female <- c(rep(0.95, 3), rep(0.95,8), rep(0.999, 6))
 deltaJ4wtct<- cbind(d_female, d_male)
 fit.t4wtct<- dlm.multivariate(y = y.manaus, Ft = Ft2, Gt = Gt2, nit = 30000, bn = 10000, thin = 1,
                           v0 = v02, s0 = s02, m0 = m02, C0 = C02, delta = deltaJ4wtct, V = V4)
 
 #The best-performing model is indicated in boldface. We have omitted in the paper
 ls1wtct <- log_pred_density(y.manaus, fit.t1wtct$mu, fit.t1wtct$V)
 ls2wtct <- log_pred_density(y.manaus, fit.t2wtct$mu, fit.t2wtct$V)
 ls3wtct <- log_pred_density(y.manaus, fit.t3wtct$mu, fit.t3wtct$V)
 ls4wtct <- log_pred_density(y.manaus, fit.t4wtct$mu, fit.t4wtct$V)
 
data.frame(
   model = c("proposta 1",
             "proposta 2",
             "proposta 3",
             "proposta 4"),
   logscore = c(ls1wtct, ls2wtct, ls3wtct,ls4wtct)
 )

##### Extrapolation and convergence between male and female population
 w= 120 ## age extrapolation
 h<- (w-last(age))/5 ###  de 5 em 5 anos 
 
 
### extrapolation for model with discount factor configuration 3 (with common term)
fit.pred.manaus3 <- predict.multivariate(y.manaus, m0, C0, h, V = fit.manaus3$V,
                                         Ft, Gt, delta=deltaJ3, ages = age)

mx.aux <- df.agregado %>%
  filter( Microrregião.IBGE == "MANAUS") 
mx.mar <-  mx.aux[, c("sex","age_start", "mx", "Dx", "Ex")]

mx.mar <- mx.aux %>%
  dplyr::select(age_start, Dx, Ex) %>%
  group_by(age_start) %>%
  summarise(
    Dx = sum(Dx),
    Ex = sum(Ex)
  ) %>%
  mutate(logmx = log(Dx/Ex))


pred.manaus3 <- predict_conv_chain(obj_pred = fit.pred.manaus3,
                              age = age,
                              h = h,
                              y_base = mx.mar %>% pull(logmx),
                              delta.aux =  0.8,
                              V = 0.1,
                              final_age = 120)

aux.names <- c("qx.fitted", "qx.lower", "qx.upper")
qxF.manaus3 = data.frame(age=seq(0,w, by=5) , rbind(qx_fitted(fit.manaus3)[[1]],
                                       setNames(pred.manaus3[[1]], aux.names)) )
qxM.manaus3 = data.frame(age=seq(0,w, by=5) , rbind(qx_fitted(fit.manaus3)[[2]],
                                                     setNames(pred.manaus3[[2]], aux.names)) )

qxF.manaus3$regiao <- "Female"
qxM.manaus3$regiao <- "Male"

out <- bind_rows(qxF.manaus3, qxM.manaus3)
out$local <- "Manaus"


ggplot() +
  geom_point(data = mx.obs, aes(x = age, y = log(mx), color = regiao), size = 1.7) +
  geom_line(data = out,aes(x = age,y = log(qx.fitted),color = regiao),linewidth = 1.2) +
  geom_ribbon(data = out, aes(x = age,ymin = log(qx.lower),ymax = log(qx.upper),fill = regiao),alpha = 0.25) +
  scale_y_continuous(expression("log rate"), limits = c(-10, 3), breaks=seq(-10,3, by=2)) +
  scale_x_continuous("Age", breaks = seq(0, 120, by = 20)) +
  scale_color_manual(values = c("Female" = "tomato","Male" = "steelblue")) +
  scale_fill_manual(values = c("Female" = "tomato","Male" = "steelblue")) +
  theme_classic(base_size = 20) +
  facet_wrap(~local)+
  theme(legend.position = c(0.85, 0.15),legend.title = element_blank(), 
        strip.background = element_rect(colour = "black", fill = "gray87"),
        panel.border = element_rect(color = "black", fill = NA)) -> pt


###### fit multivariate model without common term
res1 <- buildFtGt(J, alpha= FALSE)

Ft2 <- res1$Ft
Gt2 <- res1$Gt

m02 = rep(0, nrow(Gt2))
C02 = diag(100, nrow(Gt2))
s02 <- diag(J)*0.01  ###priori do petris vaga
v02 <- 5


d_male   <- c(rep(0.90, 4), rep(0.95,7), rep(0.999, 6))
d_female <- c(rep(0.90, 3), rep(0.95,8), rep(0.999, 6))
deltaJ3usual<- cbind(d_female, d_male)
fit.t3usual<- dlm.multivariate(y = y.manaus, Ft = Ft2, Gt = Gt2, nit = 30000, bn = 10000, thin = 1,
                          v0 = v02, s0 = s02, m0 = m02, C0 = C02, delta = deltaJ3usual, V = V4)

w= 120 ## idade da extrapolacao 
h<- (w-last(age))/5 ###  de 5 em 5 anos 

### predict
fit.pred.manaus3usual <- predict.multivariate(y.manaus, m02, C02, h, V = fit.t3usual$V,
                                         Ft2, Gt2, delta=deltaJ3usual, ages = age)

### convergence between populations
mx.aux <- df.agregado %>%
  filter( Microrregião.IBGE == "MANAUS") 
mx.mar <-  mx.aux[, c("sex","age_start", "mx", "Dx", "Ex")]

mx.mar <- mx.aux %>%
  dplyr::select(age_start, Dx, Ex) %>%
  group_by(age_start) %>%
  summarise(
    Dx = sum(Dx),
    Ex = sum(Ex)
  ) %>%
  mutate(logmx = log(Dx/Ex))


pred.manaus3usual <- predict_conv_chain(obj_pred = fit.pred.manaus3usual,
                                   age = age,
                                   h = h,
                                   y_base = mx.mar %>% pull(logmx),
                                   delta.aux =  0.9,
                                   V = 0.1,
                                   final_age = 120)

aux.names <- c("qx.fitted", "qx.lower", "qx.upper")
qxF.manaus3usual = data.frame(age=seq(0,w, by=5) , rbind(qx_fitted(fit.t3usual)[[1]],
                                                    setNames(pred.manaus3usual[[1]], aux.names)) )
qxM.manaus3usual = data.frame(age=seq(0,w, by=5) , rbind(qx_fitted(fit.t3usual)[[2]],
                                                    setNames(pred.manaus3usual[[2]], aux.names)) )

qxF.manaus3usual$regiao <- "Female"
qxM.manaus3usual$regiao <- "Male"

outusual <- bind_rows(qxF.manaus3usual, qxM.manaus3usual)
outusual$local <- "Manaus"



ggplot() +
  geom_point(data = mx.obs, aes(x = age, y = log(mx), color = regiao), size = 1.7) +
  geom_line(data = outusual,aes(x = age,y = log(qx.fitted),color = regiao),linewidth = 1.2) +
  geom_ribbon(data = outusual, aes(x = age,ymin = log(qx.lower),ymax = log(qx.upper),fill = regiao),alpha = 0.25) +
  #  geom_hline(yintercept = 1,linetype = "dashed",color = "black") +
  scale_y_continuous(expression("log rate"), limits = c(-10, 3), breaks=seq(-10,3, by=2)) +
  scale_x_continuous("Age", breaks = seq(0, 120, by = 20)) +
  scale_color_manual(values = c("Female" = "tomato","Male" = "steelblue")) +
  scale_fill_manual(values = c("Female" = "tomato","Male" = "steelblue")) +
  theme_classic(base_size = 20) +
  facet_wrap(~local)+
  theme(legend.position = c(0.85, 0.15),legend.title = element_blank(), 
        strip.background = element_rect(colour = "black", fill = "gray87"),
        panel.border = element_rect(color = "black", fill = NA)) -> pt


########================================================================================== 
###4.1.2 Borrowing strength across sexes in sparse populations
########================================================================================== 

#We next analyse the Juruá and Japurá micro-regions, which have substantially smaller
# populations and sparser mortality data than the Manaus region.


#### A. Modelling JURUA region
## population 1 - female jurua
## population 2 - male jurua
## population 3 - male manaus

## setting log mortality rate and priors
 yM.manaus<- as.vector(df.m[,8]) ;  yF.manaus<- as.vector(df.f[,8]) ; 
 yM.jurua<- as.vector(df.m[,5]) ; yF.jurua<- as.vector(df.f[,5]) ;
 
 y.jurua<-  data.frame(yF.jurua, yM.jurua, yM.manaus) #, yF.manaus)
 y.jurua<- as.matrix(y.jurua)
 
 V5= cov(y.jurua)
 V5
 J=3# 4
 age<- seq(0,80, by=5)
 
 res <- buildFtGt(J, alpha= TRUE)
 
 Ft <- res$Ft
 Gt <- res$Gt
 
 m0 = rep(0, nrow(Gt))
 C0 = diag(100, nrow(Gt))
 s0 <- diag(J)*0.01  ###priori do petris vaga
 v0 <- 5
 
#### discount matrix
 d_maleM <- c(rep(0.96, 4), rep(0.96,4), rep(0.999, 9))
 d_maleJ <- c(rep(0.96, 8), rep(0.999, 9))
 d_femaleJ <- c(rep(0.96, 5),rep(0.98, 5),rep(0.995, 7))
 d_alpha <- c(rep(0.98, 10), rep(0.999, 7))
 deltaJ5<- cbind(d_femaleJ, d_maleJ, d_maleM,  d_alpha)

 ### fit multivariate model with common term
 fit.joint.jurua<- dlm.multivariate(y = y.jurua, Ft = Ft, Gt = Gt, nit = 30000, bn = 10000, thin = 1,
                           v0 = v0, s0 = s0, m0 = m0, C0 = C0, delta = deltaJ5, V = V5)
 
### fitted probability for each population
 qxF.jurua <- data.frame(age=seq(0,80, by=5) ,qx_fitted(fit.joint.jurua)[[1]])
 qxM.jurua <- data.frame(age=seq(0,80, by=5) , qx_fitted(fit.joint.jurua)[[2]])
 qxM.manaus <- data.frame(age=seq(0,80, by=5) , qx_fitted(fit.joint.jurua)[[3]])

#### building the comparison 
 qxF.jurua$regiao <- "Female"
 qxM.jurua$regiao <- "Male"
 out.jurua <- bind_rows(qxF.jurua, qxM.jurua)
 out.jurua$local <- "Juruá (pooling with 3 populations)"
 mx.obs.jurua <- df.agregado %>%
   filter( Microrregião.IBGE == "JURUA") 
 mx.obs.jurua <-  mx.obs.jurua[, c("sex","age_start", "mx")]
 colnames(mx.obs.jurua)<- c("regiao", "age", "mx")
 mx.manaus.male <- df.agregado %>%
   filter( Microrregião.IBGE == "MANAUS", sex=="Male") 
 mx.manaus.male <-  mx.manaus.male[, c("sex","age_start", "mx")]
 colnames(mx.manaus.male)<- c("regiao", "age", "mx")
 
 
 mx.fitWithAlpha<- data.frame(mxFhat= qxF.japura[,2], mxMhat= qxM.japura[,2] ,mxM2hat=qxM.manaus[,2])
 mx.fit.multivariateWithAlpha<- as.data.frame(mx.fitWithAlpha)
 
g.juruaJoint<- ggplot() +
  geom_point(data = mx.obs.jurua, aes(x = age,y = log(mx),color = regiao),size = 1.7) +
  geom_line(data = out.jurua,aes(x = age,y = log(qx.fitted), color = regiao),linewidth = 1.2) +
  geom_ribbon(data = out.jurua,aes(x = age,ymin = log(qx.lower),ymax =log(qx.upper),fill = regiao),
              alpha = 0.25,show.legend = FALSE) +
  geom_line(data = qxM.manaus,aes(x = age, y = log(qx.fitted),color = "Manaus (male)"),linetype = "dashed",linewidth = 1) +
  scale_color_manual( values = c("Female" = "tomato","Male" = "steelblue","Manaus (male)" = "black")) +
  scale_fill_manual( values = c("Female" = "tomato","Male" = "steelblue") ) +
  scale_y_continuous( expression("log rate"), limits=c(-10,0)) +
  scale_x_continuous("Age",breaks = seq(0, 110, by = 10)) +
  labs(color = NULL) +
  theme_classic(base_size = 20) +
  facet_wrap(~local) +
  theme( legend.position = c(0.70, 0.18),  legend.title = element_blank(),strip.background = element_rect(
    colour = "black",fill = "gray87" ),panel.border = element_rect(color = "black",fill = NA))
g.juruaJoint


library(gridExtra)
#-----
### fit multivariate model without common term
## multivariado sem termo comum 
V5= cov(y.jurua)
V5
deltaJ5noalpha<- cbind(d_femaleJ, d_maleJ, d_maleM)
J=3#
res4 <- buildFtGt(J, alpha= FALSE)
Ft4 <- res4$Ft
Gt4 <- res4$Gt
m04 = rep(0, nrow(Gt4))
C04 = diag(100, nrow(Gt4))
s04 <- diag(J)*0.01  ###priori do petris vaga
v04 <- 5


fit.joint.jurua.without.alpha<- dlm.multivariate(y = y.jurua, Ft = Ft4, Gt = Gt4, nit = 30000, bn = 10000, thin = 1,
                                                 v0 = v04, s0 = s04, m0 = m04, C0 = C04, delta = deltaJ5noalpha, V = V5)


qxF.jurua <- data.frame(age=seq(0,80, by=5) ,qx_fitted(fit.joint.jurua.without.alpha)[[1]])
qxM.jurua <- data.frame(age=seq(0,80, by=5) , qx_fitted(fit.joint.jurua.without.alpha)[[2]])
qxM.manaus <- data.frame(age=seq(0,80, by=5) , qx_fitted(fit.joint.jurua.without.alpha)[[3]])

mx.fit.without.alpha<- data.frame(mxFhat= qxF.japura[,2], mxMhat= qxM.japura[,2] ,mxM2hat=qxM.manaus[,2])
mx.fit.without.alpha<- as.data.frame(mx.fit.without.alpha)

#-----
#### fit univariate model
y.jurua1<-  data.frame(yF.jurua, yM.jurua) #, yF.manaus)
y.jurua1<- as.matrix(y.jurua1)


Vunivariate<- diag(2) 
deltaJunivariate<- cbind(d_femaleJ, d_maleJ)
J=2# 4
res3 <- buildFtGt(J, alpha= FALSE)
Ft3 <- res3$Ft
Gt3 <- res3$Gt
m03 = rep(0, nrow(Gt3))
C03 = diag(100, nrow(Gt3))
s03 <- diag(J)*0.01  ###priori do petris vaga
v03 <- 5

fit.joint.jurua.univariate<- dlm.multivariate(y = y.jurua1, Ft = Ft3, Gt = Gt3, nit = 30000, bn = 10000, thin = 1,
                                   v0 = v03, s0 = s03, m0 = m03, C0 = C03, delta = deltaJunivariate, V = Vunivariate)
### para o calculo das metricas
qxF.jurua <- data.frame(age=seq(0,80, by=5) ,qx_fitted(fit.joint.jurua.univariate)[[1]])
qxM.jurua <- data.frame(age=seq(0,80, by=5) , qx_fitted(fit.joint.jurua.univariate)[[2]])

mx.fit.univariate<- data.frame(mxFhat= qxF.jurua[,2], mxMhat= qxM.jurua[,2])
mx.fit.univariate<- as.data.frame(mx.fit.univariate)

####### model comparison
#LPPD
ls.univariate <- log_pred_density(y.jurua1, fit.joint.jurua.univariate$mu,fit.joint.jurua.univariate$V)
ls.multivariate <- log_pred_density(y.jurua, fit.joint.jurua.without.alpha$mu,fit.joint.jurua.without.alpha$V)
ls.multivariateWithAlpha <- log_pred_density(y.jurua, fit.joint.jurua $mu, fit.joint.jurua $V)


 data.frame(
  model = c("model 1",
            "model 2",
            "model 3"),
  logscore = round(c(ls.univariate, ls.multivariate, ls.multivariateWithAlpha),3)
)

 ### weigthed metrics
 ExM.jurua <- extract_var(data = df.agregado,var = "Ex",regions = "JURUA",sex = "Male")
 ExF.jurua <- extract_var(data = df.agregado,var = "Ex",regions = "JURUA",sex = "Female")
 ExM.manaus <- extract_var(data = df.agregado,var = "Ex",regions = "MANAUS",sex = "Male")
 
 mxM.jurua <- extract_var(data = df.agregado,var = "mx",regions = "JURUA",sex = "Male")
 mxF.jurua <- extract_var(data = df.agregado,var = "mx",regions = "JURUA",sex = "Female")
 mxM.manaus <- extract_var(data = df.agregado,var = "mx",regions = "MANAUS",sex = "Male")
 
 Ex.all<- data.frame(ExF= ExF.jurua, ExM= ExM.jurua, ExM2= ExM.manaus)
 mx.all<- data.frame(mxF= mxF.jurua, mxM= mxM.jurua,mxM2= mxM.manaus)
 
 
 multivariate.with.alpha<-weighted_metrics(
   obs  = log(mx.all),
   fit  = log(mx.fit.multivariateWithAlpha),
   expo = Ex.all
 )
 multivariate.with.alpha

 multivariate.without.alpha<-
   weighted_metrics(
     obs  = log(mx.all),
     fit  = log(mx.fit.without.alpha),
     expo = Ex.all
   )
 
 multivariate.without.alpha

 univariate<-weighted_metrics(
   obs  = log(mx.all)[,-3],
   fit  = log(mx.fit.univariate),
   expo = Ex.all[,-3]
 )
 
 univariate

##--------------------------------------------------------------- 
#### B. Modelling JAPURA region

## population 1 - female japura
## population 2 - male japura
## population 3 - male manaus

 mx.obs.japura <- df.agregado %>%
  filter( Microrregião.IBGE == "JAPURA") 
mx.obs.japura <-  mx.obs.japura[, c("sex","age_start", "mx")]
colnames(mx.obs.japura)<- c("regiao", "age", "mx")

mx.manaus.male <- df.agregado %>%
  filter( Microrregião.IBGE == "MANAUS", sex=="Male") 
mx.manaus.male <-  mx.manaus.male[, c("sex","age_start", "mx")]
colnames(mx.manaus.male)<- c("regiao", "age", "mx")

yM.manaus<- as.vector(df.m[,8]) ;  yF.manaus<- as.vector(df.f[,8]) ; 
yM.japura<- as.vector(df.m[,3]) ; yF.japura<- as.vector(df.f[,3]) ;

y.japura<-  data.frame(yF.japura, yM.japura, yM.manaus) #, yF.manaus)
y.japura<- as.matrix(y.japura)

V6= cov(y.japura)
V6
J= 3# 4
age<- seq(0,80, by=5)

res <- buildFtGt(J, alpha= TRUE)

Ft <- res$Ft
Gt <- res$Gt

m0 = rep(0, nrow(Gt))
C0 = diag(100, nrow(Gt))
s0 <- diag(J)*0.01  ###priori do petris vaga
v0 <- 15


### discount factor matrix
d_maleM <- c(rep(0.96, 4), rep(0.96,4), rep(0.999, 9))
d_maleJ <- c(rep(0.96, 8), rep(0.996,3),rep(0.998,2), rep(0.999, 4))#c(rep(0.95, 6),rep(0.999, 11))
d_femaleJ <- c(rep(0.98, 5),rep(0.995, 5), rep(0.997, 5),rep(0.999, 2))# c(rep(0.95, 6), rep(0.999,11))
d_alpha <- c(rep(0.985, 10), rep(0.999, 7))
deltaJ6<- cbind(d_femaleJ, d_maleJ, d_maleM, d_alpha)

## fit multivariate model with common term
fit.joint.japura<- dlm.multivariate(y = y.japura, Ft = Ft, Gt = Gt, nit = 30000, bn = 10000, thin = 1,
                                   v0 = v0, s0 = s0, m0 = m0, C0 = C0, delta = deltaJ6, V = V6)


qxF.japura <- data.frame(age=seq(0,80, by=5) ,qx_fitted(fit.joint.japura)[[1]])
qxM.japura <- data.frame(age=seq(0,80, by=5) , qx_fitted(fit.joint.japura)[[2]])
qxM.manaus <- data.frame(age=seq(0,80, by=5) , qx_fitted(fit.joint.japura)[[3]])


qxF.japura$regiao <- "Female"
qxM.japura$regiao <- "Male"

out.japura <- bind_rows(qxF.japura, qxM.japura)
out.japura$local <- "Japurá (pooling with 3 populations)"

mx.fitWithAlpha<- data.frame(mxFhat= qxF.japura[,2], mxMhat= qxM.japura[,2] ,mxM2hat=qxM.manaus[,2])
mx.fit.multivariateWithAlpha<- as.data.frame(mx.fitWithAlpha)



g.japuraJoint<- ggplot() +
  geom_point(data = mx.obs.japura, aes(x = age,y = log(mx),color = regiao),size = 1.7) +
  # geom_point(data = mx.manaus.male, aes(x = age,y = 1-exp(-mx)),size = 1.7) +
  geom_line(data = out.japura,aes(x = age,y = log(qx.fitted), color = regiao),linewidth = 1.2) +
  geom_ribbon(data = out.japura,aes(x = age,ymin = log(qx.lower),ymax = log(qx.upper),fill = regiao),
              alpha = 0.25,show.legend = FALSE) +
  geom_line(data = qxM.manaus,aes(x = age, y = log(qx.fitted),color = "Manaus (male)"),linetype = "dashed",linewidth = 1) +
  # geom_line(data = qxF.manaus,aes(x = age, y = 1-exp(-qx.fitted),color = "Manaus"),linetype = "dashed",linewidth = 1) +
  scale_color_manual( values = c("Female" = "tomato","Male" = "steelblue","Manaus (male)" = "black")) +
  scale_fill_manual( values = c("Female" = "tomato","Male" = "steelblue") ) +
  scale_y_continuous( expression("log rate"),limits=c(-11,0), breaks=c(-10, -7.5, -5, -2.5,0)) +
  scale_x_continuous("Age",breaks = seq(0, 110, by = 10)) +
  labs(color = NULL) +
  theme_classic(base_size = 20) +
  facet_wrap(~local) +
  theme( legend.position = c(0.70, 0.18),  legend.title = element_blank(),strip.background = element_rect(
    colour = "black",fill = "gray87" ),panel.border = element_rect(color = "black",fill = NA))
g.japuraJoint


## fit multivariate model without common term
V5= cov(y.japura)
V5
deltaJ5noalpha<- cbind(d_femaleJ, d_maleJ, d_maleM)
J=3#
res4 <- buildFtGt(J, alpha= FALSE)
Ft4 <- res4$Ft
Gt4 <- res4$Gt
m04 = rep(0, nrow(Gt4))
C04 = diag(100, nrow(Gt4))
s04 <- diag(J)*0.01  ###priori do petris vaga
v04 <- 5


fit.joint.japura.without.alpha<- dlm.multivariate(y = y.japura, Ft = Ft4, Gt = Gt4, nit = 30000, bn = 10000, thin = 1,
                                                  v0 = v04, s0 = s04, m0 = m04, C0 = C04, delta = deltaJ5noalpha, V = V5)


qxF.japura <- data.frame(age=seq(0,80, by=5) ,qx_fitted(fit.joint.japura.without.alpha)[[1]])
qxM.japura <- data.frame(age=seq(0,80, by=5) , qx_fitted(fit.joint.japura.without.alpha)[[2]])
qxM.manaus <- data.frame(age=seq(0,80, by=5) , qx_fitted(fit.joint.japura.without.alpha)[[3]])

mx.fit.without.alpha<- data.frame(mxFhat= qxF.japura[,2], mxMhat= qxM.japura[,2] ,mxM2hat=qxM.manaus[,2])
mx.fit.without.alpha<- as.data.frame(mx.fit.without.alpha)


## fit univariate model
y.japura1<-  data.frame(yF.japura, yM.japura)
y.japura1<- as.matrix(y.japura1)

Vunivariate<- diag(2)
deltaJunivariate<- cbind(d_femaleJ, d_maleJ)
J=2# 4
res3 <- buildFtGt(J, alpha= FALSE)
Ft3 <- res3$Ft
Gt3 <- res3$Gt
m03 = rep(0, nrow(Gt3))
C03 = diag(100, nrow(Gt3))
s03 <- diag(J)*0.01 
v03 <- 5

fit.joint.japura.univariate<- dlm.multivariate(y = y.japura1, Ft = Ft3, Gt = Gt3, nit = 30000, bn = 10000, thin = 1,
                                              v0 = v03, s0 = s03, m0 = m03, C0 = C03, delta = deltaJunivariate, V = Vunivariate)

qxF.japura <- data.frame(age=seq(0,80, by=5) ,qx_fitted(fit.joint.japura.univariate)[[1]])
qxM.japura <- data.frame(age=seq(0,80, by=5) , qx_fitted(fit.joint.japura.univariate)[[2]])

mx.fit.univariate<- data.frame(mxFhat= qxF.japura[,2], mxMhat= qxM.japura[,2])
mx.fit.univariate<- as.data.frame(mx.fit.univariate)


####### model comparison
#LPPD
ls.univariate <- log_pred_density(y.japura1, fit.joint.japura.univariate$mu,fit.joint.japura.univariate$V)
ls.multivariate <- log_pred_density(y.japura, fit.joint.japura.without.alpha$mu,fit.joint.japura.without.alpha$V)
ls.multivariateWithAlpha <- log_pred_density(y.japura, fit.joint.japura $mu, fit.joint.japura$V)


data.frame(
  model = c("model 1",
            "model 2",
            "model 3"),
  logscore = round(c(ls.univariate, ls.multivariate, ls.multivariateWithAlpha),3)
)


### weigthed metrics
ExM.japura <- extract_var(data = df.agregado,var = "Ex",regions = "JAPURA",sex = "Male")
ExF.japura <- extract_var(data = df.agregado,var = "Ex",regions = "JAPURA",sex = "Female")
ExM.manaus <- extract_var(data = df.agregado,var = "Ex",regions = "MANAUS",sex = "Male")

mxM.japura <- extract_var(data = df.agregado,var = "mx",regions = "JAPURA",sex = "Male")
mxF.japura <- extract_var(data = df.agregado,var = "mx",regions = "JAPURA",sex = "Female")
mxM.manaus <- extract_var(data = df.agregado,var = "mx",regions = "MANAUS",sex = "Male")

Ex.all<- data.frame(ExF= ExF.japura, ExM= ExM.japura, ExM2= ExM.manaus)
mx.all<- data.frame(mxF= mxF.japura, mxM= mxM.japura,mxM2= mxM.manaus)


multivariate.with.alpha<-weighted_metrics(
  obs  = log(mx.all),
  fit  = log(mx.fit.multivariateWithAlpha),
  expo = Ex.all
)
multivariate.with.alpha

multivariate.without.alpha<-
  weighted_metrics(
  obs  = log(mx.all),
  fit  = log(mx.fit.without.alpha),
  expo = Ex.all
)

multivariate.without.alpha

univariate<-weighted_metrics(
  obs  = log(mx.all)[,-3],
  fit  = log(mx.fit.univariate),
  expo = Ex.all[,-3]
)
univariate


########================================================================================== 
## 4.1.3 Pooling information across microregions
#We investigate the joint modelling of male mortality across the Manaus, Rio Negro, and Japurá micro-regions in Amazonas.
#These regions exhibit substantial demographic heterogeneity and markedly different exposure sizes, ranging from the 
#comparatively large and more stable population of Manaus to the sparse mortality data observed in Japurá. 
#This application assesses whether the proposed multivariate framework can pool information across geographically
#related populations while preserving region-specific mortality dynamics.

### model comparison with different discount matrix configuration
## configuration 1
d.rionegro<- c(rep(0.90, 3), rep(0.95,6), rep(0.999, 8))
d.japura<- c(rep(0.90, 3), rep(0.99,5), rep(0.999, 8))
d.manaus<-  c(rep(0.90, 3), rep(0.95,8), rep(0.999, 6))
d.gamma<-   c(rep(0.95, 3), rep(0.99,3), rep(0.99, 11)) 
deltaJ1<-  cbind(d.rionegro, d.japura, d.manaus,d.gamma)

fit.region.all1<- dlm.multivariate(y = y.m, Ft = Ft, Gt = Gt, nit = 30000, bn = 10000, thin = 1,
                                   v0 = v0, s0 = s0, m0 = m0, C0 = C0, delta = deltaJ1, V = V)

### configuration 2
d.gamma<-   rep(0.95, 17)
deltaJ2<-  cbind(d.rionegro, d.japura, d.manaus,d.gamma)
fit.region.all2<- dlm.multivariate(y = y.m, Ft = Ft, Gt = Gt, nit = 30000, bn = 10000, thin = 1,
                                   v0 = v0, s0 = s0, m0 = m0, C0 = C0, delta = deltaJ2, V = V)

### configuration 3: the best model 
d.rionegro<- c(rep(0.90, 3), rep(0.95,6), rep(0.999, 8))
d.japura<- c(rep(0.90, 3), rep(0.99,7), rep(0.999, 7))
d.manaus<-  c(rep(0.90, 3), rep(0.95,8), rep(0.999, 6))
d.gamma<-   rep(0.999, 17)
deltaJ3<-  cbind(d.rionegro, d.japura, d.manaus,d.gamma)
fit.region.all3<- dlm.multivariate(y = y.m, Ft = Ft, Gt = Gt, nit = 30000, bn = 10000, thin = 1,
                                   v0 = v0, s0 = s0, m0 = m0, C0 = C0, delta = deltaJ3, V = V)

###configuration 4
d.gamma<-   c(rep(0.95, 10), rep(0.99, 7)) 
deltaJ4<-  cbind(d.rionegro, d.japura, d.manaus,d.gamma)
fit.region.all4<- dlm.multivariate(y = y.m, Ft = Ft, Gt = Gt, nit = 30000, bn = 10000, thin = 1,
                                   v0 = v0, s0 = s0, m0 = m0, C0 = C0, delta = deltaJ4, V = V)

ls1 <- log_pred_density(y.m, fit.region.all1$mu, fit.region.all1$V)
ls2 <- log_pred_density(y.m, fit.region.all2$mu, fit.region.all2$V)
ls3 <- log_pred_density(y.m, fit.region.all3$mu, fit.region.all3$V)
ls4 <- log_pred_density(y.m, fit.region.all4$mu, fit.region.all4$V)

data.frame(
  model = c("prop 1",
            "prop 2",
            "prop 3",
            "prop 4"),
  logscore = c(ls1, ls2, ls3,ls4)
)


### conditional on best discount matrix configuration 3:
fit.region.all<- fit.region.all3
qxaux <- qx_fitted(fit.region.all)

qx_list <- setNames(
  lapply(qxaux, function(q){
    df <- as.data.frame(q)
    df$age <- age
    df
  }),
  c("Rio Negro","Japurá","Manaus")
)
qx.all <- bind_rows(qx_list, .id = "regiao")


### compute metrics
mx.fitWithAlpha<- data.frame(qx.all$qx.fitted[qx.all$regiao=="Rio Negro"],
                             qx.all$qx.fitted[qx.all$regiao=="Japurá"],
                             qx.all$qx.fitted[qx.all$regiao=="Manaus"])
mx.fit.multivariateWithAlpha<- as.data.frame(mx.fitWithAlpha)


### Figure 7 :multivariate (second row) SSM fits with 95\% predictive intervals for Japurá,
#Manaus and Rio Negro microregions during 2022--2024.
#Dots represent crude mortality rates and solid lines represent fitted mortality curves.
 ggplot(NULL, aes(x = age)) +
  geom_point(data = mx.obs.m, aes(x = age, y = log(mx), color = regiao), size = 1.7) +
  geom_line(data = qx.all, aes(y = log(qx.fitted), color = regiao), linewidth = 0.8) +
  geom_ribbon(data = qx.all, aes(ymin = log(qx.lower), ymax = log(qx.upper), fill = regiao), alpha = 0.25, color = NA) +
  scale_y_continuous(expression("log rate"),limits=c(-10,0), breaks=c(-10, -7.5, -5, -2.5,0)) +
  scale_x_continuous("Age", breaks = seq(0, 80, by = 10)) +
  facet_wrap(~regiao) +
  scale_color_manual(values = c("Rio Negro" = "magenta4", "Japurá" = "darkgreen","Manaus" = "steelblue"
  )) +
  scale_fill_manual(values = c("Rio Negro" = "magenta4","Japurá" = "darkgreen", "Manaus" = "steelblue"
  )) +
  theme_classic(base_size = 20) +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    strip.background = element_rect(colour = "black", fill = "gray87"),
    panel.border = element_rect(color = "black", fill = NA)
  ) 

 ### Figure 8: Comparison of the fitted male mortality curves for the Manaus, Rio Negro, and Japurá micro-regions between
 #ages 0 and 30 years during 2022--2024.
 ggplot(subset(qx.all, age <= 30),aes(x = age, y = log(qx.fitted), color = regiao)) +
   geom_line(linewidth = 1.2) +
   facet_wrap(~"comparison of fitted curves - age < 30") +
   scale_y_continuous(expression("log rate"),limits=c(-10,-2.5), breaks=c(-10, -7.5, -5, -2.5,0)) +
   scale_x_continuous("Age",breaks = seq(0, 30, by = 5)) +
   scale_color_manual(values = c("Rio Negro" = "magenta4", "Japurá" = "darkgreen","Manaus" = "steelblue"
   )) +
   theme_classic(base_size = 20) +
   theme(
     legend.position = c(.80, 0.2),
     legend.title = element_blank(),
     strip.background = element_rect(colour = "black", fill = "gray87"),
     panel.border = element_rect(color = "black", fill = NA)
   ) 
 
 #### ### univariate modelling 
 Gt1<-  matrix(c(1,0,1,1), 2)
 Ft1<- matrix(c(1,0), nrow = 1)
 d.rionegro<- c(rep(0.90, 3), rep(0.95,6), rep(0.999, 8))
 d.japura<- c(rep(0.90, 3), rep(0.99,7), rep(0.999, 7))
 d.manaus<-  c(rep(0.90, 3), rep(0.95,8), rep(0.999, 6))
 d.gamma<-   rep(0.999, 17)
 
 ## manaus 
 fit.manaus<- BayesMortalityPlus::dlm(y.m[,3], Ft= Ft1, Gt= Gt1, delta =d.manaus,
                                      prior = list(m0 = rep(0, nrow(Gt1)), C0 = diag(100, nrow(Gt1))), ages=age)
 ## rio negro 
 fit.rionegro<- BayesMortalityPlus::dlm(y.m[,1], Ft= Ft1, Gt= Gt1, delta =d.rionegro,
                                        prior = list(m0 = rep(0, nrow(Gt1)), C0 = diag(100, nrow(Gt1))), ages=age)
 ## japura
 fit.japura<- BayesMortalityPlus::dlm(y.m[,2], Ft= Ft1, Gt= Gt1, delta =d.japura,
                                      prior = list(m0 = rep(0, nrow(Gt1)), C0 = diag(100, nrow(Gt1))), ages=age)
 
 ### fazer os graficos de comparacao...
 qxaux1 <- fitted(fit.japura)  %>%
   mutate(regiao = "Japurá") 
 qxaux2 <- fitted(fit.manaus) %>%
   mutate(regiao = "Manaus")
 qxaux3 <- fitted(fit.rionegro) %>%
   mutate(regiao = "Rio Negro")
 ### para colocar no face_wrap 
 qxaux1$local<- "Japurá (without pooling)"
 qxaux2$local<- "Manaus (without pooling)"
 qxaux3$local<- "Rio Negro (without pooling)"
 
 qx.all.univariate <- dplyr::bind_rows(qxaux1, qxaux2, qxaux3)
 qx.all.univariate
 
 mx.obs.m$local <- paste0(mx.obs.m$regiao,
                          " (without pooling)") 
 
 ### compute metrics
 mx.fit.univariate<- data.frame(qxaux3$qx.fitted, qxaux1$qx.fitted,qxaux2$qx.fitted)
 mx.fit.univariate<- as.data.frame(mx.fit.univariate)
 
 
 ### Figure 7: univariate (first row) SSM fits with 95\% predictive intervals for Japurá,
 #Manaus and Rio Negro microregions during 2022--2024.
 #Dots represent crude mortality rates and solid lines represent fitted mortality curves.
 ggplot(NULL, aes(x = age)) +
   geom_point(data = mx.obs.m, aes(x = age, y = log(mx), color = regiao), size = 1.7) +
   geom_line(data = qx.all.univariate, aes(y = log(qx.fitted), color = regiao), linewidth = 0.8) +
   geom_ribbon(data = qx.all.univariate, aes(ymin = log(qx.lower), ymax = log(qx.upper), fill = regiao), alpha = 0.25, color = NA) +
   scale_y_continuous(expression("log rate"),limits=c(-10,0), breaks=c(-10, -7.5, -5, -2.5,0)) +
   scale_x_continuous("Age", breaks = seq(0, 80, by = 10)) +
   facet_wrap(~local) +
   scale_color_manual(values = c("Rio Negro" = "magenta4", "Japurá" = "darkgreen","Manaus" = "steelblue"
   )) +
   scale_fill_manual(values = c("Rio Negro" = "magenta4","Japurá" = "darkgreen", "Manaus" = "steelblue"
   )) +
   theme_classic(base_size = 20) +
   theme(
     legend.position = "none",
     legend.title = element_blank(),
     strip.background = element_rect(colour = "black", fill = "gray87"),
     panel.border = element_rect(color = "black", fill = NA)
   ) 
 
 
 
 ####### Table 3: Model comparison based on the log posterior predictive density (LPPD),
 #the exposure-weighted mean absolute error (WMAE), and the exposure-weighted root mean squared error (WRMSE).
 #The best-performing model is highlighted in boldface.
 ### 
 
 d.rionegro<- c(rep(0.90, 3), rep(0.95,6), rep(0.999, 8))
 d.japura<- c(rep(0.90, 3), rep(0.99,7), rep(0.999, 7))
 d.manaus<-  c(rep(0.90, 3), rep(0.95,8), rep(0.999, 6))
 
 ### caso univariado
 Vunivariate<- diag(3)
 deltaJunivariate<- cbind(d.rionegro, d.japura, d.manaus)
 J=3
 res3 <- buildFtGt(J, alpha= FALSE)
 Ft3 <- res3$Ft
 Gt3 <- res3$Gt
 m03 = rep(0, nrow(Gt3))
 C03 = diag(100, nrow(Gt3))
 s03 <- diag(J)*0.01  ###priori do petris vaga
 v03 <- 5
 
 fit.joint.all.univariate<- dlm.multivariate(y = y.m, Ft = Ft3, Gt = Gt3, nit = 30000, bn = 10000, thin = 1,
                                             v0 = v03, s0 = s03, m0 = m03, C0 = C03, delta = deltaJunivariate, V = Vunivariate)
 
 
 
 qxaux <- qx_fitted(fit.joint.all.univariate)
 
 qx_list <- setNames(
   lapply(qxaux, function(q){
     df <- as.data.frame(q)
     df$age <- age
     df
   }),
   c("Rio Negro","Japurá","Manaus")
 )
 
 # Data frame final
 qx.all <- bind_rows(qx_list, .id = "regiao")
 qx.all
 
 
 
 
 ### para as metricas ponderadas
 mx.fit.univariate<- data.frame(qx.all$qx.fitted[qx.all$regiao=="Rio Negro"],
                                qx.all$qx.fitted[qx.all$regiao=="Japurá"],
                                qx.all$qx.fitted[qx.all$regiao=="Manaus"])
 mx.fit.univariate<- as.data.frame(mx.fit.univariate)
 
 ## multivariate modelling without common term
 V5= cov(y.m)
 V5
 deltaJ5noalpha<- cbind(d.rionegro, d.japura, d.manaus)
 J=3#
 res4 <- buildFtGt(J, alpha= FALSE)
 Ft4 <- res4$Ft
 Gt4 <- res4$Gt
 m04 = rep(0, nrow(Gt4))
 C04 = diag(100, nrow(Gt4))
 s04 <- diag(J)*0.01  ###priori do petris vaga
 v04 <- 5
 
 
 fit.joint.all.without.alpha<- dlm.multivariate(y = y.m, Ft = Ft4, Gt = Gt4, nit = 30000, bn = 10000, thin = 1,
                                                v0 = v04, s0 = s04, m0 = m04, C0 = C04, delta = deltaJ5noalpha, V = V5)
 
 
 qxaux <- qx_fitted(fit.joint.all.without.alpha)
 
 qx_list <- setNames(
   lapply(qxaux, function(q){
     df <- as.data.frame(q)
     df$age <- age
     df
   }),
   c("Rio Negro","Japurá","Manaus")
 )
 
 # Data frame final
 qx.all <- bind_rows(qx_list, .id = "regiao")
 qx.all
 
 
 
 
 ### compute metrics
 mx.fit.without.alpha<- data.frame(qx.all$qx.fitted[qx.all$regiao=="Rio Negro"],
                                   qx.all$qx.fitted[qx.all$regiao=="Japurá"],
                                   qx.all$qx.fitted[qx.all$regiao=="Manaus"])
 mx.fit.without.alpha<- as.data.frame(mx.fit.without.alpha)
 
 ### model comparison
 
 #LPPD
 ls.univariate <- log_pred_density(y.m, fit.joint.all.univariate$mu,fit.joint.all.univariate$V)
 ls.multivariate <- log_pred_density(y.m, fit.joint.all.without.alpha$mu,fit.joint.all.without.alpha$V)
 ls.multivariateWithAlpha <- log_pred_density(y.m,fit.region.all$mu, fit.region.all$V)
 
 
 data.frame(
   model = c("model 1",
             "model 2",
             "model 3"),
   logscore = round(c(ls.univariate, ls.multivariate, ls.multivariateWithAlpha),3)
 )

 
 ### weigthed metrics
 ExM. <- extract_var(data = df.agregado,var = "Ex",regions = c("RIO NEGRO","JAPURA", "MANAUS"),sex = "Male")
 mxM. <- extract_var(data = df.agregado,var = "mx",regions = c("RIO NEGRO","JAPURA", "MANAUS"),sex = "Male")
 
 Ex.all<- data.frame(ExM.$`RIO NEGRO`, ExM.$JAPURA,ExM.$MANAUS)
 mx.all<- data.frame(mxM.$`RIO NEGRO`, mxM.$JAPURA, mxM.$MANAUS)
 
 
 multivariate.with.alpha<-weighted_metrics(
   obs  = log(mx.all),
   fit  = log(mx.fit.multivariateWithAlpha),
   expo = Ex.all
 )
 multivariate.with.alpha> multivariate.with.alpha

 
 multivariate.without.alpha<-
   weighted_metrics(
     obs  = log(mx.all),
     fit  = log(mx.fit.without.alpha),
     expo = Ex.all
   )

 
 multivariate.without.alpha
 
 univariate<-weighted_metrics(
   obs  = log(mx.all),
   fit  = log(mx.fit.univariate),
   expo = Ex.all
 )
 
 univariate

 
 
 
 
 
 