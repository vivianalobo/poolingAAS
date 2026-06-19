######################################################
#### authors: figueiredo, lobo, alves and fonseca
#### output: data modelling 
####  modelling multiple populations
### data: amazonas, 2022-2024
######################################################

### multivariate dlm (ffbs, gibbs, qx) e data (log.mx, Dx, Ex)
source("amazoniadata.R")
source("multivariatefunc.R")

#### aditional functions
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
make_df <- function(qxaux, nome_delta){
  qx_list <- setNames(
    lapply(qxaux, function(q){
      df <- as.data.frame(q)
      df$age <- age
      df
    }),
    c("RIO NEGRO", "JAPURA", "MANAUS")
  )
  
  df <- dplyr::bind_rows(qx_list, .id = "regiao")
  df$delta <- nome_delta
  return(df)
}




##--------------------------------------------------------------------------------
########## male data amazon
### filtering male 
df.m<- df.agregado %>%
  filter(sex == "Male")
df.m<- df.m %>%
dplyr::select(Microrregião.IBGE, age_start, log_mx) %>%  # seleciona apenas as colunas necessárias
  pivot_wider(
    names_from = Microrregião.IBGE,  # cada microrregião vira uma coluna
    values_from = log_mx             # valores das células
  ) %>%
  arrange(age_start)  # ordena por idade

### mortality rate
nomes_regioes <- c("RIO NEGRO","JAPURA","MANAUS")
mx.obs.m <- df.agregado %>%
  filter(sex == "Male",  Microrregião.IBGE %in% nomes_regioes) 
mx.obs.m <-  mx.obs.m[, c(regiao="Microrregião.IBGE", "age_start", "mx")] 
colnames(mx.obs.m)<- c("regiao", "age", "mx")
mx.obs.m<- mx.obs.m %>%
  mutate(regiao = dplyr::recode(regiao,
                                "RIO NEGRO" = "Rio Negro",
                                "JAPURA" = "Japurá",
                                "MANAUS" = "Manaus"))

### Rio Negro, Japura e Manaus - lograte
y.m<- as.matrix(df.m[,c(2,3,8)])

## compute covariance matrix - male 
y_clean <- y.m[complete.cases(y.m), ]
V <- cov(y_clean)
V <- make_pd(V)
V <- V + diag(1e-3, nrow(V)) 
V

## setting age
age<- seq(0,80, by=5)

##------------------------------------------------------------------------
### run multivarite model w/ commun term
# dimensão do estado: (mu_1,...,mu_J, beta_1,...,beta_J, gamma)
J <- 3
res <- buildFtGt(J, alpha= TRUE)
Ft <- res$Ft
Gt <- res$Gt
dim(Gt) 
dim(Ft)  

m0 = rep(0, nrow(Gt))
C0 = diag(10, nrow(Gt))
s0 <- diag(J)*0.01  ###priori do petris vaga
v0 <- 5

###### section 3.4 Age-varying smoothness
##  Figure~\ref{fig:figagevarying} illustrates this mechanism for {male and female} mortality SSM fit in Manaus region,
#comparing three alternative specifications for discount factors.
#The first two scenarios apply a single discount factor for all ages {and populations}, with $\delta= 0.75$ and $\delta = 0.999$, respectively.
#A lower discount factor allows the model to respond rapidly to local fluctuations, producing a more irregular fit, 
#{in this case, illustrating potential over-fitting}, whereas a value very close to one enforces strong smoothness and results in an almost
#linear trajectory across ages. The third scenario adopts age-varying discount factors, assigning higher values of $\delta^{(j)}_x$ for
#age ranges in which mortality is expected to evolve smoothly, and lower values to ages characterised by sharper changes. 
#This approach ensures that the fitted curve remains smooth over adult ages while retaining sufficient flexibility at younger ages, 
#where mortality rates change more abruptly and a more realistic uncertainty structure. 

d.rionegro= d.japura= d.manaus= d.gamma= rep(0.75, nrow(y.m))
deltaJ1<-  cbind(d.rionegro, d.japura, d.manaus,d.gamma)
dim(deltaJ1)

d.rionegro= d.japura= d.manaus= d.gamma= rep(0.999, nrow(y.m)) 
deltaJ2<-  cbind(d.rionegro, d.japura, d.manaus,d.gamma)
dim(deltaJ2)

d.rionegro= d.japura= d.manaus= d.gamma= c(rep(0.85, 3), rep(0.95,8), rep(0.999, 6))
deltaJ3<-  cbind(d.rionegro, d.japura, d.manaus,d.gamma)
dim(deltaJ3)


fit.region.m1<- dlm.multivariate(y = y.m, Ft = Ft, Gt = Gt, nit = 30000, bn = 10000, thin = 1,
                                   v0 = v0, s0 = s0, m0 = m0, C0 = C0, delta = deltaJ1, V = V)

fit.region.m2<- dlm.multivariate(y = y.m, Ft = Ft, Gt = Gt, nit = 30000, bn = 10000, thin = 1,
                                   v0 = v0, s0 = s0, m0 = m0, C0 = C0, delta = deltaJ2, V = V)

fit.region.m3<- dlm.multivariate(y = y.m, Ft = Ft, Gt = Gt, nit = 30000, bn = 10000, thin = 1,
                                   v0 = v0, s0 = s0, m0 = m0, C0 = C0, delta = deltaJ3, V = V)


qxaux1 <- qx_fitted(fit.region.m1)
qxaux2 <- qx_fitted(fit.region.m2)
qxaux3 <- qx_fitted(fit.region.m3)

make_df <- function(qxaux, nome_delta){
  qx_list <- setNames(
    lapply(qxaux, function(q){
      df <- as.data.frame(q)
      df$age <- age
      df
    }),
    c("RIO NEGRO", "JAPURA", "MANAUS")
  )
  
  df <- dplyr::bind_rows(qx_list, .id = "regiao")
  df$delta <- nome_delta
  return(df)
}

df3 <- make_df(qxaux3, "different per age")
df2 <- make_df(qxaux2, "0.999")
df1 <- make_df(qxaux1, "0.75")
qx.all <- dplyr::bind_rows(df1, df2, df3)
qx.all

#save.image("outputMale3Populations.RData")

### olhando os tres ajustes
pdf("Fig1PoolingMale.pdf", width=12, height=5)
ggplot(NULL, aes(x = age)) +
  geom_point(data = mx.obs.m, aes(x = age, y = mx, color = regiao), size = 1.7) +
  geom_line(data = qx.all, aes(y = 1-exp(-qx.fitted), color = regiao), linewidth = 0.8) +
  geom_ribbon(data = qx.all, aes(ymin = 1-exp(-qx.lower), ymax = 1-exp(-qx.upper), fill = regiao), alpha = 0.25, color = NA) +
  scale_y_continuous(expression(q[x]), trans = 'log10',
                     limits= c(1e-5, 1),
                     breaks = c(1e-5,1e-4,1e-3, 1e-2, 1e-1, 1),
                     labels = scales::label_scientific()) +
  scale_x_continuous("Age", breaks = seq(0, 80, by = 10)) +
  facet_wrap(~delta) +
  scale_color_manual(values = c("RIO NEGRO" = "magenta4", "JAPURA" = "darkgreen","MANAUS" = "steelblue"
  )) +
  scale_fill_manual(values = c("RIO NEGRO" = "magenta4","JAPURA" = "darkgreen", "MANAUS" = "steelblue"
  )) +
  theme_classic(base_size = 20) +
  theme(
    legend.position = c(0.9, 0.15),
    legend.title = element_blank(),
    strip.background = element_rect(colour = "black", fill = "gray87"),
    panel.border = element_rect(color = "black", fill = NA)
  ) 


plot.popJ <- function(rate, fit, population, colorreg){ 
  
  g1 <- ggplot(NULL, aes(x = age)) +
    
    geom_point(data = rate %>% dplyr::filter(regiao == population),aes(y = mx, color = regiao),size = 1.7) +
    geom_line(data = fit %>% dplyr::filter(regiao == population),aes(y = 1 - exp(-qx.fitted), color = regiao),linewidth = 0.8) +
    geom_ribbon(data = fit %>% dplyr::filter(regiao == population),aes(ymin = 1 - exp(-qx.lower),
        ymax = 1 - exp(-qx.upper),fill = regiao),alpha = 0.25,color = NA) +
    facet_wrap(~delta) +
    scale_y_continuous(expression(q[x]),trans = "log10",limits = c(1e-5, 1),
      breaks = c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1),labels = scales::label_scientific()) +
     scale_x_continuous("Age",breaks = seq(0, 80, by = 10)) +
     scale_color_manual(values = setNames(colorreg, population)) +
     scale_fill_manual(values = setNames(colorreg, population)) +
     theme_classic(base_size = 20) +
    theme(legend.position = c(0.9, 0.15),legend.title = element_blank(),strip.background = element_rect(
        colour = "black",fill = "gray87"),panel.border = element_rect(color = "black",fill = NA))
  
  return(g1)
}

### manaus 
#pdf("Fig2aPoolingMaleManaus.pdf", width=12, height=5)
p1<- plot.popJ(rate = mx.obs.m,fit = qx.all,population = "MANAUS",colorreg = "steelblue")
p1
graphics.off()
### japura 
pdf("Fig2aPoolingMaleJapura.pdf", width=12, height=5)
p2<- plot.popJ(rate = mx.obs.m,fit = qx.all,population = "JAPURA",colorreg = "darkgreen")
p2
graphics.off()
### rio negro 
pdf("Fig2aPoolingMaleRioNegro.pdf", width=12, height=5)
p3<- plot.popJ(rate = mx.obs.m,fit = qx.all,population = "RIO NEGRO",colorreg = "magenta4")
p3
graphics.off()

library(patchwork)
p1 / p2 / p3


# descontos altos para alpha_x implicam que o padrão comum evolui lentamente;
# favorecendo forte borrowing strength; e produzindo curvas mais coerentes entre populações.
# descontos menores para alpha_x permitem que a estrutura compartilhada varie mais rapidamente;
# reduzindo a influência comum; e permitindo maior heterogeneidade entre populações.
# 
# Então o desconto do termo comum controla, na prática:
#   - o grau de persistência da informação compartilhada;
#   - e a intensidade do borrowing strength ao longo das idades.


#-------------------------------------------------------------------------
#### mudando o desconto para todas as três populacoes e fator comum - na ordem
##### E COMPARANDO OS MODELOS VIA LOG PREDICITVE DENSITY... QUAL O MELHOR EM TERMO DE DESCONTO?
# d.rionegro<- c(rep(0.90, 3), rep(0.95,6), rep(0.999, 8))
# d.japura<- c(rep(0.90, 3), rep(0.99,7), rep(0.999, 7))
# d.manaus<-  c(rep(0.90, 3), rep(0.95,8), rep(0.999, 6))

d.rionegro<- c(rep(0.90, 3), rep(0.95,6), rep(0.999, 8))
d.japura<- c(rep(0.90, 3), rep(0.99,5), rep(0.999, 8))
d.manaus<-  c(rep(0.90, 3), rep(0.95,8), rep(0.999, 6))
d.gamma<-   c(rep(0.95, 3), rep(0.99,3), rep(0.99, 11)) 
deltaJ1<-  cbind(d.rionegro, d.japura, d.manaus,d.gamma)
dim(deltaJ1)
fit.region.all1<- dlm.multivariate(y = y.m, Ft = Ft, Gt = Gt, nit = 30000, bn = 10000, thin = 1,
                                     v0 = v0, s0 = s0, m0 = m0, C0 = C0, delta = deltaJ1, V = V)

d.gamma<-   rep(0.95, 17)

deltaJ2<-  cbind(d.rionegro, d.japura, d.manaus,d.gamma)
dim(deltaJ2)
fit.region.all2<- dlm.multivariate(y = y.m, Ft = Ft, Gt = Gt, nit = 30000, bn = 10000, thin = 1,
                                   v0 = v0, s0 = s0, m0 = m0, C0 = C0, delta = deltaJ2, V = V)

### melhor modelo
d.rionegro<- c(rep(0.90, 3), rep(0.95,6), rep(0.999, 8))
d.japura<- c(rep(0.90, 3), rep(0.99,7), rep(0.999, 7))
d.manaus<-  c(rep(0.90, 3), rep(0.95,8), rep(0.999, 6))
d.gamma<-   rep(0.999, 17)

deltaJ3<-  cbind(d.rionegro, d.japura, d.manaus,d.gamma)
dim(deltaJ3)
fit.region.all3<- dlm.multivariate(y = y.m, Ft = Ft, Gt = Gt, nit = 30000, bn = 10000, thin = 1,
                                   v0 = v0, s0 = s0, m0 = m0, C0 = C0, delta = deltaJ3, V = V)


d.gamma<-   c(rep(0.95, 10), rep(0.99, 7)) 

deltaJ4<-  cbind(d.rionegro, d.japura, d.manaus,d.gamma)
dim(deltaJ4)
fit.region.all4<- dlm.multivariate(y = y.m, Ft = Ft, Gt = Gt, nit = 30000, bn = 10000, thin = 1,
                                   v0 = v0, s0 = s0, m0 = m0, C0 = C0, delta = deltaJ4, V = V)


### comparando diversos descontos.
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
# model   logscore
# 1 prop 1   6.032237
# 2 prop 2 -13.334267
# 3 prop 3   6.762527
# 4 prop 4 -14.309257

# saveRDS(fit.region.all3, "outModelAmazon3pop.rds")

### sem termo comum para deltaJ1

### com termo comum para deltaJ1



### A PARTIR DA ESCOLHA DO MELHOR MODELO
### leitura dos dados acima e funcao Gt Ft
fit.region.all3<- readRDS("outModelAmazon3pop.rds")

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

# Data frame final
qx.all <- bind_rows(qx_list, .id = "regiao")
qx.all




### para as metricas ponderadas
mx.fitWithAlpha<- data.frame(qx.all$qx.fitted[qx.all$regiao=="Rio Negro"],
                             qx.all$qx.fitted[qx.all$regiao=="Japurá"],
                             qx.all$qx.fitted[qx.all$regiao=="Manaus"])
mx.fit.multivariateWithAlpha<- as.data.frame(mx.fitWithAlpha)

g1 <- ggplot(NULL, aes(x = age)) +
  geom_point(data = mx.obs.m, aes(x = age, y = mx, color = regiao), size = 1.7) +
  geom_line(data = qx.all, aes(y = 1-exp(-qx.fitted), color = regiao), linewidth = 0.8) +
  geom_ribbon(data = qx.all, aes(ymin = 1-exp(-qx.lower), ymax = 1-exp(-qx.upper), fill = regiao), alpha = 0.25, color = NA) +
  scale_y_continuous(expression(q[x]), trans = 'log10',
                     limits= c(1e-4, 1),
                     breaks = c(1e-4,1e-3, 1e-2, 1e-1, 1),
                     labels = scales::label_scientific()) +
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

pdf("FigPoolingMaleAmazon3pop.pdf", width=12, height=5)
### rio negro, japura, manaus
g1
graphics.off()


g1lograte <- ggplot(NULL, aes(x = age)) +
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

pdf("FigPoolingMaleAmazon3popLogScale.pdf", width=12, height=5)
### rio negro, japura, manaus
g1lograte
graphics.off()


pdf("FigPoolingMaleAmazon3popAge30.pdf", width=6, height=6)
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
graphics.off()

### univariate modelling - usando o mesmo fator de desconto do multivariato (neste caso nao tem para o termo comum)
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

### para as metricas ponderadas
mx.fit.univariate<- data.frame(qxaux3$qx.fitted, qxaux1$qx.fitted,qxaux2$qx.fitted)
mx.fit.univariate<- as.data.frame(mx.fit.univariate)



g2 <- ggplot(NULL, aes(x = age)) +
 geom_point(data = mx.obs.m, aes(x = age, y = mx, color = regiao), size = 1.7) +
  geom_line(data = qx.all.univariate, aes(y = 1-exp(-qx.fitted), color = regiao), linewidth = 0.8) +
  geom_ribbon(data = qx.all.univariate, aes(ymin = 1-exp(-qx.lower), ymax = 1-exp(-qx.upper), fill = regiao), alpha = 0.25, color = NA) +
  scale_y_continuous(expression(q[x]), trans = 'log10',
                     limits= c(1e-4, 1),
                     breaks = c(1e-4,1e-3, 1e-2, 1e-1, 1),
                     labels = scales::label_scientific()) +
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

pdf("FigWithoutPoolingMaleAmazon3pop.pdf", width=12, height=5)
g2
graphics.off()


g2lograte <- ggplot(NULL, aes(x = age)) +
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


pdf("FigWithoutPoolingMaleAmazon3popLogScale.pdf", width=12, height=5)
g2lograte
graphics.off()



### comparacao de modelos univariado, multivariado, multivariado com termo comum


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

## multivariado sem termo comum 
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




### para as metricas ponderadas
mx.fit.without.alpha<- data.frame(qx.all$qx.fitted[qx.all$regiao=="Rio Negro"],
                             qx.all$qx.fitted[qx.all$regiao=="Japurá"],
                             qx.all$qx.fitted[qx.all$regiao=="Manaus"])
mx.fit.without.alpha<- as.data.frame(mx.fit.without.alpha)

### comparando os modelos 

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
# model logscore
# 1 model 1  -14.023
# 2 model 2  -13.981
# 3 model 3    6.763




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
# Population       WLMSE       WMAE
# 1 mxM...RIO.NEGRO. 0.023723469 0.12533417
# 2      mxM..JAPURA 0.240097663 0.39906818
# 3      mxM..MANAUS 0.003503916 0.04810075
# 4            Total 0.006732307 0.05483414

multivariate.without.alpha<-
  weighted_metrics(
    obs  = log(mx.all),
    fit  = log(mx.fit.without.alpha),
    expo = Ex.all
)
# multivariate.without.alpha
# Population     WLMSE      WMAE
# 1 mxM...RIO.NEGRO. 0.1110968 0.2753643
# 2      mxM..JAPURA 0.2196518 0.3586652
# 3      mxM..MANAUS 0.1111267 0.2875995
# 4            Total 0.1122246 0.2878156

multivariate.without.alpha

univariate<-weighted_metrics(
  obs  = log(mx.all),
  fit  = log(mx.fit.univariate),
  expo = Ex.all
)

univariate
# univariate
# Population     WLMSE      WMAE
# 1 mxM...RIO.NEGRO. 0.1118165 0.2758832
# 2      mxM..JAPURA 0.2162579 0.3561333
# 3      mxM..MANAUS 0.1110609 0.2873513
# 4            Total 0.1121574 0.2875759


















### old 
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
################# outro ajuste
## pooling para japura pegando todas as outras regioes (sem manaus) e emprestando informacao para japura

df.comp <- bind_rows(df_manaus, df_sem_manaus) %>%
  arrange(sex, grupo, age_start)
df.comp

# pdf("FigAggSemManaus.pdf",  width=10, height=6)
ggplot(df.comp, aes(x = age, y = 1-exp(-mx),color = grupo, linetype = grupo)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  facet_wrap(~sex) +
 # scale_y_continuous(trans = "log10",  labels = scales::comma,name = expression(m[x])) +
  scale_y_continuous(
    expression(q[x]),
    trans = 'log10',
    limits = c(1e-4, 1),
    breaks = c(1e-4,1e-3, 1e-2, 1e-1, 1),
    labels = scales::label_scientific()
  ) +
  scale_x_continuous("Idade", breaks = seq(0, 90, by = 10)) +
  scale_color_manual(values = c("MANAUS" = "black","Others" = "gray50")) +
  scale_linetype_manual(values = c("MANAUS" = "solid","Others" = "dashed")) +
  theme_classic(base_size = 20) +
  theme(legend.position = c(0.9, 0.15),
        legend.title = element_blank(),legend.text = element_text(size = 12),
        strip.background = element_rect(colour = "black", fill = "gray87"),
        strip.text = element_text(size = 10),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
# graphics.off()



### japura y.m[,2]
aux1<- df_sem_manaus %>%
  filter(sex=="Male")
df.others<- as.vector(aux1$log_mx) 
y.m 
y.m.Jothers<- data.frame(y.m, OTHERS= df.others) 
y.m.Jothers<- as.matrix(y.m.Jothers)
head(y.m.Jothers)
is.matrix(y.m.Jothers)

y.m.J<- y.m.Jothers[,c(2,4)]
V1<- cov(y.m.J)
V1

### run multivarite model w/ commun term
J=2
res <- build_Ft_Gt(J)
Ft <- res$Ft
Gt <- res$Gt

m0 = rep(0, nrow(Gt))
C0 = diag(100, nrow(Gt))
s0 <- diag(J)*0.01  ###priori do petris vaga
v0 <- 5
d3<- c(rep(0.85, 3), rep(0.95,8), rep(0.999, 6))

fit.region.m4<- dlm.multivariate(y = y.m.J, Ft = Ft, Gt = Gt, nit = 30000, bn = 10000, thin = 1,
                                 v0 = v0, s0 = s0, m0 = m0, C0 = C0, delta = d3, V = V1)

age<- seq(0,80, by=5)
qxaux <- qx_fitted(fit.region.m4)

qx_list <- setNames(
  lapply(qxaux, function(q){
    df <- as.data.frame(q)
    df$age <- age
    df
  }),
  c("JAPURA", "baseline")
)


qx_list

qx.all <- dplyr::bind_rows(qx_list, .id = "regiao")
head(qx.all)

### ATENCAO AQUI!
### adicionar os pontos aqui.... do baseline e japura..
### pensando tambem em aidionar a curva de manaus para ver o nivel...
#pdf("Fig2MotivationFemaleRioNegro.pdf",  width=10, height=6) 
ggplot() +
  geom_line(data = qx.all,aes(x = age, y = 1-exp(-qx.fitted), color= regiao), linewidth = 1.5) +
  geom_ribbon(data = qx.all, aes(x = age, ymin = 1-exp(-qx.lower), ymax = 1-exp(-qx.upper),fill = regiao), alpha = 0.2) +
  scale_y_continuous(expression(q[x]),trans = 'log10',
                     limits= c(1e-6,1),
                     breaks = c(1e-6, 1e-5, 1e-4, 1e-3,1e-2,  1e-1, 1),
                     labels = scales::label_scientific())+
  scale_x_continuous("Age", breaks = seq(0, 80, by = 10)) +
  scale_color_manual(values = c("JAPURA" = "darkgreen","baseline"= "darkblue")) +
  scale_fill_manual(values = c("JAPURA" = "darkgreen", "baseline"= "darkblue")) +
  facet_wrap(~regiao, ncol = 2) +
  theme_classic(base_size = 20) +
  theme(legend.position = c(0.9, 0.15),
        legend.title = element_blank(),
        strip.background = element_rect(colour = "black", fill = "gray87"),
        panel.border = element_rect(color = "black", fill = NA)
  )
#graphics.off()
ggplot(qx.all, aes(x = age)) +
  
  # intervalo de credibilidade
  geom_ribbon(aes(
    ymin = 1 - exp(-qx.lower),
    ymax = 1 - exp(-qx.upper),
    fill = regiao
  ), alpha = 0.2) +
  
  # curva principal
  geom_line(aes(
    y = 1 - exp(-qx.fitted),
    color = regiao
  ), linewidth = 1.5) +
  
  scale_y_continuous(
    expression(q[x]),
    trans = 'log10',
    limits = c(1e-6, 1),
    breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1),
    labels = scales::label_scientific()
  ) +
  
  scale_x_continuous("Age", breaks = seq(0, 80, by = 10)) +
  
  scale_color_manual(values = c(
    "JAPURA" = "darkgreen",
    "baseline" = "darkblue"
  )) +
  
  scale_fill_manual(values = c(
    "JAPURA" = "darkgreen",
    "baseline" = "darkblue"
  )) +
  
  theme_classic(base_size = 20) +
  
  theme(
    legend.position = c(0.85, 0.2),
    legend.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA)
  )



#### old



### constroi Ft e Gt para qualquer J com termo comum
# build_Ft_Gt <- function(J){
#   m <- 2*J + 1
#   Gt <- matrix(0, nrow = m, ncol = m)
#   for (j in 1:J) {
#     Gt[j, j] <- 1
#     Gt[j, J + j] <- 1
#   }
#   for (j in 1:J) {
#     Gt[J + j, J + j] <- 1
#   }
#   Gt[m, m] <- 1
#   
#   Ft <- matrix(0, nrow = J, ncol = m)
#   for (j in 1:J) {
#     Ft[j, j] <- 1        # mu_j
#     Ft[j, m] <- 1        # gamma
#   }
#   return(list(Ft = Ft, Gt = Gt))
# }


