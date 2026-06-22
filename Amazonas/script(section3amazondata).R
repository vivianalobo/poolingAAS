######################################################
#### authors: figueiredo, lobo, alves and fonseca
#### output: age-varying smoothness
### data: amazonas, 2022-2024
######################################################

#### section3: Amazonas mortality data

### run amazon dataset and auxiliary functions 
source("amazoniadata.R")
source("multivariatefunc.R")
source("predict_dlm_at.R")
source("predict_bivar2_att_chain.R")
source("funcmetrics.R")
source("auxfunc.R")



#--------------------------------------------------------------------------------------------------------------------------------------
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


### delta = 0.75
d_male  = d_female=  rep(0.75, length(age))
d_alpha = rep(0.95, length(age))
deltaJ1<- cbind(d_female, d_male, d_alpha)
deltaJ1<- cbind(d_female, d_male, d_alpha)
### delta= 0.99
d_male  = d_female= d_alpha = rep(0.999, length(age))
deltaJ2<- cbind(d_female, d_male, d_alpha)
### delta different per age
d_male   <- c(rep(0.90, 4), rep(0.95,7), rep(0.999, 6))
d_female <- c(rep(0.90, 3), rep(0.95,8), rep(0.999, 6))
d_alpha  <- c(rep(0.95, 5), rep(0.999, 12))

deltaJ3<- cbind(d_female, d_male, d_alpha)

fit.manaus1<- dlm.multivariate(y = y.manaus, Ft = Ft, Gt = Gt, nit = 30000, bn = 10000, thin = 1,
                               v0 = v0, s0 = s0, m0 = m0, C0 = C0, delta = deltaJ1, V = V4)
fit.manaus2<- dlm.multivariate(y = y.manaus, Ft = Ft, Gt = Gt, nit = 30000, bn = 10000, thin = 1,
                               v0 = v0, s0 = s0, m0 = m0, C0 = C0, delta = deltaJ2, V = V4)
fit.manaus3<- dlm.multivariate(y = y.manaus, Ft = Ft, Gt = Gt, nit = 30000, bn = 10000, thin = 1,
                               v0 = v0, s0 = s0, m0 = m0, C0 = C0, delta = deltaJ3, V = V4)



qxaux1 <- qx_fitted(fit.manaus1)
qxaux2 <- qx_fitted(fit.manaus2)
qxaux3 <- qx_fitted(fit.manaus3)

make_df <- function(qxaux, nome_delta){
  qx_list <- setNames(
    lapply(qxaux, function(q){
      df <- as.data.frame(q)
      df$age <- age
      df
    }),
    c("Female", "Male")
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
names(qx.all)




#### plotting in log-scale
ggplot(qx.all, aes(x = age)) +
  geom_point(data = mx.obs, aes(x = age, y = log(mx), color = regiao), size = 1.7) +
  geom_ribbon(aes( ymin = log(qx.lower), ymax = log(qx.upper), fill = regiao), alpha = 0.2) +
  geom_line(aes(y = log(qx.fitted),color = regiao), linewidth = 1.2) +
  scale_y_continuous(expression("log rate")) +
  scale_x_continuous("Age", breaks = seq(0, 80, by = 10)) +
  scale_color_manual(values = c("Female" = "tomato","Male" = "steelblue")) +
  scale_fill_manual(values = c("Female" = "tomato","Male" = "steelblue")) +
  theme_classic(base_size = 20) +
  facet_wrap(~delta)+
  theme(legend.position = c(0.9, 0.15),legend.title = element_blank(), 
        strip.background = element_rect(colour = "black", fill = "gray87"),
        panel.border = element_rect(color = "black", fill = NA))
