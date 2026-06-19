########################################################
############ data read - amazonas dataset 
############       pooling model        
####### authors: figueiredo, lobo, fonseca and alves         
########################################################
library(dplyr)
library(ggplot2)
library(BayesMortalityPlus)
library(patchwork)
#library(tidyverse)
require(stringr)
library(prodlim)
library(purrr)
library(tidyr)
library(scales)
require(gridExtra)
library(MASS)
library(stringi)


grad <- seq(0,110,20)
grad[(grad/5)%%2 != 0] <- ""
point <- format_format(big.mark = " ", decimal.mark = ".", scientific = FALSE)

prep_counts <- function(file, varname, sex, year){
  
  read.table(file,
             sep = ";",
             header = TRUE,
             quote = "\"",
             na.strings = c("-", ""),
             check.names = FALSE) %>%
    
    filter(`Microrregião IBGE` != "Total") %>%
    mutate(across(-`Microrregião IBGE`, as.numeric)) %>%
    mutate(sex = sex, year = year) %>%
    
    pivot_longer(
      cols = -c(`Microrregião IBGE`, sex, year),
      names_to = "age_group",
      values_to = varname
    )
}


# Dx
dx_all <- bind_rows(
  prep_counts("obt_F_2022.txt", "Dx", "Female", 2022),
  prep_counts("obt_M_2022.txt", "Dx", "Male",   2022),
  
  prep_counts("obt_F_2023.txt", "Dx", "Female", 2023),
  prep_counts("obt_M_2023.txt", "Dx", "Male",   2023),
  
  prep_counts("obt_F_2024.txt", "Dx", "Female", 2024),
  prep_counts("obt_M_2024.txt", "Dx", "Male",   2024)
  
)

print(dx_all, n=100)

dx_04 <- dx_all %>%
  filter(age_group %in% c(
    "0 a 6 dias",
    "7 a 27 dias",
    "28 a 364 dias",
    "1 a 4 anos"
  )) %>%
  group_by(`Microrregião IBGE`, sex, year) %>%
  summarise(Dx = sum(Dx, na.rm = TRUE), .groups = "drop") %>%
  mutate(age_group = "0 a 4 anos")


dx_rest <- dx_all %>%
  filter(!age_group %in% c(
    "0 a 6 dias",
    "7 a 27 dias",
    "28 a 364 dias",
    "1 a 4 anos"
  ))


dx_all <- bind_rows(dx_rest, dx_04) %>%
  mutate(age_group = factor(age_group,
                            levels = c(
                              "0 a 4 anos",
                              "5 a 9 anos",
                              "10 a 14 anos",
                              "15 a 19 anos",
                              "20 a 24 anos",
                              "25 a 29 anos",
                              "30 a 34 anos",
                              "35 a 39 anos",
                              "40 a 44 anos",
                              "45 a 49 anos",
                              "50 a 54 anos",
                              "55 a 59 anos",
                              "60 a 64 anos",
                              "65 a 69 anos",
                              "70 a 74 anos",
                              "75 a 79 anos",
                              "80 anos ou mais",
                              "Total")
  )) %>%
  arrange(`Microrregião IBGE`, sex, year, age_group)

dx_all <- dx_all %>%
  filter(age_group != "Total")

print(dx_all, n=100)


# Ex
ex_all <- bind_rows(
  prep_counts("exp_F_2022.txt", "Ex", "Female", 2022),
  prep_counts("exp_M_2022.txt", "Ex", "Male",   2022),
  
  prep_counts("exp_F_2023.txt", "Ex", "Female", 2023),
  prep_counts("exp_M_2023.txt", "Ex", "Male",   2023),
  
  prep_counts("exp_F_2024.txt", "Ex", "Female", 2024),
  prep_counts("exp_M_2024.txt", "Ex", "Male",   2024)
)

print(ex_all, n=100)

ex_all <- ex_all %>%
  filter(age_group != "Total")

ex_all <- ex_all %>%
  mutate(age_group = gsub("^De\\s*", "", age_group))
ex_all

df_counts <- left_join(dx_all, ex_all,
                       by = c("Microrregião IBGE", "sex", "year", "age_group"))

print(df_counts, n=100)


df_counts <- df_counts %>%
  mutate(
    age_start = as.numeric(str_extract(age_group, "\\d+")),
    age = case_when(
      str_detect(age_group, "80") ~ 85,
      TRUE ~ age_start + 2.5
    )
  )
df_counts

######## 

df_long  <- df_counts %>%
  mutate(
    mx = Dx / Ex,
    mx = ifelse(Ex == 0, NA, mx),   # só trata divisão por zero
    log_mx = log(mx),
    zero = is.na(mx)
  ) 


df_long

df_long <- df_long %>%
  rename(Microrregião.IBGE = `Microrregião IBGE`)

print(df_long, n=100)



df_agregado <- df_long %>%
  # mutate(Dx = ifelse(is.na(Dx), 0, Dx)) %>%
  group_by(Microrregião.IBGE, sex, age_group, age_start, age) %>%
  summarise(
    Dx = sum(Dx, na.rm = TRUE),
    Ex = sum(Ex, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    mx = Dx / Ex,
    log_mx = log(mx)
  )


df_agregado <- df_agregado %>%
  arrange(Microrregião.IBGE, sex, age_start) %>%
  mutate(Microrregião.IBGE = str_remove(Microrregião.IBGE, "^\\d+\\s+"))

df_aux <- df_long %>%
  mutate(Dx = ifelse(is.na(Dx), 0, Dx)) %>%
  group_by(Microrregião.IBGE, sex, age_group, age_start, age) %>%
  summarise(
    Dx = sum(Dx),
    Ex = sum(Ex),
    .groups = "drop"
  ) %>%
  mutate(
    mx = Dx / Ex,
    log_mx = log(mx)
  )

df_aux <- df_aux %>%
  arrange(Microrregião.IBGE, sex, age_start) %>%
  mutate(Microrregião.IBGE = str_remove(Microrregião.IBGE, "^\\d+\\s+"))


df.agregado<- df_aux %>%
  mutate(mx = (Dx + 0.5) / Ex, log_mx = log(mx))


print(df.agregado, n=100)
