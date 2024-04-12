rm(list=ls())
par(mfrow=c(1,1))

setwd("/analysis")


library("ggplot2")
library("purrr")
library("tidyr")
library("stringr")
library("reshape2")
library("insect")
library("seqinr")
library("adegenet")
library("tools")
library("latexpdf")
library("viridis")
library("MASS")
library("gridExtra")
library("multcomp")
library("lattice")
library("car")
library("Matrix")
library("glmmTMB")
library("DHARMa")



A1_vs_A2 <- read.csv("A1_vs_A2.csv")
A1_vs_A3 <- read.csv("A1_vs_A3.csv")
A1_vs_A4 <- read.csv("A1_vs_A4.csv")

A1_vs_B1 <- read.csv("A1_vs_B1.csv")
A1_vs_B2 <- read.csv("A1_vs_B2.csv")
A1_vs_B3 <- read.csv("A1_vs_B3.csv")
A1_vs_B4 <- read.csv("A1_vs_B4.csv")

A1_vs_C1 <- read.csv("A1_vs_C1.csv")
A1_vs_C2 <- read.csv("A1_vs_C2.csv")
A1_vs_C3 <- read.csv("A1_vs_C3.csv")
A1_vs_C4 <- read.csv("A1_vs_C4.csv")

A1_vs_D1 <- read.csv("A1_vs_D1.csv")
A1_vs_D2 <- read.csv("A1_vs_D2.csv")
A1_vs_D3 <- read.csv("A1_vs_D3.csv")
A1_vs_D4 <- read.csv("A1_vs_D4.csv")


all_variants <- read.csv("all_variants_HPV16.csv")



df_list <- list(all_variants, A1_vs_A2, A1_vs_A3, A1_vs_A4, A1_vs_B1, A1_vs_B2, A1_vs_B3, A1_vs_B4, A1_vs_C1, A1_vs_C2, A1_vs_C3, A1_vs_C4, A1_vs_D1, A1_vs_D2, A1_vs_D3, A1_vs_D4)
df_list2 <- lapply(df_list, function(x) x[, names(x) != "X"])
df_list3 <- lapply(df_list2, function(x) x[, names(x) != "diff"])





df_variants_all<-Reduce(function(x, y) merge(x, y, by.x = c("reference", "position"), 
                            by.y = c("HPV16_A1", "pos"), all.x = TRUE), 
       df_list3)


df_variants_all <- df_variants_all %>%
  mutate(across(starts_with("HPV16_"), as.character),
         minor_base = as.character(minor_base))



df_check <- df_variants_all %>%
  mutate(across(starts_with("HPV16_"), ~ ifelse(. == minor_base, 1, 0)))


df_check_no_na <- replace(df_check, is.na(df_check), 0)

df_check_no_na_select<- df_check_no_na %>% dplyr::select(position, minor_base, HPV16_A2, HPV16_A3)


df_sum <- df_check_no_na  %>%
  group_by(sample_id2) %>%
  summarise(across(starts_with("HPV16_"), sum))
                   


## Same with major variants

df_variants_all_major<-Reduce(function(x, y) merge(x, y, by.x = c("reference", "position"), 
                                             by.y = c("HPV16_A1", "pos"), all.x = TRUE), 
                        df_list3)


df_variants_all_major <- df_variants_all_major %>%
  mutate(across(starts_with("HPV16_"), as.character),
         major_base = as.character(major_base))


df_check_major <- df_variants_all_major %>%
  mutate(across(starts_with("HPV16_"), ~ ifelse(. == major_base, 1, 0)))


df_check_no_na_major <- replace(df_check_major, is.na(df_check_major), 0)

df_sum_major <- df_check_no_na_major  %>%
  group_by(sample_id2) %>%
  summarise(across(starts_with("HPV16_"), sum))
