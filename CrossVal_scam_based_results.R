

library(tidyverse)
library(reshape2)
library(mgcv)
library(gridExtra)
library(openxlsx)
library(datetimeutils)
library(mapdata)
library(lubridate)
library(marmap)  
library(ggeffects)
library(corrplot)
library(plotmo)
library(stringr)
library(sjPlot)
library(lattice)
library(scam)
library(mboost)
library(gamboostLSS)

###selection methods
library(CAST)
library(caret)
library(subselect)


load("~/hake_data_SeaWise.Rdata")
source("~/function_scam_CrossVal_selection.R")

### remove NAs
dat = hake_data_SeaWise %>% dplyr::ungroup() %>% dplyr::filter(!is.na(Phyc_surface),
                                                               !is.na(roughness),
                                                               !is.na(S_mean_100m),
                                                               !is.na(T_surface),
                                                               !is.na(slope))
dat$survey = as.factor(dat$Survey.y)
dat = dat %>% filter(!is.na(HaulDur),HaulDur>0)


###############################
### run cross validation ####
##########################

### JUVENILES 
dat$y = dat$occurrence_juv
basef = y ~ survey #s(depth, k = 10, bs = "cv", m = 2, sp = 1e-04) + survey 
variables = c("T_deep","Mix_layer","O2_deep","Phyc_mean_100m","S_mean_100m","Nppv_mean_100m") ## provisional variables

k_cross = 10
p_train = 0.8
t.test.alpha = 0.05
N_test_positive=2
k=8

CrossVal.scam(basef=basef,vars = variables,k=k,dat = dat,k_cross=k_cross,
              p_train = p_train,t.test.alpha=t.test.alpha,sp=0.00001,N_test_positive=N_test_positive)



## run final model ## 
dat$y = dat$occurrence_juv

## only environ variables
formu_juv_occur_Env =y ~ survey + s(Phyc_mean_100m, m = 2, bs = "cv", k = 8) + 
  s(S_mean_100m, m = 2, bs = "cv", k = 8) + s(T_deep, m = 2, bs = "cv", k = 8) + 
  s(Nppv_mean_100m, m = 2, bs = "cv", k = 8) + s(Mix_layer, m = 2, bs = "cv", k = 8) + 
  s(O2_deep, m = 2, bs = "cv", k = 8)
model_juv_occur_Env <- scam(formu_juv_occur_Env, family=binomial(link="logit"), data=dat,sp=rep(0.00001,3))
summary(model_juv_occur_Env)
plotmo(model_juv_occur_Env)



### ADULTS

dat$y = dat$occurrence_adu
basef = y ~ survey # s(depth, k = 10, bs = "cv", m = 2, sp = 1e-04) + survey 
variables = c("T_deep","Mix_layer", "O2_deep","Phyc_mean_100m","S_mean_100m","Nppv_mean_100m")  ## provisional variables

k_cross = 10
p_train = 0.8
t.test.alpha = 0.05
N_test_positive=2
k=8

CrossVal.scam(basef=basef,vars = variables,k=k,dat = dat,k_cross=k_cross,
              p_train = p_train,t.test.alpha=t.test.alpha,sp=0.00001,N_test_positive=N_test_positive)





## run final model ## 
dat$y = dat$occurrence_adu

## only environ variables
formu_adu_occur_Env = y ~ survey + s(T_deep, m = 2, bs = "cv", k = 8) + 
  s(O2_deep,m = 2, bs = "cv", k = 8) + s(Phyc_mean_100m, m = 2, bs = "cv", k = 8)
model_adu_occur_Env <- scam(formu_adu_occur_Env, family=binomial(link="logit"), data=dat,sp=rep(0.00001,3))
summary(model_adu_occur_Env)
plotmo(model_adu_occur_Env)




