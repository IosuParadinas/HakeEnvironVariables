

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


load("C:/Users/jparadinas/Dropbox/AZTI/SeaWise/HAKE/data/hake_data_SeaWise.Rdata")

### remove NAs
dat = hake_data_SeaWise %>% dplyr::ungroup() %>% dplyr::filter(!is.na(Phyc_surface),
                                                               !is.na(roughness),
                                                               !is.na(S_mean_100m),
                                                               !is.na(T_surface),
                                                               !is.na(slope))
dat$survey = as.factor(dat$Survey.y)
dat = dat %>% filter(!is.na(HaulDur),HaulDur>0)


################
### example ###
##############

source("C:/Users/ip30/Dropbox/AZTI/SeaWise/SC_GAM/function_scam_CrossVal_selection.R")

dat$y = dat$occurrence_adu

basef = y ~ survey #s(depth, k = 10, bs = "cv", m = 2, sp = 1e-04) + survey 

variables = c("T_deep","Mix_layer",
              "O2_deep","Phyc_mean_100m","S_mean_100m","Nppv_mean_100m")


k_cross = 10
p_train = 0.8
t.test.alpha = 0.05
N_test_positive=2
k=8

CrossVal.scam(basef=basef,vars = variables,k=k,dat = dat,k_cross=k_cross,
              p_train = p_train,t.test.alpha=t.test.alpha,sp=0.00001,N_test_positive=N_test_positive)


##########################
### cross validation ####
#######################

###################
### JUVENILES ###
###############

#### occurrence 

## only environ variables

formu_juv_occur_Env =y ~ survey + s(Phyc_mean_100m, m = 2, bs = "cv", k = 8) + 
  s(S_mean_100m, m = 2, bs = "cv", k = 8) + s(T_deep, m = 2, bs = "cv", k = 8) + 
  s(Nppv_mean_100m, m = 2, bs = "cv", k = 8) + s(Mix_layer, m = 2, bs = "cv", k = 8) + 
  s(O2_deep, m = 2, bs = "cv", k = 8)
model_juv_occur_Env <- scam(formu_juv_occur_Env, family=binomial(link="logit"), data=dat,sp=rep(0.00001,3))
summary(model_juv_occur_Env)
plotmo(model_juv_occur_Env)

## inlcuding depth + survey
dat$y = dat$occurrence_juv

formu_juv_occur = y ~ survey+s(depth,m=2,bs='cv',k=10)+s(sun_alt,m=2,bs='cv',k=10)+s(O2_deep,m=2,bs='cv',k=10)
model_juv_occur <- scam(formu_juv_occur, family=binomial(link="logit"), data=dat,sp=rep(0.00001,3))
summary(model_juv_occur)
plotmo(model_juv_occur)

### abundance
# 

################
### ADULTS ###
############

#### occurrence 
dat$y = dat$occurrence_adu

## only environ variables
formu_adu_occur_Env = y ~ survey + s(T_deep, m = 2, bs = "cv", k = 8) + 
  s(O2_deep,m = 2, bs = "cv", k = 8) + s(Phyc_mean_100m, m = 2, bs = "cv", k = 8)
model_adu_occur_Env <- scam(formu_adu_occur_Env, family=binomial(link="logit"), data=dat,sp=rep(0.00001,3))
summary(model_adu_occur_Env)
plotmo(model_adu_occur_Env)

## inlcuding depth + survey
dat$y = dat$occurrence_adu
formu_adu_occur = y ~ s(O2_deep, k = 10, bs = "cv", m = 2, sp = 1e-04) + survey 
model_adu_occur <- scam(formu_adu_occur, family=binomial(link="logit"), data=dat,sp=rep(0.00001,3))
summary(model_adu_occur)
plotmo(model_adu_occur)


# $smod$null
# y ~ survey
# 
# $smod$depth
# [1] "y ~ survey+s(depth,m=2,bs='cv',k=8)"
# 
# $smod$T_mean_100m
# [1] "y ~ survey+s(depth,m=2,bs='cv',k=8)+s(T_mean_100m,m=2,bs='cv',k=8)"
# 
# $smod$sun_alt
# [1] "y ~ survey+s(depth,m=2,bs='cv',k=8)+s(T_mean_100m,m=2,bs='cv',k=8)+s(sun_alt,m=2,bs='cv',k=8)"
# 
# $smod$O2_deep
# [1] "y ~ survey+s(depth,m=2,bs='cv',k=8)+s(T_mean_100m,m=2,bs='cv',k=8)+s(sun_alt,m=2,bs='cv',k=8)+s(O2_deep,m=2,bs='cv',k=8)"
# 
# $smod$roughness
# [1] "y ~ survey+s(depth,m=2,bs='cv',k=8)+s(T_mean_100m,m=2,bs='cv',k=8)+s(sun_alt,m=2,bs='cv',k=8)+s(O2_deep,m=2,bs='cv',k=8)+s(roughness,m=2,bs='cv',k=8)"
# 
# 
# $svars
# [1] "depth"       "T_mean_100m" "sun_alt"     "O2_deep"     "roughness"  
# 
# $final_model
# 
# Family: binomial 
# Link function: logit 
# 
# Formula:
#   y ~ survey + s(depth, m = 2, bs = "cv", k = 8) + s(T_mean_100m, 
#                                                      m = 2, bs = "cv", k = 8) + s(sun_alt, m = 2, bs = "cv", k = 8) + 
#   s(O2_deep, m = 2, bs = "cv", k = 8) + s(roughness, m = 2, 
#                                           bs = "cv", k = 8)
# <environment: 0x000001a81dfddf18>
#   
#   Estimated degrees of freedom:
#   0 0 6 0 0  total = 6 
# 
# UBRE score: 10.189
# rank: 6/45




### abundance
# 

