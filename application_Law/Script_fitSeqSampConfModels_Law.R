###########################################################################
##### Example of a modeling analysis using the dynConfiR package    #######
###########################################################################

# Sebastian Hellmann, 15.05.2024

# 1) Read in data and remove possible outliers
# 2) Fit several dynamical confidence models
# 3) Visualize Information Criteria for the different models
# 4) Predict rating and rt distribution for the four best-fitting models
# 5) Visualization of observations and predictions in different aggregations

###   Preamble and imports    ####
rm(list = ls())
# # insert path by hand:
# script_path <- paste0("C:/Users/PPA859/Documents", # insert right path, here
#                       "/application_Law")
# or use RStudio to find script file path
script_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_path)

library(BayesFactor) # to compute BFs to check whether participants
                     # have a bias or are at chance level accuracy
library(tidyverse)
library(dynConfiR)   # for fitting and predicting the models
library(ggh4x)       # for nested_wrap (creates better strip labels)
require(ggpubr)      # for arranging the IC-plots

## 1. Read in data and remove possible outliers             ----

Data <-read.delim("data_Law_unpub.csv", sep = ",")
#Data <- read_delim("data_Law_unpub.csv", delim=",", )
Data <- Data %>% 
  filter(!isCalibTrial) %>%
  select(Subj_idx, Stimulus, Response, Confidence, RT_decConf,
                        coh_level, coherence) %>%
  rename_all(tolower) %>%
  rename(participant=subj_idx,
         rt=rt_decconf)
head(Data)
dim(Data)

### We decided to allow for biases:
Data %>% group_by(participant) %>%
  summarise(pstimleft = sum(stimulus==1)/length(response),
            prespleft = sum(response==1)/length(response),
            AboveChance = as.numeric(try(extractBF(proportionBF(y=sum(response==1), N=length(response), p=1/2))$bf))) %>%
  mutate(decBias= ifelse(AboveChance>3, 1, ifelse(AboveChance<.33, -1, 0)))


## Preprocessing:
# Remove very high and low response times (assumed to be due to lapses) 
cutoffs <- Data %>%group_by(participant) %>%
  reframe(maxrt = mean(rt)+4*sd(rt),
          minrt=median(rt)-1*sd(rt))
ggplot(Data, aes(x=rt))+
  geom_density(bw=0.02)+
  geom_vline(data=cutoffs, aes(xintercept=maxrt), color="red")+
  geom_vline(data=cutoffs, aes(xintercept=minrt), color="red")+
  facet_wrap(.~participant, scales = "free")
Nrowges <- nrow(Data)
Data <- Data %>%
  #filter(rt<10) %>%
  group_by(participant) %>%
  filter(rt < mean(rt)+4*sd(rt) & rt >median(rt)-sd(rt)) %>%
  ungroup()
Nrowfiltered <- nrow(Data)
1-Nrowfiltered/Nrowges

## Remove participants performing at chance level
BadParts <- Data %>% group_by(participant) %>%
  summarise(accuracy=sum(response==stimulus)/length(response),
            AboveChance = as.numeric(try(extractBF(proportionBF(y=sum(response==stimulus), N=length(response), p=1/2))$bf))) %>%
  mutate(goodPart= ifelse(AboveChance>3, 1, ifelse(AboveChance<.33, -1, 0)))
BadParts
## --> remove participant 9
BadParts[BadParts$goodPart==-1,]
subset(Data, participant%in% BadParts[BadParts$goodPart==-1,"participant"]) %>% 
  group_by(participant, coh_level) %>%
  summarise(acc = sum(response==stimulus)/n(),
            coh = mean(coherence))
Data <- filter(Data, !participant%in% BadParts[BadParts$goodPart==-1,"participant"])


## 2. Fit the four models to the three experimental data   ----

### a) Start a new model fit (may take some days...)
### Rather skip to b)...

# parallel::detectCores()
# parfits <- fitRTConfModels(Data, models=c("dynaViTE", "dynWEV",
#                                           "2DSD", "PCRMt", "IRMt"), nRatings=4,
#                 restr_tau = "simult_conf", logging = TRUE,
#                 opts=list(nAttempts=4, nRestarts=4), parallel="both", n.cores = c(5, 4),
#                 condition = "coh_level", rating="confidence")
# save(file="parfits.RData", parfits)

### b) Load already fitted parameters         
load("parfits.RData")

#__________________________________________________________----

# 3. Visualize Information Criteria for the different models  
source("IC_analysis_1.R")


# 4) Predict rating and rt distribution for the models

# # ##### Only run, when no precitions had been computed before     ######
# prediction_ConfDist <- predictConfModels(parfits, simult_conf=TRUE, parallel = TRUE)
# prediction_RTConfDist <- predictRTModels(parfits, simult_conf=TRUE, maxrt=6,
#                                          scaled=TRUE, DistConf = prediction_ConfDist,
#                                          parallel=TRUE)
# 
# save(file="fitsandpredictions_1.RData",
#      Data, parfits, prediction_ConfDist, prediction_RTConfDist)
load("fitsandpredictions_1.RData")
options(digits=2)
print(head(prediction_ConfDist), row.names=FALSE)
# Sanity check: Sum of probabiltities should be 1 per stimulus & condition
sum(prediction_ConfDist$p)  
length(unique(Data$participant))*5*2*5 # 15 sbjs x 5 conditions x 2 stimuli x 5 models
# 2(stimulus)*5(conditions)*15(participants)*4(models) = 600


print(head(prediction_RTConfDist), row.names=FALSE)


# 5) Visualization of observations and predictions in different aggregations
Data <- Data %>% mutate(condition = coh_level, 
                        rating = confidence, 
                        correct= as.numeric(stimulus==response))

#### Discrete rating distribution 
source("plotscripts/plot_rating_dist_1.R")

#### Accuracy and mean confidence
source("plotscripts/plot_accuracies_meanconf_1.R")
source("plotscripts/plot_meanconf_by_response_1.R")


#### Response time quantiles against confidence ratings and coherence level

### Aggregate predicted response time distribution across participants 
# (We compute the quantiles over the pooles empirical RT distributions for all
# observations over all participants. To get an equivalent compuation of the predicted 
# RT quantiles, we have to form a weighted sum of the predicted RT density across 
# participants with the weight equal to the number of trials for the respective 
# condition and participant.)
Ns_part <- Data %>% group_by(participant) %>% 
  summarise(N=n(), MinRT = min(rt), .groups = "drop")  %>%
  select(participant, N)
model_levels <-  c("PCRMt", "IRMt", "dynaViTE", "dynWEV", "2DSD")
Preds_RTdens_corr_cond_rating <-  prediction_RTConfDist %>% 
  left_join(Ns_part, by="participant") %>%
  ungroup() %>%
  select(-participant) %>%
  group_by(rating, condition, model, correct, rt) %>%
  summarise(dens = sum(dens*N)/nrow(Data), densscaled = sum(N*densscaled)/nrow(Data), .groups = "drop") %>%
  mutate(model=factor(model, levels=model_levels, ordered=TRUE))
source("plotscripts/plot_RTQuants_conf_1.R")
source("plotscripts/plot_RTQuants_cond_1.R")





#####################################################################
#####################################################################
##########   Exploration with fixed parameters       ################
#####################################################################
#####################################################################

# ### a) Start a new model fit (may take some days...)
# ### Rather skip to b)...
# dir.create("fitting_fixed_noiseparameters")
# setwd("fitting_fixed_noiseparameters")
# parfits_fixed <- fitRTConfModels(Data, models=c("dynaViTE", "dynWEV", "2DSD"),
#                                           nRatings=4,
#                            restr_tau = "simult_conf", logging = TRUE,
#                            fixed = list(sz=0, st0=0, svis=1),
#                            opts=list(nAttempts=4, nRestarts=4), parallel="both", n.cores = c(5, 4),
#                            condition = "coh_level", rating="confidence")
# 
# setwd("..")
# parfits_fixed[, setdiff(names(parfits), names(parfits_fixed))] <- NA
# parfits_fixed <- mutate(parfits_fixed, model = paste0(model, " (fixed)"))
# save(file="parfits_fixed_and_free.RData", parfits, parfits_fixed)

### b) Load already fitted parameters         
load("parfits_fixed_and_free.RData")

### c) Combine the fits for the different models  

allfits <- rbind(parfits, parfits_fixed)
parfits$model
parfits_fixed$model
allfits$model

#__________________________________________________________----
#__________________________________________________________----
unique(prediction_ConfDist$model)


# 3. Visualize Information Criteria for the different models  
source("IC_analysis_2.R")


# 4) Predict rating and rt distribution for the new models
# ##### Only run, when no precitions had been computed before     ######
# prediction_ConfDist_fixed <- predictConfModels(parfits_fixed, simult_conf=TRUE, parallel = TRUE)
# prediction_RTConfDist_fixed <- predictRTModels(parfits_fixed, simult_conf=TRUE, maxrt=6,
#                                          scaled=TRUE, DistConf = prediction_ConfDist_fixed,
#                                          parallel=TRUE)
# save(file="fitsandpredictions_allfits_2.RData",
#      Data, allfits, parfits, parfits_fixed,
#      prediction_ConfDist, prediction_RTConfDist,
#      prediction_ConfDist_fixed, prediction_RTConfDist_fixed)
load("fitsandpredictions_allfits_2.RData")

head(prediction_ConfDist_fixed)
# Sanity check: Sum of probabiltities should be 1 per stimulus & condition
sum(prediction_ConfDist_fixed$p)  # 15 sbjs x 2 stimuli x 5 conditions x 3 models
# 2(stimulus)*5(conditions)*15(participants)*4(models) = 600

head(prediction_RTConfDist)


# 5) Visualization of observations and predictions in different aggregations
#### Discrete rating distribution 
source("plotscripts/plot_rating_dist_2.R")

#### Accuracy and mean confidence
source("plotscripts/plot_accuracies_meanconf_2.R")

#### Response time quantiles against confidence ratings and coherence level
### Aggregate predicted response time distribution across participants 
# (We compute the quantiles over the pooles empirical RT distributions for all
# observations over all participants. To get an equivalent compuation of the predicted 
# RT quantiles, we have to form a weighted sum of the predicted RT density across 
# participants with the weight equal to the number of trials for the respective 
# condition and participant.)
Ns_part <- Data %>% group_by(participant) %>% 
  summarise(N=n(), MinRT = min(rt), .groups = "drop")  %>%
  select(participant, N)
model_levels_2 <- c("dynaViTE\n(fixed)","dynaViTE","dynWEV\n(fixed)","dynWEV", "2DSD\n(fixed)", "2DSD")
model_colors_2<-c(  "#F0E442",      "#D55E00" ,     "#56B4E9",    "#0072B2","#000000",  "#CC79A7")
Preds_RTdens_corr_cond_rating <-  
  rbind(prediction_RTConfDist, prediction_RTConfDist_fixed) %>%
  filter(!grepl("RM", model)) %>%
  mutate(model= str_replace(model, " ", "\n")) %>%
  mutate(model=factor(model, levels=model_levels_2, ordered = TRUE)) %>%
  left_join(Ns_part, by="participant") %>%
  ungroup() %>%
  select(-participant) %>%
  group_by(rating, condition, model, correct, rt) %>%
  summarise(dens = sum(dens*N)/nrow(Data), densscaled = sum(N*densscaled)/nrow(Data), .groups = "drop")
source("plotscripts/plot_RTQuants_conf_2.R")
source("plotscripts/plot_RTQuants_cond_2.R")




